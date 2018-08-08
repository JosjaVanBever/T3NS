#include <stdlib.h>
#include <stdio.h>

#include "tensorproducts.h"
#include "macros.h"
#include "debug.h"
#include "hamiltonian.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

/* Builds a naive sectors list, just a direct product of the different ranges possible from the 
 * other symmsectors */
static void build_all_sectors( struct symsecs * const res, const struct symsecs * const sectors1, 
    const struct symsecs * const sectors2 );

/* ============================================================================================ */

void find_goodqnumbersectors( int ****dimarray, int ****qnumbersarray, int *total, 
    const struct symsecs symarr[] )
{
  /** Loop over bond 1 and 2, tensorproduct them to form bond 3 and then look at the ones that 
   * actually exist in bond 3. First do it for the first resulting symsec,
   * after that, do it for all the rest.
   */
  int sym1, sym2, i;
  int prevsym[ 2 ][ bookie.nr_symmetries ];
  int min_irrep[ bookie.nr_symmetries ];
  int nr_irreps[ bookie.nr_symmetries ];
  int step     [ bookie.nr_symmetries ];
  int max_irrep[ bookie.nr_symmetries ];

  *dimarray      = safe_malloc( symarr[ 0 ].nr_symsec, int** );
  *qnumbersarray = safe_malloc( symarr[ 0 ].nr_symsec, int** );
  *total = 0;

  for( i = 0 ; i < bookie.nr_symmetries ; i++ )
  {
    prevsym[ 0 ][ i ] = symarr[ 0 ].irreps[ 0 * bookie.nr_symmetries + i ];
    prevsym[ 1 ][ i ] = symarr[ 1 ].irreps[ 0 * bookie.nr_symmetries + i ];
    tensprod_irrep( &min_irrep[ i ], &nr_irreps[ i ], &step[ i ], prevsym[ 0 ][ i ],
        prevsym[ 1 ][ i ], 1, bookie.sgs[ i ] );
    max_irrep[ i ] = min_irrep[ i ] + step[ i ] * ( nr_irreps[ i ] - 1 );
  }

  for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; sym1++ )
  {
    (*dimarray)[ sym1 ]      = NULL;
    (*qnumbersarray)[ sym1 ] = NULL;
    if( symarr[ 0 ].dims[ sym1 ] == 0 )
      continue;

    (*dimarray)[ sym1 ]      = safe_malloc( symarr[ 1 ].nr_symsec, int* );
    (*qnumbersarray)[ sym1 ] = safe_malloc( symarr[ 1 ].nr_symsec, int* );

    for( i = 0 ; i < bookie.nr_symmetries ;i++ )
    {
      if( symarr[ 0 ].irreps[ sym1 * bookie.nr_symmetries + i ] != prevsym[ 0 ][ i ] )
      {
        prevsym[ 0 ][ i ] = symarr[ 0 ].irreps[ sym1 * bookie.nr_symmetries + i ];
        tensprod_irrep( &min_irrep[ i ], &nr_irreps[ i ], &step[ i ], prevsym[ 0 ][ i ],
            prevsym[ 1 ][ i ], 1, bookie.sgs[ i ] );
        max_irrep[ i ] = min_irrep[ i ] + step[ i ] * ( nr_irreps[ i ] - 1 );
      }
    }

    for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; sym2++ )
    {
      int irrep[ bookie.nr_symmetries ];
      int dim         = symarr[ 0 ].dims[ sym1 ] * symarr[ 1 ].dims[ sym2 ];
      int totalirreps = 1;
      int count       = -1;
      int curr        = 0;
      (*dimarray)[ sym1 ][ sym2 ]      = NULL;
      (*qnumbersarray)[ sym1 ][ sym2 ] = NULL;
      if( symarr[ 1 ].dims[ sym2 ] == 0 )
        continue;

      for( i = 0 ; i < bookie.nr_symmetries ;i++ )
      {
        if( symarr[ 1 ].irreps[ sym2 * bookie.nr_symmetries + i ] != prevsym[ 1 ][ i ] )
        {
          prevsym[ 1 ][ i ] = symarr[ 1 ].irreps[ sym2 * bookie.nr_symmetries + i ];
          tensprod_irrep( &min_irrep[ i ], &nr_irreps[ i ], &step[ i ], prevsym[ 0 ][ i ],
              prevsym[ 1 ][ i ], 1, bookie.sgs[ i ] );
          max_irrep[ i ] = min_irrep[ i ] + step[ i ] * ( nr_irreps[ i ] - 1 );
        }

        irrep[ i ] = min_irrep[ i ];
        totalirreps *= nr_irreps[ i ];
      }

      (*dimarray)[ sym1 ][ sym2 ]           = safe_malloc( totalirreps, int );
      (*qnumbersarray)[ sym1 ][ sym2 ]      = safe_malloc( 1 + totalirreps, int );
      (*qnumbersarray)[ sym1 ][ sym2 ][ 0 ] = totalirreps;

      while( ++count < totalirreps )
      {
        int ind = search_symmsec( irrep, &symarr[ 2 ] );
        if( ind != -1 && symarr[ 2 ].dims[ ind ] )
        {
          (*total)++;
          (*dimarray)[ sym1 ][ sym2 ][ curr ] = dim * symarr[ 2 ].dims[ ind ];
          (*qnumbersarray)[ sym1 ][ sym2 ][ curr + 1 ] = ind;
          ++curr;
        }

        for( i = 0 ; i < bookie.nr_symmetries ; i++ )
        {
          if( ( irrep[ i ] += step[ i ] ) > max_irrep[ i ] )
            irrep[ i ] = min_irrep[ i ];
          else
            break;
        }
      }
      assert( i == bookie.nr_symmetries );
      assert( irrep[ bookie.nr_symmetries - 1 ]  == min_irrep[ bookie.nr_symmetries - 1 ] );

      if( curr == 0 )
        safe_free( (*dimarray)[ sym1 ][ sym2 ] );
      else
        (*dimarray)[ sym1 ][ sym2 ] = realloc( (*dimarray)[ sym1 ][ sym2 ],  curr * sizeof( int ));

      (*qnumbersarray)[ sym1 ][ sym2 ] 
        = realloc( (*qnumbersarray)[ sym1 ][ sym2 ],  ( curr + 1 ) * sizeof( int ));
      (*qnumbersarray)[ sym1 ][ sym2 ][ 0 ] = curr;

      if( ( curr != 0 && (*dimarray)[ sym1 ][ sym2 ] == NULL ) ||
          (*qnumbersarray)[ sym1 ][ sym2 ] == NULL )
      {
        fprintf( stderr, "%s@%s: realloc failed: curr = %d\n", __FILE__, __func__, curr );
        exit( EXIT_FAILURE );
      }
    }
  }
}

void destroy_dim_and_qnumbersarray( int ****dimarray, int ****qnumbersarray, const struct symsecs 
    symarr[] )
{
  int i;

  for( i = 0 ; i < symarr[ 0 ].nr_symsec ; ++i )
  {
    int j;

    for( j = 0 ; j < symarr[ 1 ].nr_symsec ; ++j )
    {
      safe_free( (*dimarray)[ i ][ j ] );
      safe_free( (*qnumbersarray)[ i ][ j ] );
    }
    safe_free( (*dimarray)[ i ] );
    safe_free( (*qnumbersarray)[ i ] );
  }
  safe_free( *dimarray );
  safe_free( *qnumbersarray );
}

void find_qnumbers_with_index_in_array( const int id, const int idnr, int *** qnumbersarray, 
    int ***dimarray, const struct symsecs symarr[], QN_TYPE **res_qnumbers, int ** res_dim, 
    int *length )
{
  int sym1, sym2, sym3;
  int counter;
  const int dim1 = symarr[ 0 ].nr_symsec;
  const int dim2 = symarr[ 0 ].nr_symsec * symarr[ 1 ].nr_symsec;
  const int NULLflag = ( res_qnumbers == NULL ) && ( res_dim == NULL );

  assert( ( res_qnumbers == NULL ) == ( res_dim == NULL ) &&
      "Both res_qnumbers and res_dim should be null-pointers or valid pointers" );

  *length = 0;
  switch ( idnr )
  {
    case 0:
      assert( id < symarr[ 0 ].nr_symsec );
      sym1 = id;
      for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; ++sym2 )
      {
        *length += qnumbersarray[ sym1 ][ sym2 ][ 0 ];
      }
      if( NULLflag )
        break;

      *res_qnumbers = safe_malloc( *length, QN_TYPE );
      *res_dim      = safe_malloc( *length, int );

      counter = 0;
      for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; ++sym2 )
      {
        QN_TYPE ind = sym1 + dim1 * sym2;
        for( sym3 = 0 ; sym3 < qnumbersarray[ sym1 ][ sym2 ][ 0 ] ; ++sym3 )
        {
          (*res_qnumbers)[ counter ] = ind + qnumbersarray[ sym1 ][ sym2 ][ 1 + sym3 ] * dim2;
          (*res_dim)[ counter ]      = dimarray[ sym1 ][ sym2 ][ sym3 ];
          ++counter;
        }
      }
      break;
    case 1:
      assert( id < symarr[ 1 ].nr_symsec );
      sym2 = id;
      for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; ++sym1 )
      {
        *length += qnumbersarray[ sym1 ][ sym2 ][ 0 ];
      }
      if( NULLflag )
        break;

      *res_qnumbers = safe_malloc( *length, QN_TYPE );
      *res_dim      = safe_malloc( *length, int );

      counter = 0;
      for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; ++sym1 )
      {
        QN_TYPE ind = sym1 + dim1 * sym2;
        for( sym3 = 0 ; sym3 < qnumbersarray[ sym1 ][ sym2 ][ 0 ] ; ++sym3 )
        {
          (*res_qnumbers)[ counter ] = ind + qnumbersarray[ sym1 ][ sym2 ][ 1 + sym3 ] * dim2;
          (*res_dim)[ counter ]      = dimarray[ sym1 ][ sym2 ][ sym3 ];
          ++counter;
        }
      }
      break;
    case 2:
      assert( id < symarr[ 2 ].nr_symsec );

      for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; ++sym1 )
        for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; ++sym2 )
          for( sym3 = 0 ; sym3 < qnumbersarray[ sym1 ][ sym2 ][ 0 ] ; ++sym3 )
          {
            *length += qnumbersarray[ sym1 ][ sym2 ][ 1 + sym3 ] == id;
          }
      if( NULLflag )
        break;

      *res_qnumbers = safe_malloc( *length, QN_TYPE );
      *res_dim      = safe_malloc( *length, int );

      counter = 0;
      for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; ++sym1 )
        for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; ++sym2 ){
          QN_TYPE ind = sym1 + dim1 * sym2;
          for( sym3 = 0 ; sym3 < qnumbersarray[ sym1 ][ sym2 ][ 0 ] ; ++sym3 )
            if( qnumbersarray[ sym1 ][ sym2 ][ 1 + sym3 ] == id )
            {
              (*res_qnumbers)[ counter ] = ind + qnumbersarray[ sym1 ][ sym2 ][ 1 + sym3 ] * dim2;
              (*res_dim)[ counter ]      = dimarray[ sym1 ][ sym2 ][ sym3 ];
              ++counter;
            }
        }
      assert( counter == *length );
      break;
    default:
      fprintf( stderr, "ERROR: Wrong idnr (%d) passed in find_qnumbers_with_index_in_array\n",idnr);
      exit( EXIT_FAILURE );
  }
}

void tensprod_symsecs( struct symsecs * const res, const struct symsecs * const sectors1, 
    const struct symsecs * const sectors2, const int sign, const char o )
{
  /* First I make a 'worst-case' res, where a rough first guess of the symmsecs that will occur 
   * is initialized in. After this, this 'worst-case' res can be simplified by kicking out all 
   * the symmsecs with D = 0.
   */
  int i;
  assert( o == 'f' || o == 'n' );

  res->nr_symsec = 0;
  res->irreps    = NULL;
  res->fcidims   = NULL;
  res->dims      = NULL;
  res->totaldims = 0;

  /* This function will give a rough first irreps array with also a lot of forbidden symmsecs. */
  build_all_sectors( res, sectors1, sectors2 );
  res->fcidims = safe_calloc( res->nr_symsec,  double );
  if( o == 'n' ) res->dims = safe_calloc( res->nr_symsec,  int );

  for( i = 0 ; i < sectors1->nr_symsec ; i++ )
  {
    int j;
    /* zero dimension symmsec */
    if( ( o == 'f' && sectors1->fcidims[ i ] < 0.5 ) || ( o == 'n' && sectors1->dims[ i ] == 0 ) )
      continue;

    for( j = 0 ; j < sectors2->nr_symsec ; j++ )
    {
      int nr_symmsecs;
      int *resultsymmsec;
      if( ( o == 'f' && sectors2->fcidims[ j ] < 0.5 ) || ( o == 'n' && sectors2->dims[ j ] == 0 ) )
        continue;

      /* for non-abelian symmetries, like SU(2), there are multiple irreps that are valid as
       * result of the tensorproduct of two irreps */
      tensprod_symmsec( &resultsymmsec, &nr_symmsecs, &sectors1->irreps[ i * bookie.nr_symmetries ], 
                        &sectors2->irreps[ j * bookie.nr_symmetries ], sign, bookie.sgs,
                        bookie.nr_symmetries );

      for( nr_symmsecs-- ; nr_symmsecs >= 0 ; nr_symmsecs-- )
      {
        int pos_symmsec = search_symmsec( resultsymmsec + bookie.nr_symmetries * nr_symmsecs, res );
        if( pos_symmsec < 0 )
          break;
        if( o == 'f' )
          res->fcidims[ pos_symmsec ] += sectors1->fcidims[ i ] * sectors2->fcidims[ j ];
        if( o == 'n' )
        {
          res->fcidims[ pos_symmsec ] = 1;
          res->dims[ pos_symmsec ] = 1;
        }
      }
      safe_free( resultsymmsec );
    }
  }

  /* now we have the 'worst-case' res. Kick out all the symmsecs with dimension 0. */
  kick_empty_symsecs( res, o );
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void build_all_sectors( struct symsecs * const res, const struct symsecs * const sectors1, 
    const struct symsecs * const sectors2 )
{
  int max_irrep[ bookie.nr_symmetries ];
  int indices[ bookie.nr_symmetries ];
  int i;
  int cnt;
  int cnt2;
  int nr_symsec = 1;
  res->irreps    = NULL;
  res->dims      = NULL;
  res->fcidims   = NULL;
  res->totaldims = 0;

  for( i = 0 ; i < bookie.nr_symmetries ; i++ )
  {
    indices[ i ] = 0;
    max_irrep[ i ] = get_max_irrep( sectors1->irreps + i, sectors1->nr_symsec, sectors2->irreps + i, 
        sectors2->nr_symsec, bookie.nr_symmetries, bookie.sgs[ i ] );

    nr_symsec *= max_irrep[ i ];
  }

  res->nr_symsec = 0;
  cnt = 0;
  while( cnt != nr_symsec )
  {
    res->nr_symsec += consistent_state( bookie.sgs, indices, bookie.nr_symmetries );

    for( i = 0 ; i < bookie.nr_symmetries ; i++ )
    {
      indices[ i ]++;
      if( indices[ i ] == max_irrep[ i ] )
        indices[ i ] = 0;
      else
        break;
    }
    cnt++;
  }
  res->irreps = safe_malloc( res->nr_symsec * bookie.nr_symmetries, int );

  cnt = 0;
  cnt2 = 0;
  while( cnt != nr_symsec )
  {
    if( consistent_state( bookie.sgs, indices, bookie.nr_symmetries ) )
    {
      for( i = 0 ; i < bookie.nr_symmetries ; i++ )
        res->irreps[ cnt2 * bookie.nr_symmetries + i ] = indices[ i ];
      cnt2++;
    }

    for( i = 0 ; i < bookie.nr_symmetries ; i++ )
    {
      indices[ i ]++;
      if( indices[ i ] == max_irrep[ i ] )
        indices[ i ] = 0;
      else
        break;
    }
    cnt++;
  }
  assert( ( i == bookie.nr_symmetries ) && ( indices[ i - 1 ] == 0 ) && "Not all symmsecs looped" );
}