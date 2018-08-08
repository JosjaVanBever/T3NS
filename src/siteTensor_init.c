#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "bookkeeper.h"
#include "network.h"
#include "lapack.h"
#include "tensorproducts.h"
#include "siteTensor.h"
#include "macros.h"
#include "debug.h"
#include "sort.h"

/* ============================================================================================== */
/* ================================= DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================== */

/* Makes the blocks out of the dimarray and qnumbersarray returned by find_goodqnumbersectors */
static void make_1sblocks( struct siteTensor * const tens, int ***dimarray, int ***qnumbersarray, 
    const struct symsecs symarr[] );

/* This makes the new internal symsecs and the qnumbers array and the beginblock in the 
 * sparseblock structure, no internal symsecs are made for tens->nrsites = 1  */
static void make_new_internalsymsecs_and_tensor( struct siteTensor * const tens, 
    struct symsecs internalsymsec[] );

/* contracts the correct tensor objects in T3NS to a new big site and destroys them */
static void contractsiteTensors( struct siteTensor * const tens, struct siteTensor * const T3NS, 
    struct symsecs internalsymsec[] );

/* changes the passed internal symsecs linked to passed siteTensor in the bookkeeper */
static void change_internals_in_bookkeeper( struct symsecs internalsymsec[], 
    struct siteTensor * const tens );

static void make_multisitetensor( struct siteTensor * tens, const struct symsecs internalsymsec[], 
    const int internalbonds[], const int nr_internal );

static void get_ss_with_internal( struct symsecs symsec[], const struct symsecs internalsymsec[],
    const int interalbonds[], const int nr_internal, const int bonds[] );

static void destroy_ss_with_internal( struct symsecs symsec[], const int internalbonds[], 
    const int nr_internal, const int bonds[] );

static void get_sites_to_use( const struct siteTensor * const tens, const int internalbonds[], 
    int sites_to_use[], const int nr_internal, int * const innersite );

/* Makes the qnumbers array in the tens struct or counts the number of blocks and saves it in the 
 * tens struct. For the first case the dimension of every block is also stored im dim_of_blocks */
static void make_qnumbers_and_dims( int ***qnumbersarray[], int ***dimarray[], struct siteTensor * 
    const tens, int ** const dim_of_blocks, const struct symsecs all_symarr[], 
    const int all_bonds[], const int internalbonds[], const int nr_internal );

/* Sorts the qnumbers array in tens, and makes the appropriate sparseBlocks struct.
 * Also destroys the passed dim_of_blocks array that was passed. */
static void sort_and_make( struct siteTensor * const tens, int ** const dim_of_blocks );

static QN_TYPE change_newtooldqnumber( QN_TYPE new, int * newtoold, int * maxdims, const int newdim,
    const int map );

static int * make_newtoold( const struct symsecs * const internalss, const int bond );

/* =================================== INIT & DESTROY ========================================== */

void init_null_siteTensor( struct siteTensor * const tens )
{
  tens->nrsites  = 0;
  tens->sites    = NULL;
  tens->nrblocks = 0;
  tens->qnumbers = NULL;
  init_null_sparseblocks( &tens->blocks );
}

void init_1siteTensor( struct siteTensor * const tens, const int site, const char o )
{
  /* One-site is the only type of siteTensor I should make out of thin air */
  int i;
  int ***dimarray      = NULL;
  int ***qnumbersarray = NULL;
  const int nrind = 3;
  int couplings[ nrind ];
  struct symsecs symarr[ nrind ];

  tens->nrsites    = 1;
  tens->sites      = safe_malloc( tens->nrsites, int );
  tens->sites[ 0 ] = site;

  siteTensor_give_couplings( tens, couplings );

  get_symsecs_arr( symarr, couplings, nrind );
  find_goodqnumbersectors( &dimarray, &qnumbersarray, &tens->nrblocks, symarr );

  make_1sblocks( tens, dimarray, qnumbersarray, symarr );

  /* Clean the symarr array. And destroy appropriate symsecs (e.g. the ones linked to a physical) */
  clean_symsecs_arr( symarr, couplings, nrind );

  /* initialization of the tel array */
  switch( o )
  {
    case 'n':
      tens->blocks.tel = NULL;
      break;
    case 'm':
      tens->blocks.tel = safe_malloc( tens->blocks.beginblock[ tens->nrblocks ], EL_TYPE );
      break;
    case '0':
      tens->blocks.tel = safe_calloc( tens->blocks.beginblock[ tens->nrblocks ], EL_TYPE );
      break;
    case 'r':
      tens->blocks.tel = safe_malloc( tens->blocks.beginblock[ tens->nrblocks ], EL_TYPE );
      srand( time( NULL ) );
      for( i = 0 ; i <  tens->blocks.beginblock[ tens->nrblocks ] ; ++i )
        tens->blocks.tel[ i ] = ( rand() * 1. ) / RAND_MAX;
      break;
    default:
      fprintf( stderr, "%s@%s: Unknown option \'%c\' was inputted.\n", __FILE__, __func__, o );
      exit( EXIT_FAILURE );
  }
}

void destroy_siteTensor( struct siteTensor * const tens )
{
  tens->nrsites = 0;
  safe_free( tens->sites );
  tens->nrblocks = 0;
  safe_free( tens->qnumbers );
  destroy_sparseblocks( &tens->blocks );
}

void makesiteTensor( struct siteTensor * const tens, struct siteTensor * const T3NS, 
    const int sitelist[] )
{
  /**
   * Thus for this, first, I need all to calculate all the internal bonds for this site object.
   * After that, I calculate for all sites, the possible symmetry sectors of the internal bonds of 
   * every site. 
   * Ie. 
   * If I have a tensor with two outer and one inner bond, I calculate all the existing
   * tensorproducts of the symsecs of the outer bonds to the inner bond.
   * If I have a tensor with more than inner bonds, I first calculate the possible inner bonds 
   * through other tensors and then check the possible tensorproducts to occur with this one.
   * So I can have 3 different cases ( exempting onesite optimization at the moment ):
   *
   * twosite optimization:
   *     * Two tenors with two outer bonds and 1 inner bond ( which is a common one ).
   *     => First calculate the possible tensorproducts of the two outer bonds to the inner bond
   *        for both sites. Than compare the two results and only keep the symsecs that occur in
   *        both.
   * threesite optimization:
   *     * two tensors with two outer bonds and 1 inner bond, one tensor with one outer and two 
   *       inner bonds, which are common bonds with the inners of the two other tensors.
   *     => First calculate the possible tensorproducts of the sites with the two outer bonds to
   *        inner bond. Than calculate the possible tensorproduct at the third (middle) site with
   *        the inner symsecs you just obtained. symsecs that can not occur you kick out of the two
   *        inner bonds.
   * foursite optimization:
   *     * three tensors with two outer bonds and 1 inner bond, one tensor with three inner bonds, 
   *     which are common bonds with the inners of the three other tensors.
   *     => First calculate the possible tensorproducts of the sites with the two outer bonds to
   *        inner bond. Than calculate the possible tensorproduct at the fourth (middle) site with
   *        the inner symsecs you just obtained. symsecs that can not occur you kick out of the two
   *        inner bonds.
   * 
   * What would be the optimal site ordering?
   *
   * While doing all these tensorproducts, you should also already calculate the dimensions linked 
   * with the sparse blocks for it.  After that you obtained the inner symsecs, you can create the 
   * siteObject. Make the different qnumbers and their dimensions. and now do the contracts
   * by looping over all the different qnumbers and matching them with the qnumbers of the original
   * tensors and than contracting the blocks.
   *
   * Now you can destroy the original tensors (maybe dont do that in this function).
   * You should also now set the symsecs of the inner bonds to dim=1. The original dimensions can
   * now be forgotten.
   */
  struct symsecs internalsymsec[ 3 ];
  int i;

  for( tens->nrsites = 0 ; tens->nrsites < 4 ; ++tens->nrsites ) 
    if( sitelist[ tens->nrsites ] == -1 ) break;
  assert( tens->nrsites <= 2 && "At this moment only two-siteoptimization" );

  tens->sites = safe_malloc( tens->nrsites, int );
  for( i = 0 ; i < tens->nrsites ; ++i ) tens->sites[ i ] = sitelist[ i ];
  assert( siteTensor_give_nr_internalbonds( tens ) <= 3 );

  /* This makes the new internal symsecs and the qnumbers array and the beginblock in the 
   * sparseblock structure, no internal symsecs are made for tens->nrsites = 1  */
  make_new_internalsymsecs_and_tensor( tens, internalsymsec );
  /* contracts the correct tensor objects in T3NS to a new big site and destroys them */
  contractsiteTensors( tens, T3NS, internalsymsec );
  /* changes the passed internal symsecs linked to passed siteTensor in the bookkeeper */
  
  change_internals_in_bookkeeper( internalsymsec, tens );
}

void deep_copy_siteTensor( struct siteTensor * const copy, const struct siteTensor * const tocopy )
{
  int i;
  copy->nrsites = tocopy->nrsites;
  copy->nrblocks = tocopy->nrblocks;

  copy->sites = safe_malloc( copy->nrsites, int );
  copy->qnumbers = safe_malloc( copy->nrsites * copy->nrblocks, QN_TYPE );
  for( i = 0 ; i < copy->nrsites ; ++i ) copy->sites[i] = tocopy->sites[i];
  for( i = 0 ; i < copy->nrsites * copy->nrblocks ; ++i ) copy->qnumbers[i] = tocopy->qnumbers[i];

  deep_copy_sparseblocks( &copy->blocks, &tocopy->blocks, tocopy->nrblocks );
}

/* ============================================================================================== */
/* ================================== DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================== */

static void make_1sblocks( struct siteTensor * const tens, int ***dimarray, int ***qnumbersarray, 
    const struct symsecs symarr[] )
{
  int sym1, sym2, sym3;
  int cnt = 0;
  int i;
  int *tempdims           = safe_malloc( tens->nrblocks, int );
  int *idx;
  QN_TYPE *tempqnumbers   = safe_malloc( tens->nrblocks, QN_TYPE );
  tens->qnumbers          = safe_malloc( tens->nrblocks, QN_TYPE );
  tens->blocks.beginblock = safe_malloc( tens->nrblocks + 1, int );
  assert( tens->nrsites == 1 && "make_1sblocks not defined for more than 1 site" );

  for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; ++sym1  )
  {
    for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; ++sym2 )
    {
      const QN_TYPE ind = sym1 + sym2 * symarr[ 0 ].nr_symsec;
      const QN_TYPE increment = symarr[ 0 ].nr_symsec * symarr[ 1 ].nr_symsec;

      for( sym3 = 0 ; sym3 < qnumbersarray[ sym1 ][ sym2 ][ 0 ] ; ++sym3 )
      {
        if( dimarray[ sym1 ][ sym2 ][ sym3 ] == 0 )
          continue;

        tempdims[ cnt ]     = dimarray[ sym1 ][ sym2 ][ sym3 ];
        tempqnumbers[ cnt ] = ind + qnumbersarray[ sym1 ][ sym2 ][ sym3 + 1 ] * increment;
        ++cnt;
      }

      safe_free( dimarray[ sym1 ][ sym2 ] );
      safe_free( qnumbersarray[ sym1 ][ sym2 ] );
    }
    safe_free( dimarray[ sym1 ] );
    safe_free( qnumbersarray[ sym1 ] );
  }
  assert( cnt == tens->nrblocks );

  safe_free( dimarray );
  safe_free( qnumbersarray );

  /* Reform leading order, and I could kick this order */
  idx = qnumbersSort( tempqnumbers, siteTensor_give_nr_of_couplings( tens ), tens->nrblocks );

  tens->blocks.beginblock[ 0 ] = 0;
  for( i = 0 ; i < tens->nrblocks ; ++i )
  {
    tens->qnumbers[ i ] = tempqnumbers[ idx[ i ] ];
    tens->blocks.beginblock[ i + 1 ] = tens->blocks.beginblock[ i ] + tempdims[ idx[ i ] ];
  }
  safe_free( tempdims );
  safe_free( tempqnumbers );
  safe_free( idx );
}

static void make_new_internalsymsecs_and_tensor( struct siteTensor * const tens, 
    struct symsecs internalsymsec[] )
{
  /* 
   * Naive and slower way.
   *
   * execute tensprod_symsecs for the formation of the new internal symsecs. 
   * For this I need nr_sites - 1 different sites.
   * The last site, I execute find_goodqnumbersectors and I kick out unused symsecs of the internal
   * symsecs. (by putting them to zero and kick_empty_symsec.)
   * All the internal symsecs are updated. And after this I execute find_goodqnumbersectors for 
   * every site. And match the things.
   */
  const int nr_internal = siteTensor_give_nr_internalbonds( tens );
  int internalbonds[ nr_internal ];
  int sites_to_use[ nr_internal ]; /* These sites you use to make the initial inner symsec */
  int innersite; /* This site you use to refine the innersymsec through find_goodqnumbersectors */
  int i;

  siteTensor_give_internalbonds( tens, internalbonds );
  get_sites_to_use( tens, internalbonds, sites_to_use, nr_internal, &innersite );

  for( i = 0 ; i < nr_internal ; ++i )
  {
    int bonds[ 3 ];
    struct symsecs symsec[ 2 ];
    int cnt = 0;
    int j;
    int sign;

    get_bonds_of_site( sites_to_use[ i ], bonds );
    sign = internalbonds[ i ] == bonds[ 2 ] ? 1 : -1;

    for( j = 0 ; j < 3 ; ++j )
      if( bonds[ j ] != internalbonds[ i ] ) get_symsecs( &symsec[ cnt++ ], bonds[ j ] );
    assert( cnt == 2 );
    tensprod_symsecs( &internalsymsec[ i ], &symsec[ 0 ], &symsec[ 1 ], sign, 'n' );

    cnt = 0;
    for( j = 0 ; j < 3 ; ++j )
      if( bonds[ j ] != internalbonds[ i ] ) clean_symsecs( &symsec[ cnt++ ], bonds[ j ] );
  }

  { /* inner site */
    int ***dimarray;
    int ***qnumbersarray;
    int total;
    int bonds[ 3 ];
    struct symsecs symsec[ 3 ];
    get_bonds_of_site( innersite, bonds );
    get_ss_with_internal( symsec, internalsymsec, internalbonds, nr_internal, bonds );

    find_goodqnumbersectors( &dimarray, &qnumbersarray, &total, symsec );

    for( i = 0 ; i < 3 ; ++i )
    {
      int j;
      for( j = 0 ; j < nr_internal ; ++j )
        if( bonds[ i ] == internalbonds[ j ] ) break;
      if( j != nr_internal )
      {
        int ss;
        for( ss = 0 ; ss < internalsymsec[ j ].nr_symsec ; ++ss )
        {
          int length;
          find_qnumbers_with_index_in_array( ss, i, qnumbersarray, dimarray, symsec, NULL, NULL, 
              &length );
          if( length == 0 )
            internalsymsec[ j ].fcidims[ ss ] = 0;
        }
        kick_empty_symsecs( &internalsymsec[ j ], 'n' );
      }
    }

    destroy_dim_and_qnumbersarray( &dimarray, &qnumbersarray, symsec );
    destroy_ss_with_internal( symsec, internalbonds, nr_internal, bonds );
  }

  make_multisitetensor( tens, internalsymsec, internalbonds, nr_internal );
}

static void make_multisitetensor( struct siteTensor * tens, const struct symsecs internalsymsec[], 
    const int internalbonds[], const int nr_internal )
{
  int ***dimarray[ tens->nrsites ];
  int ***qnumbersarray[ tens->nrsites ];
  int *dim_of_blocks;
  int i;
  struct symsecs all_symarr[ 3 * tens->nrsites ];
  int all_bonds[ 3 * tens->nrsites ];

  for( i = 0 ; i < tens->nrsites ; ++i )
  {
    int * bonds = &all_bonds[ 3 * i ];
    struct symsecs * symarr = &all_symarr[ 3 * i ];

    int total;
    get_bonds_of_site( tens->sites[ i ], bonds );
    get_ss_with_internal( symarr, internalsymsec, internalbonds, nr_internal, bonds );

    find_goodqnumbersectors( &dimarray[ i ], &qnumbersarray[ i ], &total, symarr );
  }

  make_qnumbers_and_dims( qnumbersarray, dimarray, tens, NULL, all_symarr, all_bonds, internalbonds,
      nr_internal );
  make_qnumbers_and_dims( qnumbersarray, dimarray, tens, &dim_of_blocks, all_symarr, all_bonds,
      internalbonds, nr_internal );

  for( i = 0 ; i < tens->nrsites ; ++i )
  {
    int * bonds = &all_bonds[ 3 * i ];
    struct symsecs * symarr = &all_symarr[ 3 * i ];
    
    destroy_dim_and_qnumbersarray( &dimarray[ i ], &qnumbersarray[ i ], symarr );
    destroy_ss_with_internal( symarr, internalbonds, nr_internal, bonds );
  }

  sort_and_make( tens, &dim_of_blocks );
}

static void change_internals_in_bookkeeper( struct symsecs internalsymsec[], 
    struct siteTensor * const tens )
{
  const int nr_internal = siteTensor_give_nr_internalbonds( tens );
  int internalbonds[ nr_internal ];
  int i;

  siteTensor_give_internalbonds( tens, internalbonds );
  for( i = 0 ; i < nr_internal ; ++i  )
  {
    struct symsecs * const psymsec = &bookie.list_of_symsecs[ internalbonds[ i ] ];
    destroy_symsecs( psymsec );
    *psymsec = internalsymsec[ i ];
  }
}

static void get_ss_with_internal( struct symsecs symsec[], const struct symsecs internalsymsec[],
    const int internalbonds[], const int nr_internal, const int bonds[] )
{
  int i;
  for( i = 0 ; i < 3 ; ++i )
  {
    int j;
    for( j = 0 ; j < nr_internal ; ++j )
      if( bonds[ i ] == internalbonds[ j ] ) break;
    if( j != nr_internal )
      symsec[ i ] = internalsymsec[ j ];
    else
      get_symsecs( &symsec[ i ], bonds[ i ] );
  }
}

static void destroy_ss_with_internal( struct symsecs symsec[], const int internalbonds[], 
    const int nr_internal, const int bonds[] )
{
  int i;
  for( i = 0 ; i < 3 ; ++i )
  {
    int j;
    for( j = 0 ; j < nr_internal ; ++j )
      if( bonds[ i ] == internalbonds[ j ] ) break;
    if( j == nr_internal )
      clean_symsecs( &symsec[ i ], bonds[ i ] );
  }
}

static void get_sites_to_use( const struct siteTensor * const tens, const int internalbonds[], 
    int sites_to_use[], const int nr_internal, int * const innersite )
{
  assert( tens->nrsites == 2 && "At this moment only for twosite tensors" );
  sites_to_use[ 0 ] = tens->sites[ 0 ];
  *innersite = tens->sites[ 1 ];
}

static void make_qnumbers_and_dims( int ***qnumbersarray[], int ***dimarray[], struct siteTensor * 
    const tens, int ** const dim_of_blocks, const struct symsecs all_symarr[], 
    const int all_bonds[], const int internalbonds[], const int nr_internal )
{
  int flag = dim_of_blocks == NULL; /* do I want to calculate the number of blocks or do I want to 
                                     * make the qnumbers? */
  int curr_nr_blocks = 0;
  int internal_site;
  int sym[ 3 ];
  int common_with_internal[ tens->nrsites ][ 2 ];
  int i;

  /* Now search the internal site. */
  for( internal_site = 0 ; internal_site < tens->nrsites ; ++internal_site )
  {
    int nr_of_internal_bonds = 0 ;
    for( i = 0 ; i < 3 ; ++i )
    {
      int j;
      for( j = 0 ; j < nr_internal ; ++j )
        if( internalbonds[ j ] == all_bonds[ internal_site * 3 + i ] )
        {
          ++nr_of_internal_bonds;
          break;
        }
    }
    if( nr_of_internal_bonds > 1 )
      break;
  }

  assert( internal_site != tens->nrsites || tens->nrsites == 2 );
  if( tens->nrsites == 2 ) internal_site = 1;

  /* make the common_with_internal array */
  for( i = 0 ; i < tens->nrsites ; ++i )
  {
    int j, k;
    if( i == internal_site )
      break;

    for( j = 0 ; j < 3 ; ++j )
    {
      for( k = 0 ; k < 3 ; ++k )
        if( all_bonds[ i * 3 + j ] == all_bonds[ internal_site * 3 + k ] )
          break;
      if( k != 3 )
        break;
    }
    common_with_internal[ i ][ 0 ] = j;
    common_with_internal[ i ][ 1 ] = k;
  }

  if( !flag )
  {
    tens->qnumbers = safe_malloc( tens->nrblocks * tens->nrsites, QN_TYPE );
    *dim_of_blocks = safe_malloc( tens->nrblocks, int );
  }

  for( sym[ 0 ] = 0 ; sym[ 0 ] < all_symarr[ internal_site * 3 ].nr_symsec ; ++sym[ 0 ] )
  {
    for( sym[ 1 ] = 0 ; sym[ 1 ] < all_symarr[ internal_site * 3 + 1 ].nr_symsec ; ++sym[ 1 ] )
    {
      for( sym[ 2 ] = 0 ; sym[ 2 ] < qnumbersarray[ internal_site ][ sym[ 0 ] ][ sym[ 1 ] ][ 0 ] ;
          ++sym[ 2 ] )
      {
        int site;
        int curr_block_size = 1;
        QN_TYPE *res_qnumbers[ tens->nrsites ];
        int *res_dim[ tens->nrsites ];
        int length[ tens->nrsites ];
        length[ internal_site ] = 1;

        for( site = 0 ; site < tens->nrsites ; ++site )
          if( site != internal_site )
          {
            QN_TYPE **p_res_qnumbers = flag ? NULL : &res_qnumbers[ site ];
            int **p_res_dim          = flag ? NULL : &res_dim[ site ];
            find_qnumbers_with_index_in_array( sym[ common_with_internal[ site ][ 1 ] ], 
                common_with_internal[ site ][ 0 ], qnumbersarray[ site ], dimarray[ site ], 
                &all_symarr[ 3 * site ], p_res_qnumbers, p_res_dim, &length[ site ] );
            curr_block_size *= length[ site ];
          }
        if( flag )
          curr_nr_blocks += curr_block_size;
        else
        {
          int indexes[ tens->nrsites ];
          int upper_bound = curr_nr_blocks + curr_block_size;
          int internal_dim = dimarray[ internal_site ][ sym[ 0 ] ][ sym[ 1 ] ][ sym[ 2 ] ];
          QN_TYPE internal_qnumber = sym[ 0 ] + sym[ 1 ] * all_symarr[ internal_site * 3 ].nr_symsec
            + qnumbersarray[ internal_site ][ sym[ 0 ] ][ sym[ 1 ] ][ sym[ 2 ] + 1 ] *
            all_symarr[ internal_site * 3 ].nr_symsec * all_symarr[ internal_site * 3+1].nr_symsec;

          for( i = 0 ; i < tens->nrsites ; ++i ) indexes[ i ] = 0;

          for( ; curr_nr_blocks < upper_bound ; ++curr_nr_blocks )
          {
            (*dim_of_blocks)[ curr_nr_blocks ] = 1;
            for( i = 0 ; i < tens->nrsites ; ++i )
            {
              if( i == internal_site )
              {
                tens->qnumbers[ curr_nr_blocks * tens->nrsites + i ] = internal_qnumber;
                (*dim_of_blocks)[ curr_nr_blocks ] *= internal_dim;
              }
              else
              {
                tens->qnumbers[ curr_nr_blocks * tens->nrsites + i ] = res_qnumbers[i][indexes[i]];
                (*dim_of_blocks)[ curr_nr_blocks ] *= res_dim[ i ][ indexes[ i ] ];
              }
            }

            for( i = 0 ; i < tens->nrsites ; ++i )
            {
              ++indexes[ i ];
              if( indexes[ i ] < length[ i ] )
                break;
              indexes[ i ] = 0;
            }
          }
          assert( i == tens->nrsites );
        }

        for( site = 0 ; site < tens->nrsites * !flag ; ++site )
        {
          if( site != internal_site )
          {
            safe_free( res_qnumbers[ site ] );
            safe_free( res_dim[ site ] );
          }
        }
      }
    }
  }

  assert( flag || curr_nr_blocks == tens->nrblocks );
  if( flag ) tens->nrblocks = curr_nr_blocks;
}

static void sort_and_make( struct siteTensor * const tens, int ** const dim_of_blocks )
{
  int * idx;
  QN_TYPE * sorted_qnumbers = safe_malloc( tens->nrblocks * tens->nrsites, QN_TYPE );
  int i;

  tens->blocks.beginblock = safe_malloc( tens->nrblocks + 1, int );
  idx = qnumbersSort( tens->qnumbers, tens->nrsites, tens->nrblocks );

  tens->blocks.beginblock[ 0 ] = 0; 
  for( i = 0 ; i < tens->nrblocks ; ++i )
  {
    int j;
    for( j = 0 ; j < tens->nrsites ; ++j )
      sorted_qnumbers[ i * tens->nrsites + j ] = tens->qnumbers[ idx[ i ] * tens->nrsites + j ];
    tens->blocks.beginblock[ i + 1 ] = (*dim_of_blocks)[ idx[ i ] ] + tens->blocks.beginblock[ i ];
  }

  safe_free( *dim_of_blocks );
  safe_free( tens->qnumbers );
  safe_free( idx );
  tens->qnumbers = sorted_qnumbers;
  tens->blocks.tel = safe_calloc( tens->blocks.beginblock[ tens->nrblocks ], EL_TYPE );
}

static void contractsiteTensors( struct siteTensor * const tens, struct siteTensor * const T3NS, 
    struct symsecs internalsymsec[] )
{
  const int nrbonds =  siteTensor_give_nr_of_couplings( tens ) * 3;
  int bonds[ nrbonds ];
  struct symsecs symarr[ nrbonds ];
  const int nr_internal = siteTensor_give_nr_internalbonds( tens );
  int maxdims[ nrbonds ];
  int internalbonds[ nr_internal ];
  int map[ 2 ];
  int site;
  int resblock;
  int * newtoold;
  
  struct siteTensor tens1 = T3NS[ tens->sites[ 0 ] ];
  struct siteTensor tens2 = T3NS[ tens->sites[ 1 ] ];

  const char NOTRANS = 'N';
  const int I_ONE = 1;
  const double ONE = 1;
  const double ZERO = 0;

  siteTensor_give_internalbonds( tens, internalbonds );
  newtoold = make_newtoold( &internalsymsec[ 0 ], internalbonds[ 0 ] );
  siteTensor_give_qnumberbonds( tens, bonds );
  get_maxdims_of_bonds( maxdims, bonds, nrbonds );
  get_symsecs_arr( symarr, bonds, nrbonds );

  assert( tens->nrsites == 2 && nr_internal == 1 && "At this moment only two-site optimization" );
  for( site = 0 ; site < tens->nrsites ; ++site )
    for( map[ site ] = 0 ; map[ site ] < 3 ; ++map[ site ] )
      if( bonds[ 3 * site + map[ site ] ] == internalbonds[ 0 ] )
        break;
  assert( map[ 0 ] == 2 && "The last bond of the first tensor is the internal one" );
  assert( map[ 1 ] < 2  && "The first or second bond of the second tensor is the internal one" );

  /* Two site optimization and furthermore last bond of the first tensor is contracted 
   * with the first or second bond of the second tensor */
  for( resblock = 0 ; resblock < tens->nrblocks ; ++resblock )
  {
    QN_TYPE oldqnumeros[ tens->nrsites ];
    int indexes[ 3 * tens->nrsites ];
    int M1, N1, M2, N2, K2;
    int i, m, k;
    int block1, block2;
    for( i = 0 ; i < tens->nrsites ; ++i )
      oldqnumeros[ i ] = tens->qnumbers[ tens->nrsites * resblock + i ];

    oldqnumeros[ 0 ] = change_newtooldqnumber( oldqnumeros[ 0 ], newtoold, &maxdims[ 0 ],
        internalsymsec[ 0 ].nr_symsec, map[ 0 ] );
    oldqnumeros[ 1 ] = change_newtooldqnumber( oldqnumeros[ 1 ], newtoold, &maxdims[ 3 ],
        internalsymsec[ 0 ].nr_symsec, map[ 1 ] );

    if( oldqnumeros[ 0 ] < 0 || oldqnumeros[ 1 ] < 0 )
      continue;

    block1 = qnumbersSearch( &oldqnumeros[ 0 ], 1, tens1.qnumbers, 1, tens1.nrblocks );
    block2 = qnumbersSearch( &oldqnumeros[ 1 ], 1, tens2.qnumbers, 1, tens2.nrblocks );

    if( block1 < 0 || block2 < 0 )
      continue;
    assert( tens1.qnumbers[ block1 ] == oldqnumeros[ 0 ] );
    assert( tens2.qnumbers[ block2 ] == oldqnumeros[ 1 ] );

    EL_TYPE *tel1   = get_tel_block( &tens1.blocks, block1 );
    EL_TYPE *tel2   = get_tel_block( &tens2.blocks, block2 );
    EL_TYPE *telres = get_tel_block( &tens->blocks, resblock );

    for( site = 0 ; site < tens->nrsites ; ++site )
    {
      //QN_TYPE curr_QN = tens->qnumbers[ resblock * tens->nrsites + site ];
      QN_TYPE curr_QN = oldqnumeros[ site ];
      int i;
      for( i = 0 ; i < 3 ; ++i )
      {
        indexes[ site * 3 + i ] = curr_QN % maxdims[ site * 3 + i ];
        curr_QN /= maxdims[ site * 3 + i ];
      }
      assert( curr_QN == 0 );
    }

    /* calculate the different dimensions of the blocks. */
    assert( indexes[ map[ 0 ] ] == indexes[ 3 + map[ 1 ] ] );
    M1 = symarr[ 0 ].dims[ indexes[ 0 ] ] * symarr[ 1 ].dims[ indexes[ 1 ] ];
    N1 = symarr[ 2 ].dims[ indexes[ 2 ] ];

    M2 = 1;
    for( i = 0 ; i < map[ 1 ] ; ++i ) M2 *= symarr[ 3 + i ].dims[ indexes[ 3 + i ] ];
    N2 = symarr[ 3 + i ].dims[ indexes[ 3 + i ] ];
    ++i;
    K2 = 1;
    for( ; i < 3 ; ++i ) K2 *= symarr[ 3 + i ].dims[ indexes[ 3 + i ] ];

    assert( N1 == N2 );
    assert( M1 * N1 == get_size_block( &tens1.blocks, block1 ) );
    assert( M2 * N2 * K2 == get_size_block( &tens2.blocks, block2 ) );
    assert( M1 * M2 * K2 == get_size_block( &tens->blocks, resblock ) );

    for( m = 0 ; m < M2 ; ++m )
      for( k = 0 ; k < K2 ; ++k )
        dgemv_( &NOTRANS, &M1, &N1, &ONE, tel1, &M1, tel2 + m + k * M2 * N2, &M2, &ZERO, 
            telres + m * M1 + k * M1 * M2, &I_ONE );
  }

  safe_free( newtoold );
  clean_symsecs_arr( symarr, bonds, nrbonds );
}

static QN_TYPE change_newtooldqnumber( QN_TYPE new, int * newtoold, int * maxdims, const int newdim,
    const int map )
{
  const int olddim = maxdims[ map ];
  int i;
  int indices[ 3 ];
  maxdims[ map ] = newdim;
  for( i = 0 ; i < 3 ; ++i )
  {
    indices[ i ] = new % maxdims[ i ];
    new /= maxdims[ i ];
  }
  assert( new == 0 );
  indices[ map ] = newtoold[ indices[ map ] ];
  if( indices[ map ] == -1 )
  {
    maxdims[ map ] = olddim;
    return -1;
  }

  maxdims[ map ] = olddim;
  for( i = 2 ; i >= 0 ; --i )
  {
    new *= maxdims[ i ];
    new += indices[ i ];
  }
  return new;
}

static int * make_newtoold( const struct symsecs * const internalss, const int bond )
{
  struct symsecs newss;
  int * result = safe_malloc( internalss->nr_symsec, int );
  int i = 0;
  int currj = 0;
  get_symsecs( &newss, bond );

  for( i = 0 ; i < internalss->nr_symsec ; ++i )
  {
    int j;
    result[ i ] = -1;
    for( j = currj ; j < newss.nr_symsec ; ++j )
    {
      int k;
      for( k = 0 ; k < bookie.nr_symmetries ; ++k )
        if( internalss->irreps[ i * bookie.nr_symmetries + k ] != 
            newss.irreps[ j * bookie.nr_symmetries + k ] )
          break;
      if( k == bookie.nr_symmetries )
      {
        result[ i ] = j;
        currj = j + 1;
        break;
      }
    }
  }
  clean_symsecs( &newss, bond );
  return result;
}