#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "hamiltonian.h"
#include "hamiltonian_qc.h"
#include "bookkeeper.h"
#include "symmetries.h"

enum hamtypes ham;

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

/** Sets the Hamiltonian to the right internal enum. **/
static int set_hamiltonian( char hamiltonian[] );

/* ============================================================================================ */

void readinteraction( char interactionstring[] )
{
  if( !set_hamiltonian( interactionstring ) )
    exit(EXIT_FAILURE);

  switch( ham )
  {
    case QC :
    case QCSU2 :
      QC_make_hamiltonian( interactionstring );
      break;
    default:
      fprintf( stderr, "ERROR : unrecognized interaction %s.\n", interactionstring );
      exit( EXIT_FAILURE );
  }
}

void get_physsymsecs( struct symsecs *res, int bond )
{
  switch( ham )
  {
    case QC :
    case QCSU2 :
      QC_get_physsymsecs( res, bond );
      break;
    default:
      fprintf( stderr, "Unrecognized Hamiltonian.\n");
  }
}

void get_hamiltoniansymsecs( struct symsecs * const res, const int bond )
{
  switch( ham )
  {
    case QC :
    case QCSU2 :
      QC_get_hamiltoniansymsecs( res, bond );
      break;
    default:
      fprintf( stderr, "%s@%s: Unrecognized Hamiltonian.\n", __FILE__, __func__ );
  }
}

int get_nr_hamsymsec( void )
{
  switch( ham )
  {
    case QC :
    case QCSU2 :
      return QC_get_nr_hamsymsec();
      break;
    default:
      fprintf( stderr, "Unrecognized Hamiltonian.\n");
      return 0;
  }
}

int get_trivialhamsymsec( void )
{
  switch( ham )
  {
    case QC :
    case QCSU2 :
      return QC_get_trivialhamsymsec( );
      break;
    default:
      fprintf( stderr, "Unrecognized Hamiltonian.\n");
      return -1;
  }
}

int give_hermhamsymsec( const int orighamsymsec )
{
  switch( ham )
  {
    case QC :
    case QCSU2 :
      return QC_give_hermhamsymsec( orighamsymsec );
      break;
    default:
      fprintf( stderr, "Unrecognized Hamiltonian.\n");
      return -1;
  }
}

int consistencynetworkinteraction( void )
{
  switch( ham )
  {
    case QC:
    case QCSU2:
      return QC_consistencynetworkinteraction();
    default:
      fprintf( stderr, "Unrecognized Hamiltonian.\n");
      exit( EXIT_FAILURE );
  }

  return 0;
}

int get_hamsymsec_site( const int siteoperator, const int site )
{
  switch( ham )
  {
    case QC :
      return QC_get_hamsymsec_site( siteoperator, site );
    case QCSU2 :
    default:
      fprintf( stderr, "%s@%s: unrecognized Hamiltonian.\n", __FILE__, __func__ );
      exit( EXIT_FAILURE );
  }
}

double get_site_element( const int siteoperator, const int braindex, const int ketindex )
{
  switch( ham )
  {
    case QC :
      return QC_get_site_element( siteoperator, braindex, ketindex );
    case QCSU2 :
    default:
      fprintf( stderr, "%s@%s: unrecognized Hamiltonian.\n", __FILE__, __func__ );
      exit( EXIT_FAILURE );
  }
}

void hamiltonian_tensor_products( int * const nr_of_prods, int ** const possible_prods, const int
    resulting_hamsymsec, const int site )
{
  switch( ham )
  {
    case QC :
      QC_hamiltonian_tensor_products( nr_of_prods, possible_prods, resulting_hamsymsec, site );
      break;
    case QCSU2 :
    default:
      fprintf( stderr, "%s@%s: unrecognized Hamiltonian.\n", __FILE__, __func__ );
      exit( EXIT_FAILURE );
  }
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static int set_hamiltonian( char hamiltonian[] )
{
  char *ext = strrchr( hamiltonian, '.' );
  if( ext )
  {
    char *extfcidump = "FCIDUMP";

    ext++;
    while( *ext && tolower( *( ext++ ) ) == tolower( *( extfcidump++ ) ) );

    /* extension is fcidump */
    if( *ext == *extfcidump ){
      enum symmetrygroup symmQC[] = { Z2, U1, U1 };
      enum symmetrygroup symmQCSU2[] = { Z2, U1, SU2 };
      int i;
      if( bookie.nr_symmetries != 3 && bookie.nr_symmetries != 4 )
      {
        fprintf( stderr, "Invalid symmetry groups for quantum chemistry were inputted!\n" );
        return 0;
      }
      if( bookie.nr_symmetries == 4 && bookie.sgs[ 3 ] < C1 )
      {
        fprintf( stderr, "Invalid symmetry groups for quantum chemistry were inputted!\n" );
        return 0;
      }
      for( i = 0 ; i < 3 ; i++ )
        if( symmQC[ i ] != bookie.sgs[ i ] )
          break;
      if( i == 3 )
      {
        ham = QC;
        return 1;
      }

      for( i = 0 ; i < 3 ; i++ )
        if( symmQCSU2[ i ] != bookie.sgs[ i ] )
          break;
      if( i == 3 )
      {
        ham = QCSU2;
        return 1;
      }
      fprintf( stderr, "Invalid symmetry groups for quantum chemistry were inputted!\n" );
      return 0;
    }
  }

  fprintf( stderr, "ERROR : Interaction %s is an unknown interaction.\n", hamiltonian );
  return 0;
}
