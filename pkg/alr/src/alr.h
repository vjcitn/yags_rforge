/* /usr16/stdevs/stdev0f/SLIBS/alr.dev/SCCS/s.alr.h */
/* alr.h 1.2 98/02/23 */

#include "stdlib.h"
#include "chanmat.h"

#define set_matrix_array( arrname, nel ) \
   if (!( arrname = (MATRIX **)malloc((unsigned)(nel * sizeof(MATRIX *))))) \
     { fprintf(stderr, \
	    "set_matrix_array (mac): out of memory for arrname nel elements\n"); \
     exit(1); \
     }

#define SMALL .00000000001


static double dmaxarg1, dmaxarg2;
#define IMAX(a,b) (dmaxarg1=(a), dmaxarg2=(b), (dmaxarg1) > (dmaxarg2) ? \
(int)(dmaxarg1) : (int)(dmaxarg2))

static double dminarg1, dminarg2;
#define IMIN(a,b) (dminarg1=(a), dminarg2=(b), (dminarg1) < (dminarg2) ? \
(int)(dminarg1) : (int)(dminarg2))

#define MAX_NUM_CLUSTS 2000
#define MAX_NUM_CLASS  10

/*---------------------------------------------------------------------------*/
int dgecoXXY_();
int split();
int nchanges();
