/* /usr16/stdevs/stdev0f/SLIBS/alr.dev/SCCS/s.chanmat.h */
/* chanmat.h 3.1 97/03/23 */

/* alr support */
#include "chanmatstruct.h"
#include "chanmatfuns.h"
#include <setjmp.h>
extern jmp_buf env;

#define errorbranch( exitlabel ) \
if ( setjmp(env) != 0 ) \
        { \
        fprintf(stderr,"chanmat error detected, returning to caller\n"); \
        goto exitlabel ; \
        }

