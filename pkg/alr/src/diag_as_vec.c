/* /usr16/stdevs/stdev0f/SLIBS/alr.dev/SCCS/s.diag_as_vec.c */
/* alr support diag_as_vec.c 3.1 97/03/23 */


#include "chanmat.h"
#include<setjmp.h>

jmp_buf env;

MATRIX *diag_as_vec(inmat)
MATRIX *inmat;
{
  int i;
  MATRIX *outmat;

  if(inmat->ncols!=inmat->nrows)
    {
      fprintf(stderr,"M+-: diag_as_vec: arg is not a square matrix. Dies.\n");
      fprintf(stderr,"\nNumber of columns = %d",inmat->ncols);
      fprintf(stderr,"\nNumber of rows    = %d\n",inmat->nrows);
      Seterr_and_terminate(DIAG_AS_VEC_ARG_BAD);
    }

  outmat= create_matrix(inmat->nrows,1,EPHEMERAL);
  for(i= 0;i<inmat->nrows;i++)
    {
      *(ELREF(outmat,i,0))=  *(ELREF(inmat,i,i));
    }
  return outmat;
}



MATRIX *matsqrt(x)
MATRIX *x;
{
  int i,j;
  MATRIX *tmp;
  tmp= matcopy(x);
  for(i= 0;i<x->ncols;i++)
    {
      for(j= 0;j<x->nrows;j++)
	{
	  MEL(tmp,i,j)= sqrt(MEL(x,i,j));
	}
    }
  if(is_ephemeral(x))destroy_matrix(x);
  return tmp;
}

MATRIX *mat1over(x)
MATRIX *x;
{
  int i,j;
  MATRIX *tmp;
  tmp= matcopy(x);
  for(i=0;i<x->ncols;i++)
    {
      for(j=0;j<x->nrows;j++)
	{
	  MEL(tmp,i,j)= 1./(MEL(x,i,j));
	}
    }
  if(is_ephemeral(x))destroy_matrix(x);
  return tmp;
}
