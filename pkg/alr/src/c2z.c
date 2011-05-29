
/* c2z 1.2 95/03/23 */
/*---------------------------------------------------------------------------*/
/* class2z.c                                                                 */
/*---------------------------------------------------------------------------*/

#include "alr.h"


void class2z( class, id, dmat, z, zid, k, n, nrz, p, q, nclust, flag)

     double *class, *id, *dmat, *z, *zid;
     int    *k, *n, *nrz, *p, *q, *nclust, *flag;

{
  MATRIX *C[MAX_NUM_CLUSTS], *ID[MAX_NUM_CLUSTS], *cin, *idin;
  MATRIX *D[MAX_NUM_CLASS][MAX_NUM_CLASS], *din, *Z, *Zid, *ID_num;
  int    i, j, a, b, one, count, ni, cc, rr;

  one=1;
  from_S( class, n, &one, cin);
  from_S( id, n, &one, idin);
  from_S( dmat, p, q, din);
  from_S( z, nrz, q, Z);
  from_S( zid, nrz, &one, Zid);

printf("[1]\n");

  count=0;
  for(i=0;i< *k;i++){
    for(j=i;j< *k;j++){
      D[i][j]=create_matrix(1,*q,EPHEMERAL);
      for(rr=0;rr<*q;rr++) MEL(D[i][j],0,rr)=MEL(din,count,rr);
      count++;
    }
  }
  destroy_matrix(din);

printf("[2]\n");

  split(cin,idin,C);
  split(idin,idin,ID);
  destroy_matrix(idin);
  destroy_matrix(cin);

  ID_num=create_matrix(*nclust,one,EPHEMERAL);
  for(i=0;i<*nclust;i++){
    MEL(ID_num,i,0)=MEL(ID[i],0,0);
    destroy_matrix(ID[i]);
  }

printf("[3]\n");

  count=0;
  for(cc=0;cc<*nclust;cc++){
    ni=C[cc]->nrows;
    if(ni==1){
      for(rr=0;rr<*q;rr++) MEL(Z,count,rr)=0.0;
      MEL(Zid,count,0)=MEL(ID_num,cc,0);
      count++; 
    }else{
      for(i=0;i<(ni-1);i++){
	for(j=(i+1);j<ni;j++){
	  a=IMIN(MEL(C[cc],i,0),MEL(C[cc],j,0));
	  b=IMAX(MEL(C[cc],i,0),MEL(C[cc],j,0));
	  for(rr=0;rr<*q;rr++) MEL(Z,count,rr)=MEL(D[a-1][b-1],0,rr);
	  MEL(Zid,count,0)=MEL(ID_num,cc,0);
	  count++;
	}
      }
    }

  }

printf("[4]\n");

  if( count != Z->nrows ) *flag=1;

    to_S( Z, z);
    to_S( Zid, zid);

    for(i=0;i<*nclust;i++) destroy_matrix(C[i]);
    destroy_matrix(ID_num);

  for(i=0;i< *k;i++){
    for(j=i;j< *k;j++){
      destroy_matrix(D[i][j]);
    }
  }

}
/*---------------------------------------------------------------------------*/


/* start of refab.c */

int oldr( n, i, j )
int n, i, j;   /* n is cluster size, i and j are zero based indices < n */
{
int P;
if (i == 0) return(j - 1);
P = .5*(i*i + i);
return(i*n - P + j - i - 1);
}

int newr( n, i, j )
int n, i, j;   /* n is cluster size, i and j are zero based indices < n */
{
int off;
if (j > i) off = j-1;
else off = j;
return( (n-1)*i + off );
}

MATRIX *arr2mat( arr, nel )
MATRIX **arr;
int nel;
{
int nrow = 0, i ,j,k, ncol, outrow;
MATRIX *out;
ncol = arr[0]->ncols;
for (i = 0 ; i < nel ; i++ )
    nrow += arr[i]->nrows;
out = create_matrix( nrow, ncol, PERMANENT );
outrow = 0;
for (i = 0 ; i < nel ; i++ )
  {
  for (j = 0 ; j < arr[i]->nrows; j++)
   {
    for (k = 0 ; k < ncol ; k++ )
    {
    MEL(out,outrow,k) = MEL(arr[i],j,k);
    }
   outrow++;
   }
  }
  return out;
}



void refab(  z, zid, nrz, q , newzp)

/* this program refabricates the output of Heagerty's class2z */
/* function to obtain the redundant "permutation-invariant" form of Z */
/* input is the Z matrix and associated discriminator "id", */
/* number of rows and columns of Z, and a pointer to space */
/* holding the invariant Z -- this space will be filled by*/
/* the new Z matrix */

     double *z, *zid, *newzp;
     int    *nrz, *q ;

{
	MATRIX *Zin, *Zidin, *nnz;
  MATRIX **Z, **NEWZ;
  int    nc , one ,newzrow, nel, thisnrows, thisncols;
int oldrow, newrow,i,j,k,m;

  one=1;
  from_S( z, nrz, q, Zin);
  from_S( zid, nrz, &one, Zidin);

  nc = nchanges(Zidin);
  set_matrix_array( Z, nc );
  set_matrix_array( NEWZ, nc );
split(Zin, Zidin, Z);

for ( i = 0 ; i < nc ; i++ )
   {
 thisnrows = Z[i]->nrows;
 thisncols = Z[i]->ncols;
NEWZ[i] = create_matrix(thisnrows*2, thisncols, PERMANENT);
 nel = (int)((1+sqrt((double)1.+8.*((double)thisnrows)))/2.);
   for (j = 0 ; j < nel ; j++)
      {
       for (k = 0 ; k < nel ; k++ )
        {
	if (j == k) continue;
	if (j < k)
	   {
	   oldrow = oldr(nel,j,k);
	   newrow = newr(nel,j,k);
	   for (m = 0; m < thisncols; m++ )
	       MEL(NEWZ[i],newrow,m) = MEL(Z[i],oldrow,m);
	   }
	else
	   {
	   oldrow = oldr(nel,k,j);
	   newrow = newr(nel,j,k);
	   for (m = 0; m < thisncols; m++ )
	       MEL(NEWZ[i],newrow,m) = MEL(Z[i],oldrow,m);
	   }
	}
       }
}
nnz = arr2mat( NEWZ, nc );
to_S(nnz,newzp)
}
