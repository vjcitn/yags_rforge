/* revision for symmetry and weighting: */
/* /usr16/stdevs/stdev0f/SLIBS/alr4/SCCS/s.alr.c */
/* 4.11 98/02/17 */
/* genesis: */
/* /usr16/stdevs/stdev0f/SLIBS/alr.dev/SCCS/s.alr.c */
/* 3.4 97/03/23 */
/* the following memory allocation redefinitions */
/* should be checked for compatibility with S */
/*
#define calloc(x,y) malloc((unsigned)(x)*(y))
*/

/* alr.c @(#) chanlib revisions: alr.c 3.4 97/03/23 */
/* alternating logistic regression */

#include "alr.h"

#define sweep luinv 

MATRIX *stack(), *matantilogit(), *star(), *matxdiagasvec(), *oneminus(); 
MATRIX *pluck_Z(), *flip(), *flipz(), *flipinds(), *luinv(), *make_Rvar();
MATRIX *make_irvar(), *matsort();

double matmaxabs();
extern int tchoose2();
int nrfun();

double estcond();

typedef struct matrix_pair_p   /* appended! */
	{
	MATRIX *X1;
	MATRIX *X2;
	MATRIX **SB;
	} MATRIX_PAIR_P;

typedef struct matrix_triplet_p
	{
	MATRIX *X1;
	MATRIX *X2;
	MATRIX *X3;
	MATRIX **SA;
	} MATRIX_TRIPLET_P;

typedef struct matwcond
	{
	MATRIX *mat;
	double condition;
	} MATWCOND;

MATWCOND *invwcond();


MATRIX_PAIR_P *get_beta_comps(), *beta_comps;
MATRIX_PAIR_P *get_beta_comps_exch();
MATRIX_PAIR_P *get_beta_comps_2cl();
MATRIX_TRIPLET_P *get_alpha_comps(), *alpha_comps, *falpha_comps;
MATRIX_TRIPLET_P *get_alpha_comps_exch(), *get_alpha_comps_2cl();
		
/* main()
{
MATWCOND *x;
MATRIX *X;
X = matread("X");
x = invwcond( X );
matdump( matmult( x->mat, X ) );
}
*/

void alr( x, y, z, z_is_mast, id, zid, zlocs, nobs, nzrows, p, q, 
	     initbeta, initalpha, cov,
	     tol, maxiter , exch, twocl, bindep, cl2p, silent, weights)

double *x, *y, *id, *initbeta, *initalpha, *cov, *tol, *z, *zid, *zlocs, *cl2p;
double *weights;
int *nzrows, *nobs, *p, *q, *maxiter, *z_is_mast, *exch, *twocl, *bindep, *silent;

#define Matdump( x )
#define Printf( x )
#define IPrintf( x )
#define sMatdump( x ) { printf(" x \n"); matdump( x ); }
#define sPrintf( x ) { printf(" x -> %lf\n", x ); }
#define sIPrintf( x ) { printf(" x -> %d\n", x ); }
{

MATRIX *xin, *yin, *idin, *zin, *beta, *alpha, *zidin, *zlin, *cl2in;
MATRIX **X, **Y, **Z, **ZL, **CL2, **W;
MATRIX *thisZ;
MATRIX *mu1, *psi2, *C, *B, *op, *del1, *del2, *nbeta;
MATRIX *varmat, *S, *del1s, *del2s, *adel1s, *adel2s, *delta;
MATRIX *D, *G, *E, *ystar, *Zeta;
MATRIX *aop, *aH1, *saH2, *alvarmat, *DtEG, *aldel, *nalpha, *fulldel;
MATRIX *robvmat, *Ni, *Win;
int one = 1, *onep, nclust, c=0, iter, ni, nc2, i=0, j=0, n;
int r0, k, l, flip, invariant, is_exchangeable, is_twocl;
double mu2, off, mui, muj, psi2ij, a, chkrt;
double gamr, mpr, tmpm, mu2ij, tmpd, toff, tmpalogit, tmpn1;
double dNdal, tmpN1, dDdalk, dndalk, xil, xjl, dmidbl, dmjdbl;
double tmpf2, tmpf2d, Zij, Zrk, Alphak, dadalk, Dn, chk;

double dadbl, dnrdbl, tmpf, td; 

onep = &one;

fprintf(stderr,"@(#) alr.c: @(#) alr.c 4.11 98/02/17, chanlib version.\n");
is_exchangeable = *exch;
is_twocl = *twocl;

from_S( y, nobs, onep, yin )
from_S( id, nobs, onep, idin )
from_S( weights, nobs, onep, Win )
from_S( x, nobs, p, xin )
from_S( initbeta, p, onep, beta )
from_S( initalpha, q, onep, alpha )

make_permanent( beta );
make_permanent( alpha );

if (!is_exchangeable && !is_twocl)
{
from_S( z, nzrows, q, zin )
make_permanent( zin );
}

nclust = nchanges( idin );
n = *nobs;

set_matrix_array( X, nclust );
set_matrix_array( Y, nclust );
set_matrix_array( W, nclust );

if (is_twocl)
{
from_S( cl2p, nobs, onep, cl2in )
set_matrix_array( CL2, nclust );
split( cl2in, idin, CL2 );
destroy_matrix( cl2in );
}

split( xin, idin, X );
split( yin, idin, Y );
split( Win, idin, W );

destroy_matrix( xin );
destroy_matrix( yin );
destroy_matrix( Win );

if (!is_exchangeable && !is_twocl)
{
if (!(*z_is_mast)) 
	{
	from_S( zid, nzrows, onep, zidin )
	set_matrix_array( Z, nclust );
	split( zin, zidin, Z );
	destroy_matrix( zin );
	destroy_matrix( zidin );
	}
else
	{
	from_S( zlocs, nobs, onep, zlin )
	set_matrix_array( ZL, nclust );
	split( zlin, idin, ZL );
	destroy_matrix( zlin );
	}
}

destroy_matrix( idin );
Ni = create_matrix( nclust, 1, PERMANENT );
for (i = 0 ; i < nclust ; i++ )
	MEL(Ni,i,0) = Y[i]->nrows;

iter = 0;
do {
	if ( iter > 0 )
		{
		for ( i = 0 ; i < nclust ; i++ )
			{
			destroy_matrix( beta_comps->SB[i] );
			if (MEL(Ni,i,0) >1) destroy_matrix( alpha_comps->SA[i] );
			}
		cfree( beta_comps->SB );
		cfree( alpha_comps->SA ); 
		destroy_matrix( beta_comps->X1 );
		destroy_matrix( beta_comps->X2 );
		destroy_matrix( alpha_comps->X1 );
		destroy_matrix( alpha_comps->X2 );
		destroy_matrix( alpha_comps->X3 );
		/*
		cfree( beta_comps );
		cfree( alpha_comps ); 
		*/
		}
	if (!(*silent)) printf("Starting beta\n");

	if ( is_exchangeable )
		beta_comps = get_beta_comps_exch( nclust, X, Y, beta, 
				    alpha, *bindep , W);
	else if ( is_twocl )
		beta_comps = get_beta_comps_2cl( nclust, X, Y, beta, 
				    alpha , CL2 , *bindep, W );
	else beta_comps = get_beta_comps( nclust, X, Y, Z, beta, 
				    alpha, ZL, zin, *z_is_mast, *bindep, W );

	if (!(*silent)) printf("Done beta\n");

	make_permanent( beta_comps->X1 );
	make_permanent( beta_comps->X2 );
	delta = matmult( sweep(beta_comps->X1), beta_comps->X2 );

	make_permanent(delta);
	nbeta =  matadd( beta, delta );
	destroy_matrix(beta);
	beta = matcopy( nbeta );

	destroy_matrix(nbeta);
	make_permanent(beta);

	Matdump(beta);

	/* CONDIT LOGIST */
	if (!(*silent)) printf("Starting alpha\n");

	flip = 0;
	if (is_exchangeable) alpha_comps = get_alpha_comps_exch( nclust, X, Y, beta,
				 alpha, flip, W );
	else if (is_twocl) alpha_comps = get_alpha_comps_2cl( nclust, X, Y, beta,
				 alpha, flip, CL2, W );
	else alpha_comps = get_alpha_comps( nclust, X, Y, Z, beta,
				 alpha, ZL, zin, *z_is_mast , flip, W );

	make_permanent( alpha_comps->X1 );
	make_permanent( alpha_comps->X2 );
	make_permanent( alpha_comps->X3 );

	aldel = matmult( sweep(alpha_comps->X1), alpha_comps->X2 );

	make_permanent(aldel);
	nalpha = matadd( alpha, aldel );
	destroy_matrix(alpha);

	alpha = matcopy( nalpha );
	make_permanent(alpha);

	Matdump( alpha );

	fulldel = stack( delta, aldel );
	destroy_matrix(delta);
	destroy_matrix(aldel);

	chk = matmaxabs( fulldel );

	if (!(*silent)) printf("%d, %lf -- iter, curdel\n", iter, chk );

	} while ( iter++ < *maxiter && ( chk > *tol ) );

	/* compat chk */


	invariant = 0;
	if (!(invariant))
		robvmat = make_Rvar( nclust, beta_comps, alpha_comps, Ni ); 
	else {
		flip = 1;
/*		if ( is_exchangeable) falpha_comps = get_alpha_comps_exch( nclust, X, Y, beta,
				alpha, flip );
		else if ( is_twocl) falpha_comps = get_alpha_comps_2cl( nclust, X, Y, beta,
				alpha, flip , CL2 );
		else falpha_comps = get_alpha_comps( nclust, X, Y, Z, beta,
				alpha, ZL, zin, *z_is_mast , flip ); */
		/* robvmat = make_rvar( nclust, beta_comps, falpha_comps ); */
/*		robvmat = make_irvar( nclust, beta_comps, alpha_comps, falpha_comps, Ni ); */
		}

	/* matdump(robvmat); */

	for ( i = 0 ; i < nclust ; i++ )
		{
		destroy_matrix( X[i] );
		destroy_matrix( Y[i] );
		if ( is_twocl ) destroy_matrix( CL2[i] );
		if (!is_exchangeable && !is_twocl)
		{
		if ( *z_is_mast ) destroy_matrix(ZL[i]);
			else destroy_matrix( Z[i] );
		}
		}

	cfree( X );
	cfree( Y );
	if (!is_exchangeable)
		{
	if ( *z_is_mast ) cfree( ZL )
	   else cfree( Z )
		}

	to_S( robvmat, cov )
	to_S( beta, initbeta )
	to_S( alpha, initalpha )

        *maxiter = iter;
}

/* pluck_Z uses 1-based location indices */

MATRIX *pluck_Z( Z, zl )  /* if this works, these computations should */
MATRIX *Z, *zl;            /* be embedded wherever Z is used */
{
int ni, nc, nzr, nic2, nn, i, j, m, this_row, k, find, sind;
MATRIX *zret;

ni = zl->nrows;
nc = zl->ncols;
if ( nc > 1 ) fprintf(stderr,"pluck_Z: zl should be column.\n");

nc = Z->ncols;
nzr = Z->nrows;

nic2 = tchoose2( ni );

zret = create_matrix( nic2, nc, EPHEMERAL );

nn = (1 + (int)sqrt((double)(1. + 4.*(double)nzr )))/2;

k = 0;
for ( i = 0 ; i < ni ; i++ )
	{
	for ( j = 0; j < ni ; j++ )
		{
		if (i == j) continue;
		find = (int) MEL( zl, i, 0 );
		sind = (int) MEL( zl, j, 0 );
		this_row = nrfun( nn, find-1, sind-1 ); /* 0 based */
		for ( m = 0 ; m < nc ; m++ )
			{
			MEL( zret, k, m ) = MEL( Z, this_row, m );
			}
		k++;
		}
	}
return zret;
}

MATRIX_PAIR_P *get_beta_comps( nclust, X, Y, Z, beta, 
				    alpha, ZL, zin, z_is_mast , bindep, W)
int nclust, z_is_mast, bindep;
MATRIX **X, **Y, **Z, *beta, *alpha, **ZL, *zin, **W;
{
int c, ni, nc2, i,j, r0, p, condcount ;
MATRIX *thisZ, *mu1, *psi2, *B, *C, *op, *del1, *del2, *varmat, 
    *del1s, *del2s, *S;
double psi2ij, mui, muj, a, chkrt, mu2;
MATRIX_PAIR_P *ans;
MATWCOND *Binvc;
double mincond = 1.;
double thisw;

p = X[0]->ncols;

ans = (MATRIX_PAIR_P *)calloc( 1, (unsigned)sizeof(struct matrix_pair_p));
ans->SB = ( MATRIX **)calloc(nclust,(unsigned)sizeof(struct matrix));
if ( ans->SB == (MATRIX **)NULL)
	{
	fprintf(stderr,"MATRIX_PAIR_P out of space for SB\n");
	exit(1);
	}

condcount = 0;

for ( c = 0 ; c < nclust ; c++ )
	{
thisw = MEL(W[c],0,0);
	if (!(z_is_mast)) thisZ = matcopy( Z[c] );
	else thisZ = pluck_Z( zin, ZL[c] );  /* should pluck */
	/* else thisZ = matcopy( zin );  *//* should pluck */

	ni = X[c]->nrows;
	nc2 = tchoose2( ni );

	mu1 = matantilogit( matmult( X[c], beta ) );
	psi2 = matexp( matmult( thisZ, alpha ) );

	make_permanent(mu1);
	make_permanent(psi2);

	if (!bindep)
	{

	B = create_matrix( ni, ni, PERMANENT );

	for ( i = 0 ; i < ni ; i++ )
		{
		for ( j = 0 ; j < ni ; j++ )
			{

			if (i == j) continue;
			r0 = nrfun( ni, i, j );  /* zerobased r */

			psi2ij = MEL( psi2, r0, 0 );
			mui = MEL( mu1, i, 0 );
			muj = MEL( mu1, j, 0 );

			a = 1. + (mui+muj)*(psi2ij-1.);

			chkrt = a*a - 4.*psi2ij*(psi2ij-1.)*mui*muj;

			if (chkrt < 0.) 
			  fprintf(stderr,"beta comps: Mu2 calcs have invalid sqrt.\n");

			if (fabs(psi2ij-1.) < SMALL )
			  fprintf(stderr,"Mu2 calcs have invalid denom, c= %d i=%d, j=%d, r0=%d.\n",c,i,j,r0);

			if ( psi2ij == 1.) mu2 = mui*muj;
			else mu2 =
			  ( a - sqrt( chkrt ) )/( 2.*(psi2ij-1.) );

			MEL(B,i,j) = mu2 - mui*muj;  /* while it is around */
			MEL(B,j,i) = MEL(B,i,j);
			}
			     
		}
	for ( i = 0 ; i < ni ; i++ )
		MEL(B,i,i) = MEL( mu1, i, 0 )*(1. - MEL( mu1, i, 0 ) );
	}
	else
           {
           B = create_matrix( ni, 1, PERMANENT );
           for ( i = 0 ; i < ni ; i++ )
                MEL(B,i,0) = 1./(MEL( mu1, i, 0 )*(1. - MEL( mu1, i, 0 ) ));
           }



	C = px1_times_pxq( 
		 px1_times_pxq( mu1, oneminus( mu1 ) ) , X[c] );

	make_permanent(C);

	if ( !bindep )
	   op = matmult( transp(C), sweep(B) );  /* code restored 5-29-2011 */
/* 5-29-2011: code dropped owing to linpack interface problem
	   {
	   Binvc = invwcond(B);
	   if (Binvc->condition < mincond) mincond = Binvc->condition;
	   if ( Binvc->condition < .03 && condcount < 6) 
	       {
	       condcount++;
	       fprintf(stderr,
			"get_beta_comps: full weighting has low condition #\n");
	       fprintf(stderr,
			"get_beta_comps: cluster %d, cond(dgeco) = %lf\n", c,
							   Binvc->condition );
	       fprintf(stderr,
			"get_beta_comps: consider independence model for beta\n");
		if ( condcount > 4 ) 
		   fprintf(stderr,"No more condition msgs this iter.\n");
	       }
	   op = matmult( transp(C), Binvc->mat );
	   }
*/
	else op = matxdiagasvec( transp(C), B );

	make_permanent(op);

	del1 = matmult( op, C );

	S = scalar_times_matrix( thisw, matsub( Y[c] , mu1 ) );

	del2 = matmult( op, S );
	ans->SB[c] = matcopy( del2 );

	make_permanent( del1 );
	make_permanent( del2 );

	if ( c == 0 )
		{
		/* varmat = create_matrix( p, p, PERMANENT ); */
		del1s = create_matrix( p, p, PERMANENT );
		del2s = create_matrix( p, 1, PERMANENT );
		}

	/* varmat = matadd( varmat, matmult( del2, transp(del2) ) ); */

	del1s = matadd( del1s, del1 );

	del2s = matadd( del2s, del2 );

	destroy_matrix(mu1);
	destroy_matrix(psi2);
	destroy_matrix(B);
	destroy_matrix(C);
	destroy_matrix(op);
	destroy_matrix(del1);
	destroy_matrix(del2);
	}
make_permanent( del1s );
make_permanent( del2s );
ans->X1 = del1s;
ans->X2 = del2s;
return ans;
}

/* as long as Z is correctly formatted, should work */

MATRIX_TRIPLET_P *get_alpha_comps( nclust, X, Y, Z, beta,
				 alpha, ZL, zin, z_is_mast , Flip , W)
int nclust, z_is_mast, Flip;
MATRIX **X, **Y, **Z, *beta, *alpha, **ZL, *zin, **W;
{
MATRIX_TRIPLET_P *alpha_comps;
int ni, nc2, r0, c, i, j, k, l, p, q, lz;
double chkrt, a, mui, muj, psi2ij, gamr, tmpm, mu2ij, tmpd, ttt,
   tmpN1, tmpf2d, tmpf2, Zrk, Alphak, dadalk, dNdal,
   Dn, mpr, xil, xjl, dmidbl, dmjdbl, dadbl, dnrdbl, tmpf,
   dDdalk, dndalk, Zij, toff, tmpalogit;
MATRIX *D, *G, *E, *ystar, *mu1, *psi2, *Zeta, *adel1s, *adel2s,
   *aH1, *saH2, *aop, *alvarmat, *thisZ, *DtEG, *td, *thisx, *thisy, *wstar;
MATRIX_TRIPLET_P *ans;
int firstnonsing = 1;
double thisw;

ans = (MATRIX_TRIPLET_P *)calloc( 1, (unsigned)sizeof(struct matrix_triplet_p));
ans->SA = ( MATRIX **)calloc(nclust,(unsigned)sizeof(struct matrix));
if ( ans->SA == (MATRIX **)NULL)
	{
	fprintf(stderr,"MATRIX_TRIPLE_P out of space for SA\n");
	exit(1);
	}

p = X[0]->ncols;

for ( c = 0 ; c < nclust ; c++ )
	{

	ni = X[c]->nrows;
	thisw = MEL(W[c],0,0);

	if ( ni == 1 ) continue;
	

	if ( Flip )
		{
		thisx = flip(X[c]);
		thisy = flip(Y[c]);
		}
	else
		{
		thisx = matcopy(X[c]);
		thisy = matcopy(Y[c]);
		}
	make_permanent( thisx );
	make_permanent( thisy );


	nc2 = tchoose2( ni );

	G = create_matrix( nc2, p, PERMANENT );
	E = create_matrix( nc2, 1, PERMANENT );

	ystar = star(thisy);
	wstar = star(W[c]);
	make_permanent( ystar );
	make_permanent( wstar );

	Zeta = create_matrix( nc2, 1, PERMANENT );

	if (z_is_mast) lz = zin->nrows;
	else lz = Z[c]->nrows;

	if (!Flip)
	{
	if (!(z_is_mast)) thisZ = matcopy( Z[c] );
			else thisZ = pluck_Z( zin, ZL[c] );  /* should pluck */
	make_permanent( thisZ );
	}
	else
	{
	if (!(z_is_mast)) thisZ = matcopy( flipz(Z[c]) );
			else thisZ = pluck_Z( flipz(zin), flipinds(ZL[c],lz) );  /* should pluck */
	make_permanent( thisZ );
	}

	if ( Flip && (c == 826) )
		{
		Matdump( thisx );
		Matdump( thisy );
		Matdump( thisZ );
		}

	q = thisZ->ncols;

	D = create_matrix( nc2, q, PERMANENT );

	mu1 = matantilogit( matmult( thisx, beta ));
	psi2 = matexp( matmult( thisZ, alpha ));


	make_permanent(mu1);

	make_permanent(psi2);

/* formerly canonical indexing: */
/*	for ( i = 0 ; i < (ni-1) ; i++ )
		{
		for ( j = (i+1) ; j < ni ; j++ )
			{  */
/* now need to be redundant */
        for ( i = 0; i < ni; i++ )
                {
                for (j = 0 ; j < ni ; j ++ )
                        {
			if (j == i) continue;
			r0 = nrfun( ni, i, j );  /* zerobased r */

/* printf("r0 = %d\n",r0); */

			psi2ij = MEL( psi2, r0, 0 );

			gamr = log( psi2ij );

			mui = MEL( mu1, i, 0 );
			muj = MEL( mu1, j, 0 );

			mpr = mui*muj;

			a = 1. + (mui+muj)*(psi2ij-1.);

			chkrt = a*a - 4.*psi2ij*(psi2ij-1.)*mpr;

			if (chkrt < 0.) 
			  fprintf(stderr,"Mu2 calcs have invalid sqrt.\n");

			if (fabs(psi2ij-1.) < SMALL )
			  fprintf(stderr,"Mu2 calcs have invalid denom. (psi2ij = %lf)\n",
psi2ij);


			tmpm = a - sqrt(chkrt) ;

			if ( psi2ij == 1.) mu2ij = mpr;

			else mu2ij =
			  ( a - sqrt( chkrt ) )/( 2.*(psi2ij-1.) );

		        tmpd = 1. - ( mui + muj ) + mu2ij ;

			toff = (mui - mu2ij)/ tmpd;

			if ( toff < SMALL )
			  fprintf(stderr,"offset calc log of %f, omitted \n",
						 toff);

			tmpalogit = gamr*MEL(thisy,j,0) + log(toff);

			if ( c < 0 )
			printf("%lf <- tmpalogit\n",tmpalogit);

			MEL(Zeta,r0,0) = exp(tmpalogit)/(1.+exp(tmpalogit));

			tmpN1 = 1./sqrt(chkrt);
			tmpf2d = (mui-mu2ij)*tmpd;
			tmpf2 = (muj-1.)/((mui-mu2ij)*tmpd);
			tmpf = (mui - mu2ij)/tmpd;

			Zij = MEL(Zeta,r0,0);

			for ( k = 0 ; k < q ; k++ )
				{
/* int ooo; */
				Zrk = MEL( thisZ, r0, k );
/* printf("Zrk = %lf\n",Zrk);
if (Zrk < 1.0) scanf("%d\n",&ooo); */
				Alphak = MEL( alpha, k, 0 );
				dadalk = (mui+muj)*Zrk*psi2ij;
				Dn = 2.*(psi2ij-1.);

				dNdal = dadalk - .5*tmpN1*
				   (2.*a*dadalk - 
					4.*mpr*(2.*psi2ij -1.)*Zrk*psi2ij );

				dDdalk = 2.*Zrk*psi2ij;


				if ( fabs(Dn) < SMALL )
					{
					fprintf(stderr,"Alpha comps: Dn small. (%lf)\n",Dn);
					dndalk = 0.;
					}
				else dndalk = (dNdal*Dn-dDdalk*tmpm)/(Dn*Dn);

				MEL( D, r0 , k ) = Zij*(1.-Zij)*
					     (Zrk*MEL(thisy,j,0)+
						dndalk*tmpf2);
				}

			/* */
			for ( l = 0 ; l < p ; l++ )
				{
				xil = MEL(thisx,i,l);
				xjl = MEL(thisx,j,l);

				dmidbl = xil*mui*(1.-mui);
				dmjdbl = xjl*muj*(1.-muj);

				dadbl = (psi2ij-1.)*(dmidbl+dmjdbl);
				dnrdbl = ( dadbl - .5*tmpN1*(2.*a*dadbl-4.*psi2ij*
						     (psi2ij-1.)*
						     (dmidbl*muj+dmjdbl*mui)))/
								(2.*psi2ij-2.);
				
				ttt = Zij*(1.-Zij)*(1./tmpf)*
						((dmidbl-dnrdbl)*(1.-muj)
						+dmjdbl*(mui-mu2ij))/(tmpd*tmpd);
				MEL(G,r0,l) = Zij*(1.-Zij)*(1./tmpf)*
						((dmidbl-dnrdbl)*(1.-muj)
						+dmjdbl*(mui-mu2ij))/(tmpd*tmpd);
				if ( 822 < 0 && c < 828 )
					{
					Printf( (double)r0 );
					Printf( (double)l );
					Printf( ttt );
					}
				}
			/* */
			MEL(E,r0,0) = 1./(Zij*(1.-Zij));
			} /* j */
		} /* i */

	if (c < 0 )
		{
		Matdump(D);
		Matdump(E);
		Matdump(G);
		}

	aop = matxdiagasvec( transp(D), E );

	make_permanent( aop );

	aH1 = matmult( aop , D );

	saH2 = matmult( aop, scalar_times_matrix( thisw, matsub( ystar, Zeta ) ) );

	ans->SA[c] = matcopy(saH2);

	td = matmult( aop, G ); 


	if ( firstnonsing )
		{
		firstnonsing = 0;
		alvarmat = create_matrix( q, q, PERMANENT );
		adel1s = create_matrix( q, q, PERMANENT );
		adel2s = create_matrix( q, 1, PERMANENT );
		DtEG = create_matrix( q, p, PERMANENT );
		}

	adel1s = matadd( adel1s, aH1 );
	adel2s = matadd( adel2s, saH2 );
	DtEG = matadd( DtEG, td );

	if ( 820 < 0 && c < 830 ) 
		{
		printf("%d\n",c);
		Matdump( DtEG );
		Matdump( D );
		Matdump( E );
		Matdump( G );
		}
	
	destroy_matrix(D);
	destroy_matrix(G);
	destroy_matrix(E);
	destroy_matrix(mu1);
	destroy_matrix(thisZ);
	destroy_matrix(psi2);
	destroy_matrix(ystar);
	destroy_matrix(wstar);
	destroy_matrix(Zeta);


	}
make_permanent( adel1s );
make_permanent( adel2s );
make_permanent( DtEG );
ans->X1 = adel1s;
ans->X2 = adel2s;
ans->X3 = DtEG;
return ans;
}

MATRIX *flip( x ) MATRIX *x;  
{
int i, j, k;
MATRIX *o;
o = create_matrix( x->nrows, x->ncols , EPHEMERAL );
j = 0;
for ( i = (x->nrows)-1 ; i > -1 ; i-- )
	{
	for ( k = 0 ; k < x->ncols ; k++ )
		{
		MEL( o, j, k ) = MEL( x, i, k );
		}
	j++;
	}
return o;
}

MATRIX *flipz( z ) MATRIX *z;
{
int n, nr, nc, i, j, zoutrow, zinrow, c;
MATRIX *zf;
nr = z->nrows;
nc = z->ncols;
if (nr == 1) return z;
zf = create_matrix( nr, nc, EPHEMERAL );
n = (int)(.5*(double)(1+sqrt((double)(1+8*nr))));
zoutrow = 0;
for ( i = n-1 ; i>0 ; i-- )
	{
	for ( j = (i-1); j > -1 ; j-- )
		{
		zinrow = nrfun(n,i,j);
		for ( c = 0 ; c < nc ; c++ )
			{
			MEL(zf, zoutrow, c) = MEL( z, zinrow, c );
			}
		zoutrow++;
		}
	}
return zf;
}

MATRIX *make_rvar( nclust, bc, ac )
int nclust;
MATRIX_PAIR_P *bc;
MATRIX_TRIPLET_P *ac;
{
MATRIX *H1, *H2, *H2A, *H2B, *H2O, *ans, *H1I;
int i, p, q;

p = (bc->X1)->ncols;
q = (ac->X1)->ncols;

H2A = create_matrix( p, p, PERMANENT );
H2B = create_matrix( q, q, PERMANENT );
H2O = create_matrix( p, q, PERMANENT );

for ( i = 0 ; i < nclust ; i++ )
	{
	make_permanent( bc->SB[i] );
	make_permanent( ac->SA[i] );
	H2A = matadd( H2A, matmult( bc->SB[i], transp( bc->SB[i] ) ) );
	H2B = matadd( H2B, matmult( ac->SA[i], transp( ac->SA[i] ) ) );
	H2O = matadd( H2O, matmult( bc->SB[i], transp( ac->SA[i] ) ) );
	destroy_matrix( bc->SB[i] );
	destroy_matrix( ac->SA[i] );
	}
make_permanent( H2A );
make_permanent( H2B );
make_permanent( H2O );

H1 = create_matrix( p+q, p+q, PERMANENT );
H2 = create_matrix( p+q, p+q, PERMANENT );

plug( bc->X1, H1, 0, 0 );
plug( ac->X2, H1, 0, p );
plug( ac->X1, H1, p, p );


H1I = luinv(H1);
make_permanent( H1I );

plug( H2A, H2, 0, 0 );
plug( H2O, H2, 0, p );
plug( transp(H2O), H2, p, 0 );
plug( H2B, H2, p, p );

ans = matmult( matmult ( H1I , H2 ), transp(H1I) );
destroy_matrix( H1I );
destroy_matrix( H2A );
destroy_matrix( H2B );
destroy_matrix( H2O );
destroy_matrix( H2 );
return ans;
}

/* should return 1-based location indices for a flipping */
/* if input locs are (1,2,4) for n=4, then should */
/* get back (1,3,4) */

MATRIX *flipinds( zi, n ) MATRIX *zi; int n;
{
MATRIX *tmp;
double n2;
int i;
make_permanent( zi );
tmp = matsub( zi, zi );
n2 = .5*(1.+sqrt(1.+8.*(double)n));
for ( i = 0 ; i < tmp->nrows ; i++ )
	{
	MEL( tmp, i, 0 ) = n2 - MEL( zi, i, 0 ) + 1.;
	}
return matsort( tmp );
}

int comp_mat(), swap_mat();

MATRIX *matsort(mat)
MATRIX *mat;
{
MATRIX *qSmat;
double *qhead;
int nelem, width;

width = 8;

qSmat = matcopy(mat);
nelem = qSmat->nrows * qSmat->ncols;

qhead = qSmat->data;

qsort( (char *) qhead, nelem, width, comp_mat );
return qSmat;
}

int comp_mat(i,j)
double *i, *j;
{
if (*i<*j) return -1;
if (*i>*j) return 1;
return 0;
}


/* main()
{
MATRIX *X, *XL;
X = matread("X");
XL = matread("XL");
matdump(pluck_Z(X,XL)); 
matdump(flipinds(X,(int)*XL->data));
}*/


/* MATRIX *make_irvar(nclust, beta_comps, alpha_comps, falpha_comps)
int nclust;
MATRIX_PAIR_P *beta_comps;
MATRIX_TRIPLET_P *alpha_comps, *falpha_comps;
{} */

MATRIX *make_Rvar( nclust, bc, ac, Ni )
int nclust;
MATRIX *Ni;
MATRIX_PAIR_P *bc;
MATRIX_TRIPLET_P *ac ;
{
MATRIX *H1, *H2, *H2A, *H2B, *H2O, *ans, *H1I;
MATRIX *fH2B, *fH2O;
int i, p, q;

p = (bc->X1)->ncols;
q = (ac->X1)->ncols;

H2A = create_matrix( p, p, PERMANENT );
H2B = create_matrix( q, q, PERMANENT );
H2O = create_matrix( p, q, PERMANENT );

for ( i = 0 ; i < nclust ; i++ )
	{

	if ( MEL(Ni,i,0)<2 ) continue; /* skip singleton */

	make_permanent( bc->SB[i] );
	make_permanent( ac->SA[i] );
	H2A = matadd( H2A, matmult( bc->SB[i], transp( bc->SB[i] ) ) );
	H2B = matadd( H2B, matmult( ac->SA[i], transp( ac->SA[i] ) ) );
	H2O = matadd( H2O, matmult( bc->SB[i], transp( ac->SA[i] ) ) );
	destroy_matrix( bc->SB[i] );
	destroy_matrix( ac->SA[i] );
	}

make_permanent( H2A );
make_permanent( H2B );
make_permanent( H2O );

H1 = create_matrix( p+q, p+q, PERMANENT );
H2 = create_matrix( p+q, p+q, PERMANENT );


plug( bc->X1, H1, 0, 0 );
plug( ac->X3, H1, p, 0 );
plug( ac->X1, H1, p, p );

matdump(H1);


H1I = luinv(H1);
make_permanent( H1I );

Matdump( H1 );
Matdump( H1I );

plug( H2A, H2, 0, 0 );
plug( H2O, H2, 0, p );
plug( transp(H2O), H2, p, 0 );
plug( H2B, H2, p, p );

Matdump( H2 );

ans = matmult( matmult ( H1I , H2 ), transp(H1I) );
destroy_matrix( H1I );
destroy_matrix( H2A );
destroy_matrix( H2B );
destroy_matrix( H2O );
destroy_matrix( H2 );
return ans;
}

MATRIX *make_irvar( nclust, bc, ac, fac, Ni )
int nclust;
MATRIX *Ni;
MATRIX_PAIR_P *bc;
MATRIX_TRIPLET_P *ac, *fac;
{
MATRIX *H1, *H2, *H2A, *H2B, *H2O, *ans, *H1I;
MATRIX *fH2B, *fH2O;
int i, p, q;

p = (bc->X1)->ncols;
q = (ac->X1)->ncols;

H2A = create_matrix( p, p, PERMANENT );
H2B = create_matrix( q, q, PERMANENT );
H2O = create_matrix( p, q, PERMANENT );
fH2B = create_matrix( q, q, PERMANENT );
fH2O = create_matrix( p, q, PERMANENT );

for ( i = 0 ; i < nclust ; i++ )
	{
			  /* NULL test not a very portable approach? */

	if ( MEL(Ni,i,0)<2 ) continue; /* skip singleton */

	make_permanent( bc->SB[i] );
	make_permanent( ac->SA[i] );
	make_permanent( fac->SA[i] );
	H2A = matadd( H2A, matmult( bc->SB[i], transp( bc->SB[i] ) ) );
	H2B = matadd( H2B, matmult( ac->SA[i], transp( ac->SA[i] ) ) );
	H2O = matadd( H2O, matmult( bc->SB[i], transp( ac->SA[i] ) ) );
	fH2B = matadd( fH2B, matmult( fac->SA[i], transp( fac->SA[i] ) ) );
	fH2O = matadd( fH2O, matmult( bc->SB[i], transp( fac->SA[i] ) ) );
	destroy_matrix( bc->SB[i] );
	destroy_matrix( ac->SA[i] );
	destroy_matrix( fac->SA[i] );
	}

H2B = scalar_times_matrix( .5,  matadd( H2B, fH2B ) );
H2O = scalar_times_matrix( .5,  matadd( H2O, fH2O ) );
make_permanent( H2A );
make_permanent( H2B );
make_permanent( H2O );

H1 = create_matrix( p+q, p+q, PERMANENT );
H2 = create_matrix( p+q, p+q, PERMANENT );

Matdump( ac->X3 );
Matdump( fac->X3 );
Matdump( ac->X1 );
Matdump( fac->X1 );

plug( bc->X1, H1, 0, 0 );
plug( scalar_times_matrix( .5,  matadd( ac->X3, fac->X3) ), H1, p, 0 );
plug( scalar_times_matrix( .5, matadd( ac->X1, fac->X1) ), H1, p, p );


H1I = luinv(H1);
make_permanent( H1I );

Matdump( H1 );
Matdump( H1I );

plug( H2A, H2, 0, 0 );
plug( H2O, H2, 0, p );
plug( transp(H2O), H2, p, 0 );
plug( H2B, H2, p, p );

Matdump( H2 );

ans = matmult( matmult ( H1I , H2 ), transp(H1I) );
destroy_matrix( H1I );
destroy_matrix( H2A );
destroy_matrix( H2B );
destroy_matrix( H2O );
/* destroy_matrix( fH2B );
destroy_matrix( fH2O ); */
destroy_matrix( H2 );
return ans;
}

MATRIX_TRIPLET_P *get_alpha_comps_exch( nclust, X, Y, beta,
				 alpha, Flip , W)
int nclust, Flip;
MATRIX **X, **Y, *beta, *alpha, **W;
{
MATRIX_TRIPLET_P *alpha_comps;
int ni, nc2, r0, c, i, j, k, l, p, q, lz;
double chkrt, a, mui, muj, psi2ij, gamr, tmpm, mu2ij, tmpd, ttt,
   tmpN1, tmpf2d, tmpf2, Zrk, Alphak, dadalk, dNdal,
   Dn, mpr, xil, xjl, dmidbl, dmjdbl, dadbl, dnrdbl, tmpf,
   dDdalk, dndalk, Zij, toff, tmpalogit, Dij, Eij, DtE, Gijl;
MATRIX *D, *G, *E, *ystar, *mu1, *psi2, *Zeta, *adel1s, *adel2s,
   *aH1, *saH2, *aop, *alvarmat, *thisZ, *DtEG, *td, *thisx, *thisy ;
MATRIX_TRIPLET_P *ans;
int firstnonsing = 1;
double thisw;

ans = (MATRIX_TRIPLET_P *)calloc( 1, (unsigned)sizeof(struct matrix_triplet_p));
ans->SA = ( MATRIX **)calloc(nclust,(unsigned)sizeof(struct matrix));
if ( ans->SA == (MATRIX **)NULL)
	{
	fprintf(stderr,"MATRIX_TRIPLE_P out of space for SA\n");
	exit(1);
	}

p = X[0]->ncols;
q = 1;

adel1s = create_matrix( q, q, PERMANENT );
adel2s = create_matrix( q, 1, PERMANENT );
DtEG = create_matrix( q, p, PERMANENT );

for ( c = 0 ; c < nclust ; c++ )
	{
thisw = MEL(W[c],0,0);
	ni = X[c]->nrows;
	if ( ni == 1 ) continue;  /* skip a singleton */

	if ( Flip )
		{
		thisx = flip(X[c]);
		thisy = flip(Y[c]);
		}
	else
		{
		thisx = matcopy(X[c]);
		thisy = matcopy(Y[c]);
		}
	make_permanent( thisx );
	make_permanent( thisy );

	ans->SA[c] = create_matrix( 1, 1, PERMANENT );



	q = 1;

	mu1 = matantilogit( matmult( thisx, beta ));
	psi2 = matexp( alpha ); /* it is "scalar" in exch case */

	make_permanent(mu1);
	make_permanent(psi2);

	for ( i = 0 ; i < ni ; i++ )
		{
		for ( j = 0 ; j < ni ; j++ )
			{
			if (i == j) continue;
			psi2ij = MEL( psi2, 0, 0 );

			gamr = log( psi2ij );

			mui = MEL( mu1, i, 0 );
			muj = MEL( mu1, j, 0 );

			mpr = mui*muj;

			a = 1. + (mui+muj)*(psi2ij-1.);

			chkrt = a*a - 4.*psi2ij*(psi2ij-1.)*mpr;

			if (chkrt < 0.) 
			  fprintf(stderr,"Mu2 calcs have invalid sqrt.\n");

			if (fabs(psi2ij-1.) < SMALL )
			  fprintf(stderr,"Mu2 calcs have invalid denom.\n");

			tmpm = a - sqrt(chkrt) ;

			if ( psi2ij == 1.) mu2ij = mpr;

			else mu2ij =
			  ( a - sqrt( chkrt ) )/( 2.*(psi2ij-1.) );

		        tmpd = 1. - ( mui + muj ) + mu2ij ;

			toff = (mui - mu2ij)/ tmpd;

			if ( toff < SMALL )
			  fprintf(stderr,"offset calc log of %f, omitted \n",
						 toff);

			tmpalogit = gamr*MEL(thisy,j,0) + log(toff);

			if ( c < 0 )
			printf("%lf <- tmpalogit\n",tmpalogit);

			Zij = exp(tmpalogit)/(1.+exp(tmpalogit));

			tmpN1 = 1./sqrt(chkrt);
			tmpf2d = (mui-mu2ij)*tmpd;
			tmpf2 = (muj-1.)/((mui-mu2ij)*tmpd);
			tmpf = (mui - mu2ij)/tmpd;

			Alphak = MEL( alpha, 0, 0 );
			dadalk = (mui+muj)*psi2ij;
			Dn = 2.*(psi2ij-1.);

			dNdal = dadalk - .5*tmpN1*
			   (2.*a*dadalk - 4.*mpr*(2.*psi2ij -1.)*psi2ij );

			dDdalk = 2.*psi2ij;

			if ( fabs(Dn) < SMALL )
				{
				fprintf(stderr,"Dn small.\n");
				dndalk = 0.;
				}
			else dndalk = (dNdal*Dn-dDdalk*tmpm)/(Dn*Dn);

			Dij = Zij*(1.-Zij)*(MEL(thisy,j,0)+
					dndalk*tmpf2);

			Eij = 1./(Zij*(1.-Zij));

			DtE = Dij*Eij;
	
			for ( l = 0 ; l < p ; l++ )
				{
				xil = MEL(thisx,i,l);
				xjl = MEL(thisx,j,l);
	
				dmidbl = xil*mui*(1.-mui);
				dmjdbl = xjl*muj*(1.-muj);
	
				dadbl = (psi2ij-1.)*(dmidbl+dmjdbl);
				dnrdbl = ( dadbl - .5*tmpN1*(2.*a*dadbl-4.*psi2ij*
						     (psi2ij-1.)*
						     (dmidbl*muj+dmjdbl*mui)))/
								(2.*psi2ij-2.);
				
				Gijl = Zij*(1.-Zij)*(1./tmpf)*
						((dmidbl-dnrdbl)*(1.-muj)
						+dmjdbl*(mui-mu2ij))/(tmpd*tmpd);
				MEL(DtEG,0,l) = MEL(DtEG,0,l) + DtE*Gijl;
				}


			MEL(adel1s,0,0) = MEL(adel1s,0,0) + DtE*Dij;
			MEL(adel2s,0,0) = MEL(adel2s,0,0) + DtE*thisw*(MEL(thisy,i,0)-Zij);
			MEL(ans->SA[c],0,0) = MEL(ans->SA[c],0,0)+DtE*thisw*(MEL(thisy,i,0)-Zij);
			     /* collapse the score contributions */
			} /* j */
		} /* i */

	}
make_permanent( adel1s );
make_permanent( adel2s );
make_permanent( DtEG );
ans->X1 = adel1s;
ans->X2 = adel2s;
ans->X3 = DtEG;
return ans;
}

MATRIX_PAIR_P *get_beta_comps_exch( nclust, X, Y, beta, 
				    alpha , bindep, W)
int nclust, bindep;
MATRIX **X, **Y, *beta, *alpha, **W;
{
int c, ni, nc2, i,j, r0, p;
MATRIX *thisZ, *mu1, *B, *C, *op, *del1, *del2, *varmat, 
    *del1s, *del2s, *S;
double psi2ij, mui, muj, a, chkrt, mu2, thisw;
MATRIX_PAIR_P *ans;
MATWCOND *Binvc;
double mincond = 1.;
int condcount = 0;

p = X[0]->ncols;

ans = (MATRIX_PAIR_P *)calloc( 1, (unsigned)sizeof(struct matrix_pair_p));
ans->SB = ( MATRIX **)calloc(nclust,(unsigned)sizeof(struct matrix));
if ( ans->SB == (MATRIX **)NULL)
	{
	fprintf(stderr,"MATRIX_PAIR_P out of space for SB\n");
	exit(1);
	}

for ( c = 0 ; c < nclust ; c++ )
	{
	thisw = MEL(W[c],0,0);

	ni = X[c]->nrows;

	mu1 = matantilogit( matmult( X[c], beta ) );

	make_permanent(mu1);

	if ( !bindep )
	   {
	   B = create_matrix( ni, ni, PERMANENT );
	

	for ( i = 0 ; i < ni ; i++ )
		{
		for ( j = 0 ; j < ni ; j++ )
			{
			if (i == j) continue;

			psi2ij = exp( MEL(alpha,0,0) );

			mui = MEL( mu1, i, 0 );
			muj = MEL( mu1, j, 0 );

			a = 1. + (mui+muj)*(psi2ij-1.);

			chkrt = a*a - 4.*psi2ij*(psi2ij-1.)*mui*muj;

			if (chkrt < 0.) 
			  fprintf(stderr,"Mu2 calcs have invalid sqrt.\n");

			if (fabs(psi2ij-1.) < SMALL )
			  fprintf(stderr,"Mu2 calcs have invalid denom.\n");

			if ( psi2ij == 1.) mu2 = mui*muj;
			else mu2 =
			  ( a - sqrt( chkrt ) )/( 2.*(psi2ij-1.) );

			MEL(B,i,j) = (mu2 - mui*muj);  /* while it is around */
			MEL(B,j,i) = MEL(B,i,j);
			}
			     
		}
	
	for ( i = 0 ; i < ni ; i++ )
		MEL(B,i,i) = MEL( mu1, i, 0 )*(1. - MEL( mu1, i, 0 ) );
	}
/* note that "B" is already inverted if "bindep" is used */
	else
	   {
	   B = create_matrix( ni, 1, PERMANENT );
           for ( i = 0 ; i < ni ; i++ )
                MEL(B,i,0) = 1./(MEL( mu1, i, 0 )*(1. - MEL( mu1, i, 0 ) ));
	   }

	C = px1_times_pxq( 
		 px1_times_pxq( mu1, oneminus( mu1 ) ) , X[c] );

	make_permanent(C);

	/* matdump( corner( B , 5, 5 ) ); */

	/* printf("[%d]B cond: %lf\n", c, estcond(B) ); */

	if ( !bindep )
	   op = matmult( transp(C), sweep(B) );  /* code restored 5-29-2011 */
	   /* op = matmult( transp(C), sweep(B) ); */
/* code dropped 5-29-2011, see invwcond commented out above 
           {
           Binvc = invwcond(B);
           if (Binvc->condition < mincond) mincond = Binvc->condition;
           if ( mincond < .03 && condcount < 6 )
               {
	       condcount++ ;
               fprintf(stderr,
                        "get_beta_comps: full weighting has low condition #\n");
               fprintf(stderr,
                        "get_beta_comps: cluster %d, cond(dgeco) = %lf\n", c,
                                                           Binvc->condition );
               fprintf(stderr,
                        "get_beta_comps: consider independence model for beta\n");
		if ( condcount > 4 ) 
		   fprintf(stderr,"No more condition msgs this iter.\n");
               }
           op = matmult( transp(C), Binvc->mat );
           }
*/

	else  /* B already inverted */
	   op = matxdiagasvec( transp(C), B );

	make_permanent(op);


	del1 =  matmult( op, C );

	S = scalar_times_matrix( thisw, matsub( Y[c] , mu1 ) );

	del2 = matmult( op, S );

	ans->SB[c] = matcopy( del2 );

	make_permanent( del1 );
	make_permanent( del2 );

	if ( c == 0 )
		{
		del1s = create_matrix( p, p, PERMANENT );
		del2s = create_matrix( p, 1, PERMANENT );
		}

	del1s = matadd( del1s, del1) ;

	del2s = matadd( del2s, del2);

	destroy_matrix(mu1);
	destroy_matrix(B);
	destroy_matrix(C);
	destroy_matrix(op);
	destroy_matrix(del1);
	destroy_matrix(del2);
	}
make_permanent( del1s );
make_permanent( del2s );
ans->X1 = del1s;
ans->X2 = del2s;
return ans;
}

MATRIX_PAIR_P *get_beta_comps_2cl( nclust, X, Y, beta, 
				    alpha , CL2, bindep, W)
int nclust, bindep;
MATRIX **X, **Y, *beta, *alpha, **CL2, **W;
{
int c, ni, nc2, i,j, r0, p;
MATRIX *thisZ, *mu1, *B, *C, *op, *del1, *del2, *varmat, 
    *del1s, *del2s, *S;
double psi2ij, mui, muj, a, chkrt, mu2;
MATRIX_PAIR_P *ans;
MATWCOND *Binvc;
double mincond = 1.;
int condcount = 0;
double thisw;

p = X[0]->ncols;

ans = (MATRIX_PAIR_P *)calloc( 1, (unsigned)sizeof(struct matrix_pair_p));
ans->SB = ( MATRIX **)calloc(nclust,(unsigned)sizeof(struct matrix));
if ( ans->SB == (MATRIX **)NULL)
	{
	fprintf(stderr,"MATRIX_PAIR_P out of space for SB\n");
	exit(1);
	}

for ( c = 0 ; c < nclust ; c++ )
	{
	printf("C%d %c",c, ((c % 10)==9) ? '\n' : ' '); 
/*	printf("CLUST %d/%d\n",c,nclust); */
	thisw = MEL(W[c],0,0);

	ni = X[c]->nrows;

	mu1 = matantilogit( matmult( X[c], beta ) );

	make_permanent(mu1);

	if (!bindep)
	{
	B = create_matrix( ni, ni, PERMANENT );

	for ( i = 0 ; i < ni ; i++ )
		{
		for ( j = 0 ; j < ni ; j++ )
			{
			if (i == j) continue;

			if ( MEL(CL2[c],i,0) != MEL(CL2[c],j,0) )
				psi2ij = exp( MEL(alpha,0,0) );
			else psi2ij = exp( MEL(alpha,0,0) + MEL(alpha,1,0) );

			mui = MEL( mu1, i, 0 );
			muj = MEL( mu1, j, 0 );

			a = 1. + (mui+muj)*(psi2ij-1.);

			chkrt = a*a - 4.*psi2ij*(psi2ij-1.)*mui*muj;

			if (chkrt < 0.) 
			  fprintf(stderr,"Mu2 calcs have invalid sqrt.\n");

			if (fabs(psi2ij-1.) < SMALL )
			  fprintf(stderr,"Mu2 calcs have invalid denom.\n");

			if ( psi2ij == 1.) mu2 = mui*muj;
			else mu2 =
			  ( a - sqrt( chkrt ) )/( 2.*(psi2ij-1.) );

			MEL(B,i,j) = mu2 - mui*muj;  /* while it is around */
			MEL(B,j,i) = MEL(B,i,j);
			}
			     
		}
	for ( i = 0 ; i < ni ; i++ )
		MEL(B,i,i) = MEL( mu1, i, 0 )*(1. - MEL( mu1, i, 0 ) );
	}
	else
	{
	B = create_matrix( ni, 1, PERMANENT );
        for ( i = 0 ; i < ni ; i++ )
                MEL(B,i,0) = 1./(MEL( mu1, i, 0 )*(1. - MEL( mu1, i, 0 ) ));
	}

	C = px1_times_pxq( 
		 px1_times_pxq( mu1, oneminus( mu1 ) ) , X[c] );

	make_permanent(C);

	if ( !bindep )
	   op = matmult( transp(C), sweep(B) );  /* restored 5-29-2011 */
/* droped owing to linpack interface problem
           {
           Binvc = invwcond(B);
           if (Binvc->condition < mincond) mincond = Binvc->condition;
           if ( mincond < .03 && condcount < 6)
               {
	       condcount++;
               fprintf(stderr,
                        "get_beta_comps: full weighting has low condition #\n");
               fprintf(stderr,
                        "get_beta_comps: cluster %d, cond(dgeco) = %lf\n", c,
                                                           Binvc->condition );
               fprintf(stderr,
                        "get_beta_comps: consider independence model for beta\n");
		if ( condcount > 4 ) 
		   fprintf(stderr,"No more condition msgs this iter.\n");
               }
           op = matmult( transp(C), Binvc->mat );
           }
*/

	else
	   op = matxdiagasvec( transp(C), B );

	make_permanent(op);

	del1 = matmult( op, C );

	S = scalar_times_matrix( thisw, matsub( Y[c] , mu1 ) );
	/* S = matsub( Y[c] , mu1 ); */

	del2 = matmult( op, S );
	ans->SB[c] = matcopy( del2 );

	make_permanent( del1 );
	make_permanent( del2 );

	if ( c == 0 )
		{
		/* varmat = create_matrix( p, p, PERMANENT ); */
		del1s = create_matrix( p, p, PERMANENT );
		del2s = create_matrix( p, 1, PERMANENT );
		}

	/* varmat = matadd( varmat, matmult( del2, transp(del2) ) ); */

	del1s = matadd( del1s, del1 );

	del2s = matadd( del2s, del2 );

	destroy_matrix(mu1);
	destroy_matrix(B);
	destroy_matrix(C);
	destroy_matrix(op);
	destroy_matrix(del1);
	destroy_matrix(del2);
	}
make_permanent( del1s );
make_permanent( del2s );
ans->X1 = del1s;
ans->X2 = del2s;
return ans;
}

MATRIX_TRIPLET_P *get_alpha_comps_2cl( nclust, X, Y, beta,
				 alpha, Flip , CL2, W)
int nclust, Flip;
MATRIX **X, **Y, *beta, *alpha, **CL2, **W;
{
MATRIX_TRIPLET_P *alpha_comps;
int ni, nc2, r0, c, i, j, k, l, p, q, lz;
double chkrt, a, mui, muj, psi2ij, gamr, tmpm, mu2ij, tmpd, ttt,
   tmpN1, tmpf2d, tmpf2, Zrk, Alphak, dadalk, dNdal,
   Dn, mpr, xil, xjl, dmidbl, dmjdbl, dadbl, dnrdbl, tmpf,
   dDdalk, dndalk, Zetaij, toff, tmpalogit, Dij, Eij, DtE, Gijl;
double Z2, dadalk2, dNdal2, dDdalk2, dndalk2, Dij2, DtE2;

MATRIX *D, *G, *E, *ystar, *mu1, *psi2, *Zeta, *adel1s, *adel2s,
   *aH1, *saH2, *aop, *alvarmat, *thisZ, *DtEG, *td, *thisx, *thisy, 
   *thiscl2;
MATRIX_TRIPLET_P *ans;
int firstnonsing = 1, samesub=0;
double thisw;

ans = (MATRIX_TRIPLET_P *)calloc( 1, (unsigned)sizeof(struct matrix_triplet_p));
ans->SA = ( MATRIX **)calloc(nclust,(unsigned)sizeof(struct matrix));
if ( ans->SA == (MATRIX **)NULL)
	{
	fprintf(stderr,"MATRIX_TRIPLE_P out of space for SA\n");
	exit(1);
	}

p = X[0]->ncols;
q = 2;

adel1s = create_matrix( q, q, PERMANENT );
adel2s = create_matrix( q, 1, PERMANENT );
DtEG = create_matrix( q, p, PERMANENT );

for ( c = 0 ; c < nclust ; c++ )
	{
	ni = X[c]->nrows;
thisw = MEL(W[c],0,0);
	if ( ni == 1 ) continue;  /* skip a singleton */

	/*printf("%d ",c); */
	printf("A%d %c",c, ((c % 10)==9) ? '\n' : ' '); 
	if ( Flip )
		{
		thisx = flip(X[c]);
		thisy = flip(Y[c]);
		thiscl2 = flip(CL2[c]);
		}
	else
		{
		thisx = matcopy(X[c]);
		thisy = matcopy(Y[c]);
		thiscl2 = matcopy(CL2[c]);
		}
	make_permanent( thisx );
	make_permanent( thisy );

	ans->SA[c] = create_matrix( q, 1, PERMANENT );



	mu1 = matantilogit( matmult( thisx, beta ));

	make_permanent(mu1);

	for ( i = 0 ; i < ni ; i++ )
		{
		for ( j = 0 ; j < ni ; j++ )
			{
			if ( i == j ) continue;

			if ( MEL(thiscl2,i,0) == MEL(thiscl2,j,0) )
				samesub = 1;
			else samesub = 0;

			if ( !samesub )
				psi2ij = exp( MEL(alpha,0,0) );
			else psi2ij = exp( MEL(alpha,0,0) + MEL(alpha,1,0) );

			gamr = log( psi2ij );

			mui = MEL( mu1, i, 0 );
			muj = MEL( mu1, j, 0 );

			mpr = mui*muj;

			a = 1. + (mui+muj)*(psi2ij-1.);

			chkrt = a*a - 4.*psi2ij*(psi2ij-1.)*mpr;

			if (chkrt < 0.) 
			  fprintf(stderr,"Mu2 calcs have invalid sqrt.\n");

			if (fabs(psi2ij-1.) < SMALL )
			  fprintf(stderr,"Mu2 calcs have invalid denom.\n");

			tmpm = a - sqrt(chkrt) ;

			if ( psi2ij == 1.) mu2ij = mpr;

			else mu2ij =
			  ( a - sqrt( chkrt ) )/( 2.*(psi2ij-1.) );

		        tmpd = 1. - ( mui + muj ) + mu2ij ;

			toff = (mui - mu2ij)/ tmpd;

			if ( toff < SMALL )
			  fprintf(stderr,"offset calc log of %f, omitted \n",
						 toff);

			tmpalogit = gamr*MEL(thisy,j,0) + log(toff);

			if ( c < 0 )
			printf("%lf <- tmpalogit\n",tmpalogit);

			Zetaij = exp(tmpalogit)/(1.+exp(tmpalogit));

			tmpN1 = 1./sqrt(chkrt);
			tmpf2d = (mui-mu2ij)*tmpd;
			tmpf2 = (muj-1.)/((mui-mu2ij)*tmpd);
			tmpf = (mui - mu2ij)/tmpd;

			if (samesub) Z2 = 1.;
			    else     Z2 = 0;

			Alphak = MEL( alpha, 0, 0 );
			dadalk = (mui+muj)*psi2ij;
			dadalk2 = (mui+muj)*psi2ij*Z2;
			Dn = 2.*(psi2ij-1.);

			dNdal = dadalk - .5*tmpN1*
			   (2.*a*dadalk - 4.*mpr*(2.*psi2ij -1.)*psi2ij );
			dNdal2 = dadalk2 - .5*tmpN1*
			   (2.*a*dadalk2 - 4.*mpr*(2.*psi2ij -1.)*psi2ij*Z2 );

			dDdalk = 2.*psi2ij;
			dDdalk2 = 2.*psi2ij*Z2;

			if ( fabs(Dn) < SMALL )
				{
				fprintf(stderr,"Dn small.\n");
				dndalk = 0.;
				}
			else dndalk = (dNdal*Dn-dDdalk*tmpm)/(Dn*Dn);

			if ( fabs(Dn) < SMALL )
				{
				fprintf(stderr,"Dn small.\n");
				dndalk2 = 0.;
				}
			else dndalk2 = (dNdal2*Dn-dDdalk2*tmpm)/(Dn*Dn);

			Dij = Zetaij*(1.-Zetaij)*(MEL(thisy,j,0)+
					dndalk*tmpf2);

			Dij2 = Zetaij*(1.-Zetaij)*(Z2*MEL(thisy,j,0)+
					dndalk2*tmpf2);


			Eij = 1./(Zetaij*(1.-Zetaij));

			DtE = Dij*Eij;
			DtE2 = Dij2*Eij;
			
			for ( l = 0 ; l < p ; l++ )
				{
				xil = MEL(thisx,i,l);
				xjl = MEL(thisx,j,l);
	
				dmidbl = xil*mui*(1.-mui);
				dmjdbl = xjl*muj*(1.-muj);
	
				dadbl = (psi2ij-1.)*(dmidbl+dmjdbl);
				dnrdbl = ( dadbl - .5*tmpN1*(2.*a*dadbl-4.*psi2ij*
						     (psi2ij-1.)*
						     (dmidbl*muj+dmjdbl*mui)))/
								(2.*psi2ij-2.);
				
				Gijl = Zetaij*(1.-Zetaij)*(1./tmpf)*
						((dmidbl-dnrdbl)*(1.-muj)
						+dmjdbl*(mui-mu2ij))/(tmpd*tmpd);
				MEL(DtEG,0,l) = MEL(DtEG,0,l) + DtE*Gijl;
				MEL(DtEG,1,l) = MEL(DtEG,1,l) + DtE2*Gijl;
				/* if ( c < 2 && i < 3 && j < 4 )
					{
					printf("DtE %lf, DtE2 %lf\n",DtE, DtE2);
					matdump(DtEG);
					} */
				}
		
			MEL(adel1s,0,0) = MEL(adel1s,0,0) + DtE*Dij;
			MEL(adel1s,0,1) = MEL(adel1s,0,1) + DtE*Dij2;
			MEL(adel1s,1,1) = MEL(adel1s,1,1) + DtE2*Dij2;
			/* it is finished below */


			MEL(adel2s,0,0) = MEL(adel2s,0,0) + DtE*thisw*(MEL(thisy,i,0)-Zetaij);
			MEL(adel2s,1,0) = MEL(adel2s,1,0) + DtE2*thisw*(MEL(thisy,i,0)-Zetaij);

			MEL(ans->SA[c],0,0) = MEL(ans->SA[c],0,0)+DtE*thisw*(MEL(thisy,i,0)-Zetaij);
			MEL(ans->SA[c],1,0) = MEL(ans->SA[c],1,0)+DtE2*thisw*(MEL(thisy,i,0)-Zetaij);
			     /* collapse the score contributions */
			} /* j */
		} /* i */
		destroy_matrix(thisx);
		destroy_matrix(thisy);
		destroy_matrix(thiscl2);
		destroy_matrix(mu1);
	}
make_permanent( adel1s );
make_permanent( adel2s );
make_permanent( DtEG );
MEL(adel1s,1,0)=MEL(adel1s,0,1);
ans->X1 = adel1s;
ans->X2 = adel2s;
ans->X3 = DtEG;
return ans;
}

MATWCOND *invwcond( x )
MATRIX *x;
{
MATWCOND *out;
MATRIX *tx, *save;
double *head, *z, cond, *det, *work;
int nr, *ipvt, job, nargs;

cond = 0.;
out = (MATWCOND *)calloc( 1, (unsigned)sizeof(struct matwcond) );
out->mat = (MATRIX *) calloc( 1, (unsigned)sizeof(struct matrix) );

save = matcopy(x);

tx = transp(x);
nr = tx->nrows;
ipvt = (int *)calloc( nr, (unsigned)sizeof(int) );
z = (double *)calloc( nr, (unsigned)sizeof(double) );
head = tx->data;

nargs = 6;
dgecoXXY_( head, &nr, &nr, ipvt, &cond, z, &nargs);

det = (double *)calloc( 2, (unsigned)sizeof(double) );
work = (double *)calloc( nr, (unsigned)sizeof(double) );

/*

job = 1;
nargs = 7;

dgediXXY_( head, &nr, &nr, ipvt, det, work, &job, &nargs );
*/

out->condition = cond;
save = sweep(save);  
out->mat = save;

return out;
}

double estcond( x )
MATRIX *x;
{
MATRIX *tx;
double *head, *z, cond;
int nr, *ipvt, nargs;

tx = transp(x);
nr = tx->nrows;
ipvt = (int *)calloc( nr, (unsigned)sizeof(int) );
z = (double *)calloc( nr, (unsigned)sizeof(double) );
head = tx->data;
nargs = 6;
dgecoXXY_( head, &nr, &nr, ipvt, &cond, z , &nargs);
destroy_matrix( tx );
return cond;
}

int P( n ) int n;
	{ return ((n*(n+1))/2); }

int nrfun( n, i, j ) int n,i,j;
	{

/* > rs
function(n, i, j)
(n - 1) * i + ifelse(j > i, j - 1, j)
*/

return (n-1)*i + ( j>i ? j-1 : j );

}

