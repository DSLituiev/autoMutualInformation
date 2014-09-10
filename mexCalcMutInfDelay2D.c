/*******************************************************************
A MEX-adapted version of a C program "minfo.c" by Eric R. Weeks

Adaptation: Dmytro S. Lituiev, March 2013, University of Zurich, Switzerland

Handles a 2D arrays, with the signals running along the first dimension.
i.e. in an array [ T x N ] -- N time series, each with T time points

**************************************************
minfo.c -- Eric R. Weeks -- started 2/28/97

does the mutual information algorithm discussed by Fraser & Swinney
(Phys Rev A 33 (1986) p1134-1140)

v01:  2-28-97: taken from shell.c (7/1/96)
			quicksort routine taken from sane.c (original version)
v02:  2-28-97: revised sorting for s[] (different than q[])
			sped up math
v03:  2-28-97: add in tau loop
	 3-01-97: fix for variable number of input; add -b option
v04:  3-01-97: take out chi2 tests for substructure
v05:  3-01-97: realize that with chi2 tests taken out, number()
			function doesn't need to be called very often.  remove
			a[] and b[][] arrays!  Much faster now.

This program is public domain, although please leave my name and
email address attached.

email: weeks@physics.emory.edu
web: http://www.physics.emory.edu/~weeks/

explanation of how to use this program:
    http://www.physics.emory.edu/~weeks/software/minfo.html

 *******************************************************************/
#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.14159265358979323846264338328
#define EE 2.71828182845904523536
#define MAXNUM 100000
#define KMAX 25
#define IS_REAL_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && \
mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_SCALAR(P) (IS_REAL_2D_FULL_DOUBLE(P) && mxGetNumberOfElements(P) == 1)

 mwIndex  s[MAXNUM],q[MAXNUM];
 mwIndex  Mx, n0, mmax;
double pop[MAXNUM];
 mwIndex  pindex[MAXNUM];
 mwIndex  pow2[KMAX];
 mwIndex  bin;
 mwIndex  globalivalue[4];
double logtwo;

 mwIndex  number( mwIndex  karray2[], mwIndex  m2);
double findisq();
double ffunct( mwIndex  kmarray[], mwIndex  m);
void saneqsort( mwIndex  p, mwIndex  r);
 mwIndex  qpartition_neurons( mwIndex  p, mwIndex  r);

// mexCalcMutInfDelay
void mexFunction( mwIndex  nlhs, mxArray *plhs[],  mwIndex  nrhs, const mxArray *prhs[])
 {
 /* Macros for the ouput and input arguments */
    #define B_OUT plhs[0]
	#define ARRAY_IN prhs[0]
	#define TAU_IN prhs[1]
	#define BIN_IN prhs[2]
	double *B,  *array, info;
	mwIndex  Nx, nn, i, j, tau, tauMax, Mx;
	if(nrhs < 1 || nrhs > 3) /* Check the number of arguments */
	    mexErrMsgTxt("Wrong number of input arguments.");
	else if(nlhs > 1)
		mexErrMsgTxt("Too many output arguments.");
		
	if(!IS_REAL_2D_FULL_DOUBLE(ARRAY_IN)) /* Check A */
		mexErrMsgTxt("A must be a real 2D full double array.");
	
	// get the input array dimensions and the pointer
	Mx = (mwIndex) mxGetM(ARRAY_IN); /* Get the dimensions */
	Nx = (mwIndex) mxGetN(ARRAY_IN);
	array = mxGetPr(ARRAY_IN);  /* Get the pointer to the data of A */
	
	if(nrhs < 2) /* If p is unspecified, set it to a default value */
	  	tauMax = min(20, Mx-5);
	else /* If P was specified, check that it is a real double scalar */
		if(!IS_REAL_SCALAR(TAU_IN) )
		        mexErrMsgTxt("tau_max must be a real double scalar.");
		else
       			tauMax = ( mwIndex ) mxGetScalar( TAU_IN ); /* Get tauMax */
	if(nrhs < 3) /* If p is unspecified, set it to a default value */
	  	bin = -1;
	else /* If P was specified, check that it is a real double scalar */
		if(!IS_REAL_SCALAR(BIN_IN) )
		    mexErrMsgTxt("tau_max must be a real double scalar.");
		else
			bin = ( mwIndex )mxGetData( BIN_IN ); /* Get bin number */


	B_OUT = mxCreateDoubleMatrix( (mwSize)tauMax+1, (mwSize)Nx, mxREAL); /* Create the output matrix */
	B = mxGetPr(B_OUT); /* Get the pointer to the data of B */
 	
	/*----- done reading in data -----*/

	logtwo = 1.0/log(2.0);

	// mexPrintf("Nx:\t%u\tM:\t%u\n\n", Nx, Mx);
	
for (nn = 0; nn<Nx; nn++) { // a loop through the second dimension

	n0 = 1;
	while ((n0+tauMax)<= Mx)  n0 *= 2;
	n0 /= 2;

	/* n0 = Mx - tauMax; */
	j = n0;
	for (i=0; i<KMAX;i++)  {
		pow2[i] = j;
		j /= 2;
	}
	// mmax = (( mwIndex )(log(((float)n0))*logtwo+0.1));
	// mexPrintf( "n0: %d mmax: %d\n", n0, mmax);

	for (i=0; i<n0;i++)  {
		pop[i] = array[i + (Mx)*nn];
		// mexPrintf("Nx*nn\t%u\tu:\t%u\tpop:\t%f\n", (Mx)*nn, i, pop[i]);
		pindex[i] = i;		
	}
	saneqsort(0, n0-1);
	/* for (i=0;i<Mx-tauMax;i++)  s[pindex[i]] = i; */
	/* WARNING!!!!!! Note that this definition is somewhat opposite
	 * of what 'makes sense' but will make number() routine faster */
	for (i=0;i<n0;i++)  s[i] = pindex[i]; // s[i + (tauMax+1)*nn] = pindex[i + (tauMax+1)*nn];  //

	for (tau=0;tau<=tauMax;tau++)  {
		/* now do tau offset for q[] array */
		for (i=tau; i<tau+n0; i++)  {
			pop[i-tau] = array[i + (Mx)*nn];
			pindex[i-tau] = i-tau;			
		}
		
		saneqsort(0, n0-1);
		
		for (i=0;i<n0;i++)		q[pindex[i]] = i;
		/* for (i=0;i<n0;i++)  fprintf(stderr,"%d ",q[i]); */
		/* fprintf(stderr,"\n"); */


		/* assume at this time that s[], q[] contain integers from
		 * 0 to 2^Nx-1.   q[] is based on the time lagged data from
		 * array[].  s[] is just array[].
		 *
		 * example: s[] = 0 4 5 3 6 1 2 7
		 *          q[] = 7 4 5 3 6 2 0 1
		 *
		 *     o.......      so the first partition is s: 0-3 | 4-7
		 *     ......o.                                q: 0-3 | 4-7
		 *     .....o..
		 *     ....o...      with 3 in LL, 1 in UL, 3 in UR, 1 in LR
		 *     ...o....      quadrants.
		 *     .o......
		 *     .......o      second partition is s: 01 | 23 | 45 | 67
		 *     ..o.....                          q: 01 | 23 | 45 | 67
		 *                   with distribution: 1001
		 *				                    0020
		 *								1100
		 *								0101
		 */

		/* now find I(S,Q) according to formula (19) */
		info = findisq();		
		// mexPrintf("tau = %u\n", tau + (tauMax+1)*nn );
        B[ tau + (tauMax+1)*nn ] = info;
		// printf("%d %f\n",tau,info);
		// fprintf(stderr,"%d %f\n",tau,info);
	}
}
	// exit(0);
}	/* END OF MAIN */



double findisq()
{
	double info;
	mwIndex  kmarray[KMAX];
	double x,y;

	kmarray[0] = 0;
	// mexPrintf("kmarray= %u\n", kmarray );
	x = ((double) n0);
	y = ffunct(kmarray,0);
	info = (1.0/x)*y - log(x)*logtwo;

	return info;
}

double ffunct( mwIndex  kmarray[], mwIndex  m)
{
	/* THIS FUNCTION CAN CALL ITSELF RECURSIVELY */
	double value;
	 mwIndex  n,j;
	 mwIndex  temparray[KMAX];

	for (j=0;j<=m;j++)  temparray[j] = kmarray[j];

	n = number(temparray,m);
	
	value = ((double) n);
	if (n<=1)  {
		value = 0.0;
	} else if (n==2)  {
		value = 4.0;
	} else if (m==bin)  {
		/* no substructure */
		value = value*log(value)*logtwo;
	} else {
		/* assume substructure exists */
		value = value*2.0;
		for (j=0;j<=3;j++)  {
			temparray[m+1] = j;
			value += ffunct(temparray,m+1);
		}
	}

	return value;
}


 mwIndex  number( mwIndex  karray2[], mwIndex  m2)
{
	/* THIS FUNCTION IS NOT RECURSIVE */
	 mwIndex  ivalue;
	 mwIndex  los,his,loq,hiq;
	register  mwIndex  i,j;

	if (m2>0)  {
		los = 0;loq = 0;
		his = n0; hiq = n0;
		for (i=1;i<=m2;i++)  {
			if (karray2[i]%2==0)  his -= pow2[i];
			else                  los += pow2[i];
			if (karray2[i]<2)     hiq -= pow2[i];
			else                  loq += pow2[i];
		}
		ivalue = 0;
		for (i=los;i<his;i++)  {
			j = q[s[i]];
			if ((j>=loq)&&(j<hiq))  ivalue++;
		}
	} else {
		ivalue = n0;
	}

	return ivalue;
}



/* quicksort, taken from sane.c by David E. Moriarty */
/* modified for this purpose */

void saneqsort(p, r)
   mwIndex  p,r;
{
   mwIndex  q;
  if (p<r) {
    q = qpartition_neurons(p,r);
    saneqsort(p, q);
    saneqsort(q+1, r);
  }
}

/** partition function for saneqsort.
     **/

 mwIndex  qpartition_neurons(p,r)
   mwIndex  p,r;
{
  register  mwIndex  i,j;
  double x,temp;
   mwIndex  tempi;

  x = pop[p];
  i = p - 1;
  j = r + 1;
  while(1) {
    do{
      --j;
    }while(pop[j] < x);
    do{
      ++i;
    }while(pop[i] > x);
    if (i < j) {
	 /* here's where the in-place swap takes place */
	 /* fortunately pindex[] stores the *original* location */
      temp = pop[i]; pop[i] = pop[j]; pop[j] = temp;
	 tempi = pindex[i]; pindex[i] = pindex[j]; pindex[j] = tempi;
    }
    else
      return j;
  }
}

