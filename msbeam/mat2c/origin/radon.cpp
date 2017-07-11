/* KOMPILIERUNG DER MEX-ROUTINE:
	 WINDOWS-VERSION

   OPTIMIZED VERSION:
     mex -O radon.cpp;

	 DEBUG VERSION
     mex -g radon.cpp;

*/

#include <math.h>
#include "mex.h"
#include "memory.h"

// Input Arguments
#define	mexF					prhs[0]				// Activity map
#define mexAlpha	    prhs[1]				// Vector of rotation angles (in DEGREES)

// Output Arguments
#define	mexGOut				plhs[0]				// (Output) radon transform

#define PI 3.14159265358979323846264

// Variables
int Nx, Ny, Nxy, NAlpha;
double *alpha, *f, *fRot, *g, *xy;

//***************************************************************************
//**********************   M A I N   R O U T I N E S   **********************
//***************************************************************************

int iround(double x)
{
    return((int)floor(x+0.5));
}

//***************************************************************************

double x2Grid(double x)
{
	return (x+0.5)*(double)(Nx-1);
}

//***************************************************************************

void MatlabData2C(int nrhs, const mxArray *prhs[])
{
    int n;

  // Check for the proper number of arguments
  if ((nrhs<2) || (nrhs>2))
    mexErrMsgTxt("Wrong number of input arguments.\n");

  // Get size of input
	Nx  = mxGetM(mexF);
	Ny  = mxGetN(mexF);
	Nxy = Nx*Ny;

	if (Nxy==0)
		mexErrMsgTxt("Activity matrix F is empty.\n");

	if (Nx!=Ny)
		mexErrMsgTxt("Activity matrix F must be square.\n");

	f  = mxGetPr(mexF);

	// Get rotation angle vector
	NAlpha = mxGetM(mexAlpha)*mxGetN(mexAlpha);
	if (NAlpha==0)
		mexErrMsgTxt("Rotation angle vector is empty.\n");

	alpha = mxGetPr(mexAlpha);

	// Get memory for other variables
	xy = (double*)mxMalloc(Nx*sizeof(double));
	g  = (double*)mxMalloc(NAlpha*Nx*sizeof(double));
	
	fRot = (double*)mxMalloc(Nx*Ny*sizeof(double));

	for (n=0; n<Nx; n++)
		xy[n] = ((double)n / (double)(Nx-1)) - 0.5;
}

//***************************************************************************

void Clip(double *px1, double *py1, double *px2, double *py2)
{
	// Clip (px1,py1)->(px2,py2) against the square [-0.5,0.5] x [-0.5,0.5]

	if (*px1<-0.5) { *py1 = *py2 + (*py1-*py2) * (-0.5-*px2) / (*px1-*px2); *px1=-0.5; }
    if (*py1<-0.5) { *px1 = *px2 + (*px1-*px2) * (-0.5-*py2) / (*py1-*py2); *py1=-0.5; }
	if (*px1> 0.5) { *py1 = *py2 + (*py1-*py2) * ( 0.5-*px2) / (*px1-*px2); *px1= 0.5; }
    if (*py1> 0.5) { *px1 = *px2 + (*px1-*px2) * ( 0.5-*py2) / (*py1-*py2); *py1= 0.5; }

	if (*px2<-0.5) { *py2 = *py1 + (*py2-*py1) * (-0.5-*px1) / (*px2-*px1); *px2=-0.5; }
    if (*py2<-0.5) { *px2 = *px1 + (*px2-*px1) * (-0.5-*py1) / (*py2-*py1); *py2=-0.5; }
	if (*px2> 0.5) { *py2 = *py1 + (*py2-*py1) * ( 0.5-*px1) / (*px2-*px1); *px2= 0.5; }
    if (*py2> 0.5) { *px2 = *px1 + (*px2-*px1) * ( 0.5-*py1) / (*py2-*py1); *py2= 0.5; }
}

//***************************************************************************

void RotateImage(double *ISrc, double *IDst, double a)
{
	int nx, ny, indSrc, indDst, NxMin, NxMax, rx, ry;
	double px1, py1, px2, py2, sina, cosa, fx, fy, px, py;

	memset(IDst, 0, Nx*Ny*sizeof(double));
	sina = sin(a);
	cosa = cos(a);

	for (ny=0; ny<Ny; ny++)
	{
		//
		// Find out from where to where we have to compute the rotation
		//

		// Set (p1,p2) to the backwards-rotated line segment (-0.5,y)->(0.5,y)
		px1 = -0.5*cosa + xy[ny]*sina;
		py1 =  0.5*sina + xy[ny]*cosa;

		px2 =  0.5*cosa + xy[ny]*sina;
		py2 = -0.5*sina + xy[ny]*cosa;

		// Clip the line segment against the square [-0.5,0.5] x [-0.5,0.5]
		Clip(&px1, &py1, &px2, &py2);
		
		// Rotate back (only x-coordinate)
		NxMin = iround(x2Grid(px1*cosa - py1*sina));
		NxMax = iround(x2Grid(px2*cosa - py2*sina));

		mxAssert((NxMin>=0) && (NxMax<Nx), NULL);

		// printf("%i:%i->%i\n", ny, NxMin, NxMax);

		indDst = ny*Nx + NxMin;
		for (nx=NxMin; nx<=NxMax; nx++, indDst++)
		{
			px = x2Grid(  xy[nx]*cosa + xy[ny]*sina);
			py = x2Grid(- xy[nx]*sina + xy[ny]*cosa);

			if (px<0) px=0; if (px>Nx-1.00001) px=Nx-1.00001;
			if (py<0) py=0; if (py>Ny-1.00001) py=Ny-1.00001;

			rx = (int)floor(px);
			ry = (int)floor(py);

			fx = px-rx;
			fy = py-ry;

			//printf("  %i %i\n", rx, ry);

			mxAssert((rx>=0) && (ry>=0) && (rx<Nx-1) && (ry<Ny-1), NULL);
			mxAssert((fx>=0) && (fx<=1) && (fy>=0) && (fy<=1), NULL);

			indSrc = ry*Nx+rx;

			//printf("%2i,%2i|", rx, ry);
			IDst[indDst] = (ISrc[indSrc   ]*(1-fx) + ISrc[indSrc   +1]*fx) * (1-fy) + 
										 (ISrc[indSrc+Nx]*(1-fx) + ISrc[indSrc+Nx+1]*fx) *    fy;
		}
		
		//printf("\n");
	}
}

//***************************************************************************

void Radon(mxArray *plhs[])
{
	int nAlpha, nx, ny, indG, ind;

	indG = 0;
	for (nAlpha=0; nAlpha<NAlpha; nAlpha++)
	{
		double a = alpha[nAlpha] * PI / 180.0;

		// Compute a rotated version of f
		RotateImage(f, fRot, a);

		// Compute the output image g by integrating over the rows of fRot
		indG = nAlpha;
		for (nx=0; nx<Nx; nx++)
		{
			g[indG] = 0;
			ind = nx;

			for (ny=0; ny<Ny; ny++, ind+=Nx)
				g[indG] += fRot[ind];

			indG += NAlpha;
		}
	}

	// Write output data
  mexGOut = mxCreateDoubleMatrix(NAlpha,Nx,mxREAL);
	memcpy(mxGetPr(mexGOut), g, NAlpha*Nx*sizeof(double));


        mxFree(xy);
        mxFree(g);
        mxFree(fRot);
}

//***************************************************************************
//***************************************************************************
//***************************************************************************

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
  // Load signals
  MatlabData2C(nrhs, prhs);

	// Main computational routine
	Radon(plhs);

  return;  
}

