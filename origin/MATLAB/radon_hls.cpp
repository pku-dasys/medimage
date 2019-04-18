#include <math.h>

#define PI 3.14159265358979323846264
#define imageHeight         512
#define imageWidth          512
#define measurementHeight   512
#define measurementWidth    180
#define edgeHeight          512
#define edgeWidth           512
#define NAlpha              180

int iround(double x)
{
    return((int)floor(x+0.5));
}

//***************************************************************************

double x2Grid(double x)
{
	return (x+0.5)*(double)(imageWidth-1);
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

void RotateImage(double ISrc[imageHeight][imageWidth], double IDst[imageHeight][imageWidth], double a)
{
	int nx, ny, indSrc, indDst, NxMin, NxMax, rx, ry;
	double px1, py1, px2, py2, sina, cosa, fx, fy, px, py;

	sina = sin(a);
	cosa = cos(a);
    
	for (ny=0; ny<imageHeight; ny++)
	{
		//
		// Find out from where to where we have to compute the rotation
		//
        
		// Set (p1,p2) to the backwards-rotated line segment (-0.5,y)->(0.5,y)
        double xy_y = ((double)ny / (double)(imageHeight-1)) - 0.5;
		px1 = -0.5*cosa + xy_y*sina;
		py1 =  0.5*sina + xy_y*cosa;
        
		px2 =  0.5*cosa + xy_y*sina;
		py2 = -0.5*sina + xy_y*cosa;
        
		// Clip the line segment against the square [-0.5,0.5] x [-0.5,0.5]
		Clip(&px1, &py1, &px2, &py2);
		
		// Rotate back (only x-coordinate)
		NxMin = iround(x2Grid(px1*cosa - py1*sina));
		NxMax = iround(x2Grid(px2*cosa - py2*sina));
        
		//mxAssert((NxMin>=0) && (NxMax<Nx), NULL);
        
		// printf("%i:%i->%i\n", ny, NxMin, NxMax);
        
        for (nx=0; nx<imageWidth; nx++)
            IDst[ny][nx] = 0;
        
		for (nx=NxMin; nx<=NxMax; nx++)
		{
            double xy_x = ((double)nx / (double)(imageWidth-1)) - 0.5;
			px = x2Grid(  xy_x*cosa + xy_y*sina);
			py = x2Grid(- xy_x*sina + xy_y*cosa);
            
			if (px<0) px=0; if (px>imageWidth-1.00001) px=imageWidth-1.00001;
			if (py<0) py=0; if (py>imageHeight-1.00001) py=imageHeight-1.00001;
            
			rx = (int)floor(px);
			ry = (int)floor(py);
            
			fx = px-rx;
			fy = py-ry;
            
			//printf("  %i %i\n", rx, ry);
            
			//mxAssert((rx>=0) && (ry>=0) && (rx<Nx-1) && (ry<Ny-1), NULL);
			//mxAssert((fx>=0) && (fx<=1) && (fy>=0) && (fy<=1), NULL);
            
			IDst[ny][nx] = (ISrc[ry][rx]*(1-fx) + ISrc[ry][rx+1]*fx) * (1-fy) +
            (ISrc[ry+1][rx]*(1-fx) + ISrc[ry+1][rx+1]*fx) *    fy;
		}
		//printf("\n");
	}
}

void A(double f[imageHeight][imageWidth], double g[measurementHeight][measurementWidth])
{
	int nAlpha, nx, ny;
    double fRot[imageHeight][imageWidth];

	for (nAlpha=0; nAlpha<NAlpha; nAlpha++)
	{
		double a = nAlpha * PI / 180.0;
        
		// Compute a rotated version of f
		RotateImage(f, fRot, a);
        
		// Compute the output image g by integrating over the rows of fRot
		for (nx=0; nx<imageWidth; nx++)
		{
			g[nx][nAlpha] = 0;
			for (ny=0; ny<imageHeight; ny++)
				g[nx][nAlpha] += fRot[ny][nx];
		}
	}

}

void AStar(double f[imageHeight][imageWidth], double g[measurementHeight][measurementWidth]){
    
    int nAlpha, nx, ny, n;
    double lines[imageHeight][imageWidth];
	for (nAlpha=0; nAlpha<NAlpha; nAlpha++)
	{
		double a = nAlpha * PI / 180.0;
        
		// Spread nAlpha-th row of g over lines
		for (nx=0; nx<imageWidth; nx++)
		{
			for (ny=0; ny<imageHeight; ny++)
				lines[ny][nx] = g[nx][nAlpha];
		}
        
		// Rotate lines BACKWARDS and add it to fOut
		RotateImage(lines, f, -a);
	}
}

/*
 % This function minimizes the Ambrosio-Tortorelli functional for a fixed
 % v in the f variable.
 %
 % The functional is approximately minimized by the method of steepest descent
 % with a fixed number of iteration steps.
 %
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Input/Variables
 % g                      l x k measurement
 % f                      n x m reconstruction
 % v                      n x m edge indicator
 % alpha                  weight for smoothnes penalty
 % A                      linear operator given as a function handle
 % AStar                  adjoint operator of A given as a function handle
 % IterImage              number of descend steps
 %
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */

void SD(double g[measurementHeight][measurementWidth], double f[imageHeight][imageWidth], double v[edgeHeight][edgeWidth], double alpha, double epsilon, int IterImage)
{
    int i, j, k, nAlpha;
    double d[imageHeight][imageWidth], tmp[measurementHeight][measurementWidth];
    double c_1, nV_KepsNablad, dSkalard, normAd;
    for (k=1;k<=IterImage;k++){
        A(f, tmp);
        for (i=0;i<measurementHeight;i++)
            for (j=0;j<measurementWidth;j++)
                tmp[i][j] = g[i][j] - tmp[i][j];
        AStar(d, tmp);
        
        //% d = - \nabla_f AT(f,v);
        //% d = A*(g - Af) + 2alpha div(v^2 \nabla f) +0.01epsilon^2 LaPlace f:
        nV_KepsNablad = 0;
        dSkalard = 0;
        for (i=0;i<imageHeight;i++)
            for(j=0;j<imageWidth;j++)
            {
                double v_center = v[i][j];
                double f_center = f[i][j];
                double v_left = j>0 ? v[i][j-1] : 1;
                double v_up = i>0 ? v[i-1][j]: 1;
                double f_left = j>0 ? f[i][j-1] : 0;
                double f_right = j<imageWidth-1 ? f[i][j+1] : 0;
                double f_up = i>0 ? f[i-1][j]:0;
                double f_down = i<imageHeight-1 ? f[i+1][j]:0;
                d[i][j] = d[i][j] + alpha * (v_center*v_center*(f_down-f_center))+v_center*v_center*(f_right-f_center)-v_up*v_up*(f_center-f_up)-v_left*v_left*(f_center-f_left)+0.01*epsilon*epsilon*(f_down+f_right+f_up+f_left-4*f_center);
                nV_KepsNablad += v_center * v_center + 0.01*epsilon*epsilon*((f_center-f_up)*(f_center-f_up)+(f_center-f_left)*(f_center-f_left));
                dSkalard += d[i][j]*d[i][j];
            }
        A(d, tmp);
        normAd = 0;
        for (i=0;i<measurementHeight;i++)
            for (j=0;j<measurementWidth;j++)
                normAd += tmp[i][j];
        c_1         = dSkalard/(normAd + alpha*nV_KepsNablad);
        for (i=0;i<imageHeight;i++)
            for (j=0;j<imageWidth;j++)
                f[i][j] += c_1 * d[i][j];
    }
}
