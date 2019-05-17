#include "c_imp.h"

#include <string.h>
#include <stdlib.h>
#include <math.h>

void swap(double *a,double *b) {
    double t = *a; *a = *b; *b = t;
}

double padded_f[imageHight+2][imageWidth+2];

void medfilt2(double f[imageHight][imageWidth]) {
    memset(padded_f, 0, sizeof(padded_f));
	int i,j;
    for (i = 1; i<=imageHight; ++i)
        for (j = 1; j<=imageWidth; ++j)
            padded_f[i][j] = f[i-1][j-1];
    for (i = 1; i<=imageHight; ++i)
        for (j = 1; j<=imageWidth; ++j) {
            double t[9];
            t[0] = padded_f[i][j];
            t[1] = padded_f[i-1][j-1];
            t[2] = padded_f[i-1][j];
            t[3] = padded_f[i-1][j+1];
            t[4] = padded_f[i][j+1];
            t[5] = padded_f[i+1][j+1];
            t[6] = padded_f[i+1][j];
            t[7] = padded_f[i+1][j-1];
            t[8] = padded_f[i][j-1];
			int p,q;
            for (p = 0; p<9; ++p)
                for (q = p+1; q<9; ++q)
                    if (t[p]>t[q]) swap(&t[p],&t[q]);
            f[i-1][j-1] = t[4];
        }
}

double sqr(double x) {
    return x*x;
}

void randn(double x[nAlpha][imageHight]) {
    double u[nAlpha*imageHight+1];
	int i,j;
    for (i = 0; i<nAlpha*imageHight+1; ++i)
        u[i] = (double)rand()/RAND_MAX;
    int k = 0;
    double square, amp, angle;
    for (i = 0; i<nAlpha; ++i)
        for (j = 0; j<imageHight; ++j) {
            if (k%2==0) {
                square = -2. * log(u[k]);
                if (square<0) square = 0;
                amp = sqrt(square);
                angle = 2. * M_PI * u[k+1];
                x[i][j] = amp*sin(angle);
            }
            else
                x[i][j] = amp*cos(angle);
            ++k;
        }
}

int iround(double x) {
    return((int)floor(x+0.5));
}

double x2Grid(double x) {
    return (x+0.5)*(double)(imageHight-1);
}

void Clip(double *px1, double *py1, double *px2, double *py2) {
    if (*px1<-0.5) { *py1 = *py2 + (*py1-*py2) * (-0.5-*px2) / (*px1-*px2); *px1=-0.5; }
    if (*py1<-0.5) { *px1 = *px2 + (*px1-*px2) * (-0.5-*py2) / (*py1-*py2); *py1=-0.5; }
    if (*px1> 0.5) { *py1 = *py2 + (*py1-*py2) * ( 0.5-*px2) / (*px1-*px2); *px1= 0.5; }
    if (*py1> 0.5) { *px1 = *px2 + (*px1-*px2) * ( 0.5-*py2) / (*py1-*py2); *py1= 0.5; }

    if (*px2<-0.5) { *py2 = *py1 + (*py2-*py1) * (-0.5-*px1) / (*px2-*px1); *px2=-0.5; }
    if (*py2<-0.5) { *px2 = *px1 + (*px2-*px1) * (-0.5-*py1) / (*py2-*py1); *py2=-0.5; }
    if (*px2> 0.5) { *py2 = *py1 + (*py2-*py1) * ( 0.5-*px1) / (*px2-*px1); *px2= 0.5; }
    if (*py2> 0.5) { *px2 = *px1 + (*px2-*px1) * ( 0.5-*py1) / (*py2-*py1); *py2= 0.5; }
}

void RotateImage(double ISrc[imageHight][imageWidth], double IDst[imageHight][imageWidth],
                   double a,double xy[imageHight]) {
    int nx, ny, NxMin, NxMax, rx, ry;
    double px1, py1, px2, py2, sina, cosa, fx, fy, px, py;

    memset(IDst, 0, imageHight*imageWidth*sizeof(double));
    sina = sin(a);
    cosa = cos(a);

    for (ny=0; ny<imageWidth; ny++) {
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

        // printf("%i:%i->%i\n", ny, NxMin, NxMax);

        for (nx=NxMin; nx<=NxMax; nx++)
        {
            px = x2Grid(  xy[nx]*cosa + xy[ny]*sina);
            py = x2Grid(- xy[nx]*sina + xy[ny]*cosa);

            if (px<0) px=0; if (px>imageHight-1.00001) px=imageHight-1.00001;
            if (py<0) py=0; if (py>imageWidth-1.00001) py=imageWidth-1.00001;

            rx = (int)floor(px);
            ry = (int)floor(py);

            fx = px-rx;
            fy = py-ry;

            IDst[nx][ny] = (ISrc[rx][ry]*(1-fx)+ISrc[rx+1][ry]*fx)*(1-fy)
                           +(ISrc[rx][ry+1]*(1-fx)+ISrc[rx+1][ry+1]*fx)*fy;
        }
        
        //printf("\n");
    }
}
