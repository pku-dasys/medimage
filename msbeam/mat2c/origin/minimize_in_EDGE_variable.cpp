#include "mex.h"

#include <cstdio>
#include <cmath>
#include <cstring>

using namespace std;

double *f, *v;
int imageHight, imageWidth, IterEdge;
double alpha, beta, epsilon, ED_TOL;

double *vupdate, *vold;
double *g, *d, *laplace_v, *d_xf, *d_yf, *n, *tempZeros, *tempOnes, *padded_f;
double nNablaD, nD, c;

// function vupdate = minimize_in_EDGE_variable(f,v,imageHight,imageWidth,alpha,beta,epsilon,ED_TOL,IterEdge)
void loadData(int nrhs, const mxArray *prhs[]) {
    f = mxGetPr(prhs[0]);
    v = mxGetPr(prhs[1]);

    imageHight = mxGetScalar(prhs[2]);
    imageWidth = mxGetScalar(prhs[3]);
    alpha = mxGetScalar(prhs[4]);
    beta = mxGetScalar(prhs[5]);
    epsilon = mxGetScalar(prhs[6]);
    ED_TOL = mxGetScalar(prhs[7]);
    IterEdge = mxGetScalar(prhs[8]);

    //mexPrintf("H=%d W=%d alpha=%f beta=%f eps=%f ED=%f Iter=%d\n",
    //          imageHight, imageWidth, alpha, beta, epsilon, ED_TOL, IterEdge);

    vupdate = (double*)mxMalloc(imageHight*imageWidth*sizeof(double));
    vold = (double*)mxMalloc(imageHight*imageWidth*sizeof(double));
}

void zeros(double *&a,int x,int y) {
    a = (double*)mxMalloc(x*y*sizeof(double));
    memset(a, 0, x*y*sizeof(double));
}

void ones(double *&a,int x,int y) {
    a = (double*)mxMalloc(x*y*sizeof(double));
    for (int i = 0; i<x; ++i)
        for (int j = 0; j<y; ++j)
            a[i*x+j] = 1.0;
}

double sqr(double x) {
    return x*x;
}

void minimize_in_EDGE_variable(mxArray *plhs[]) {
    // vupdate = v;
    for (int i = 0; i<imageHight; ++i)
        for (int j = 0; j<imageWidth; ++j)
            vupdate[i*imageHight+j] = v[i*imageHight+j];

    // g                 = zeros(imageHight,imageWidth); % gradient \nabla_v AT(f,v);
    zeros(g, imageHight, imageWidth);
    zeros(d, imageHight, imageWidth); // WARNING

    // laplace_v          = zeros(imageHight,imageWidth);
    zeros(laplace_v, imageHight, imageWidth);

    // d_xf              = zeros(imageHight,imageWidth);
    zeros(d_xf, imageHight, imageWidth);

    // d_yf              = zeros(imageHight,imageWidth);
    zeros(d_yf, imageHight, imageWidth);

    // n                 = zeros(imageHight,imageWidth); % auxillary variable
    zeros(n, imageHight, imageWidth);

    // tempZeros = zeros(imageHight+2, imageWidth+2);
    zeros(tempZeros, imageHight+2, imageWidth+2);

    // tempOnes  = ones(imageHight+2, imageWidth+2);
    ones(tempOnes, imageHight+2, imageWidth+2);

    // padded_f  = zeros(imageHight+2, imageWidth+2);
    zeros(padded_f, imageHight+2, imageWidth+2);

    // padded_f(2:imageHight+1,2:imageWidth+1)  = f;
    for (int i = 1; i<=imageHight; ++i)
        for (int j = 1; j<=imageWidth; ++j) {
            padded_f[i*(imageHight+2)+j] = f[(i-1)*imageHight+(j-1)];
        }

    for (int k = 1; k<=IterEdge; ++k) {
        memcpy(vold, vupdate, imageHight*imageWidth*sizeof(double));
        for (int i = 1; i<=imageHight; ++i)
            for (int j = 1; j<=imageWidth; ++j)
                tempOnes[i*(imageHight+2)+j] = vold[(i-1)*imageHight+(j-1)];
        for (int i = 1; i<=imageHight; ++i)
            for (int j = 1; j<=imageWidth; ++j) {
                d_xf[(i-1)*imageHight+(j-1)] = 0.5*(padded_f[(i+1)*(imageHight+2)+j]-padded_f[(i-1)*(imageHight+2)+j]);
                d_yf[(i-1)*imageHight+(j-1)] = 0.5*(padded_f[i*(imageHight+2)+(j+1)]-padded_f[i*(imageHight+2)+(j-1)]);
                n[(i-1)*imageHight+(j-1)] = sqr(d_xf[(i-1)*imageHight+(j-1)])+sqr(d_yf[(i-1)*imageHight+(j-1)]);
                laplace_v[(i-1)*imageHight+(j-1)] = 0.25*(tempOnes[(i+1)*(imageHight+2)+j]
                                                          +tempOnes[i*(imageHight+2)+(j+1)]
                                                          -4*tempOnes[i*(imageHight+2)+j]
                                                          +tempOnes[(i-1)*(imageHight+2)+j]
                                                          +tempOnes[i*(imageHight+2)+(j-1)]);
                g[(i-1)*imageHight+(j-1)] = alpha*tempOnes[i*(imageHight+2)+j]*n[(i-1)*imageHight+(j-1)]
                                            -beta*epsilon*(laplace_v[(i-1)*imageHight+(j-1)])
                                            +beta*(tempOnes[i*(imageHight+2)+j]-1)/(4*epsilon);
            }

        for (int i = 0; i<imageHight; ++i)
            for (int j = 0; j<imageWidth; ++j)
                d[i*imageHight+j] = -g[i*imageHight+j];
        
        for (int i = 1; i<=imageHight; ++i)
            for (int j = 1; j<=imageWidth; ++j)
                tempZeros[i*(imageHight+2)+j] = d[(i-1)*imageHight+(j-1)];
        nNablaD = 0;
        for (int i = 1; i<=imageHight; ++i)
            for (int j = 1; j<=imageWidth; ++j)
                nNablaD = nNablaD+sqr(0.5*(tempZeros[(i+1)*(imageHight+2)+j]-tempZeros[(i-1)*(imageHight+2)+j]))
                                 +sqr(0.5*(tempZeros[i*(imageHight+2)+(j+1)]-tempZeros[i*(imageHight+2)+(j-1)]));

        nD = 0;
        for (int i = 0; i<imageHight; ++i)
            for (int j = 0; j<imageWidth; ++j)
                nD += sqr(d[i*imageHight+j]);
        // c = nD/(alpha*sum(sum(n.*d.^2)) + beta*epsilon*nNablaD + (beta/(4*epsilon))*nD);
        double tmp = 0;
        for (int i = 0; i<imageHight; ++i)
            for (int j = 0; j<imageWidth; ++j)
                tmp += n[i*imageHight+j]*sqr(d[i*imageHight+j]);
        c = nD/(alpha*tmp+beta*epsilon*nNablaD+(beta/(4*epsilon))*nD);
        //mexPrintf("It=%d nD=%f c=%f\n",k,(double)nD,(double)c);

        double maxCD = 0;
        for (int i = 0; i<imageHight; ++i)
            for (int j = 0; j<imageWidth; ++j) {
                vupdate[i*imageHight+j] = vold[i*imageHight+j]+c*d[i*imageHight+j];
                double CD = fabs(c*d[i*imageHight+j]);
                if (CD>maxCD) maxCD = CD;
            }
        
        for (int i = 0; i<imageHight; ++i)
            for (int j = 0; j<imageWidth; ++j) {
                if (vupdate[i*imageHight+j]>1) vupdate[i*imageHight+j] = 1;
                if (vupdate[i*imageHight+j]<0) vupdate[i*imageHight+j] = 0;
            }

        if (k==IterEdge) {
            mexPrintf("Min Edge: maxIt reached\n");
            break;
        }
        //mexPrintf("maxCD=%f\n",maxCD);
        if (maxCD<ED_TOL) {
            mexPrintf("Min Edge: ED_Tol reached, %d Iterations \n", k);
            break;
        }
    }

    
    plhs[0] = mxCreateDoubleMatrix(imageHight,imageWidth,mxREAL);
    memcpy(mxGetPr(plhs[0]), vupdate, imageHight*imageWidth*sizeof(double));
}

//***************************************************************************
//***************************************************************************
//***************************************************************************

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    loadData(nrhs, prhs);
    minimize_in_EDGE_variable(plhs);
    return;
}

