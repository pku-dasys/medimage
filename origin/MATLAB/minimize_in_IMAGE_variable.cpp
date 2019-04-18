#include "mex.h"

#include <cstdio>
#include <cmath>
#include <cstring>

using namespace std;

double *g, *f, *v;
int gM, gN;
int imageHight, imageWidth, IterImage;
double alpha, epsilon, IM_TOL;

double *tempZeros, *padded_v, *Af, *Ad, *r, *d, *p, *q, *tempF, *myF;
double AT_f_old, nV_KepsNablap, c_1, AT_f_new, pSkalard, normAp, AdSkalarAp, c_2;

// function fupdate = minimize_in_IMAGE_variable(g,f,v,imageHight,imageWidth,alpha,epsilon,IM_TOL,IterImage)
void loadData(int nrhs, const mxArray *prhs[]) {
    g = mxGetPr(prhs[0]);

    gM = mxGetM(prhs[0]);
    gN = mxGetN(prhs[0]);
        
    f = mxGetPr(prhs[1]);
    v = mxGetPr(prhs[2]);

    imageHight = mxGetScalar(prhs[3]);
    imageWidth = mxGetScalar(prhs[4]);
    alpha = mxGetScalar(prhs[5]);
    epsilon = mxGetScalar(prhs[6]);
    IM_TOL = mxGetScalar(prhs[7]);
    IterImage = mxGetScalar(prhs[8]);
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

void opA(double *f,double *Af) {
    mxArray *IN[1];
    mxArray *OUT[1];
    IN[0] = mxCreateDoubleMatrix(imageHight, imageWidth, mxREAL);
    memcpy(mxGetPr(IN[0]), f, imageHight*imageWidth*sizeof(double));
    mexCallMATLAB(1,OUT,1,IN,"Anm");
    memcpy(Af, mxGetPr(OUT[0]), gM*gN*sizeof(double));
}

void opAS(double *f,double *ASf) {
    mxArray *IN[1];
    mxArray *OUT[1];
    IN[0] = mxCreateDoubleMatrix(gM, gN, mxREAL);
    memcpy(mxGetPr(IN[0]), f, gM*gN*sizeof(double));
    mexCallMATLAB(1,OUT,1,IN,"AStern");
    memcpy(ASf, mxGetPr(OUT[0]), imageHight*imageHight*sizeof(double));
}

void domedfilt2(double *f, double *mf) {
    mxArray *IN[1];
    mxArray *OUT[1];
    IN[0] = mxCreateDoubleMatrix(imageHight, imageWidth, mxREAL);
    memcpy(mxGetPr(IN[0]), f, imageHight*imageHight*sizeof(double));
    mexCallMATLAB(1,OUT,1,IN,"medfilt2");
    memcpy(mf, mxGetPr(OUT[0]), imageHight*imageHight*sizeof(double));
}

void minimize_in_IMAGE_variable(mxArray *plhs[]) {
    myF = (double*)mxMalloc(imageHight*imageWidth*sizeof(double));
    memcpy(myF, f, imageHight*imageWidth*sizeof(double));

    zeros(tempZeros, imageHight+2, imageWidth+2);

    ones(padded_v, imageHight+2, imageWidth+2);
    for (int i = 1; i<=imageHight; ++i)
        for (int j = 1; j<=imageWidth; ++j) {
            padded_v[i*(imageHight+2)+j] = v[(i-1)*imageHight+(j-1)];
        }

    Af = (double*)mxMalloc(gM*gN*sizeof(double));
    r = (double*)mxMalloc(gM*gN*sizeof(double));
    // r = g-A(f);
    opA(myF,Af);
    for (int i = 0; i<gM; ++i)
        for (int j = 0; j<gN; ++j)
            r[j*gM+i] = g[j*gM+i]-Af[j*gM+i];
            
    for (int i = 1; i<=imageHight; ++i)
        for (int j = 1; j<=imageWidth; ++j) {
            tempZeros[i*(imageHight+2)+j] = myF[(i-1)*imageHight+(j-1)];
        }
        
    d = (double*)mxMalloc(imageHight*imageWidth*sizeof(double));
    opAS(r,d);
    
    FILE *fid = fopen("my_d.txt","w");
    for (int i = 0; i<imageHight; ++i) {
        for (int j = 0; j<imageWidth; ++j) {
            fprintf(fid,"%f ",d[i*imageHight+j]);
        }
        fprintf(fid,"\n");
    }
    fclose(fid);
    
    for (int i = 1; i<=imageHight; ++i)
        for (int j = 1; j<=imageWidth; ++j) {
            d[(i-1)*imageHight+(j-1)] += alpha*(
                +sqr(padded_v[i*(imageHight+2)+j])*(tempZeros[(i+1)*(imageHight+2)+j]-tempZeros[i*(imageHight+2)+j])
                +sqr(padded_v[i*(imageHight+2)+j])*(tempZeros[i*(imageHight+2)+(j+1)]-tempZeros[i*(imageHight+2)+j])
                -sqr(padded_v[(i-1)*(imageHight+2)+j])*(tempZeros[i*(imageHight+2)+j]-tempZeros[(i-1)*(imageHight+2)+j])
                -sqr(padded_v[i*(imageHight+2)+(j-1)])*(tempZeros[i*(imageHight+2)+j]-tempZeros[i*(imageHight+2)+(j-1)])
                +0.01*sqr(epsilon)*(
                    tempZeros[(i+1)*(imageHight+2)+j]
                    +tempZeros[i*(imageHight+2)+(j+1)]
                    -4*tempZeros[i*(imageHight+2)+j]
                    +tempZeros[(i-1)*(imageHight+2)+j]
                    +tempZeros[i*(imageHight+2)+(j-1)]
                )
            );
        }
    AT_f_old = 0;
    for (int i = 0; i<gM; ++i)
        for (int j = 0; j<gN; ++j)
            AT_f_old += sqr(Af[i*gM+j]-g[i*gM+j]);
    mexPrintf("AT_f_old 1: %f\n",AT_f_old);
    
    for (int i = 1; i<=imageHight; ++i)
        for (int j = 1; j<=imageWidth; ++j) {
            tempZeros[i*(imageHight+2)+j] = myF[(i-1)*imageHight+(j-1)];
        }
    for (int i = 2; i<=imageHight-1; ++i)
        for (int j = 2; j<=imageWidth-1; ++j)
            AT_f_old += alpha*(sqr(padded_v[i*(imageHight+2)+j])+0.01*sqr(epsilon))*(
                sqr(tempZeros[i*(imageHight+2)+j]-tempZeros[(i-1)*(imageHight+2)+j])
                +sqr(tempZeros[i*(imageHight+2)+j]-tempZeros[i*(imageHight+2)+(j-1)])
            );
    mexPrintf("AT_f_old 2: %f\n",AT_f_old);
    
    p = (double*)mxMalloc(imageHight*imageWidth*sizeof(double));
    memcpy(p, d, imageHight*imageWidth*sizeof(double));
       
    q = (double*)mxMalloc(gM*gN*sizeof(double));
    
    tempF = (double*)mxMalloc(imageHight*imageWidth*sizeof(double));
    Ad = (double*)mxMalloc(gM*gN*sizeof(double));
    
    for (int k = 1; k<=IterImage; ++k) {
        opA(p,q);
        
        for (int i = 1; i<=imageHight; ++i)
            for (int j = 1; j<=imageWidth; ++j) {
                tempZeros[i*(imageHight+2)+j] = p[(i-1)*imageHight+(j-1)];
            }
        nV_KepsNablap = 0;
        for (int i = 1; i<=imageHight; ++i)
            for (int j = 1; j<=imageWidth; ++j) {
                nV_KepsNablap += (sqr(padded_v[i*(imageHight+2)+j])+0.01*sqr(epsilon))*(
                    sqr(tempZeros[i*(imageHight+2)+j]-tempZeros[(i-1)*(imageHight+2)+j])
                    +sqr(tempZeros[i*(imageHight+2)+j]-tempZeros[i*(imageHight+2)+(j-1)])
                );
            }
        mexPrintf("It=%d, nV_KepsNablap=%f\n",k,nV_KepsNablap);
        pSkalard = 0;
        for (int i = 0; i<imageHight; ++i)
            for (int j = 0; j<imageWidth; ++j)
                pSkalard += p[i*imageHight+j]*d[i*imageHight+j];
        normAp = 0;
        for (int i = 0; i<gM; ++i)
            for (int j = 0; j<gN; ++j)
                normAp += sqr(q[i*gM+j]);
        c_1 = pSkalard/(normAp+alpha*nV_KepsNablap);
        
        memcpy(tempF, myF, imageHight*imageWidth*sizeof(double));
        double maxCP = 0;
        for (int i = 0; i<imageHight; ++i)
            for (int j = 0; j<imageWidth; ++j) {
                myF[i*imageHight+j] += c_1*p[i*imageHight+j];
                double CP = fabs(c_1*p[i*imageHight+j]);
                if (maxCP<CP) maxCP = CP;
            }
        
        if (k%(int(IterImage/2+0.0001))==0) {
            domedfilt2(f, f);
        }
        
        if (k==IterImage) {
            mexPrintf("Min Image: maxIt reached\n");
            break;
        }
        if (maxCP<IM_TOL) {
            mexPrintf("Min Image: IM_Tol reached, %i Iterations \n",k);
            break;
        }
        
        opA(myF,Af);
        AT_f_new = 0;
        for (int i = 0; i<gM; ++i)
            for (int j = 0; j<gN; ++j)
                AT_f_new += sqr(Af[i*gM+j]-g[i*gM+j]);
                
        for (int i = 1; i<=imageHight; ++i)
            for (int j = 1; j<=imageWidth; ++j) {
                tempZeros[i*(imageHight+2)+j] = myF[(i-1)*imageHight+(j-1)];
            }
        for (int i = 2; i<=imageHight-1; ++i)
            for (int j = 2; j<=imageWidth-1; ++j)
                AT_f_new += alpha*(sqr(padded_v[i*(imageHight+2)+j])+0.01*sqr(epsilon))*(
                    sqr(tempZeros[i*(imageHight+2)+j]-tempZeros[(i-1)*(imageHight+2)+j])
                    +sqr(tempZeros[i*(imageHight+2)+j]-tempZeros[i*(imageHight+2)+(j-1)])
                );
        if (1.1*AT_f_old<AT_f_new && 1<k) {
            mexPrintf("Min Image: Instable, %i Iterations \n", k);
            memcpy(myF, tempF, imageHight*imageWidth*sizeof(double));
            break;
        }
        AT_f_old = AT_f_new;
        for (int i = 0; i<gM; ++i)
            for (int j = 0; j<gN; ++j)
                r[i*gM+j] = g[i*gM+j]-c_1*q[i*gM+j];
                
        for (int i = 1; i<=imageHight; ++i)
            for (int j = 1; j<=imageWidth; ++j) {
                tempZeros[i*(imageHight+2)+j] = myF[(i-1)*imageHight+(j-1)];
            }
        opAS(r,d);
        for (int i = 1; i<=imageHight; ++i)
            for (int j = 1; j<=imageWidth; ++j) {
                d[(i-1)*imageHight+(j-1)] += alpha*(
                    +sqr(padded_v[i*(imageHight+2)+j])*(tempZeros[(i+1)*(imageHight+2)+j]-tempZeros[i*(imageHight+2)+j])
                    +sqr(padded_v[i*(imageHight+2)+j])*(tempZeros[i*(imageHight+2)+(j+1)]-tempZeros[i*(imageHight+2)+j])
                    -sqr(padded_v[(i-1)*(imageHight+2)+j])*(tempZeros[i*(imageHight+2)+j]-tempZeros[(i-1)*(imageHight+2)+j])
                    -sqr(padded_v[i*(imageHight+2)+(j-1)])*(tempZeros[i*(imageHight+2)+j]-tempZeros[i*(imageHight+2)+(j-1)])
                    +0.01*sqr(epsilon)*(
                        tempZeros[(i+1)*(imageHight+2)+j]
                        +tempZeros[i*(imageHight+2)+(j+1)]
                        -4*tempZeros[i*(imageHight+2)+j]
                        +tempZeros[(i-1)*(imageHight+2)+j]
                        +tempZeros[i*(imageHight+2)+(j-1)]
                    )
                );
            }
        opA(d, Ad);
        AdSkalarAp = 0;
        for (int i = 0; i<gM; ++i)
            for (int j = 0; j<gN; ++j)
                AdSkalarAp += Ad[i*gM+j]*q[i*gM+j];
        c_2 = -AdSkalarAp/normAp;
        
        for (int i = 0; i<imageHight; ++i)
            for (int j = 0; j<imageWidth; ++j) {
                p[i*imageHight+j] = d[i*imageHight+j]+c_2*p[i*imageHight+j];
            }
    }
    
    plhs[0] = mxCreateDoubleMatrix(imageHight,imageWidth,mxREAL);
    memcpy(mxGetPr(plhs[0]), myF, imageHight*imageWidth*sizeof(double));
}

//***************************************************************************
//***************************************************************************
//***************************************************************************

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    loadData(nrhs, prhs);
    minimize_in_IMAGE_variable(plhs);
    return;
}

