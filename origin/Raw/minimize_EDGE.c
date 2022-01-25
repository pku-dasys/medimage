#include "minimize_EDGE.h"
#include "c_imp.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

void minimize_in_EDGE_variable(double out[imageHight][imageWidth], double f[imageHight][imageWidth],
                                 double v[imageHight][imageWidth],double alpha,double beta,double epsilon,
                                 double ED_TOL,int IterEdge) {
    double nNablaD, nD, c;
	int i,j,k;
	
    for (i = 0; i<imageHight; ++i)
        for (j = 0; j<imageWidth; ++j)
            update[i][j] = v[i][j];

    memset(g,0,sizeof(g));
    memset(d,0,sizeof(d));
    memset(laplace_v,0,sizeof(laplace_v));
    memset(d_xf,0,sizeof(d_xf));
    memset(d_yf,0,sizeof(d_yf));
    memset(n,0,sizeof(n));
    memset(tempZeros,0,sizeof(tempZeros));
    for (i = 0; i<imageHight+2; ++i)
        for (j = 0; j<imageWidth+2; ++j)
            tempOnes[i][j] = 1.0;
    memset(padded,0,sizeof(padded));

    // padded(2:imageHight+1,2:imageWidth+1)  = f;
    for (i = 1; i<=imageHight; ++i)
        for (j = 1; j<=imageWidth; ++j) {
            padded[i][j] = f[i-1][j-1];
        }

    for (k = 1; k<=IterEdge; ++k) {
        for (i = 0; i<imageHight; ++i)
            for (j = 0; j<imageWidth; ++j)
                temp[i][j] = update[i][j];
        for (i = 1; i<=imageHight; ++i)
            for (j = 1; j<=imageWidth; ++j)
                tempOnes[i][j] = temp[i-1][j-1];
        for (i = 1; i<=imageHight; ++i)
            for (j = 1; j<=imageWidth; ++j) {
                d_xf[i-1][j-1] = 0.5*(padded[i+1][j]-padded[i-1][j]);
                d_yf[i-1][j-1] = 0.5*(padded[i][j+1]-padded[i][j-1]);
                n[i-1][j-1] = sqr(d_xf[i-1][j-1])+sqr(d_yf[i-1][j-1]);
                laplace_v[i-1][j-1] = 0.25*(tempOnes[i+1][j]+tempOnes[i][j+1]
                                            -4*tempOnes[i][j]
                                            +tempOnes[i-1][j]+tempOnes[i][j-1]);
                g[i-1][j-1] = alpha*tempOnes[i][j]*n[i-1][j-1]
                              -beta*epsilon*laplace_v[i-1][j-1]
                              +beta*(tempOnes[i][j]-1)/(4*epsilon);
            }

        for (i = 0; i<imageHight; ++i)
            for (j = 0; j<imageWidth; ++j)
                d[i][j] = -g[i][j];
        
        for (i = 1; i<=imageHight; ++i)
            for (j = 1; j<=imageWidth; ++j)
                tempZeros[i][j] = d[i-1][j-1];
        nNablaD = 0;
        for (i = 1; i<=imageHight; ++i)
            for (j = 1; j<=imageWidth; ++j)
                nNablaD = nNablaD+sqr(0.5*(tempZeros[i+1][j]-tempZeros[i-1][j]))
                                 +sqr(0.5*(tempZeros[i][j+1]-tempZeros[i][j-1]));

        nD = 0;
        for (i = 0; i<imageHight; ++i)
            for (j = 0; j<imageWidth; ++j)
                nD += sqr(d[i][j]);
        // c = nD/(alpha*sum(sum(n.*d.^2)) + beta*epsilon*nNablaD + (beta/(4*epsilon))*nD);
        double tmp = 0;
        for (i = 0; i<imageHight; ++i)
            for (j = 0; j<imageWidth; ++j)
                tmp += n[i][j]*sqr(d[i][j]);
        c = nD/(alpha*tmp+beta*epsilon*nNablaD+(beta/(4*epsilon))*nD);
        printf("minimize_EDGE It=%d nD=%f c=%f\n",k,(double)nD,(double)c);

        double maxCD = 0;
        for (i = 0; i<imageHight; ++i)
            for (j = 0; j<imageWidth; ++j) {
                update[i][j] = temp[i][j]+c*d[i][j];
                double CD = fabs(c*d[i][j]);
                if (CD>maxCD) maxCD = CD;
            }
        
        for (i = 0; i<imageHight; ++i)
            for (j = 0; j<imageWidth; ++j) {
                if (update[i][j]>1) update[i][j] = 1;
                if (update[i][j]<0) update[i][j] = 0;
            }

        if (k==IterEdge) {
            printf("Min Edge: maxIt reached\n");
            break;
        }
        //mexPrintf("maxCD=%f\n",maxCD);
        if (maxCD<ED_TOL) {
            printf("Min Edge: ED_Tol reached, %d Iterations \n", k);
            break;
        }
    }

    for (i = 0; i<imageHight; ++i)
        for (j = 0; j<imageWidth; ++j)
            out[i][j] = update[i][j];
}

