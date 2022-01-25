#include "minimize_IMAGE.h"
#include "c_imp.h"
#include "A.h"
#include "AStar.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

void minimize_in_IMAGE_variable(double out[imageHight][imageWidth], double g[nAlpha][imageHight],
                                  double f[imageHight][imageWidth],double v[imageHight][imageWidth],
                                  double alpha,double epsilon,double IM_TOL,int IterImage) {    
    double AT_f_old, nV_KepsNablap, c_1, AT_f_new, pSkalard, normAp, AdSkalarAp, c_2;
	int i,j,k;

    for (i = 0; i<imageHight; ++i)
        for (j = 0; j<imageWidth; ++j)
            update[i][j] = f[i][j];

    memset(tempZeros, 0, sizeof(tempZeros));
    
    for (i = 0; i<imageHight+2; ++i)
        for (j = 0; j<imageWidth+2; ++j)
            padded[i][j] = 1.0;
    for (i = 1; i<=imageHight; ++i)
        for (j = 1; j<=imageWidth; ++j) {
            padded[i][j] = v[i-1][j-1];
        }

    A(Af, update);
    for (i = 0; i<nAlpha; ++i)
        for (j = 0; j<imageHight; ++j)
            r[i][j] = g[i][j]-Af[i][j];

    for (i = 1; i<=imageHight; ++i)
        for (j = 1; j<imageWidth; ++j)
            tempZeros[i][j] = update[i-1][j-1];
    
    AStar(d, r);
    
    for (i = 1; i<=imageHight; ++i)
        for (j = 1; j<=imageWidth; ++j) {
            d[i-1][j-1] += alpha*(
                +sqr(padded[i][j])*(tempZeros[i+1][j]-tempZeros[i][j])
                +sqr(padded[i][j])*(tempZeros[i][j+1]-tempZeros[i][j])
                -sqr(padded[i-1][j])*(tempZeros[i][j]-tempZeros[i-1][j])
                -sqr(padded[i][j-1])*(tempZeros[i][j]-tempZeros[i][j-1])
                +0.01*sqr(epsilon)*(
                    tempZeros[i+1][j]+tempZeros[i][j+1]
                    -4*tempZeros[i][j]
                    +tempZeros[i-1][j]+tempZeros[i][j-1]
                )
            );
        }
    AT_f_old = 0;
    for (i = 0; i<nAlpha; ++i)
        for (j = 0; j<imageHight; ++j)
            AT_f_old += sqr(Af[i][j]-g[i][j]);
    
    for (i = 1; i<=imageHight; ++i)
        for (j = 1; j<=imageWidth; ++j) {
            tempZeros[i][j] = update[i-1][j-1];
        }
    for (i = 2; i<=imageHight-1; ++i)
        for (j = 2; j<=imageWidth-1; ++j)
            AT_f_old += alpha*(sqr(padded[i][j])+0.01*sqr(epsilon))*(
                sqr(tempZeros[i][j]-tempZeros[i-1][j])+sqr(tempZeros[i][j]-tempZeros[i][j-1])
            );
    
    for (i = 0; i<imageHight; ++i)
        for (j = 0; j<imageWidth; ++j)
            p[i][j] = d[i][j];

    for (k = 1; k<=IterImage; ++k) {
        A(q, p);
        
        for (i = 1; i<=imageHight; ++i)
            for (j = 1; j<=imageWidth; ++j) {
                tempZeros[i][j] = p[i-1][j-1];
            }
        nV_KepsNablap = 0;
        for (i = 1; i<=imageHight; ++i)
            for (j = 1; j<=imageWidth; ++j) {
                nV_KepsNablap += (sqr(padded[i][j])+0.01*sqr(epsilon))*(
                    sqr(tempZeros[i][j]-tempZeros[i-1][j])+sqr(tempZeros[i][j]-tempZeros[i][j-1])
                );
            }
        pSkalard = 0;
        for (i = 0; i<imageHight; ++i)
            for (j = 0; j<imageWidth; ++j)
                pSkalard += p[i][j]*d[i][j];
        normAp = 0;
        for (i = 0; i<nAlpha; ++i)
            for (j = 0; j<imageHight; ++j)
                normAp += sqr(q[i][j]);
        c_1 = pSkalard/(normAp+alpha*nV_KepsNablap);
        
        for (i = 0; i<imageHight; ++i)
            for (j = 0; j<imageWidth; ++j)
                temp[i][j] = update[i][j];
        double maxCP = 0;
        for (i = 0; i<imageHight; ++i)
            for (j = 0; j<imageWidth; ++j) {
                update[i][j] += c_1*p[i][j];
                double CP = fabs(c_1*p[i][j]);
                if (maxCP<CP) maxCP = CP;
            }
        
        if (k%((int)(IterImage/2+0.0001))==0) {
            medfilt2(update);
        }
        
        if (k==IterImage) {
            printf("Min Image: maxIt reached\n");
            break;
        }
        if (maxCP<IM_TOL) {
            printf("Min Image: IM_Tol reached, %i Iterations \n",k);
            break;
        }
        
        A(Af, update);
        AT_f_new = 0;
        for (i = 0; i<nAlpha; ++i)
            for (j = 0; j<imageHight; ++j)
                AT_f_new += sqr(Af[i][j]-g[i][j]);
                
        for (i = 1; i<=imageHight; ++i)
            for (j = 1; j<=imageWidth; ++j) {
                tempZeros[i][j] = update[i-1][j-1];
            }
        for (i = 2; i<=imageHight-1; ++i)
            for (j = 2; j<=imageWidth-1; ++j)
                AT_f_new += alpha*(sqr(padded[i][j])+0.01*sqr(epsilon))*(
                    sqr(tempZeros[i][j]-tempZeros[i-1][j])
                    +sqr(tempZeros[i][j]-tempZeros[i][j-1])
                );
        if (1.1*AT_f_old<AT_f_new && 1<k) {
            printf("Min Image: Instable, %i Iterations \n", k);
            for (i = 0; i<imageHight; ++i)
                for (j = 0; j<imageWidth; ++j)
                    update[i][j] = temp[i][j];
            break;
        }
        AT_f_old = AT_f_new;
        for (i = 0; i<nAlpha; ++i)
            for (j = 0; j<imageWidth; ++j)
                r[i][j] -= c_1*q[i][j];
                
        for (i = 1; i<=imageHight; ++i)
            for (j = 1; j<=imageWidth; ++j) {
                tempZeros[i][j] = update[i-1][j-1];
            }
        AStar(d, r);
        for (i = 1; i<=imageHight; ++i)
            for (j = 1; j<=imageWidth; ++j) {
                d[i-1][j-1] += alpha*(
                    +sqr(padded[i][j])*(tempZeros[i+1][j]-tempZeros[i][j])
                    +sqr(padded[i][j])*(tempZeros[i][j+1]-tempZeros[i][j])
                    -sqr(padded[i-1][j])*(tempZeros[i][j]-tempZeros[i-1][j])
                    -sqr(padded[i][j-1])*(tempZeros[i][j]-tempZeros[i][j-1])
                    +0.01*sqr(epsilon)*(
                        tempZeros[i+1][j]+tempZeros[i][j+1]
                        -4*tempZeros[i][j]
                        +tempZeros[i-1][j]+tempZeros[i][j-1]
                    )
                );
            }
        A(Ad, d);
        AdSkalarAp = 0;
        for (i = 0; i<nAlpha; ++i)
            for (j = 0; j<imageHight; ++j)
                AdSkalarAp += Ad[i][j]*q[i][j];
        c_2 = -AdSkalarAp/normAp;
        
        for (i = 0; i<imageHight; ++i)
            for (j = 0; j<imageWidth; ++j) {
                p[i][j] = d[i][j]+c_2*p[i][j];
            }
        printf("It=%d, AT_f_old=%f, nV_KepsNablap=%f, pSkalard=%f, normAp=%f, c_1=%f, c_2=%f\n",k,AT_f_old,nV_KepsNablap,pSkalard,normAp,c_1,c_2);
    }
    
    for (i = 0; i<imageHight; ++i)
        for (j = 0; j<imageWidth; ++j)
            out[i][j] = update[i][j];
}
