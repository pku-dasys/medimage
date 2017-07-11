#include <stdio.h>
#include <string.h>
#include <math.h>

#include "config.h"
#include "proj.h"

double init_f[imageHight][imageWidth];
double f[imageHight][imageWidth],v[imageHight][imageWidth];
double gDelta[nAlpha][projWidth];

/* variable used by minimize_*.c */
double update[imageHight][imageWidth], d[imageHight][imageWidth];
double g[imageHight][imageWidth], laplace_v[imageHight][imageWidth];
double d_xf[imageHight][imageWidth], d_yf[imageHight][imageWidth];
double n[imageHight][imageWidth];
double tempZeros[imageHight+2][imageWidth+2], tempOnes[imageHight+2][imageWidth+2];
double padded[imageHight+2][imageWidth+2];

double p[imageHight][imageWidth];
double temp[imageHight][imageWidth];
double Af[nAlpha][projWidth], r[nAlpha][projWidth], q[nAlpha][projWidth], Ad[nAlpha][projWidth];
double padded_f[imageHight+2][imageWidth+2];

void minimize_in_IMAGE_variable(double out[imageHight][imageWidth], double g[nAlpha][projWidth],
                                  double f[imageHight][imageWidth],double v[imageHight][imageWidth],
                                  double alpha,double epsilon,int IterImage);
void minimize_in_EDGE_variable(double out[imageHight][imageWidth], double f[imageHight][imageWidth],
                                 double v[imageHight][imageWidth],double alpha,double beta,double epsilon,
                                 int IterEdge);
void medfilt2(double f[imageHight][imageWidth]);

double sqr(double x) { return x*x; }

void MumfordShah(double alpha,double beta,double epsilon) {
    memset(f,0,sizeof(f));
    int i,j;
    for (i = 0; i<imageHight; ++i)
        for (j = 0; j<imageWidth; ++j)
            v[i][j] = 1.0;
    
    printf("Begin alternative minimization...\n");
    for (i = 1; i<=_IterAlt; ++i) {
        printf("Iteration %d...\n",i);
        minimize_in_IMAGE_variable(f, gDelta,f,v,alpha,epsilon,_IterImage);
        minimize_in_EDGE_variable(v, f,v,alpha,beta,epsilon,_IterEdge);
    }
    printf("Leaving MumfordShah_reg...\n");
}

int main() {
    FILE *fin = fopen(PHANTOM_FILE, "r"), *fout;
    int i,j;

    for (i = 0; i<imageHight; ++i)
        for (j = 0; j<imageWidth; ++j)
            fscanf(fin,"%lf",&init_f[i][j]);
    fclose(fin);

    A(gDelta, init_f);
    MumfordShah(ALPHA,BETA,EPSILON);

    fout = fopen("I_result.txt","w");
    for (i = 0; i<imageHight; ++i) {
        for (j = 0; j<imageWidth; ++j)
            fprintf(fout,"%f ",f[i][j]);
        fprintf(fout,"\n");
    }
    fclose(fout);
    
    fout = fopen("E_result.txt","w");
    for (i = 0; i<imageHight; ++i) {
        for (j = 0; j<imageWidth; ++j)
            fprintf(fout,"%f ",v[i][j]);
        fprintf(fout,"\n");
    }
    fclose(fout);
    return 0;
}

void minimize_in_IMAGE_variable(double out[imageHight][imageWidth], double g[nAlpha][projWidth],
                                  double f[imageHight][imageWidth],double v[imageHight][imageWidth],
                                  double alpha,double epsilon,int IterImage) {    
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
        for (j = 0; j<projWidth; ++j)
            r[i][j] = g[i][j]-Af[i][j];

    for (i = 1; i<=imageHight; ++i)
        for (j = 1; j<=imageWidth; ++j)
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
        for (j = 0; j<projWidth; ++j)
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
            for (j = 0; j<projWidth; ++j)
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
        if (maxCP<1e-5) {
            printf("Min Image: IM_Tol reached, %i Iterations \n",k);
            break;
        }
        
        A(Af, update);
        AT_f_new = 0;
        for (i = 0; i<nAlpha; ++i)
            for (j = 0; j<projWidth; ++j)
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
            for (j = 0; j<projWidth; ++j)
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
            for (j = 0; j<projWidth; ++j)
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

void minimize_in_EDGE_variable(double out[imageHight][imageWidth], double f[imageHight][imageWidth],
                                 double v[imageHight][imageWidth],double alpha,double beta,double epsilon,
                                 int IterEdge) {
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
        if (maxCD<1e-5) {
            printf("Min Edge: ED_Tol reached, %d Iterations \n", k);
            break;
        }
    }

    for (i = 0; i<imageHight; ++i)
        for (j = 0; j<imageWidth; ++j)
            out[i][j] = update[i][j];
}

void swap(double *a,double *b) {
    double t = *a; *a = *b; *b = t;
}

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