#include "header.h"
#include "A.h"
#include "AStar.h"
#include "c_imp.h"
#include "minimize_IMAGE.h"
#include "minimize_EDGE.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

void read_files() {
    // f = phantom(512);
    FILE *fin = fopen(PHANTOM_FILE, "r");
	int i,j;
    for (i = 0; i<imageHight; ++i)
        for (j = 0; j<imageWidth; ++j) {
            fscanf(fin,"%lf",&init_f[i][j]);
        }
    fclose(fin);
}

void MumfordShah_reg(double g[nAlpha][projWidth],double alpha,double beta,double epsilon) {
    memset(f,0,sizeof(f));
	int i,j;
    for (i = 0; i<imageHight; ++i)
        for (j = 0; j<imageWidth; ++j)
            v[i][j] = 1.0;
    
    printf("Begin alternative minimization...\n");
    for (i = 1; i<=_IterAlt; ++i) {
        printf("Iteration %d...\n",i);
        minimize_in_IMAGE_variable(f, g,f,v,alpha,epsilon,_IM_TOL,_IterImage);
        minimize_in_EDGE_variable(v, f,v,alpha,beta,epsilon,_ED_TOL,_IterEdge);
    }
    printf("Leaving MumfordShah_reg...\n");
}

void hls_point() {
    double g[nAlpha][projWidth];

    //g = A(init_f);
    A(g, init_f);

    int i,j;
    
    double maxG = 0;
    for (i = 0; i<nAlpha; ++i)
        for (j = 0; j<projWidth; ++j) {
            double tmp = fabs(g[i][j]);
            if (maxG<tmp) maxG = tmp;
        }
    
    //double c = 0.05;
    //double n[nAlpha][imageHight];
    //randn(n);
    for (i = 0; i<nAlpha; ++i)
        for (j = 0; j<projWidth; ++j)
            //gDelta[i][j] = g[i][j]+c*maxG*n[i][j];
            gDelta[i][j] = g[i][j];
    MumfordShah_reg(gDelta,ALPHA,BETA,0.00001);
}

void output_files() {
    FILE *fout;
    int i,j;
	
    fout = fopen("I_result.txt","w");
    for (i = 0; i<imageHight; ++i) {
        for (j = 0; j<imageWidth; ++j) {
            fprintf(fout,"%f ",f[i][j]);
        }
        fprintf(fout,"\n");
    }
    fclose(fout);
    
    fout = fopen("E_result.txt","w");
    for (i = 0; i<imageHight; ++i) {
        for (j = 0; j<imageWidth; ++j) {
            fprintf(fout,"%f ",v[i][j]);
        }
        fprintf(fout,"\n");
    }
    fclose(fout);
}

int main() {
    read_files();

    hls_point();
    
    output_files();
    return 0;
}
