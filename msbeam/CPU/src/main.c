#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "proj.h"
#include "proto.h"
#include "utility.h"

float g[NPROJ*NRAY];
float filter_g[NPROJ*NRAY];

float f[IMGSIZE*IMGSIZE];
float v[IMGSIZE*IMGSIZE];

float lambda;

void MSbeam();

void write_phantom(float *img, const char *file_name) {
    FILE *fin = fopen(file_name,"w");
    int i,j;
    for (i = 0; i<NPROJ; ++i) {
        for (j = 0; j<NRAY; ++j)
            fprintf(fin,"%f ",img[i*NRAY+j]);
        fprintf(fin,"\n");
    }
    fclose(fin);
}

void filter(float *a) {
    int i,j;
    for (i = 0; i<IMGSIZE; ++i)
        for (j = 0; j<IMGSIZE; ++j)
            if (64<=i && i<256 && 64<=j && j<256);
            else if (a[i*IMGSIZE+j]<1);
            else a[i*IMGSIZE+j] = log(a[i*IMGSIZE+j])/log(2)*32.;
}

int main() {
    read_phantom(f, "phantom.dat");
    
    normalize(f);
    write_file(f, "std.dat");
    
    A(g, f);
    write_phantom(g, "g.dat");
    
    MSbeam();
    
    normalize(f);
    
    write_file(f, "img.dat");
    
    write_file(v, "edge.dat");
	
    /* */
    filter(f);
    normalize(v);
    filter(v);
    write_file(f, "img_filter.dat");
    write_file(v, "edge_filter.dat");
    /* */
    
    
    read_phantom(f, "phantom.dat");
    
    normalize(f);
    
	int i,j;
    for (i = 0; i<IMGSIZE; ++i)
        for (j = 0; j<IMGSIZE; ++j)
            if (64<=i && i<256 && 64<=j && j<256);
            else f[i*IMGSIZE+j] = 0;
                
    A(filter_g, f);
    
    for (i = 0; i<NPROJ*NRAY; ++i)
        if (filter_g[i]==0) {
            if (g[i]<1);
            else filter_g[i] = log(g[i])/log(2)*32.;
        }
    write_phantom(filter_g, "g_filter.dat");
    
    
    return 0;
}

void minIMAGE(int np,int nr,int *line,float *weight,int numb,float snorm) {
    int i;

    float d[IMGSIZE*2];
    
    float Af = 0.;
    for (i = 0; i<numb; ++i) {
        int ind = line[i];
        Af += f[ind]*weight[i];
    }
    Af -= g[np*NRAY+nr];
    /*if (nr==383) printf("%f\n",Af);*/
    
    /* calculate div (v^2 \grad f) and Af-g*/
    for (i = 0; i<numb; ++i) {
        int ind = line[i];
        int x = ind/IMGSIZE, y = ind%IMGSIZE;

        float tmp = 0.;
        float lap = 0.;
        
        if (x+1<IMGSIZE) tmp += sqr(v[ind])*(f[ind+IMGSIZE]-f[ind]);
        else             tmp += sqr(v[ind])*(   0          -f[ind]);
        
        if (y+1<IMGSIZE) tmp += sqr(v[ind])*(f[ind+1]-f[ind]);
        else             tmp += sqr(v[ind])*(   0    -f[ind]);
        
        if (x-1>=0)      tmp -= sqr(v[ind-IMGSIZE])*(f[ind]-f[ind-IMGSIZE]);
        else             tmp -=                     (f[ind]-0        );
        
        if (y-1>=0)      tmp -= sqr(v[ind-1])*(f[ind]-f[ind-1]);
        else             tmp -=               (f[ind]-0       );
        
        if (x+1<IMGSIZE) lap += f[ind+IMGSIZE];
        if (y+1<IMGSIZE) lap += f[ind+1];
        if (x-1>=0)      lap += f[ind-IMGSIZE];
        if (y-1>=0)      lap += f[ind-1];
        lap -= 4*f[ind];

        d[i] = -Af*weight[i]+ALPHA*(tmp+sqr(EPSILON)*lap);
    }
        
    for (i = 0; i<numb; ++i) {
        int ind = line[i];
        f[ind] += lambda*d[i];
        if (f[ind]<0) f[ind] = 0;
        if (f[ind]>255) f[ind] = 255;
    }
}

void minEDGE(int np,int nr,int *line,float *weight,int numb) {
    int i;

    float d[IMGSIZE*2];
    
    for (i = 0; i<numb; ++i) {
        int ind = line[i];
        int x = ind/IMGSIZE, y = ind%IMGSIZE;
        
        float a = 0.;
        float b = 0.;
        float c = 0.;

        if (x-1>=0)      a += sqr(f[ind]-f[ind-IMGSIZE]);
        else             a += sqr(f[ind]-0        );
        
        if (y-1>=0)      a += sqr(f[ind]-f[ind-1]);
        else             a += sqr(f[ind]-0       );
        
        a *= v[ind];
        
        b = v[ind]-1;
                
        if (x+1<IMGSIZE) c += v[ind+IMGSIZE];
        if (y+1<IMGSIZE) c += v[ind+1];
        if (x-1>=0)      c += v[ind-IMGSIZE];
        if (y-1>=0)      c += v[ind-1];
        c -= 4*v[ind];
        
        d[i] = -ALPHA*a-BETA/(4*EPSILON)*b+BETA*EPSILON*c;
    }
    
    for (i = 0; i<numb; ++i) {
        int ind = line[i];
        v[ind] += lambda*d[i];
        if (v[ind]<0) v[ind] = 0;
        if (v[ind]>1) v[ind] = 1;
    }
}

float AT_g[NPROJ*NRAY];

float AT,G,H;

void compute_AT(int flag) {
    A(AT_g, f);
    
    int i;
    float part_1 = 0., part_2 = 0., part_3 = 0.;
    
    for (i = 0; i<NPROJ*NRAY; ++i)
        part_1 += sqr(AT_g[i]-g[i]);
        
    for (i = 0; i<IMGSIZE*IMGSIZE; ++i) {
        float tmp = 0;
        float x = i/IMGSIZE, y = i%IMGSIZE;
        
        if (x-1 >= 0) tmp += sqr(f[i] - f[i-IMGSIZE]);
        else          tmp += sqr(f[i] - 0);
        if (y-1 >= 0) tmp += sqr(f[i] - f[i-1]);
        else          tmp += sqr(f[i] - 0);
        
        part_2 += tmp*sqr(v[i]);
    }
    part_2 *= ALPHA;
    
    for (i = 0; i<IMGSIZE*IMGSIZE; ++i) {
        float tmp = 0.;
        float x = i/IMGSIZE, y = i%IMGSIZE;
        
        if (x-1 >= 0) tmp += sqr(v[i] - v[i-IMGSIZE]);
        else          tmp += sqr(v[i] - 1);
        if (y-1 >= 0) tmp += sqr(v[i] - v[i-1]);
        else          tmp += sqr(v[i] - 1);
        
        part_3 += tmp*EPSILON + sqr(1-v[i])/(4*EPSILON);
    }
    part_3 *= BETA;
    
    
    if (flag==1) {
        AT = part_1+part_2+part_3;
        G = part_1+part_2;
        H = part_2+part_3;
    }
    printf("AT(f,v) = %f (%.2f%%),\tG(f) = %f (%.2f%%),\tH(v) = %f (%.2f%%)\n",part_1+part_2+part_3,(part_1+part_2+part_3)/AT*100.,
        part_1+part_2,(part_1+part_2)/G*100.,part_2+part_3,(part_2+part_3)/H*100.);
}

void min_wrapper(int np,int nr) {
    int line[IMGSIZE*2];
    float weight[IMGSIZE*2];
    int numb;
    float snorm;
    
    wray(np, nr, line, weight, &numb, &snorm);
    
    minIMAGE(np,nr,line,weight,numb,snorm);
    minEDGE(np,nr,line,weight,numb); 
}

void MSbeam() {
    int i,j;
    for (i = 0; i<IMGSIZE*IMGSIZE; ++i) {
        f[i] = 0.;
        v[i] = 1.;
    }
    
    lambda = 0.001;
        
    printf("Begin MSbeam minimization ...\n");
    for (i = 1; i<=ALL_ITER; ++i) {
        printf("Iteration %d ...\n", i);
        printf("lambda = %f\n",lambda);
        int np = 0,nr = 0, total_rays = NPROJ*NRAY;
        for (j = 0; j<total_rays; ++j) {
            min_wrapper(np,nr);
            pick(&np,&nr);
        }
        compute_AT(i);
        /*
        float minval = 1e30, maxval = -1e30;
        for (j = 0; j<IMGSIZE*IMGSIZE; ++j) {
            if (minval>f[j]) minval = f[j];
            if (maxval<f[j]) maxval = f[j];
        }
        printf("min f_value = %f, max f_value = %f\n",minval,maxval);
        */
        lambda = lambda/(1 + 500 * lambda);
    }
    printf("MSbeam minimization done.\n");
}
