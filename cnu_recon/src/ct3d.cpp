#include "utility.h"
#include "main.h"
#include "ct3d.h"
#include "tracing.h"
#include <cstdio>
#include <boost/lexical_cast.hpp>
#include <cstdlib>
#include <cstring>
#include <iostream>
using namespace std;
using namespace boost;

float *f;
ushort *g;
float lambda;

int max_numb;

void load_data(ushort* data, int w, int h, int p, const char* data_name) {
    FILE* input = fopen(data_name, "rb");

    if (!input) {
        fprintf(stderr, "cannot open data file\n");
        exit(1);
    }
    printf("loading ... ");
    fflush(stdout);
    fseek(input, 1024, SEEK_SET);//ignore header
    for (int i = 0; i<w*h*p; i += 2048) {
        if (!fread(data+i, sizeof(ushort), 2048, input)) {
            fprintf(stderr, "unexpected ending\n");
            exit(1);
        }
    }
    fclose(input);
    printf("done\n");
    fflush(stdout);
}
void load_part(ushort* data, int proj, int x, int y, int ys, int ye, const char* data_name){
    //proj*x*[ys,ye), ye is not included
    FILE* input = fopen(data_name, "rb");
    if(!input){
        fprintf(stderr, "cannot open data file\n");
        exit(1);
    }
    printf("loading...");
    fflush(stdout);
    fseek(input, 1024, SEEK_SET);//ignore header
    fseek(input, ys*sizeof(ushort), SEEK_CUR);
    int i, s = proj*x;
    for(i = 0; i < s; i++)
    {
        if(!fread(data, sizeof(ushort), ye-ys, input)) {
            fprintf(stderr, "unexpected ending\n");
            exit(1);
        }
        data += y;
        fseek(input, (y-(ye-ys))*sizeof(ushort), SEEK_CUR);
    }
    fclose(input);
    printf("done\n");
    fflush(stdout);
}



void minIMAGE(float Af, int64 *line, float *weight, int numb)
{
    for (int i = 0; i<numb; i++) {
        Af += f[line[i]] * weight[i];
    }
    for (int i = 0; i<numb; i++) {
        int64 ind = line[i];
        f[ind] += lambda * (-Af * weight[i]);
        if(f[ind]<0) f[ind] = 0;
    }
}

float *cos_table;
float *sin_table;

void compute(int alpha, int detectorX, int detectorY)
/*alpha is the rotate angle of light source
(or the object, then result will be a mirror image),
from the view of +Z to -Z, counterclockwise(+X->+Y->-X->-Y)
is positive. Initial angle is +Y
the parameters should be int, as it represents the relative
coordinate of grid in the detector

here, detectorY actually represent Z-axis
*/
{
    const float sina = sin_table[alpha], cosa = cos_table[alpha];
    float srcX,srcY,srcZ = 0.0;
    float dstX,dstY,dstZ;
    float oridstX,oridstY;

    srcX = (-SOD * sina) / vx;
    srcY = ( SOD * cosa) / vy;
    oridstX = detectorX - NDX*0.5 + 0.5;
    oridstY = (SOD - SDD) / vy;
    dstX = oridstX * cosa - oridstY * sina;
    dstY = oridstX * sina + oridstY * cosa;
    dstZ = detectorY - NDY*0.5 + 0.5;
    
    srcX += NX/2.0;
    srcY += NY/2.0;
    srcZ += NZ/2.0;
    dstX += NX/2.0;
    dstY += NY/2.0;
    dstZ += NZ/2.0;

    printf("%f %f %f     %f %f %f\n",srcX, srcY, srcZ, dstX, dstY, dstZ);

    int64 *ind = new int64[MAX_RAYLEN];
    float *wgt = new float[MAX_RAYLEN];
    int numb;

    float Af = -array_3d_sino(g, alpha, detectorX, detectorY);

    
    {
    forward_proj(srcX, srcY, srcZ,
                 dstX, dstY, dstZ,
                 ind, wgt, numb);
    }
    
    
    if (numb>0) {
        printf("%d %d %d\n",alpha,detectorX, detectorY);
        printf("%d\n",numb);
        for (int i = 0; i<numb; ++i) {
            int z,x,y;
            int64 tmp = ind[i];
            z = tmp/(NX*NY);
            tmp %= NX*NY;
            x = tmp/NY;
            y = tmp%NY;
            cout<<ind[i]<<endl;
            printf("%lld,%f  -  %d %d %d\n",tmp,wgt[i],z, x, y);
        }
        printf("\n");
        fflush(stdout);
        exit(0);
    }
    
    if (max_numb<numb) max_numb = max_numb;
    
    {
    minIMAGE(Af, ind, wgt, numb);
    }
    
    delete[] ind;
    delete[] wgt;
}

void wrapper() {
/*
    for(int k = 0; k < NPROJ; ++k) {
        for(int i = 0; i < NDX; i++) {
            for(int j = 0; j < NDY; j++) {
                compute(k, i, j);
            }
        }
    }
    */
   compute(0,512,512);
}

void ct3d(float *_image_data,/* float *_edge_data,*/ ushort *_sino_data) {
//deprecated: edge data
    f = _image_data;
//    v = _edge_data;
    g = _sino_data;

    //load_data(g, NPROJ, NDX, NDY, RAW_DATA_FILE.c_str());
    
    memset(f, 0, NZ*NX*NY*sizeof(float));
//    memset(v, 0, NZ*NX*NY*sizeof(float));

    cos_table = new float[NPROJ];
    sin_table = new float[NPROJ];

    for(int i = 0; i < NPROJ; i++) {
        const float alpha = (float)i/NPROJ*2*PI;
        sin_table[i] = sin(alpha);
        cos_table[i] = cos(alpha);
    }

    max_numb = 0;

    for(int iters = 0; iters < ITERATIONS; ++iters) {
        lambda = 1.0/(1000.0+iters*50.0);

        printf("iter = %d, lambda = %lf\n", iters, lambda);
        fflush(stdout);

        wrapper();

        //write_data_3d(f, NZ, NX, NY, OUTPUT_DIR+"/"+image_filename+lexical_cast<string>(iters));
    }
    puts("iter end");
    
    printf("max_numb = %d\n",max_numb);

    delete [] sin_table;
    delete [] cos_table;
}
