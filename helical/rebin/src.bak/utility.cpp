#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <sys/types.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <sys/times.h>

#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include "main.h"
#include "utility.h"

using namespace boost;
using namespace std;

double array_2d_raw (double *data, int SIZE0, int SIZE1,int x, int y){
  return *(data + x*SIZE1 + y);
}

double array_2d_sino (double *data, int x, int y){
  return array_2d_raw(data, NPROJ, NDETECTOR, x, y);
}

double array_2d_img (double *data, int x, int y) {
  return array_2d_raw(data, NX, NY,x, y);
}

double sqr(double x) { return x*x; }

double get_img_addr(double x,double y) {
    return x*NY+y;
}

void write_data_2d(double *data,int X,int Y,string output_filename) {
    FILE *output = fopen(output_filename.c_str(),"w");
    cout<<output_filename<<endl;
    if (output == NULL) printf("fopen file failed!!!\n");
    for (int x = 0; x<X; ++x) {
        for (int y = 0; y<Y; ++y) {
            //out<<array_3d_raw(data,Z,X,Y,z,x,y)<<' ';
            fprintf(output,"%lf ",array_2d_raw(data,X,Y,x,y));
        }
        fprintf(output,"\n");
            //out_slice<<endl;
    }
    fprintf(output,"\n");
        //out_slice.close();
    fclose(output);
}

double array_3d_raw(double *data,int SIZE0,int SIZE1,int SIZE2,int x,int y,int z) {
    return *(data + z + SIZE2 * (y + x * SIZE1));
}

void write_data_3d(double *data,int Z,int X,int Y,string output_filename) {
    FILE *output = fopen(output_filename.c_str(),"w");
    for (int z = 0; z<Z; ++z) {
        for (int x = 0; x<X; ++x) {
            for (int y = 0; y<Y; ++y) {
                //out<<array_3d_raw(data,Z,X,Y,z,x,y)<<' ';
                fprintf(output,"%lf ",array_3d_raw(data,Z,X,Y,z,x,y));
            }
            fprintf(output,"\n");
            //out_slice<<endl;
        }
        fprintf(output,"\n");
        //out_slice.close();
    }
    fclose(output);
}
