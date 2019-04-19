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

double array_3d_raw(double *data,int SIZE0,int SIZE1,int SIZE2,int x,int y,int z) {
    return *(data + z + SIZE2 * (y + x * SIZE1));
}

double array_3d_img(double *data,int z,int x,int y) {
    return array_3d_raw(data,NZ,NX,NY,z,x,y);
}

double array_3d_sino(double *data,int z,int x,int y) {
    return array_3d_raw(data,NPROJ,NDETECTOR,NCHANNEL,z,x,y);
}


// int64_t timer_us( void )
// {
// #ifdef WIN32
//     struct _timeb tb;
//     _ftime(&tb);
//     return ((int64_t)tb.time * (1000) + (int64_t)tb.millitm) * (1000);
// #else
//     struct timeval tv_date;
//     gettimeofday( &tv_date, NULL );
//     return ( (int64_t) tv_date.tv_sec * 1000000 + (int64_t) tv_date.tv_usec );
// #endif
// }

// int64_t timer_s( void )
// {
// #ifdef WIN32
//     struct _timeb tb;
//     _ftime(&tb);
//     return ((int64_t)tb.time);
// #else
//     struct timeval tv_date;
//     gettimeofday( &tv_date, NULL );
//     return ( (int64_t) tv_date.tv_sec);
// #endif
// }


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

double sqr(double x) { return x*x; }

double get_img_addr(double x,double y,double z) {
    return z*NX*NY+x*NY+y;
}