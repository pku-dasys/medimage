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

float array_3d_raw(float *data,int SIZE0,int SIZE1,int SIZE2,int x,int y,int z) {
    return *(data + z + SIZE2 * (y + x * SIZE1));
}

float array_3d_img(float *data,int z,int x,int y) {
    return array_3d_raw(data,NZ,NX,NY,z,x,y);
}

ushort array_3d_raw_sino(ushort *data,int SIZE0,int SIZE1,int SIZE2,int x,int y,int z) {
    return *(data + z + SIZE2 * (y + x * SIZE1));
}

ushort array_3d_sino(ushort *data,int z,int x,int y) {
    return array_3d_raw_sino(data,NPROJ,NDX,NDY,z,x,y);
}

/*
int64_t timer_us( void )
{
#ifdef WIN32
    struct _timeb tb;
    _ftime(&tb);
    return ((int64_t)tb.time * (1000) + (int64_t)tb.millitm) * (1000);
#else
    struct timeval tv_date;
    gettimeofday( &tv_date, NULL );
    return ( (int64_t) tv_date.tv_sec * 1000000 + (int64_t) tv_date.tv_usec );
#endif
}

int64_t timer_s( void )
{
#ifdef WIN32
    struct _timeb tb;
    _ftime(&tb);
    return ((int64_t)tb.time);
#else
    struct timeval tv_date;
    gettimeofday( &tv_date, NULL );
    return ( (int64_t) tv_date.tv_sec);
#endif
}
*/

void write_data_3d(float *data,int Z,int X,int Y,string output_filename) {
    for (int z = 0; z<Z; ++z) {
        string slice_name = OUTPUT_DIR+"/"+output_filename+lexical_cast<string>(z);
        FILE *output = fopen(slice_name.c_str(),"w");
        for (int x = 0; x<X; ++x) {
            for (int y = 0; y<Y; ++y) {
                //out<<array_3d_raw(data,Z,X,Y,z,x,y)<<' ';
                fprintf(output,"%lf ",array_3d_raw(data,Z,X,Y,z,x,y));
            }
            fprintf(output,"\n");
            //out_slice<<endl;
        }
        fprintf(output,"\n");
        fclose(output);
        //out_slice.close();
    }
}

float sqr(float x) { return x*x; }

//float get_img_addr(float x,float y,float z) {
//    return z*NX*NY+x*NY+y;
//}
int64 get_img_addr(int64 x,int64 y,int64 z) {
    return z*NX*NY+x*NY+y;
}
