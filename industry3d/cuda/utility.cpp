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

#include "ct3d.h"
#include "utility.h"

using namespace std;

bool equals_draw(float x,int y) {
    return fabs(x-y)<1.25;
}

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

float sqr(float x) { return x*x; }


// u is the axis of the rotation
// R is the rotating matrix
// theta
void rotate_axis_3d(float &x,float &y,float &z,float ux,float uy,float uz,float theta) {
    float c = 1/sqrt(sqr(ux)+sqr(uy)+sqr(uz));
    ux *= c;
    uy *= c;
    uz *= c;

    float X[3][1] = {{x},{y},{z}};
    float costheta = cos(theta), sintheta = sin(theta);
    float R[3][3] = {
        {costheta+sqr(ux)*(1-costheta), ux*uy*(1-costheta)-uz*sintheta, ux*uz*(1-costheta)+uy*sintheta},
        {uy*ux*(1-costheta)+uz*sintheta, costheta+sqr(uy)*(1-costheta), uy*uz*(1-costheta)-ux*sintheta},
        {uz*ux*(1-costheta)-uy*sintheta, uz*uy*(1-costheta)+ux*sintheta, costheta+sqr(uz)*(1-costheta)}
    };

    float d[3][1] = {};
    for (int k = 0; k<3; ++k)
        for (int i = 0; i<3; ++i)
            for (int j = 0; j<1; ++j)
                d[i][j] += R[i][k]*X[k][j];
    
    x = d[0][0];
    y = d[1][0];
    z = d[2][0];
}

void rotate_2d(float &x,float &y,float theta) {
    float X[2][1] = {{x},{y}};
    float costheta = cos(theta), sintheta = sin(theta);
    float R[2][2] = {
        {costheta, -sintheta},
        {sintheta, costheta}
    };

    float d[2][1] = {};
    for (int k = 0; k<2; ++k)
        for (int i = 0; i<2; ++i)
            for (int j = 0; j<1; ++j)
                d[i][j] += R[i][k]*X[k][j];
    
    x = d[0][0];
    y = d[1][0];
}
