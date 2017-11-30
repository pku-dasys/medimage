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
