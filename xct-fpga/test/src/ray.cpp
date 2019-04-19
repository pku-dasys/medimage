#include "config.h"
#ifdef __SDSCC__
#include "sds_lib.h"
#endif

using namespace std;

#define sqr(x) ((x)*(x))

void ray(Data lambda_f, Data lambda_v, Data ALPHA, Data BETA, Data EPSILON,
        Data minus_g, Number numb, Weight wgt[MAX_ELE_RAY],
        Data f0[MAX_ELE_RAY], Data f1[MAX_ELE_RAY], Data f2[MAX_ELE_RAY], Data f3[MAX_ELE_RAY], Data f4[MAX_ELE_RAY],
        Data v0[MAX_ELE_RAY], Data v1[MAX_ELE_RAY], Data v2[MAX_ELE_RAY], Data v3[MAX_ELE_RAY], Data v4[MAX_ELE_RAY],
        Data commit_f[MAX_ELE_RAY], Data commit_v[MAX_ELE_RAY]) {

    Data Afg = minus_g;

    #define LAT 16
    Data sum[MAX_ELE_RAY + LAT];

    for (int i = 0; i < LAT; ++i) {
#pragma HLS unroll
        sum[i] = 0;
    }

    accum_Afg:
    for (int i = 0; i<numb; ++i) {
#pragma HLS pipeline
        sum[i + LAT] = sum[i] + f0[i] * wgt[i];
    }

    for (int i = 0; i < LAT; ++i) {
#pragma HLS unroll
        Afg += sum[numb + i];
    }

    update_fv:
    for (int i = 0; i<numb; ++i) {
#pragma HLS pipeline
        Data tmp = 0.;
        Data lap = 0.;

        // 0  1  2  3  4
        // x  u  d  l  r

        tmp += sqr(v0[i])*(f2[i]-f0[i]);
        tmp += sqr(v0[i])*(f4[i]-f0[i]);
        tmp -= sqr(v1[i])*(f0[i]-f1[i]);
        tmp -= sqr(v3[i])*(f0[i]-f3[i]);

        lap = f1[i]+f2[i]+f3[i]+f4[i]-4*f0[i];

        Data d_f = -Afg*wgt[i]+ALPHA*(tmp+sqr(EPSILON)*lap);

        Data a = 0.;
        Data b = 0.;
        Data c = 0.;

        a = v0[i]*(sqr(f0[i]-f1[i])+sqr(f0[i]-f3[i]));

        b = v0[i]-1;

        c = v1[i]+v2[i]+v3[i]+v4[i]-4*v0[i];

        Data d_v = -ALPHA*a-BETA/(4*EPSILON)*b+BETA*EPSILON*c;

        Data new_f = f0[i]+lambda_f*d_f;
        Data new_v = v0[i]+lambda_v*d_v;

        if (new_f<0) new_f = 0;

        if (new_v<0) new_v = 0;
        if (new_v>1) new_v = 1;

        commit_f[i] = new_f;
        commit_v[i] = new_v;
    }
}
