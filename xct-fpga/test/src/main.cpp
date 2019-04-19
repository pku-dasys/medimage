#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

#include "config.h"
#include "ray.h"

#ifdef __SDSCC__
#include <stdlib.h>
#include "sds_lib.h"
#define malloc(x) (sds_alloc(x))
#define free(x) (sds_free(x))
#endif

using namespace std;

class perf_counter {
public:
     uint64_t tot, cnt, calls;
     perf_counter() : tot(0), cnt(0), calls(0) {};
     inline void reset() { tot = cnt = calls = 0; }
     inline void start() {
#ifdef __SDSCC__
         cnt = sds_clock_counter();
#else
         cnt = clock();
#endif
         calls++;
     };
     inline void stop() {
#ifdef __SDSCC__
         tot += (sds_clock_counter() - cnt);
#else
         tot += (clock() - cnt);
#endif
     };
     inline uint64_t avg_cpu_cycles() { return ((tot+(calls>>1)) / calls); };
};

#define RAYS 11776

void load_tracing_file(Index *ind, Weight *wgt, Number *numb_set)
{
    ifstream fconfig("config.txt");
    for (int i = 0; i<RAYS; ++i) {
        int nproj, ndetector, nchannel;
        int numb;
        float table_pos, angle;
        fconfig>>nproj>>ndetector>>nchannel>>numb>>table_pos>>angle;
        numb_set[i] = numb;
    }
    fconfig.close();

    ifstream ind_in("ind.txt");
    for (int i = 0; i<RAYS; ++i) {
        for (int j = 0; j<numb_set[i]; ++j)
            ind_in >> ind[i*MAX_ELE_RAY+j];
    }
    ind_in.close();

    ifstream wgt_in("wgt.txt");
    for (int i = 0; i<RAYS; ++i) {
        for (int j = 0; j<numb_set[i]; ++j)
            wgt_in >> wgt[i*MAX_ELE_RAY+j];
    }
    wgt_in.close();
}

int main() {
    Index *ind;
    Number *numb_set;
    Weight *wgt;
    Data *f, *v;

    Index *b_line;
    Weight *b_wgt;
    Data *commit_f, *commit_v;

    Data *f0, *f1, *f2, *f3, *f4;
    Data *v0, *v1, *v2, *v3, *v4;

    ind = (Index *)malloc(RAYS * MAX_ELE_RAY * sizeof(Index));
    numb_set = (Number *)malloc(RAYS * sizeof(Number));
    wgt = (Weight *)malloc(RAYS * MAX_ELE_RAY * sizeof(Weight));
    f = (Data *)malloc(NZ * NX * NY * sizeof(Data));
    v = (Data *)malloc(NZ * NX * NY * sizeof(Data));

    b_line = (Index *)malloc(MAX_ELE_RAY * sizeof(Index));
    b_wgt = (Weight *)malloc(MAX_ELE_RAY * sizeof(Weight));

    commit_f = (Data *)malloc(MAX_ELE_RAY * sizeof(Data));
    commit_v = (Data *)malloc(MAX_ELE_RAY * sizeof(Data));

    f0 = (Data *)malloc(MAX_ELE_RAY * sizeof(Data));
    f1 = (Data *)malloc(MAX_ELE_RAY * sizeof(Data));
    f2 = (Data *)malloc(MAX_ELE_RAY * sizeof(Data));
    f3 = (Data *)malloc(MAX_ELE_RAY * sizeof(Data));
    f4 = (Data *)malloc(MAX_ELE_RAY * sizeof(Data));

    v0 = (Data *)malloc(MAX_ELE_RAY * sizeof(Data));
    v1 = (Data *)malloc(MAX_ELE_RAY * sizeof(Data));
    v2 = (Data *)malloc(MAX_ELE_RAY * sizeof(Data));
    v3 = (Data *)malloc(MAX_ELE_RAY * sizeof(Data));
    v4 = (Data *)malloc(MAX_ELE_RAY * sizeof(Data));

    for (int i = 0; i<NZ * NX * NY; ++i) {
        f[i] = 0.0;
        v[i] = 1.0;
    }

    load_tracing_file(ind, wgt, numb_set);

    perf_counter hw_ctr, commit_ctr, buf_ctr;

    for (int i = 0; i<RAYS; ++i) {
        Number numb = numb_set[i];

        Data lambda_f = 0.01, lambda_v = 0.001;

        Data minus_g = 100;

        buf_ctr.start();
        for(int k = 0; k < numb; k++) {
            Index temp_ind = ind[i * MAX_ELE_RAY + k];
            Weight temp_wgt = wgt[i * MAX_ELE_RAY + k];

            b_line[k] = temp_ind;
            b_wgt[k] = temp_wgt;

            Index indxy = temp_ind%(NX*NY);
            Index x = indxy/NY, y = indxy%NY;

            Index z = temp_ind/(NX*NY);
            Index xy = temp_ind%(NX*NY);
            Index new_ind = (z-STARTZ)*NX*NY+xy;

            f0[k] = f[new_ind];
            if (x>=1) f1[k] = f[new_ind - NY];
            else        f1[k] = 0;
            if (x+1<NX) f2[k] = f[new_ind + NY];
            else        f2[k] = 0;
            if (y>=1) f3[k] = f[new_ind - 1];
            else        f3[k] = 0;
            if (y+1<NY) f4[k] = f[new_ind + 1];
            else        f4[k] = 0;

            v0[k] = v[new_ind];
            if (x>=1) v1[k] = v[new_ind - NY];
            else        v1[k] = 1;
            if (x+1<NX) v2[k] = v[new_ind + NY];
            else        v2[k] = 1;
            if (y>=1) v3[k] = v[new_ind - 1];
            else        v3[k] = 1;
            if (y+1<NY) v4[k] = v[new_ind + 1];
            else        v4[k] = 1;
        }
        buf_ctr.stop();

        hw_ctr.start();
        ray(lambda_f, lambda_v, ALPHA, BETA, EPSILON,
            minus_g, numb, b_wgt,
            f0, f1, f2, f3, f4,
            v0, v1, v2, v3, v4,
            commit_f, commit_v);
        hw_ctr.stop();

        commit_ctr.start();
        for (int k = 0; k<numb; ++k) {
            Index ind = b_line[k];

            Index z = ind/(NX*NY);
            Index xy = ind%(NX*NY);

            Index new_ind = (z-STARTZ)*NX*NY+xy;

            f[new_ind] = commit_f[k];
            v[new_ind] = commit_v[k];
        }
        commit_ctr.stop();

        if (i%(RAYS/10)==0) {
            int t = i/(RAYS/10);
            cout << t*10 << "% percent finished" << endl;
        }
    }

    cout << "finished" << endl;

    uint64_t hw_cycles = hw_ctr.avg_cpu_cycles();
    uint64_t commit_cycles = commit_ctr.avg_cpu_cycles();
    uint64_t buf_cycles = buf_ctr.avg_cpu_cycles();
    std::cout << "Average number of CPU cycles running in hardware: "
              << hw_cycles << std::endl;
    std::cout << "Average number of CPU cycles running in commit: "
              << commit_cycles << std::endl;
    std::cout << "Average number of CPU cycles running in buf: "
              << buf_cycles << std::endl;

    std::cout << "CPU Freq: " << sds_clock_frequency() << std::endl;

    ofstream fout("fout.txt");
    for (int i = 0; i<NZ; ++i) {
        for (int j = 0; j<NX; ++j) {
            for (int k = 0; k<NY; ++k)
                fout << f[i*NX*NY+j*NY+k] << ' ';
            fout<<endl;
        }
        fout<<endl;
    }
    fout.close();

    ofstream vout("vout.txt");
    for (int i = 0; i<NZ; ++i) {
        for (int j = 0; j<NX; ++j) {
            for (int k = 0; k<NY; ++k)
                vout << v[i*NX*NY+j*NY+k]<< ' ';
            vout<<endl;
        }
        vout<<endl;
    }
    vout.close();

    free(ind);
    free(numb_set);
    free(wgt);
    free(f);
    free(v);

    free(b_line);
    free(b_wgt);

    free(commit_f);
    free(commit_v);
}
