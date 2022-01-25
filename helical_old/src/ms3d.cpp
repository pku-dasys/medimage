#include "main.h"
#include "utility.h"

#include "tracing.h"
#include "ms3d.h"

#include <omp.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include <boost/lexical_cast.hpp>

#include <fstream>
#include <algorithm>

using namespace std;
using namespace boost;

/*
double table_pos[NPROJ];
double source_cos[NPROJ];
double source_sin[NPROJ];
double detector_cos[NDETECTOR];
double detector_sin[NDETECTOR];
*/

double *table_pos;
double *source_cos;
double *source_sin;
double *detector_cos;
double *detector_sin;

double lambda;

double *f, *v, *g;

int max_numb;

void load_table(double *table_pos,double *source_cos,double *source_sin,
                double *detector_cos,double *detector_sin) {
    FILE *input = fopen(table_filename.c_str(),"r");

    if (input==NULL) {
        printf("cannot open table file\n");
        exit(1);
    }

    printf("loading table information ...\n");

    double pos,angle;
    int line_number = 0;
    while (fscanf(input,"%lf,%lf\n",&pos,&angle)==2) {
        table_pos[line_number] = pos/1000.0;
        source_cos[line_number] = cos(angle/180.0*M_PI);
        source_sin[line_number] = sin(angle/180.0*M_PI);
        ++line_number;
    }

    fclose(input);

    for (int i = 1; i<line_number; ++i) {
        table_pos[i] -= table_pos[0];
        // do we need this?
        //table_pos[i] += 5.0;
    }

    double length_per_det_rad = LENGTH_PER_DET * M_PI / 180.;
    for (int j = -HALFDET; j<HALFDET; ++j) {
        detector_cos[j + HALFDET] = cos( (j+OFF_CENTER)*length_per_det_rad + M_PI);
        detector_sin[j + HALFDET] = sin( (j+OFF_CENTER)*length_per_det_rad + M_PI);
    }

    printf("total projections: %d\n",line_number);

    assert( line_number == NPROJ );

}

double AT_sum;

void minIMAGE(double Af,int *line,double *weight,int numb) {
    int i;

    double d[MAX_ELE_RAY];
    for (i = 0; i<numb; ++i) {
        Af += f[line[i]]*weight[i];
    }
    // calculate div (v^2 \grad f) and Af-g

/* 
    for (i = 0; i<numb; ++i) {
        int ind = line[i];

        int ind2d = line[i]%(NX*NY);

        int x = ind2d/NY, y = ind2d%NY;


        double _f_ind = f[ind];
        double _v_ind = v[ind];

        Af += _f_ind*weight[i];

        double tmp = 0.;
        double lap = -_f_ind;

        if (x+1<NX) {
            tmp += sqr(_v_ind)*(f[ind+NY]-_f_ind);
            lap += f[ind+NY];
        }
        else
            tmp += sqr(_v_ind)*(    0    -_f_ind);

        if (y+1<NY) {
            tmp += sqr(_v_ind)*(f[ind+1]-_f_ind);
            lap += f[ind+1];
        }
        else
            tmp += sqr(_v_ind)*(   0    -_f_ind);

        if (x-1>=0) {
            tmp -= sqr(v[ind-NY])*(_f_ind-f[ind-NY]);
            lap += f[ind-NY];
        }
        else
            tmp -=                (_f_ind-0        );

        if (y-1>=0) {
            tmp -= sqr(v[ind-1])*(_f_ind-f[ind-1]);
            lap += f[ind-1];
        }
        else
            tmp -=               (_f_ind-0       );



        //double tmp_d = -Af_minus_g*weight[i]+ALPHA*(tmp+sqr(EPSILON)*lap);
        d[i] = +ALPHA*(tmp+sqr(EPSILON)*lap);
        //f[ind] += lambda*tmp_d;
        //if (f[ind]<0) f[ind] = 0;
    }
   
    */
    AT_sum += sqr(Af);
   
    for (i = 0; i<numb; ++i) {
        int ind = line[i];
        
        f[ind] += lambda*(-Af*weight[i]);
        //f[ind] += (-Af*weight[i]);
        //f[ind] += lambda*(-Af*weight[i]+d[i]);
        //printf("%lf %lf %lf\n",lambda*(-Af*weight[i]+d[i]),  lambda*-Af*weight[i],    lambda*d[i]);
        //f[ind] += (-Af*weight[i]+d[i]);
        if (f[ind]<0) f[ind] = 0;
        //if (f[ind]>75) f[ind] = 75;
    }

}

void minEDGE(int *line,double *weight,int numb) {
    int i;

    double d[MAX_ELE_RAY];

    for (i = 0; i<numb; ++i) {
        int ind = line[i];
        int ind2d = line[i]%(NX*NY);

        int x = ind2d/NY, y = ind2d%NY;

        double a = 0.;
        double b = 0.;
        double c = 0.;

        double _f_ind = f[ind];
        double _v_ind = v[ind];

        if (x-1>=0)      a += sqr(_f_ind-f[ind-NY]);
        else             a += sqr(_f_ind-0        );

        if (y-1>=0)      a += sqr(_f_ind-f[ind-1]);
        else             a += sqr(_f_ind-0       );

        a *= _v_ind;

        b = _v_ind-1;

        if (x+1<NX) c += v[ind+NY];
        if (y+1<NY) c += v[ind+1];
        if (x-1>=0) c += v[ind-NY];
        if (y-1>=0) c += v[ind-1];
        c -= 4*_v_ind;

        d[i] = -ALPHA*a-BETA/(4*EPSILON)*b+BETA*EPSILON*c;
        //double tmp_d = -ALPHA*a-BETA/(4*EPSILON)*b+BETA*EPSILON*c;
        //v[ind] += lambda*tmp_d;
        //if (v[ind]<0) v[ind] = 0;
        //if (v[ind]>1) v[ind] = 1;

    }

    for (i = 0; i<numb; ++i) {
        int ind = line[i];
        v[ind] += lambda*d[i];
        if (v[ind]<0) v[ind] = 0;
        if (v[ind]>1) v[ind] = 1;
    }

}

void compute(int nproj, int ndetector, int nchannel) {
    // ndetector from [-half, half]

    double source_x = SOURCE_TO_ISO*source_cos[nproj] + NX/2.0;
    double source_y = SOURCE_TO_ISO*source_sin[nproj] + NX/2.0;
    double source_z = table_pos[nproj];

    double detector_x = SOURCE_TO_DET *
                        (detector_cos[ndetector] * source_cos[nproj]
                        -detector_sin[ndetector] * source_sin[nproj])
                        + source_x;
    double detector_y = SOURCE_TO_DET *
                        (detector_sin[ndetector] * source_cos[nproj]
                        +detector_cos[ndetector] * source_sin[nproj])
                        + source_y;

    double detector_z = (nchannel - HALFDETZ + 0.5) * LENGTH_PER_DET_Z + source_z;



    int *ind = new int[MAX_ELE_RAY];
    double *wgt = new double[MAX_ELE_RAY];

    //int tid = omp_get_thread_num();
    double Af = -array_3d_sino(g,nproj,ndetector,nchannel);
    //printf("TID = %d, Af = %d\n",tid,Af);
    int numb;
    
    //printf("%lf\n",Af);

    {
    //int64_t it_start = timer_us();
    forward_proj(source_x, source_y, source_z,
                 detector_x, detector_y, detector_z,
                 ind,wgt,numb);
    //int64_t it_end = timer_us();
    //printf("forward_proj time = %ld us\n",it_end-it_start);
    }
    

    {
    //int64_t it_start = timer_us();
    minIMAGE(Af,ind,wgt,numb);
    //int64_t it_end = timer_us();
    //printf("minIMAGE time = %ld us\n",it_end-it_start);
    }

    {
    //int64_t it_start = timer_us();
    //minEDGE(ind,wgt,numb);
    //int64_t it_end = timer_us();
    //printf("minEDGE time = %ld us\n",it_end-it_start);
    }

    delete [] ind;
    delete [] wgt;
}

void wrapper() {
    /*
    int tid = omp_get_thread_num();
    for (int k = 0; k<NPROJ; ++k) {
        for (int offset = 0; offset<INTERVAL; ++offset) {
            for (int i = tid * INTERVAL + offset; i<NDETECTOR; i += INTERVAL*THREAD_NUMB) {
                for (int j = 0; j<NCHANNEL; ++j) {
                    compute(k,i,j);
                }
            }
        }
    }
    */
    int tid = omp_get_thread_num();
    for (int k = 0; k<NPROJ; ++k) {
        for (int i = 0; i<NDETECTOR; ++i) {
            compute(k,i,tid);
        }
    }
    
}

const int GROUP = 16;

double *hist_val;
double *new_hist_val;

int *group_index;

double sum_val;
int count_val;
double *sub_sum;
int *sub_count;

void center_wrapper() {
    int tid = omp_get_thread_num();
    int interval = (NZ*NX*NY)/THREAD_NUMB;
    int st = tid*interval;
    int ed = (tid+1)*interval;
    for (int i = st; i<ed; ++i) {
        double min_dis = 1e20, pos = -1;
        for (int k = 0; k<GROUP; ++k) {
            double tmp = fabs(hist_val[k]-f[i]);
            if (min_dis>tmp) {
                min_dis = tmp;
                pos = k;
            }
        }
        group_index[i] = pos;
    }
}

void count_wrapper(int group_id) {
    int tid = omp_get_thread_num();
    int interval = (NZ*NX*NY)/THREAD_NUMB;
    int st = tid*interval;
    int ed = (tid+1)*interval;
    for (int i = st; i<ed; ++i) {
        if (group_index[i]==group_id) {
            sub_sum[tid] += f[i];
            ++sub_count[tid];
        }
    }
}

void dohist_wrapper() {
    int tid = omp_get_thread_num();
    int interval = (NZ*NX*NY)/THREAD_NUMB;
    int st = tid*interval;
    int ed = (tid+1)*interval;
    for (int i = st; i<ed; ++i) {
        f[i] += lambda*500*(hist_val[group_index[i]]-f[i]);
    }
}

void print_group(double *x) {
    for (int i = 0; i<GROUP; ++i) {
        printf("%f \t",x[i]);
    }
    printf("\n");
}

bool first_time = true;

void compute_group() {
    if (first_time) {
    for (int i = 0; i<GROUP; ++i)
        hist_val[i] = f[rand()%(NZ*NX*NY)];
    sort(hist_val, hist_val+GROUP);
    first_time = false;
    }
    
    print_group(hist_val);
    
    int iters = 0;
    while (true) {
        {
#pragma omp parallel num_threads(THREAD_NUMB)
            center_wrapper();
        }
        for (int i = 0; i<GROUP; ++i) {
            
            for (int k = 0; k<THREAD_NUMB; ++k) {
                sub_sum[k] = 0;
                sub_count[k] = 0;
            }
            sum_val = 0;
            count_val = 0;
            {
    #pragma omp parallel num_threads(THREAD_NUMB)
                count_wrapper(i);
            }
            for (int k = 0; k<THREAD_NUMB; ++k) {
                sum_val += sub_sum[k];
                count_val += sub_count[k];
            }
            
            new_hist_val[i] = sum_val/(count_val+1e-6);
            
        }
        
        double diff = 0;
        for (int i = 0; i<GROUP; ++i) {
            diff += fabs(new_hist_val[i]-hist_val[i]);
        }
        printf("iteration %d, diff = %f\n",iters,diff);
        fflush(stdout);
        
        for (int i = 0; i<GROUP; ++i) {
            hist_val[i] = new_hist_val[i];
        }
        print_group(hist_val);
        
        if (diff<0.1)
            break;
        
        ++iters;
    }
    {
#pragma omp parallel num_threads(THREAD_NUMB)
        dohist_wrapper();
    }
}

void ms3d(double *_image_data,double *_edge_data,double *_sino_data) {



    f = _image_data;
    v = _edge_data;
    g = _sino_data;
    
    hist_val = new double[GROUP];
    new_hist_val = new double[GROUP];
    group_index = new int[NZ*NX*NY];
    sub_sum = new double[THREAD_NUMB];
    sub_count = new int[THREAD_NUMB];
    
    ifstream sino_in(sino_filename.c_str());
    for (int k = 0; k<NPROJ; ++k) {
        for (int j = NCHANNEL-1; j>=0; --j)
            for (int i = 0; i<NDETECTOR; ++i)
                sino_in >> g[k*NDETECTOR*NCHANNEL+i*NCHANNEL+j];
    }
    sino_in.close();
    //write_data_3d(f,NZ,NX,NY,"slices/"+image_filename+"_"+lexical_cast<string>(iters));

    table_pos = new double[NPROJ];
    source_cos = new double[NPROJ];
    source_sin = new double[NPROJ];
    detector_cos = new double[NDETECTOR];
    detector_sin = new double[NDETECTOR];

    load_table(table_pos,
               source_cos, source_sin,
               detector_cos, detector_sin);

    
    for (int i = 0; i<NZ*NX*NY; ++i) *(f+i) = 0.;
    for (int i = 0; i<NZ*NX*NY; ++i) *(v+i) = 1.;
    {
        printf("initialize initial f and v ...\n");
        fflush(stdout);
        ifstream f_in(init_f_filename);
        double min_f = 1e20;
        for (int k = PROJ_OFFSET; k<PROJ_OFFSET+INIT_PROJ; ++k) {
            for (int i = 0; i<NX; ++i)
                for (int j = 0; j<NY; ++j) {
                    f_in >> f[k*NX*NY+i*NY+j];
                    if (min_f>f[k*NX*NY+i*NY+j])
                        min_f = f[k*NX*NY+i*NY+j];
                }
        }
        for (int k = PROJ_OFFSET; k<PROJ_OFFSET+INIT_PROJ; ++k) {
            for (int i = 0; i<NX; ++i)
                for (int j = 0; j<NY; ++j) {
                    f[k*NX*NY+i*NY+j] -= min_f;
                }
        }
        f_in.close();
        /*
        ifstream v_in(init_v_filename);
        for (int k = PROJ_OFFSET; k<PROJ_OFFSET+INIT_PROJ; ++k) {
            for (int i = 0; i<NX; ++i)
                for (int j = 0; j<NY; ++j) {
                    v_in >> v[k*NX*NY+i*NY+j];
                    v[k*NX*NY+i*NY+j] = 1-v[k*NX*NY+i*NY+j];
                }
        }
        v_in.close();
        */
        printf("initialization complete.\n");
        fflush(stdout);
    }

    //lambda = 5e-3;
    int64_t start = timer_s();

    printf("using %d threads ...\n",THREAD_NUMB);
    fflush(stdout);
    
    omp_set_dynamic(0);
    omp_set_num_threads(THREAD_NUMB);
    
    for (int iters = 0; iters<MS_ITERATIONS; ++iters) {

        lambda = 1.0/(1000.0+iters*50.0);
        
        printf("iteration = %d, lambda = %lf\n",iters,lambda);
        AT_sum = 0;
        fflush(stdout);
        int64_t it_start = timer_s();
        {
#pragma omp parallel num_threads(THREAD_NUMB)
        wrapper();
        }
        printf("AT_sum = %E\n",AT_sum);
        
        //compute_group();
        
        int64_t it_end = timer_s();
        printf("iteration time = %ld s\n",it_end-it_start);
        fflush(stdout);

        /*if (iters%5==0)*/ {
            write_data_3d(f,NZ,NX,NY,OUTPUT_DIR+"/"+image_filename+"_"+lexical_cast<string>(iters));
            //write_data_3d(v,NZ,NX,NY,OUTPUT_DIR+"/"+edge_filename+"_"+lexical_cast<string>(iters));
        }

        //lambda = lambda/(1 + 100 * lambda);
    }
    int64_t end = timer_s();
    printf("total iteration time = %ld s\n",end-start);
    fflush(stdout);

    delete [] table_pos;
    delete [] source_cos;
    delete [] source_sin;
    delete [] detector_cos;
    delete [] detector_sin;
    
    delete [] hist_val;
    delete [] new_hist_val;
    delete [] group_index;
    delete [] sub_sum;
    delete [] sub_count;
}