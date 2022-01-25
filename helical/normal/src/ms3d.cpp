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
#include <ctime>
#include <chrono>
#include <sys/time.h>

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

double lambda_f;
double lambda_v;

double *f, *v, *g;

double tracing_time;
double compute_Af_time;
double min_Image_time;
double min_Edge_time;
double fetch_time;
double commit_time;


int max_numb;
int INNER_ITERATIONS = 1;

long int timer_us(void)
{
#ifdef WIN32
    struct _timeb tb;
    _ftime(&tb);
    return ((int64_t)tb.time * (1000) + (int64_t)tb.millitm) * (1000);
#else
    struct timeval tv_date;
    gettimeofday(&tv_date, NULL);
    return ((long int)tv_date.tv_sec * 1000000 + (long int)tv_date.tv_usec);
#endif
}

void load_table(double *table_pos, double *source_cos, double *source_sin,
               double *detector_cos, double *detector_sin)
{
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
        //table_pos[line_number] = pos;
        source_cos[line_number] = cos(angle/180.0*M_PI);
        source_sin[line_number] = sin(angle/180.0*M_PI);
        ++line_number;
    }

    printf("loading table information finished...\n");
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
    //NPROJ = line_number;
    assert( line_number == NPROJ );

}

double AT_sum;
double Afg, alpha_reg, beta_reg;


void compute_Af(int *line, double *weight, int numb, double *loc_f, double &Af){
    //Af = -g;
    //compute_Af based local_f
    int i;
    for (i = 0; i < numb; ++i)
    {
        //Af += f[line[i]] * weight[i];
        Af += f[i*5+0] * weight[i];
    }
}

/*
void compute_Af(int *line, double *weight, int numb, double &Af)
{
    //Af = -g;
    int i;
    for (i = 0; i < numb; ++i)
    {
        Af += f[line[i]] * weight[i];
        //Af += f[i * 5 + 0] * weight[i];
    }
}
*/

/*
void minIMAGE(double Af, int *line, double *weight, int numb)
{
    int i;

    double d[MAX_ELE_RAY];
    
    for (i = 0; i<numb; ++i) {
        int ind = line[i];
        int ind2d = line[i]%(NX*NY);

        int z=line[i]/(NX*NY);

        int x = ind2d/NY, y = ind2d%NY;


        double _f_ind = f[ind];
        double _v_ind = v[ind];

        //Af += _f_ind*weight[i];

        double tmp = 0.;
        double lap = -6*_f_ind;

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
            tmp -= sqr(_v_ind)*(_f_ind-f[ind-NY]);
            lap += f[ind-NY];
        }
        else
            tmp -=                (_f_ind-0        );

        if (y-1>=0) {
            tmp -= sqr(_v_ind)*(_f_ind-f[ind-1]);
            lap += f[ind-1];
        }
        else
            tmp -=               (_f_ind-0       );
        

        //double tmp_d = -Af*weight[i]+ALPHA*(tmp+sqr(EPSILON)*lap);
        d[i]=-Af*weight[i]+ALPHA*(tmp+sqr(EPSILON)*lap);
        //f[ind] += lambda*tmp_d;
        //if (f[ind]<0) f[ind] = 0;
    }


    for (i = 0; i<numb; ++i) {
        int ind = line[i];
        f[ind] += lambda_f*d[i];
        if (f[ind]<0) f[ind] = 0;
    }
}
*/




void minIMAGE(double Af, int *line, double *weight, int numb, double* loc_f, double *loc_v){
    //min_edge based on the local_f and local_v
    int i;
    double d[MAX_ELE_RAY];
    for (i = 0; i < numb; ++i)
    {
        int ind = line[i];
        int ind2d = line[i] % (NX * NY);
        int x = ind2d / NY, y = ind2d % NY;

        double _f_ind = loc_f[i*5 + 0];
        double _v_ind = loc_v[i*5 + 0];

        //Af += _f_ind*weight[i];

        double tmp = 0.;
        double lap = -6 * _f_ind;

        if (x + 1 < NX)
        {
            tmp += sqr(_v_ind) * (loc_f[i*5 + 3] - _f_ind);
            lap += loc_f[i*5 + 3];
        }
        else
            tmp += sqr(_v_ind) * (0 - _f_ind);

        if (y + 1 < NY)
        {
            tmp += sqr(_v_ind) * (loc_f[i*5 +1] - _f_ind);
            lap += loc_f[i*5 + 1];
        }
        else
            tmp += sqr(_v_ind) * (0 - _f_ind);

        
        if (x - 1 >= 0)
        {
            tmp -= sqr(_v_ind) * (_f_ind - loc_f[i*5 + 4]);
            lap += loc_f[i*5 + 4];
        }
        else
            tmp -= (_f_ind - 0);

        if (y - 1 >= 0)
        {
            tmp -= sqr(_v_ind) * (_f_ind - loc_f[i*5 +2]);
            lap += loc_f[i*5 +2];
        }
        else
            tmp -= (_f_ind - 0);

        
        //double tmp_d = -Af*weight[i]+ALPHA*(tmp+sqr(EPSILON)*lap);
        d[i] = -Af * weight[i] + ALPHA * (tmp + sqr(EPSILON) * lap);
        //f[ind] += lambda*tmp_d;
        //if (f[ind]<0) f[ind] = 0;
    }

    for(int i = 0; i<numb; i++){
        loc_f[i*5 + 0] += lambda_f*d[i];
    }

}




/*
void minEDGE(int *line, double *weight, int numb)
{
    int i;

    double d[MAX_ELE_RAY];

    for (i = 0; i<numb; ++i) {
        int ind = line[i];
        int ind2d = line[i]%(NX*NY);

        int x = ind2d/NY, y = ind2d%NY;
        int z = line[i] / (NX * NY);

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

        if (x+1<NX) c += v[ind+NY] - _v_ind;
        if (y+1<NY) c += v[ind+1] - _v_ind;
        if (x-1>=0) c += v[ind-NY] - _v_ind;
        if (y-1>=0) c += v[ind-1] - _v_ind;

        d[i] = -ALPHA*a-BETA/(4*EPSILON)*b+BETA*EPSILON*c;

    }

    for (i = 0; i<numb; ++i) {
        int ind = line[i];
        v[ind] += lambda_v*d[i];
        if (v[ind]<0) v[ind] = 0;
        if (v[ind]>1) v[ind] = 1;
    }

}
*/


void minEDGE(int *line, double *weight, int numb,double* loc_f, double *loc_v)
{
    // min EDge based the local_f and local_v
    int i;

    double d[MAX_ELE_RAY];

    for (i = 0; i < numb; ++i)
    {
        int ind = line[i];
        int ind2d = line[i] % (NX * NY);
        int x = ind2d / NY, y = ind2d % NY;
        double a = 0.;
        double b = 0.;
        double c = 0.;

        double _f_ind = loc_f[i*5 + 0];
        double _v_ind = loc_v[i*5 + 0];

        if (x - 1 >= 0)
            a += sqr(_f_ind - loc_f[i*5 + 4]);
        else
            a += sqr(_f_ind - 0);

        if (y - 1 >= 0)
            a += sqr(_f_ind - loc_f[i*5 + 2]);
        else
            a += sqr(_f_ind - 0);

        a *= _v_ind;

        b = _v_ind - 1;

        if (x + 1 < NX)
            c += loc_v[i*5 + 3] - _v_ind;
        if (y + 1 < NY)
            c += loc_v[i*5 + 1] - _v_ind;
        if (x - 1 >= 0)
            c += loc_v[i*5 + 4]- _v_ind;
        if (y - 1 >= 0)
            c += loc_v[i*5 +2] - _v_ind;
        

        d[i] = -ALPHA * a - BETA / (4 * EPSILON) * b + BETA * EPSILON * c;
        //double tmp_d = -ALPHA*a-BETA/(4*EPSILON)*b+BETA*EPSILON*c;
        //v[ind] += lambda*tmp_d;
        //if (v[ind]<0) v[ind] = 0;
        //if (v[ind]>1) v[ind] = 1;
    }

    for (i = 0; i < numb; ++i)
    {
        int ind = line[i];
        loc_v[i*5 + 0] += lambda_v * d[i];
        if (loc_v[i * 5 + 0] < 0)
            loc_v[i * 5 + 0] = 0;
        if (loc_v[i * 5 + 0] > 1)
            loc_v[i * 5 + 0] = 1;
    }
}


void compute(int nproj, int ndetector, int nchannel,int iters) {
    // ndetector from [-half, half]

    double* loc_f = new double[MAX_ELE_RAY*5];
    double* loc_v = new double[MAX_ELE_RAY*5];

    double source_x = SOURCE_TO_ISO*source_cos[nproj] + NX*vx/2.0;
    double source_y = SOURCE_TO_ISO*source_sin[nproj] + NX*vy/2.0;
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
    //printf("%f %f %f %f %f %f\n",source_x,source_y,source_z,detector_x,detector_y,detector_z);

    int *ind = new int[MAX_ELE_RAY];
    double *wgt = new double[MAX_ELE_RAY];

    //int tid = omp_get_thread_num();
    double Af = -array_3d_sino(g,nproj,ndetector,nchannel);
    //printf("TID = %d, Af = %d\n",tid,Af);
    int numb;

    //auto start_forward = std::chrono::system_clock::now();
    long int start_forward = timer_us();
    forward_proj(source_x, source_y, source_z,
                 detector_x, detector_y, detector_z,
                 ind,wgt,numb);
    long int end_forward = timer_us();
    tracing_time += double(end_forward - start_forward);
    //auto end_forward = std::chrono::system_clock::now();
    //std::chrono::duration<double> forward_incre = (end_forward - start_forward);
    //tracing_time += forward_incre.count();

    
    //read memory
    //auto start_fetch = std::chrono::system_clock::now();
    long int start_fetch = timer_us();
    for (int i = 0; i < numb; i++)
    {
        int temp_ind = ind[i];
        loc_f[i * 5 + 0] = f[temp_ind];
        loc_f[i * 5 + 1] = f[temp_ind + 1];
        loc_f[i * 5 + 2] = f[temp_ind - 1];
        loc_f[i * 5 + 3] = f[temp_ind + NY];
        loc_f[i * 5 + 4] = f[temp_ind - NY];

        loc_v[i * 5 + 0] = v[temp_ind];
        loc_v[i * 5 + 1] = v[temp_ind + 1];
        loc_v[i * 5 + 2] = v[temp_ind - 1];
        loc_v[i * 5 + 3] = v[temp_ind + NY];
        loc_v[i * 5 + 4] = v[temp_ind - NY];
    }
    long int end_fetch = timer_us();
    fetch_time += (end_fetch - start_fetch);
    //auto end_fetch = std::chrono::system_clock::now();
    //std::chrono::duration<double> incre_fetch = (end_fetch - start_fetch);
    //fetch_time += incre_fetch.count();
    
    //calc_Af
    //auto start_calcAf = std::chrono::system_clock::now();
    long int calc_Af_start = timer_us();
    compute_Af(ind,wgt,numb,loc_f,Af);
    //compute_Af(ind, wgt, numb, Af);
    long int calc_Af_end = timer_us();
    compute_Af_time += double(calc_Af_end - calc_Af_start);
    //auto end_calcAf = std::chrono::system_clock::now();
    //std::chrono::duration<double> incre_calcAf = (end_calcAf - start_calcAf);
    //compute_Af_time += incre_calcAf.count();

    //double init_lambda_f = 1000.0;
    //double init_lambda_v = 10000.0;
    lambda_f = 1/2000;
    lambda_v = 1/10000;

    //min_Image
    //auto start_minImage = std::chrono::system_clock::now();
    long int start_minImage = timer_us();
    minIMAGE(Af, ind, wgt, numb, loc_f, loc_v);
    //minIMAGE(Af, ind, wgt, numb);
    long int end_minImage = timer_us();
    min_Image_time += (end_minImage - start_minImage);
    //auto end_minImage = std::chrono::system_clock::now();
    //std::chrono::duration<double> incre_minImg = (end_minImage - start_minImage);
    //min_Image_time += incre_minImg.count();

    //min_Edge
    //auto start_minEdge = std::chrono::system_clock::now();
    long int start_minEdge = timer_us();
    minEDGE(ind, wgt, numb,loc_f,loc_v);
    //minEDGE(ind, wgt, numb);
    long int end_minEdge = timer_us();
    min_Edge_time += (end_minEdge - start_minEdge);
    //auto end_minEdge = std::chrono::system_clock::now();
    //std::chrono::duration<double> incre_minEdge = (end_minEdge - start_minEdge);
    //min_Edge_time += incre_minEdge.count();


    
    //commit memory
    //auto start_commit = std::chrono::system_clock::now();
    long int start_commit = timer_us();
    for(int i = 0; i<numb; i++){
        int temp_ind = ind[i];
        f[temp_ind] = loc_f[i*5 + 0];
        v[temp_ind] = loc_v[i*5 + 0];
    }
    long int end_commit = timer_us();
    commit_time += (end_commit - start_commit);
    //auto end_commit = std::chrono::system_clock::now();
    //std::chrono::duration<double> incre_commit = (end_commit - start_commit);
    //commit_time += incre_commit.count();
    

    double total_iterations = MS_ITERATIONS * INNER_ITERATIONS;

    /*
    for (int loop = 0; loop<INNER_ITERATIONS; loop++){
        double iter_percent = iters * INNER_ITERATIONS + loop;
        lambda_f = 1 / (init_lambda_f + 4 * init_lambda_f * (iter_percent / total_iterations));
        minIMAGE(Af, ind, wgt, numb);
    }
    */
    
    
    /*
    for (int loop = 0; loop < INNER_ITERATIONS; loop++){
        double iter_percent = iters * INNER_ITERATIONS + loop;
        lambda_v = 1 / (init_lambda_v + 4 * init_lambda_v * (iter_percent / total_iterations));
        minEDGE(ind, wgt, numb);
    }
    */
    
    delete [] loc_f;
    delete [] loc_v;
    delete [] ind;
    delete [] wgt;
}

/*
void wrapper() {
    int tid = omp_get_thread_num();
    for (int k = 0; k<NPROJ; ++k) {
        for (int i = 0; i<NDETECTOR; ++i) {
            compute(k,i,tid);
        }
    }

}
*/
/*
void wrapper(int iters) {
    int tid = omp_get_thread_num();
    // initiallize the sino data
    for (int k = 0; k<NPROJ; ++k) {
        for (int i = 0; i<NDETECTOR; ++i) {
            for(int j =0; j<NCHANNEL/THREAD_NUMB; j++){
              compute(k,i,j*THREAD_NUMB+tid,iters);
            }
        }
    }
}
*/

void wrapper(int iters)
{
    tracing_time = 0.0;
    compute_Af_time = 0.0;
    min_Image_time = 0.0;
    min_Edge_time = 0.0;
    fetch_time = 0.0;
    commit_time = 0.0;
    NPROJ = 4000;
    for (int k = 0; k < NPROJ; ++k)
    {
        for (int i = 0; i < NDETECTOR; ++i)
        {
            for (int j = 0; j < NCHANNEL; j++)
            {
                compute(k, i, j, iters);
            }
        }
    }

    double ave_tracing_time = tracing_time/(NPROJ*NDETECTOR*NCHANNEL);
    double ave_compute_Af_time = compute_Af_time / (NPROJ * NDETECTOR * NCHANNEL);
    double ave_min_Image_time = min_Image_time / (NPROJ * NDETECTOR * NCHANNEL);
    double ave_min_Edge_time = min_Edge_time / (NPROJ * NDETECTOR * NCHANNEL);
    double ave_fetch_time = fetch_time/(NPROJ * NDETECTOR * NCHANNEL);
    double ave_commit_time = commit_time/(NPROJ * NDETECTOR *NCHANNEL);
    printf("tracing time is: %f, calc_Af time is: %f, min_Image_time is: %f, min_Edge_time is: %f\n",
           ave_tracing_time, ave_compute_Af_time, ave_min_Image_time, ave_min_Edge_time);
    printf("fetch time is:%f, commit_time is: %f\n", ave_fetch_time, ave_commit_time);
}

void ms3d(double *_image_data,double *_edge_data,double *_sino_data) {

    printf("alpha is :%lf, beta is: %lf, epsilon is:%lf\n",ALPHA,BETA,EPSILON);
    f = _image_data;
    v = _edge_data;
    g = _sino_data;

    tracing_time = 0.0;
    compute_Af_time = 0.0;
    min_Image_time = 0.0;
    min_Edge_time = 0.0;


    ifstream sino_in(sino_filename.c_str());
    for (int k = 0; k<NPROJ; ++k) {
        for (int j = NCHANNEL-1; j>=0; --j)
            for (int i = 0; i<NDETECTOR; ++i)
                sino_in >> g[k*NDETECTOR*NCHANNEL+i*NCHANNEL+j];
    }
    /*
    for (int k = 0; k<NPROJ; ++k) {
        for (int i = 0; i<NDETECTOR; ++i)
            for (int j = 0; j < NCHANNEL; j++)
                sino_in >> g[k*NDETECTOR*NCHANNEL+i*NCHANNEL+j];
    }
    */
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
    // int64_t start = timer_s();

    printf("using %d threads ...\n",THREAD_NUMB);
    fflush(stdout);

    omp_set_dynamic(0);
    omp_set_num_threads(THREAD_NUMB);

    for (int iters = 0; iters<MS_ITERATIONS; ++iters) {

        lambda_f = 1.0/(1000.0+iters*500.0);
        lambda_v = 1.0 / (5000.0 + iters * 2000.0);

        printf("iteration = %d, lambda_f = %lf, lambda_v = %lf\n", iters, lambda_f, lambda_v);
        AT_sum = 0;
        fflush(stdout);
        time_t it_start,it_end;
        it_start = time(NULL);
        // int64_t it_start = timer_s();
        {
        wrapper(iters);
        }

        //compute_AT();
        //printf("AT_sum = %E,Afg = %E, alpha_reg = %E, beta_reg =%E\n", AT_sum, Afg, alpha_reg, beta_reg);
        it_end = time(NULL);
        //compute_group();

        // int64_t it_end = timer_s();
        printf("iteration time = %ld s\n",it_end-it_start);
        fflush(stdout);

        /*if (iters == MS_ITERATIONS-1)*/{
            string alpha;string beta; string det_z;
            stringstream ss1,ss2,ss3;
            ss1<<ALPHA;ss1>>alpha;
            ss2<<BETA;  ss2>>beta;
            ss3<<LENGTH_PER_DET_Z; ss3>>det_z;

            write_data_3d(f,NZ,NX,NY,OUTPUT_DIR+"/"+image_filename+"_"+lexical_cast<string>(iters)+"_alpha_"+alpha+"_beta_"+beta+"_detz_"+det_z);
            write_data_3d(v,NZ,NX,NY,OUTPUT_DIR+"/"+edge_filename+"_"+lexical_cast<string>(iters)+"_alpha_"+alpha+"_beta_"+beta+"_detz_"+det_z);
        }

        
    }
    
    fflush(stdout);

    delete [] table_pos;
    delete [] source_cos;
    delete [] source_sin;
    delete [] detector_cos;
    delete [] detector_sin;
}
