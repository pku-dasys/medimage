#include "main.h"
#include "ms2d.h"
#include "utility.h"
#include "tracing.h"
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <algorithm>
#include <omp.h>

using namespace std;
using namespace boost;

double *table_pos;
double *source_cos;
double *source_sin;
double *detector_cos;
double *detector_sin;

double lambda_f,lambda_v;

double *f, *v, *g;
int max_numb;
double AT_sum = 0;

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


    double length_per_det_rad = LENGTH_PER_DET * M_PI / 180.;
    for (int j = -HALFDET; j<HALFDET; ++j) {
        detector_cos[j + HALFDET] = cos( (j+OFF_CENTER)*length_per_det_rad + M_PI);
        detector_sin[j + HALFDET] = sin( (j+OFF_CENTER)*length_per_det_rad + M_PI);
    }
    NPROJ = line_number;
    printf("total projections: %d\n",line_number);
    printf("the number of line number is:%d\n",line_number);
    //assert( line_number == NPROJ );

}

void minIMAGE(double Af, int *line, double *weight, int numb){

  int i;
  double d[MAX_ELE_RAY];

  for(i = 0; i < numb; i++){
    Af += f[line[i]] * weight[i];
  }
  // calculate af -g


  for (i =0; i< numb; i++){
    int ind = line[i];
    int x = ind/NY;
    int y = ind%NY;

    double _f_ind = f[ind];
    double _v_ind = v[ind];

    double tmp = 0.;
    double lap = (-_f_ind)*4;

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

    //d[i]=-Af*weight[i]+ALPHA*(tmp+sqr(EPSILON)*lap);
    d[i]=-Af*weight[i]+ALPHA*(tmp);
  }
  //calculate the gradient of f
  AT_sum += sqr(Af);
  for (i = 0;i <numb; i++){
    int ind = line[i];
    f[ind] += lambda_f*d[i];
  }
}

void minEDGE(int *line,double *weight,int numb) {
    int i;

    double d[MAX_ELE_RAY];

    for (i = 0; i<numb; ++i) {
        int ind = line[i];
        int x = ind/NY, y = ind%NY;

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
        //AT_sum = AT_sum + ALPHA*_v_ind*a + BETA*EPSILON*sqr(c) + BETA*sqr(b)/(4*EPSILON);
    }
    for (i = 0; i<numb; ++i) {
        int ind = line[i];
        v[ind] += lambda_v*d[i];
        if (v[ind]<0) v[ind] = 0;
        if (v[ind]>1) v[ind] = 1;
    }

}

void compute(int nproj, int ndetector){
  double source_x = SOURCE_TO_ISO*source_cos[nproj] + NX*vx/2.0;
  double source_y = SOURCE_TO_ISO*source_sin[nproj] + NX*vy/2.0;

  double detector_x = SOURCE_TO_DET *
                      (detector_cos[ndetector] * source_cos[nproj]
                      -detector_sin[ndetector] * source_sin[nproj])
                      + source_x;
  double detector_y = SOURCE_TO_DET *
                      (detector_sin[ndetector] * source_cos[nproj]
                      +detector_cos[ndetector] * source_sin[nproj])
                      + source_y;

  int *ind = new int[MAX_ELE_RAY];
  double *wgt = new double[MAX_ELE_RAY];
  double Af = -array_2d_sino(g,nproj,ndetector);
  int numb;

  forward_proj(source_x, source_y,
               detector_x, detector_y,
               ind,wgt,numb);

  minIMAGE(Af,ind,wgt,numb);
  minEDGE(ind,wgt,numb);

  delete [] ind;
  delete [] wgt;
}

void compute_AT_sum( ){

  for(int x = 0; x < NX; x++)
    for(int y = 0; y < NY; y++){
      int ind = x*NY + y;

      double a = 0.;
      double b = 0.;
      double c = 0.;

      double _f_ind = f[ind];
      double _v_ind = v[ind];

      if (x-1>=0)      a += sqr(_f_ind-f[ind-NY]);
      else             a += sqr(_f_ind-0        );

      if (y-1>=0)      a += sqr(_f_ind-f[ind-1]);
      else             a += sqr(_f_ind-0       );


      b = _v_ind-1;

      if (x+1<NX) c += v[ind+NY];
      if (y+1<NY) c += v[ind+1];
      if (x-1>=0) c += v[ind-NY];
      if (y-1>=0) c += v[ind-1];
      c -= 4*_v_ind;

      AT_sum = AT_sum + ALPHA*sqr(_v_ind)*a + BETA*EPSILON*sqr(c) + BETA*sqr(b)/(4*EPSILON);
    }

}

void wrapper_thread(){
  int tid = omp_get_thread_num();
  for (int i =0; i < NPROJ ; i++){
    for (int j = 0; j< thread_offset; j++){
      compute(i,j*THREAD_NUMB+tid);
    }

  }
}

void wrapper(){
  for (int i =0; i < NPROJ ; i++){
    for (int j = 0; j< NDETECTOR; j++)
      compute(i,j);
  }
}

void ms2d(double *_image_data,double *_edge_data,double *_sino_data) {

    printf("alpha is :%lf, beta is: %lf, epsilon is:%lf\n",ALPHA,BETA,EPSILON);
    f = _image_data;
    v = _edge_data;
    g = _sino_data;

    table_pos = new double[NPROJ];
    source_cos = new double[NPROJ];
    source_sin = new double[NPROJ];
    detector_cos = new double[NDETECTOR];
    detector_sin = new double[NDETECTOR];

    load_table(table_pos,
               source_cos, source_sin,
               detector_cos, detector_sin);

    //for (int i = 0; i<NPROJ*NDETECTOR; ++i) *(g+i) = 0.;

    ifstream sino_in(sino_filename.c_str());
    for (int k = 0; k<NPROJ; ++k) {
            for (int i = 0; i<NDETECTOR; ++i)
                sino_in>>g[k*NDETECTOR+i];
    }
    sino_in.close();
    for (int i = 0; i<NX*NY; ++i) *(f+i) = 0.;
    for (int i = 0; i<NX*NY; ++i) *(v+i) = 1.;
    /*
    {
      printf("initialize initial f and v ...\n");
      fflush(stdout);
      ifstream f_in("124");
      double min_f = 1e20;
      for (int i = 181; i<NX-181; ++i)
          for (int j = 181; j<NY-181; ++j) {
              f_in >> f[i*NY+j];
              if (min_f>f[i*NY+j])
                  min_f = f[i*NY+j];
          }
      f_in.close();
    }
    */

    //write_data_3d(f,NZ,NX,NY,"slices/"+image_filename+"_"+lexical_cast<string>(iters));
    printf("using %d threads ...\n",THREAD_NUMB);
    fflush(stdout);
    AT_sum = 0;
    omp_set_dynamic(0);
    omp_set_num_threads(THREAD_NUMB);
    for (int iters = 0; iters<MS_ITERATIONS; ++iters) {

        lambda_f = 1.0/(2000.0+iters*1000.0);
        lambda_v = 1.0 /(2000.0 + iters * 1000.0);

        printf("iteration = %d, lambda_f = %lf, lambda_v = %f\n",iters,lambda_f,lambda_v);
        AT_sum = 0;
        fflush(stdout);


        {
#pragma omp parallel num_threads(THREAD_NUMB)
        wrapper_thread();
        //wrapper();
        }
        //compute_AT_sum();

        printf("AT_sum = %E\n",AT_sum);

        //compute_group();


        fflush(stdout);

        if (iters >= 4) {
            //write_data_2d(f,NX,NY,OUTPUT_DIR+"/"+image_filename+"_"+lexical_cast<string>(iters)+"_"+lexical_cast<string>(ALPHA)+"_"+lexical_cast<string>(BETA));
            //write_data_2d(f,NX,NY,OUTPUT_DIR+"/"+image_filename+"_"+lexical_cast<string>(iters));

            //write_data_3d(v,NZ,NX,NY,OUTPUT_DIR+"/"+edge_filename+"_"+lexical_cast<string>(iters));
        }
        printf("write image data finished!\n");

        //lambda = lambda/(1 + 100 * lambda);
    }


    fflush(stdout);

    delete [] table_pos;
    delete [] source_cos;
    delete [] source_sin;
    delete [] detector_cos;
    delete [] detector_sin;
    //return AT_sum;

}
