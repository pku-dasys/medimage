#include "main.h"
#include "ms2d.h"
#include "utility.h"

#include <fstream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>
#include <cstdlib>
#include <boost/foreach.hpp>
#include <omp.h>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include<sstream>

using namespace std;
using namespace boost::property_tree;

string table_filename;
string sino_filename;
string image_filename;
string edge_filename;
string init_f_filename;
string init_v_filename;
string rebined_folder;

string OUTPUT_DIR;

int NX;
int NY;
int NZ;
double vx;
double vy;
int thread_offset;
int NDETECTOR;
int NPROJ;
int MS_ITERATIONS;

double SOURCE_TO_ISO;
double SOURCE_TO_DET;

double OFF_CENTER;

double LENGTH_PER_DET;

int THREAD_NUMB;

float ALPHA ;
float BETA ;
float EPSILON;

void parse_config(int argc,char **argv) {
    string config_filename;
    if (argc==2) {
        //THREAD_NUMB = atoi(argv[1]);
        //omp_set_num_threads(THREAD_NUMB);
        config_filename = string(argv[1]);
    }
    else {
        printf("Usage: ./ms3d CONFIG_FILE_NAME\n");
        exit(1);
    }

    ptree root;
    try {
        read_json(config_filename, root);
    }
    catch (ptree_error &e) {
        printf("Could not read from the JSON file.\n");
        exit(1);
    }

    try {
        //table_filename = root.get<string>("table_filename");
        //sino_filename = root.get<string>("sino_filename");
        image_filename = root.get<string>("image_filename");
        edge_filename = root.get<string>("edge_filename");
        init_f_filename = root.get<string>("init_f_filename");
        init_v_filename = root.get<string>("init_v_filename");
        rebined_folder = root.get<string>("rebined_folder");


        NX = root.get<int>("NX");
        NY = root.get<int>("NY");
        NZ = root.get<int>("NZ");

        NDETECTOR = root.get<int>("NDETECTOR");
        NPROJ = root.get<int>("NPROJ");

        SOURCE_TO_ISO = root.get<double>("SOURCE_TO_ISO");
        SOURCE_TO_DET = root.get<double>("SOURCE_TO_DET");

        OFF_CENTER = root.get<double>("OFF_CENTER");
        LENGTH_PER_DET = root.get<double>("LENGTH_PER_DET");

        THREAD_NUMB = root.get<int>("THREAD_NUMB");
        MS_ITERATIONS = root.get<int>("MS_ITERATIONS");

        ALPHA = root.get<double>("ALPHA");
        BETA = root.get<double>("BETA");
        EPSILON=root.get<double>("EPSILON");

        vx=root.get<double>("vy");
        vy=root.get<double>("vx");

        OUTPUT_DIR = root.get<string>("OUTPUT_DIR");
        thread_offset = NDETECTOR/THREAD_NUMB;
    }
    catch (ptree_error &e) {
        printf("JSON file corrupted.\n");
        exit(1);
    }

}

void print_options() {

    printf("MS3D reconstruction flow begins!\n");
    printf("================================\n");

    printf("Here are the options:\n");

    cout<<"Rebinned projection folder:"<<rebined_folder<<endl;
    printf("Number of projections:%d\n",NPROJ);
    printf("Image size   NX=%d   NY=%d\n",NX,NY);
    printf("vx is %f, vy is:%f\n",vx,vy);
    printf("Number of detector row %d\n",NDETECTOR);

    printf("%d iterations with %d threads.\n",MS_ITERATIONS, THREAD_NUMB);

    printf("Source to ISO and detectors: %.1lf, %.1lf\n",SOURCE_TO_ISO,SOURCE_TO_DET);

    printf("Length per detector: %lf\n",LENGTH_PER_DET);

    printf("Results in the directory %s\n",OUTPUT_DIR.c_str());
    printf("Alpha is: %lf,Beta is: %lf\n",ALPHA,BETA);
    printf("Thread offset is: %d\n",thread_offset);
    fflush(stdout);
}


int main(int argc,char **argv) {

    parse_config(argc,argv);
    //NZ = 260;
    print_options();

    double *image_data = new double[NX*NY*NZ];
    double *edge_data = new double[NX*NY*NZ];

    memset(image_data, 0, (NX*NY*NZ)*sizeof(double));
    memset(edge_data, 0, (NX*NY*NZ)*sizeof(double));
    int temp_index = 0;
    MS_ITERATIONS = 10;

    for(int i = 1; i<= NZ; i++){
      NPROJ = 900;
      double *sino_data = new double[NPROJ*NDETECTOR];
      double *temp_image_data = new double[NX*NY];
      double *temp_edge_data = new double[NX*NY];

      memset(temp_image_data, 0, NX*NY*sizeof(double));
      memset(temp_edge_data, 0, (NX*NY)*sizeof(double));
      memset(sino_data, 0, (NPROJ*NDETECTOR)*sizeof(double));
      printf("================================\n");
      printf("******************start reconstruct %d slice\n",i);

      string index_str;
      stringstream ss;
      ss<<i;
      ss>> index_str;

      table_filename = rebined_folder + index_str + "_pos";
      sino_filename = rebined_folder + index_str + "_raw";
      //image_filename =  index_str + "_f";
      //edge_filename = index_str + "_v";

      printf("***************start ms2d function!\n");

      ms2d(temp_image_data,temp_edge_data,sino_data);

      printf("transfer data started!!!!!!\n");

      for(int j = 0; j < NX *NY && temp_index < NZ*NX*NY; j++){
        image_data[temp_index] = temp_image_data[j];
        edge_data[temp_index] = temp_edge_data[j];
        temp_index++;
      }

      delete [] temp_edge_data;
      delete [] temp_image_data;
      delete [] sino_data;
    }

    write_data_3d(image_data,NZ,NX,NY,OUTPUT_DIR+"/"+image_filename);
    write_data_3d(edge_data, NZ, NX, NY, OUTPUT_DIR + "/" + edge_filename);
    delete [] image_data;
    delete [] edge_data;

    return 0;
}
