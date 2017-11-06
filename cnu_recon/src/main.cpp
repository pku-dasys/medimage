#include "ct3d.h"
#include "utility.h"
#include "main.h"
#include <cstdio>

#include <boost/foreach.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

using namespace std;
using namespace boost::property_tree;

string RAW_DATA_FILE;

string OUTPUT_DIR;

int NX, NY, NZ;
int NPROJ, NDX, NDY;
int ITERATIONS;

float SOD;
float SDD;

int THREAD_NUMB;

void parse_config(int argc, char** argv) {
    string config_filename;
    if (argc==2) {
        config_filename = string(argv[1]);
    }
    else {
        printf("Usage: ./ct3d CONFIG_FILE_NAME\n");
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
        RAW_DATA_FILE = root.get<string>("RAW_DATA_FILE");

        OUTPUT_DIR = root.get<string>("OUTPUT_DIR");

        NX = root.get<int>("NX");
        NY = root.get<int>("NY");
        NZ = root.get<int>("NZ");
        
        NDX = root.get<int>("NDX");
        NDY = root.get<int>("NDY");
        NPROJ = root.get<int>("NPROJ");
        
        ITERATIONS = root.get<int>("ITERATIONS");

        SOD = root.get<float>("SOD");
        SDD = root.get<float>("SDD");

        THREAD_NUMB = root.get<int>("THREAD_NUMB");
    }
    catch (ptree_error &e) {
        printf("JSON file corrupted.\n");
        exit(1);
    }
}
/*
void print_options() {
    printf("CT3D reconstruction flow begins!\n");
    printf("================================\n");
    printf("Here are the options:\n");

    printf("Image size   NZ=%d   NX=%d   NY=%d\n",NZ,NX,NY);
    printf("Number of angles: %d\n",NPROJ);
    printf("Number of detector row and channel: %d, %d\n",NDX,NDY);

    printf("%d iterations with %d threads.\n",CT_ITERATIONS, THREAD_NUMB);

    printf("Source to ISO and detectors: %.1lf, %.1lf\n",SOD,SDD);

    printf("Length per detector length: %lf\n",LENGTH_PER_DET);

    printf("Results in the directory %s\n",OUTPUT_DIR.c_str());

    fflush(stdout);
}
*/
int main(int argc, char** argv) {

    parse_config(argc, argv);

    //print_options();

    float *image_data = new float[NZ*NX*NY];
    //float *edge_data = new float[NZ*NX*NY];
    ushort *sino_data = new ushort[NPROJ*NDX*NDY];

    memset(image_data, 0, sizeof(NZ*NX*NY)*sizeof(float));
    //memset(edge_data, 0, sizeof(NZ*NX*NY)*sizeof(float));
    memset(sino_data, 0, sizeof(NPROJ*NDX*NDY)*sizeof(ushort));

    ct3d(image_data,/* edge_data,*/ sino_data);

    write_data_3d(image_data, NZ, NX, NY, "CNU");
    //write_data_3d(edge_data, NZ, NX, NY, edge_filename);

    delete [] image_data;
    //delete [] edge_data;
    delete [] sino_data;

    return 0;
}
