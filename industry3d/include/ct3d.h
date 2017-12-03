#pragma once

#include <string>

using namespace std;

#ifndef PI
#define PI 3.1415926535897932385
#endif

typedef unsigned short ushort;

typedef ushort sino_type;
typedef float img_type;
typedef float edge_type;

class Parameter {
public:
    int NX, NY, NZ;
    int NPROJ, NDX, NDY;
    int NDY_THICK, NDY_OFFSET;
    float SAMPLESIZE, PIXELSIZE;

    int MAX_RAYLEN;

    int ITERATIONS;

    float SOD;
    float SDD;

    // [-1024, 1024] detectors
    //float LENGTH_PER_DET;

    static constexpr float ALPHA = 0;
    static constexpr float BETA = 0;
    static constexpr float EPSILON = 0.001;

    int THREAD_NUMB;

    float vx,vy,vz;

    int HALFDET;

    float HALFSIZE;

    string BEAM;

    string RAW_DATA_FILE;

    string OUTPUT_DIR;

    void parse_config(int argc, char** argv);
    
    void print_options();

    void derive();
};

// the output of the program
class CTOutput {
public:
    // pointer to the reconstructed image
    img_type *img;
    edge_type *edge;
    Parameter args;

    void allocate_img(int64_t size);
    void allocate_edge(int64_t size);

    img_type& img_data(int z,int x,int y);
    edge_type& edge_data(int z,int x,int y);

    void write_img(const string &output_dir);
    void write_edge(const string &output_dir);

    float minIMAGE(float Af, int64_t *line, float *weight, int numb, float lambda);
    void minEDGE(int64_t *line, float *weight, int numb, float lambda);

    CTOutput(const Parameter &_args);
    ~CTOutput();
};

// the input of the program
class CTInput {
public:
    // pointer to the sinogram
    sino_type *sino;
    Parameter args;

    void allocate_sino(int64_t size);

    sino_type sino_data(int p,int x,int y) const;

    void read_sino(const string &raw_data_file);

    CTInput(const Parameter &_args);
    ~CTInput();
};