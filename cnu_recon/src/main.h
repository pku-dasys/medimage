#pragma once

#include <string>
#include <cmath>
#include <iostream>
#include <fstream>

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

using namespace std;
using namespace boost;

#ifndef PI
#define PI 3.1415926535897932385
#endif

typedef unsigned short ushort;
typedef long long int64;

typedef ushort sino_type;
typedef float img_type;

struct Parameter {
    int NX, NY, NZ;
    int NPROJ, NDX, NDY;
    int NDY_THICK, NDY_OFFSET;

    const int MAX_RAYLEN = 4096;

    int ITERATIONS;

    float SOD;
    float SDD;

    // [-1024, 1024] detectors
    float LENGTH_PER_DET;

    const float ALPHA = 0;
    const float BETA = 0;
    const float EPSILON = 0.001;

    int THREAD_NUMB;

    float vx,vy,vz;

    int HALFDET;

    float HALFSIZE;

    string BEAM;

    void derive() {
        HALFDET = NDX/2;
        HALFSIZE = NX/2.0;
    }
};

// the output of the program
class CTOutput {
public:
    // pointer to the reconstructed image
    img_type *img;
    Parameter args;

    void allocate_img(int size) {
        img = new img_type[size]{};
    }

    img_type& img_data(int z,int x,int y) {
        return img[z*args.NX*args.NY+x*args.NY+y];
    }

    void write_img(const string &output_dir) {
        cout<< "Start writing images ..." <<endl;
        for (int z = 0; z<args.NZ; ++z) {
            string slice_output = output_dir+"/"+lexical_cast<string>(z);
            ofstream fou(slice_output);
            fou.precision(6);
            for (int x = 0; x<args.NX; ++x) {
                for (int y = 0; y<args.NY; ++y) {
                    fou<< img_data(z, x, y) << ' ';
                }
                fou<<endl;
            }
            fou.close();
        }
    }

    void minIMAGE(float Af, int64 *line, float *weight, int numb, float lambda) {
        for (int i = 0; i<numb; i++) {
            Af += img[line[i]] * weight[i];
        }
        for (int i = 0; i<numb; i++) {
            int64 ind = line[i];
            img_type tmp = img[ind] + lambda * (-Af * weight[i]);
            if (tmp<0) tmp = 0;
            img[ind] = tmp;
        }
    }

    CTOutput(const Parameter &_args) : args(_args) {
    }
    ~CTOutput() {
        if (img) delete [] img;
    }
};

// the input of the program
class CTInput {
public:
    // pointer to the sinogram
    sino_type *sino;
    Parameter args;

    void allocate_sino(int64 size) {
        sino = new sino_type[size]{};
    }

    sino_type sino_data(int p,int x,int y) const {
        return sino[p*args.NDX*args.NDY+x*args.NDY+y];
    }

    void read_sino(const string &raw_data_file) {
        ifstream fin(raw_data_file, ios::in | ios::binary);
        fin.seekg(1024);
        int64 size = (int64)args.NPROJ * args.NDX * args.NDY;
        allocate_sino(size);
        //cout << "read sinogram from binary file ..." << endl;
        fin.read((char*)sino, sizeof(sino_type) * size);
        //cout << "reading finished." << endl;
        fin.close();
    }

    CTInput(const Parameter &_args) : args(_args) {
    }
    ~CTInput() {
        if (sino) delete [] sino;
    }
};