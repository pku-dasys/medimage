#pragma once

#include <string>
#include <cmath>
#include <fstream>

#include <boost/lexical_cast.hpp>

using namespace std;

#ifndef PI
#define PI 3.1415926535897932385
#endif

typedef unsigned short ushort;
typedef long long int64;

typedef ushort sino_type;
typedef float img_type;

// the output of the program
class CTOutput {
public:
    // pointer to the reconstructed image
    img_type *img;

    void allocate_img(int size) {
        img = new img_type[size]{};
    }

    img_type& img_data(int Z,int X,int Y,int z,int x,int y) {
        return img[z*X*Y+x*Y+y];
    }

    void write_img(const Parameter &args, const string &output_dir) {
        for (int z = 0; z<args.NZ; ++z) {
            string slice_output = output_dir+"/"+output_filename+lexical_cast<string>(z);
            ofstream fou(slice_output);
            fou.precision(6);
            for (int x = 0; x<args.NX; ++x) {
                for (int y = 0; y<args.NY; ++y) {
                    fou<< img_data(args.NZ, args.NX, args.NY, z, x, y);
                }
                fout<<endl;
            }
            fou.close();
        }
    }

    CTOutput() {}
    ~CTOutput() {
        if (img) delete [] img;
    }
};

struct Parameter {
    int NX, NY, NZ;
    int NPROJ, NDX, NDY;
    int NDY_THICK, NDY_OFFSET;

    const int MAX_RAYLEN = 40960;

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

    void derive() {
        HALFDET = NDETECTORX/2;
        HALFSIZE = NX/2.0;
    }
};

// the input of the program
class CTInput {
public:
    // pointer to the sinogram
    sino_type *sino;


    void allocate_sino(int64 size) {
        sino = new sino_type[size]{};
    }

    void read_sino(const Parameter &args, const string &raw_data_file) {
        ifstream fin(raw_data_file, ios::in | ios::binary);
        fin.seekg(1024);
        int64 size = (int64)args.NPROJ * args.NDX * args.NDY;
        allocate_sino(size);
        cout << "read sinogram from binary file ..." << endl;
        fin.read((char*)sino, sizeof(sino_type) * size);
        cout << "reading finished." << endl;
        fin.close()
    }

    CTInput() {}    
    ~CTInput() {
        if (sino) delete [] sino;
    }
};