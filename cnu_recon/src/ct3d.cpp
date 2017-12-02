#include "ct3d.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include <boost/lexical_cast.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>

using namespace std;

/////////// Class  Parameter

void Parameter::parse_config(int argc, char** argv) {
    string config_filename;
    if (argc==2) {
        config_filename = string(argv[1]);
    }
    else {
        printf("Usage: ./ct3d CONFIG_FILE_NAME\n");
        exit(1);
    }

    boost::property_tree::ptree root;
    try {
        read_json(config_filename, root);
    }
    catch (boost::property_tree::ptree_error &e) {
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
        SAMPLESIZE = root.get<float>("SAMPLESIZE");
        PIXELSIZE = root.get<float>("PIXELSIZE");

        BEAM = root.get<string>("BEAM");
    }
    catch (boost::property_tree::ptree_error &e) {
        printf("JSON file corrupted.\n");
        exit(1);
    }
    derive();
}

void Parameter::derive() {
    // derive other parameters
    HALFDET = NDX/2;
    HALFSIZE = NX/2.0;
    MAX_RAYLEN = NX*3;
}

void Parameter::print_options() {
    cout << "CT3D reconstruction flow begins!" << endl;
    cout << "================================" << endl;
    cout << "Here are the options:" << endl;
    cout << "Image size: NZ=" << NZ << " NX=" << NX <<" NY="<< NY<< endl;
    cout << "Number of angles: "<< NPROJ <<endl;
    cout << "Number of detector row and channel: " << NDX << ", " << NDY << endl;
    cout << ITERATIONS << " iterations with " << THREAD_NUMB << " threads." << endl;
    cout << "Source to ISO and detectors: " << SOD << ", " << SDD << endl;
//    cout << "Length per detector: " << LENGTH_PER_DET << endl;
    cout << "Sample Size: " << SAMPLESIZE << endl;
    cout << "Detector Pixel Size: " << PIXELSIZE << endl;
    cout<<endl;
}

/////////// Class  CTOutput

void CTOutput::allocate_img(int size) {
    img = new img_type[size]{};
}

img_type& CTOutput::img_data(int z,int x,int y) {
    return img[z*args.NX*args.NY+x*args.NY+y];
}

void CTOutput::write_img(const string &output_dir) {
    cout<< "Start writing images ..." <<endl;
    boost::filesystem::path dir(output_dir);
    if (!boost::filesystem::exists(dir))
    boost::filesystem::create_directory(dir);
    for (int z = 0; z<args.NZ; ++z) {
        string slice_output = output_dir+"/"+boost::lexical_cast<string>(z);
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

void CTOutput::minIMAGE(float Af, int64_t *line, float *weight, int numb, float lambda) {
    for (int i = 0; i<numb; i++) {
        Af += img[line[i]] * weight[i];
    }
    for (int i = 0; i<numb; i++) {
        int64_t ind = line[i];
        img_type tmp = img[ind] + lambda * (-Af * weight[i]);
        if (tmp<0) tmp = 0;
        img[ind] = tmp;
    }
}

CTOutput::CTOutput(const Parameter &_args) : args(_args) {
}
CTOutput::~CTOutput() {
    if (img) delete [] img;
}

/////////// Class  CTInput

void CTInput::allocate_sino(int64_t size) {
    sino = new sino_type[size]{};
}

sino_type CTInput::sino_data(int p,int x,int y) const {
    return sino[p*args.NDX*args.NDY+x*args.NDY+y];
}

void CTInput::read_sino(const string &raw_data_file) {
    ifstream fin(raw_data_file, ios::in | ios::binary);
    fin.seekg(1024);
    int64_t size = (int64_t)args.NPROJ * args.NDX * args.NDY;
    allocate_sino(size);
    //cout << "read sinogram from binary file ..." << endl;
    fin.read((char*)sino, sizeof(sino_type) * size);
    //cout << "reading finished." << endl;
    fin.close();
}

CTInput::CTInput(const Parameter &_args) : args(_args) {
}
CTInput::~CTInput() {
    if (sino) delete [] sino;
}
