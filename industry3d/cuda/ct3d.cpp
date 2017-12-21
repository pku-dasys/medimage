#include "ct3d.h"
#include "utility.h"

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
        //PRETRACING_FILE = root.get<string>("PRETRACING_FILE");

        OUTPUT_DIR = root.get<string>("OUTPUT_DIR");

        NX = root.get<int>("NX");
        NY = root.get<int>("NY");
        NZ = root.get<int>("NZ");

        ALPHA = root.get<float>("ALPHA");
        BETA = root.get<float>("BETA");
        EPSILON = root.get<float>("EPSILON");

        LAMBDA_IMG = root.get<float>("LAMBDA_IMG");
        LAMBDA_EDGE = root.get<float>("LAMBDA_EDGE");
        AMPLIFIER = root.get<int>("AMPLIFIER");

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

void CTOutput::allocate() {
    int64_t size = (args.NX+2)*(args.NY+2)*(args.NZ+2);
    img = new img_type[size];
    edge = new edge_type[size];
    for (int64_t i = 0; i<size; ++i) {
        img[i] = 0.;
        edge[i] = 1.0;
    }
}

img_type& CTOutput::img_data(int z,int x,int y) {
    return img[(z+1)*(args.NX+2)*(args.NY+2)+(x+1)*(args.NY+2)+(y+1)];
}

edge_type& CTOutput::edge_data(int z,int x,int y) {
    return edge[(z+1)*(args.NX+2)*(args.NY+2)+(x+1)*(args.NY+2)+(y+1)];
}

void CTOutput::write_img(const string &output_dir,int iteration) {
    cout<< "Start writing images ... ";
    boost::filesystem::path dir(output_dir);
    if (!boost::filesystem::exists(dir))
        boost::filesystem::create_directory(dir);
    for (int z = 0; z<args.NZ; ++z) {
        string slice_output;
        if (iteration==-1)
            slice_output = output_dir+"/i_"+boost::lexical_cast<string>(z);
        else
            slice_output = output_dir+"/"+boost::lexical_cast<string>(iteration)+"_i_"+boost::lexical_cast<string>(z);
        if (iteration!=-1 && z!=args.NZ/2) continue;
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
    cout<< "done." <<endl;
}

void CTOutput::write_edge(const string &output_dir,int iteration) {
    cout<< "Start writing edges ... ";
    boost::filesystem::path dir(output_dir);
    if (!boost::filesystem::exists(dir))
        boost::filesystem::create_directory(dir);
    for (int z = 0; z<args.NZ; ++z) {
        string slice_output;
        if (iteration==-1)
            slice_output = output_dir+"/e_"+boost::lexical_cast<string>(z);
        else
            slice_output = output_dir+"/"+boost::lexical_cast<string>(iteration)+"_e_"+boost::lexical_cast<string>(z);
        if (iteration!=-1 && z!=args.NZ/2) continue;
        ofstream fou(slice_output);
        fou.precision(6);
        for (int x = 0; x<args.NX; ++x) {
            for (int y = 0; y<args.NY; ++y) {
                fou<< edge_data(z, x, y) << ' ';
            }
            fou<<endl;
        }
        fou.close();
    }
    cout<< "done." <<endl;
}

CTOutput::CTOutput(const Parameter &_args) : args(_args) {
    img = nullptr;
    edge = nullptr;
}
CTOutput::~CTOutput() {
    if (img) delete [] img;
    if (edge) delete [] edge;
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

void CTInput::read_pretracing(const string &pretracing_file) {
    ifstream fin(pretracing_file, ios::in);

    back_front = new int[args.NZ];
    back_rear = new int[args.NZ];

    for (int i = 0; i<args.NZ; ++i)
        fin >> back_front[i] >> back_rear[i];

    fin.close();

}

CTInput::CTInput(const Parameter &_args) : args(_args) {
    sino = nullptr;
    back_front = nullptr;
    back_rear = nullptr;
}
CTInput::~CTInput() {
    if (sino) delete [] sino;
    if (back_front) delete [] back_front;
    if (back_rear) delete [] back_rear;
}
