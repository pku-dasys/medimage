#include "ct3d.h"
#include "utility.h"
#include "main.h"

#include <iostream>

#include <boost/foreach.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

using namespace std;
using namespace boost::property_tree;

string RAW_DATA_FILE;

string OUTPUT_DIR;

void parse_config(Parameter &args, int argc, char** argv) {
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

        args.NX = root.get<int>("NX");
        args.NY = root.get<int>("NY");
        args.NZ = root.get<int>("NZ");
        
        args.NDX = root.get<int>("NDX");
        args.NDY = root.get<int>("NDY");
        args.NPROJ = root.get<int>("NPROJ");
        args.NDY_THICK = root.get<int>("NDY_THICK");
        args.NDY_OFFSET = root.get<int>("NDY_OFFSET");
        
        args.ITERATIONS = root.get<int>("ITERATIONS");

        args.SOD = root.get<float>("SOD");
        args.SDD = root.get<float>("SDD");

        args.THREAD_NUMB = root.get<int>("THREAD_NUMB");

        args.BEAM = root.get<string>("BEAM");
    }
    catch (ptree_error &e) {
        printf("JSON file corrupted.\n");
        exit(1);
    }

    // derive other parameters
    args.derive();
}

void print_options(const Parameter &args) {
    cout << "CT3D reconstruction flow begins!" << endl;
    cout << "================================" << endl;
    cout << "Here are the options:" << endl;
    cout << "Image size: NZ=" << args.NZ << " NX=" << args.NX <<" NY="<<args.NY<< endl;
    cout << "Number of angles: "<< args.NPROJ <<endl;
    cout << "Number of detector row and channel: " << args.NDX << ", " << args.NDY << endl;
    cout << args.ITERATIONS << " iterations with " << args.THREAD_NUMB << " threads." << endl;
    cout << "Source to ISO and detectors: " << args.SOD << ", " << args.SDD << "%.1lf" << endl;
    cout << "Length per detector: " << args.LENGTH_PER_DET << endl;
    cout<<endl;
}

int main(int argc, char** argv) {
    Parameter args;

    parse_config(args, argc, argv);
    print_options(args);

    CTInput in = CTInput(args);
    CTOutput out = CTOutput(args);

    in.read_sino(RAW_DATA_FILE);
    ct3d(args,in,out);
    out.write_img(OUTPUT_DIR);

    return 0;
}
