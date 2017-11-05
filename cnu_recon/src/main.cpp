#include "config.h"
#include "ct3d.h"
#include "utility.h"
#include "main.h"
#include <cstdio>

#include <boost/foreach.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>


using namespace std;
using namespace boost::property_tree;

string table_filename;
string sino_filename;
string image_filename;
string edge_filename;
string init_f_filename;
string init_v_filename;

string OUTPUT_DIR;

int INIT_ANGLE;

int NX, NY, NZ;
int NANGLE, NDETECTORX, NDETECTORZ;
int CT_ITERATIONS;

double vz;

double SOURCE_TO_ISO;
double SOURCE_TO_DET;

double LENGTH_PER_DET;

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
		table_filename = root.get<string>("table_filename");
		sino_filename = root.get<string>("sino_filename");
		image_filename = root.get<string>("image_filename");
//		edge_filename = root.get<string>("edge_filename");
		init_f_filename = root.get<string>("init_f_filename");
		init_v_filename = root.get<string>("init_v_filename");

		OUTPUT_DIR = root.get<string>("OUTPUT_DIR");

		INIT_ANGLE = root.get<int>("INIT_ANGLE");

		NX = root.get<int>("NX");
		NY = root.get<int>("NY");
		NZ = root.get<int>("NZ");

		NDETECTORX = root.get<int>("NDETECTORX");
		NDETECTORZ = root.get<int>("NDETECTORZ");
		NANGLE = root.get<int>("NANGLE");
		CT_ITERATIONS = root.get<int>("CT_ITERATIONS");

		SOURCE_TO_ISO = root.get<double>("SOURCE_TO_ISO");
		SOURCE_TO_DET = root.get<double>("SOURCE_TO_DET");

		LENGTH_PER_DET = root.get<double>("LENGTH_PER_DET");

        vz = root.get<double>("vz");

		THREAD_NUMB = root.get<int>("THREAD_NUMB");
	}
	catch (ptree_error &e) {
		printf("JSON file corrupted.\n");
		exit(1);
	}
}

void print_options() {
	printf("CT3D reconstruction flow begins!\n");
	printf("================================\n");
	printf("Here are the options:\n");

	printf("Image size   NZ=%d   NX=%d   NY=%d\n",NZ,NX,NY);
	printf("Number of angles: %d\n",NANGLE);
	printf("Number of detector row and channel: %d, %d\n",NDETECTORX,NDETECTORZ);

	printf("%d iterations with %d threads.\n",CT_ITERATIONS, THREAD_NUMB);

	printf("Source to ISO and detectors: %.1lf, %.1lf\n",SOURCE_TO_ISO,SOURCE_TO_DET);

	printf("Length per detector length: %lf\n",LENGTH_PER_DET);

	printf("Results in the directory %s\n",OUTPUT_DIR.c_str());

	fflush(stdout);
}

int main(int argc, char** argv) {

	parse_config(argc, argv);

	print_options();

	double *image_data = new double[NZ*NX*NY];
	double *edge_data = new double[NZ*NX*NY];
	double *sino_data = new double[NANGLE*NDETECTORX*NDETECTORZ];

	memset(image_data, 0, sizeof(NZ*NX*NY)*sizeof(double));
//	memset(edge_data, 0, sizeof(NZ*NX*NY)*sizeof(double));
	memset(sino_data, 0, sizeof(NANGLE*NDETECTORX*NDETECTORZ)*sizeof(double));

	ct3d(image_data,/* edge_data,*/ sino_data);

	write_data_3d(image_data, NZ, NX, NY, image_filename);
	write_data_3d(edge_data, NZ, NX, NY, edge_filename);

	delete [] image_data;
	delete [] edge_data;
	delete [] sino_data;

	return 0;
}
