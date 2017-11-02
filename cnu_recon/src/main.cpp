#include "config.h"
#include "ct3d.h"
#include <cstdio>

using namespace std;

string table_filename;
string sino_filename;
string image_filename;
string edge_filename;
string init_f_filename;
string init_v_filename;

int PROJ_OFFSET;
int INIT_PROJ;

int NX, NY, NZ;
int NPROJ, NDETECTOR, NCHANNEL;
int NPROJ_TURN;
int MS_ITERATIONS;

double vz;

double SOURCE_TO_ISO;
double SOURCE_TO_DET;

double OFF_CENTER;

double LENGTH_PER_DET;
double LENGTH_PER_DET_Z;

double PITCH_VALUE;

int THREAD_NUMB;

void parse_config(int argc, char** argv) {
	string config_filename;
	if (argc==2) {
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
		table_filename = root.get<string>("table_filename");
		sino_filename = root.get<string>("sino_filename");
		image_filename = root.get<string>("image_filename");
		edge_filename = root.get<string>("edge_filename");
		init_f_filename = root.get<string>("init_f_filename");
		init_v_filename = root.get<string>("init_v_filename");

		OUTPUT_DIR = root.get<string>("OUTPUT_DIR");

		PROJ_OFFSET = root.get<int>("PROJ_OFFSET");
		INIT_PROJ = root.get<int>("INIT_PROJ");

		NX = root.get<int>("NX");
		NY = root.get<int>("NY");

		NDETECTOR = root.get<int>("NDETECTOR");
		NCHANNEL = root.get<int>("NCHANNEL");
		NPROJ = root.get<int>("NPROJ");
		NPROJ_TURN = root.get<int>("NPROJ_TURN");
		MS_ITERATIONS = root.get<int>("MS_ITERATIONS");

		SOURCE_TO_ISO = root.get<double>("SOURCE_TO_ISO");
		SOURCE_TO_DET = root.get<double>("SOURCE_TO_DET");

		OFF_CENTER = root.get<double>("OFF_CENTER");
		LENGTH_PER_DET = root.get<double>("LENGTH_PER_DET");
		LENGTH_PER_DET_Z = root.get<double>("LENGTH_PER_DET_Z");

		PITCH_VALUE = root.get<double>("PITCH_VALUE");

		vz = root.get<double>("vz");

		THREAD_NUMB = root.get<int>("THREAD_NUMB");
	}
	catch (ptree_error &e) {
		printf("JSON file corrupted.\n");
		exit(1);
	}

	NZ = (NPROJ*1.0/NPROJ_TURN)*PITCH_VALUE/vz;
}

void print_options() {
	printf("MS3D reconstruction flow begins!\n");
	printf("================================\n");
	printf("Here are the options:\n");

	printf("Image size   NZ=%d   NX=%d   NY=%d\n",NZ,NX,NY);
	printf("Number of projections in all and one cycle: %d, %d\n",NPROJ,NPROJ_TURN);
	printf("Number of detector row and channel: %d, %d\n",NDETECTOR,NCHANNEL);
	printf("Slice thickness: %lf\n",vz);

	printf("%d iterations with %d threads.\n",MS_ITERATIONS, THREAD_NUMB);

	printf("Source to ISO and detectors: %.1lf, %.1lf\n",SOURCE_TO_ISO,SOURCE_TO_DET);

	printf("Length per detector and z-axis length: %lf, %lf\n",LENGTH_PER_DET, LENGTH_PER_DET_Z);

	printf("Results in the directory %s\n",OUTPUT_DIR.c_str());

	fflush(stdout);
}

int main(int argc, char** argv) {

	parse_config(argc, argv);

	print_options();

	double *image_data = new double[NZ*NX*NY];
	double *edge_data = new double[NZ*NX*NY];
	double *sino_data = new double[NPROJ*NDETECTOR*NCHANNEL];

	memset(image_data, 0, sizeof(NZ*NX*NY)*sizeof(double));
	memset(edge_data, 0, sizeof(NZ*NX*NY)*sizeof(double));
	memset(sino_data, 0, sizeof(NPROJ*NDETECTOR*NCHANNEL)*sizeof(double));

	ct3d(image_data, edge_data, sino_data);

	write_data_3d(image_data, NZ, NX, NY, image_filename);
	write_data_3d(edge_data, NZ, NX, NY, edge_filename);

	delete [] image_data;
	delete [] edge_data;
	delete [] sino_data;

	return 0;
}
