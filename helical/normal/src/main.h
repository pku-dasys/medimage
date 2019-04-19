#pragma once

#include <string>

using namespace std;

extern string table_filename;
extern string sino_filename;
extern string image_filename;
extern string edge_filename;
extern string init_f_filename;
extern string init_v_filename;

extern string OUTPUT_DIR;

extern int NX;
extern int NY;
extern int NZ;
extern int NDETECTOR;
extern int NCHANNEL;
extern int NPROJ;
extern int NPROJ_TURN;

extern int PROJ_OFFSET;
extern int INIT_PROJ;

#define MAX_ELE_RAY 2560

extern int MS_ITERATIONS;

extern double SOURCE_TO_ISO;
extern double SOURCE_TO_DET;

extern double OFF_CENTER;

extern double PITCH_VALUE;

// [-25, 25] in 736 detectors
extern double LENGTH_PER_DET;
extern double LENGTH_PER_DET_Z;

//#define vx 0.5859
//#define vy 0.5859
// slice thickness
extern double vz;
extern double vx;
extern double vy;

#define HALFDET (NDETECTOR/2)
#define HALFSIZE (NX / 2.0)
#define HALFDETZ (NCHANNEL/2)

extern float ALPHA ;
extern float BETA ;
extern float EPSILON;

extern float *alpha;
extern float *beta;

extern int THREAD_NUMB;

const int INTERVAL = 16;
