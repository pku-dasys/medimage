#pragma once

#include <string>

using namespace std;

extern string table_filename;
extern string sino_filename;
extern string image_filename;
extern string edge_filename;
extern string init_f_filename;
extern string init_v_filename;
extern string rebined_folder;

extern string OUTPUT_DIR;

extern int NX;
extern int NY;
extern int NZ;
extern int NDETECTOR;

extern int NPROJ;


#define MAX_ELE_RAY 2560

extern int MS_ITERATIONS;

extern double SOURCE_TO_ISO;
extern double SOURCE_TO_DET;

extern double OFF_CENTER;

// [-25, 25] in 736 detectors
extern double LENGTH_PER_DET;

//#define vx 1.0
//#define vy 1.0
// slice thickness
extern double vx;
extern double vy;

extern int thread_offset;

#define HALFDET (NDETECTOR/2)
#define HALFSIZE (NX / 2.0)


extern float ALPHA ;
extern float BETA ;
extern float EPSILON;
extern int THREAD_NUMB;
const int INTERVAL = 16;
