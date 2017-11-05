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

extern int NX, NY, NZ;
extern int NANGLE, NDETECTORX, NDETECTORZ;
extern int NPROJ_TURN;

extern int INIT_ANGLE;

#define MAX_ELE_RAY 1280

extern int CT_ITERATIONS;

extern double SOURCE_TO_ISO;
extern double SOURCE_TO_DET;

// [-1024, 1024] detectors
extern double LENGTH_PER_DET;

#define HALFDET (NDETECTORX / 2)
#define HALFSIZE (NX / 2.0)

#define vx 1.00
#define vy 1.00
// slice thickness
extern double vz;

const float ALPHA = 0;
const float BETA = 0;
const float EPSILON = 0.001;

extern int THREAD_NUMB;

const int INTERVAL = 16;
#ifndef PI
#define PI 3.1415926535897932385
#endif
