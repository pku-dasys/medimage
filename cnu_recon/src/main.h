#pragma once

#include <string>
#include <cmath>

using namespace std;

extern string RAW_DATA_FILE;

extern string OUTPUT_DIR;

extern int NX, NY, NZ;
extern int NPROJ, NDX, NDY;
extern int NDY_THICK, NDY_OFFSET;

#define MAX_RAYLEN 40960

extern int ITERATIONS;

extern float SOD;
extern float SDD;

// [-1024, 1024] detectors
extern float LENGTH_PER_DET;

#define HALFDET (NDETECTORX / 2)
#define HALFSIZE (NX / 2.0)

//#define vx 1.00
//#define vy 1.00
//#define vz 1.00
#define vx 0.200
#define vy 0.200
#define vz 0.200

const float ALPHA = 0;
const float BETA = 0;
const float EPSILON = 0.001;

extern int THREAD_NUMB;

const int INTERVAL = 16;

#ifndef PI
#define PI 3.1415926535897932385
#endif


typedef unsigned short ushort;
typedef long long int64;

