#ifndef _CONFIG_H_
#define _CONFIG_H_

#define STARTZ 0

#define NPROJ 14524
#define NDETECTOR 736
#define NCHANNEL 16
#define MAX_ELE_RAY 2048
#define NX 874
#define NY 874
#define NZ 32

#ifdef __SDSCC__
#include <ap_fixed.h>

typedef unsigned int Addr;
typedef ap_fixed<16,8> Data;
typedef Addr Index;
typedef ap_int<16> Number;
typedef ap_ufixed<8,2> Weight;
#else
typedef int Addr;
typedef float Data;
typedef Addr Index;
typedef int Number;
typedef float Weight;
#endif

const Data ALPHA = 0.1;
const Data BETA = 0.01;
const Data EPSILON = 0.001;

#endif
