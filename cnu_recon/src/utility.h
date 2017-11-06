#pragma once

#include "main.h"

#include <string>

using namespace std;

float array_3d_img(float *data,int z,int x,int y);

ushort array_3d_sino(ushort *data,int z,int x,int y);

//int64_t timer_us();

//int64_t timer_s();

void write_data_3d(float *data,int X,int Y,int Z,string output_filename);

float sqr(float x);

//float get_img_addr(float x,float y,float z);
int64 get_img_addr(int64 x,int64 y,int64 z);
