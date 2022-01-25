#pragma once

#include <string>

using namespace std;

double array_3d_img(double *data,int z,int x,int y);

double array_3d_sino(double *data,int z,int x,int y);

int64_t timer_us();

int64_t timer_s();

void write_data_3d(double *data,int X,int Y,int Z,string output_filename);

double sqr(double x);

double get_img_addr(double x,double y,double z);