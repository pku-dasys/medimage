#pragma once
#include<string>
using namespace std;

double array_2d_raw (double *data, int SIZE0, int SIZE1,int x, int y);
double array_2d_sino (double *data, int x, int y);
double array_2d_img (double *data, int x, int y);
double sqr(double x);
double get_img_addr(double x,double y);

void write_data_2d(double *data,int X,int Y,string output_filename);
double array_3d_raw(double *data,int SIZE0,int SIZE1,int SIZE2,int x,int y,int z);
void write_data_3d(double *data,int X,int Y,int Z,string output_filename);
