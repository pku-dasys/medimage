#pragma once

#include <string>

using namespace std;

int64_t timer_us();

int64_t timer_s();

float sqr(float x);

void rotate_2d(float &x,float &y,float theta);
void rotate_axis_3d(float &x,float &y,float &z,float ux,float uy,float uz,float theta);

bool equals_draw(float x,int y);