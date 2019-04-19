#include "main.h"
#include "utility.h"
#include "tracing.h"

#include <cmath>

using namespace std;

#define LAMBDA_X(i, x_s, x_d, L) (L*((double)i-x_s)/(x_d-x_s))
#define LAMBDA_Y(j, y_s, y_d, L) (L*((double)j-y_s)/(y_d-y_s))
#define LAMBDA_Z(k, z_s, z_d, L) (L*((double)k-z_s)/(z_d-z_s))

#ifndef MAX
#define MAX( x, y ) ( ((x) > (y)) ? (x) : (y) )
#endif
#ifndef MIN
#define MIN( x, y ) ( ((x) < (y)) ? (x) : (y) )
#endif

#define ABS_VALUE(x) ( (x < 0) ? -(x) : (x) )

void forward_proj(double sx,double sy,double dx,double dy,int *ind, double *wgt, int &numb){

  int index =0;
  double ray_x, ray_y;
  double len_x, len_y;
  double absvalue_x, absvalue_y;
  double lambda_x, lambda_y;

  double L;
  double lambda_min = 0.0;
  double lambda_max;
  double lambda0, lambdaN;
  double temp;

  int signx, signy;
  int v_x, v_y;

  ray_x = dx - sx;
  ray_y = dy - sy;
  //ray vector

  L = sqrt( sqr(ray_x)+ sqr(ray_y));
  lambda_max = L;

  //the the direction of increment in x, y and z
  signx = (ray_x > 0) ? 1 : -1;
  signy = (ray_y > 0) ? 1 : -1;

  //increment in x, y and z direction
  absvalue_x = fabs(ray_x);
  absvalue_y = fabs(ray_y);

  //get x=1 Lx Ly Lz
  len_x = (absvalue_x > 1.e-4) ? (L / absvalue_x) * vx : 1.e6;
  len_y = (absvalue_y > 1.e-4) ? (L / absvalue_y) * vy : 1.e6;

  //initialize the values
  numb = 0;

  //get the entry and exit point between Ray & Image
  //distance between source and entry point
  double tempx, tempy;

  tempx = NX * vx;
  tempy = NY * vy;

  lambda0 = LAMBDA_X(0, sx, dx, L);
  lambdaN = LAMBDA_X(tempx, sx, dx, L);
  temp    = MIN(lambda0, lambdaN);
  lambda_min = MAX(lambda_min, temp);
  if (lambda_min == temp)
      index = 1;
  temp    = MAX(lambda0, lambdaN);
  lambda_max = MIN(lambda_max, temp);  // start x plane

  lambda0 = LAMBDA_Y(0, sy, dy, L);
  lambdaN = LAMBDA_Y(tempy, sy, dy, L);
  temp    = MIN(lambda0, lambdaN);
  lambda_min = MAX(lambda_min, temp);
  if (lambda_min == temp)
      index = 2;
  temp    = MAX(lambda0, lambdaN);
  lambda_max = MIN(lambda_max, temp);   // start y plane

  if (lambda_min >= lambda_max)
      return;

  lambda0 = lambda_min;   // lambda = lambda_min

  if (index == 1)
  {
      if (signx == 1)
          v_x = 0;
      else
          v_x = NX - 1;
      lambda_x = lambda0 + len_x;

      v_y = (sy + lambda0 * ray_y / L) / vy;
      tempy = v_y * vy;
      lambda_y = (absvalue_y < 1.e-4) ? 1.e6 : LAMBDA_Y(tempy + (signy > 0) * vy, sy, dy, L);
  }

  else
  {
      if (signy == 1)
          v_y = 0;
      else
          v_y = NY - 1;
      lambda_y = lambda0 + len_y;

      v_x = (sx + lambda0 * ray_x / L) / vx;
      tempx = v_x * vx;
      lambda_x = (absvalue_x < 1.e-4) ? 1.e6 : LAMBDA_X(tempx + (signx > 0) * vx, sx, dx, L);
  }


  while (lambda0 < lambda_max - 5.e-2)
  {
      if (lambda_x <= lambda_y)
      {
          //(*sino)  += (lambda_x - lambda0) * Data(imageDataPtr, NX, NY, NZ, v_x, v_y, v_z);

          ind[numb] = get_img_addr(v_x, v_y);
          wgt[numb] = lambda_x - lambda0;
          //Af += f[ind[numb]]*wgt[numb];
          ++numb;

          lambda0   = lambda_x;
          lambda_x += len_x;
          v_x      += signx;
      }
      else
      {
          //(*sino)  += (lambda_y - lambda0) * Data(imageDataPtr, NX, NY, NZ, v_x, v_y, v_z);

          ind[numb] = get_img_addr(v_x, v_y);
          wgt[numb] = lambda_y - lambda0;
          //Af += f[ind[numb]]*wgt[numb];
          ++numb;

          lambda0   = lambda_y;
          lambda_y += len_y;
          v_y      += signy;
      }
  }
}
