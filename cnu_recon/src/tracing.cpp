#include "main.h"
#include "tracing.h"
#include "utility.h"

#include <cmath>

using namespace std;

#define LAMBDA_X(i, x_s, x_d, L) (L*((float)i-x_s)/(x_d-x_s))
#define LAMBDA_Y(j, y_s, y_d, L) (L*((float)j-y_s)/(y_d-y_s))
#define LAMBDA_Z(k, z_s, z_d, L) (L*((float)k-z_s)/(z_d-z_s))

#ifndef MAX
#define MAX( x, y ) ( ((x) > (y)) ? (x) : (y) )
#endif
#ifndef MIN
#define MIN( x, y ) ( ((x) < (y)) ? (x) : (y) )
#endif
#ifndef MAX3
#define MAX3(x,y,z) MAX(MAX(x,y),z)
#endif
#ifndef MIN3
#define MIN3(x,y,z) MIN(MIN(x,y),z)
#endif
#define ABS_VALUE(x) ( (x < 0) ? -(x) : (x) )

int64 get_img_addr(const Parameter &args,int x,int y,int z) {
    return (int64)z*args.NX*args.NY+x*args.NY+y;
}

void forward_proj(const Parameter &args,
                  float sx,float sy,float sz,
                  float dx,float dy,float dz,
                  int64 *ind,float *wgt,int &numb) {

    int index = 0;

    float ray_x, ray_y, ray_z;
    float len_x, len_y, len_z;
    float absvalue_x, absvalue_y, absvalue_z;
    float lambda_x, lambda_y, lambda_z;
    float L;
    float lambda_min = 0.0;
    float lambda_max;
    float lambda0, lambdaN;
    float temp;
    int signx, signy, signz;
    int v_x, v_y, v_z;
    
    // ray vector
    ray_x = dx - sx;
    ray_y = dy - sy;
    ray_z = dz - sz;

    // distance
    L = sqrt( sqr(ray_x)+sqr(ray_y)+sqr(ray_z) );
    lambda_max = L;
    
    //the the direction of increment in x, y and z
    signx = (ray_x > 0) ? 1 : -1;
    signy = (ray_y > 0) ? 1 : -1;
    signz = (ray_z > 0) ? 1 : -1;

    //increment in x, y and z direction
    absvalue_x = fabs(ray_x);
    absvalue_y = fabs(ray_y);
    absvalue_z = fabs(ray_z);
    
    //get x=1 Lx Ly Lz
    len_x = (absvalue_x > 1.e-4) ? (L / absvalue_x) : 1.e6;
    len_y = (absvalue_y > 1.e-4) ? (L / absvalue_y) : 1.e6;
    len_z = (absvalue_z > 1.e-4) ? (L / absvalue_z) : 1.e6;

    //initialize the values
    numb = 0;

    //get the entry and exit point between Ray & Image
    //distance between source and entry point
    float tempx, tempy, tempz;

    //tempx = NX * vx;
    //tempy = NY * vy;
    //tempz = NZ * vz;

    tempx = args.NX;
    tempy = args.NY;
    tempz = args.NZ;
    
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

    lambda0 = LAMBDA_Z(0, sz, dz, L);
    lambdaN = LAMBDA_Z(tempz, sz, dz, L);
    temp    = MIN(lambda0, lambdaN);
    lambda_min = MAX(lambda_min, temp);
    if (lambda_min == temp)
        index = 3;
    temp    = MAX(lambda0, lambdaN);
    lambda_max = MIN(lambda_max, temp);  //  start z plane

    if (lambda_min >= lambda_max)
        return;

    lambda0 = lambda_min;   // lambda = lambda_min
    if (index == 1)
    {
        if (signx == 1)
            v_x = 0;
        else
            v_x = args.NX - 1;
        lambda_x = lambda0 + len_x;

        v_y = (sy + lambda0 * ray_y / L) ;
        tempy = v_y ;
        lambda_y = (absvalue_y < 1.e-4) ? 1.e6 : LAMBDA_Y(tempy + (signy > 0) , sy, dy, L);

        v_z = (sz + lambda0 * ray_z / L) ;
        tempz = v_z ;
        lambda_z = (absvalue_z < 1.e-4) ? 1.e6 : LAMBDA_Z(tempz + (signz > 0) , sz, dz, L);
    }
    else if (index == 2)
    {
        if (signy == 1)
            v_y = 0;
        else
            v_y = args.NY - 1;
        lambda_y = lambda0 + len_y;

        v_x = (sx + lambda0 * ray_x / L) ;
        tempx = v_x ;
        lambda_x = (absvalue_x < 1.e-4) ? 1.e6 : LAMBDA_X(tempx + (signx > 0) , sx, dx, L);

        v_z = (sz + lambda0 * ray_z / L) ;
        tempz = v_z ;
        lambda_z = (absvalue_z < 1.e-4) ? 1.e6 : LAMBDA_Z(tempz + (signz > 0) , sz, dz, L);
    }
    else  // if (index == 3)
    {
        if (signz == 1)
            v_z = 0;
        else
            v_z = args.NZ - 1;
        lambda_z = lambda0 + len_z;

        v_x = (sx + lambda0 * ray_x / L) ;
        tempx = v_x ;
        lambda_x = (absvalue_x < 1.e-4) ? 1.e6 : LAMBDA_X(tempx + (signx > 0) , sx, dx, L);

        v_y = (sy + lambda0 * ray_y / L) ;
        tempy = v_y ;
        lambda_y = (absvalue_y < 1.e-4) ? 1.e6 : LAMBDA_Y(tempy + (signy > 0) , sy, dy, L);
    }
    //  printf("V %d %d %d \n",v_x,v_y,v_z);

    //finale! output part
    while (lambda0 < lambda_max - 5.e-2)
    {
        if (lambda_x <= lambda_y && lambda_x <= lambda_z)
        {
            //(*sino)  += (lambda_x - lambda0) * Data(imageDataPtr, NX, NY, NZ, v_x, v_y, v_z);
            
            ind[numb] = get_img_addr(args, v_x, v_y, v_z);
            wgt[numb] = lambda_x - lambda0;
            //Af += f[ind[numb]]*wgt[numb];
            ++numb;
            
            lambda0   = lambda_x;
            lambda_x += len_x;
            v_x      += signx;
        }
        else if (lambda_y <= lambda_z)
        {
            //(*sino)  += (lambda_y - lambda0) * Data(imageDataPtr, NX, NY, NZ, v_x, v_y, v_z);
            
            ind[numb] = get_img_addr(args, v_x, v_y, v_z);
            wgt[numb] = lambda_y - lambda0;
            //Af += f[ind[numb]]*wgt[numb];
            ++numb;
            
            lambda0   = lambda_y;
            lambda_y += len_y;
            v_y      += signy;
        }
        else
        {
            //(*sino)  += (lambda_z - lambda0) * Data(imageDataPtr, NX, NY, NZ, v_x, v_y, v_z);
            
            ind[numb] = get_img_addr(args, v_x, v_y, v_z);
            wgt[numb] = lambda_z - lambda0;
            //Af += f[ind[numb]]*wgt[numb];
            ++numb;
            
            lambda0   = lambda_z;
            lambda_z += len_z;
            v_z      += signz;
        }
        //  printf("V %d %d %d\n",v_x,v_y,v_z);
        //  printf("Lambda0 %f Lambda_Max %f Diff %f\n",lambda0,lambda_max,lambda_max-lambda0);
    }
}
