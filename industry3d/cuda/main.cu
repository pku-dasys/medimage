#include "ct3d.h"
#include "utility.h"
#include "tracing.h"

#include <iostream>

#include <omp.h>

#include <cmath>

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
using namespace std;
using boost::format;


class Debug{
	
	public:
	static const long size;
	
	int *h_length;
	int *d_length;

	Debug(){
	
		h_length=(int*)malloc(size*sizeof(int));
    	cudaError_t err = cudaSuccess;
		err =cudaMalloc((void**)&d_length,size*sizeof(int));
		if (err != cudaSuccess){
			fprintf(stderr, "Failed to allocate debug vector (error code %s)!\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE);
		}

	}
	
	void transferToHost(){
    	cudaError_t err = cudaSuccess;
		err=cudaMemcpy(h_length, d_length, sizeof(int)*size, cudaMemcpyDeviceToHost);
		if (err != cudaSuccess){
			fprintf(stderr, "Failed to copy debug vector (error code %s)!\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE);
		}

	}

	void print(){
			for (int j=0; j<128*128; ++j)
				for (int k=0; k<384; ++k)
					printf("%d ", h_length[j*384+k]);
	}


};

long const Debug::size = 128*128*384;


int max_numb;
float ave_numb;

__device__ int d_max_numb;
__device__ float d_ave_numb;

__device__ float cud_sqr(float x) {return x*x;}


__device__ void cud_rotate_2d(float &x,float &y,float theta) {
    float X[2];
	X[0]=x;
	X[1]=y;

    float costheta = __cosf(theta);
	float sintheta = __sinf(theta);
    float R[2][2];
	R[0][0]=costheta;
	R[0][1]=-sintheta;
	R[1][0]=sintheta;
	R[1][1]=costheta;

    float d[2];
	d[0]=0.f;
	d[1]=1.f;
    for (int k = 0; k<2; ++k)
        for (int i = 0; i<2; ++i)
                d[i] += R[i][k]*X[k];
    
    x = d[0];
    y = d[1];
}


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

__device__ int64_t cud_get_img_addr(const Parameter &args,int x,int y,int z) {
    return (int64_t)z*args.NX*args.NY+x*args.NY+y;
}

__device__ void cud_forward_proj(int NX, int NY, int NZ,
                  float sx,float sy,float sz,
                  float dx,float dy,float dz,
                  int64_t *ind,float *wgt,int &numb) {

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
    L = sqrt( cud_sqr(ray_x)+cud_sqr(ray_y)+cud_sqr(ray_z) );
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

    tempx = NX;
    tempy = NY;
    tempz = NZ;
    
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
            v_x = NX - 1;
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
            v_y = NY - 1;
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
            v_z = NZ - 1;
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
            
    		ind[numb] =  (int64_t)v_z*NX*NY+v_x*NY+v_y;
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
            
    		ind[numb] =  (int64_t)v_z*NX*NY+v_x*NY+v_y;
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
            
    		ind[numb] =  (int64_t)v_z*NX*NY+v_x*NY+v_y;
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

__global__ void kernel(float LAMBDA_IMG, float LAMBDA_EDGE, float *f, float *v, float *g, int NX, int NY, int NZ, int NPROJ, int NDX, int NDY, int HALFDET, float PIXELSIZE, float SOD, float SDD, const int MAX_RAYLEN, float SAMPLESIZE, float HALFSIZE, float ALPHA, float BETA, float EPSILON,
		int *length) {


	int detectorX=blockIdx.x%NDX;
	//int alpha=blockIdx.y%NPROJ;
	int alpha=64;
	int detectorY=threadIdx.x%NDY;

    float srcX,srcY,srcZ;
    float dstX,dstY,dstZ;

	//cud_parallel

    srcX = detectorY+0.5 - HALFDET;
    srcX *= PIXELSIZE;
    srcY = -SOD;

    srcZ = HALFDET - detectorX - 0.5;
    srcZ *= PIXELSIZE;

    dstX = detectorY+0.5 - HALFDET;
    dstX *= PIXELSIZE;
    dstY = SDD-SOD;

    dstZ = HALFDET - detectorX - 0.5;
    dstZ *= PIXELSIZE;

    float theta = (float)alpha/NPROJ*2*PI;

    cud_rotate_2d(srcX, srcY, theta);
    cud_rotate_2d(dstX, dstY, theta);

    srcX /= SAMPLESIZE;
    srcY /= SAMPLESIZE;
    srcZ/= SAMPLESIZE;
    dstX /= SAMPLESIZE;
    dstY /= SAMPLESIZE;
    dstZ /= SAMPLESIZE;

    srcX += HALFSIZE;
    srcY += HALFSIZE;
    srcZ += HALFSIZE;
    dstX += HALFSIZE;
    dstY += HALFSIZE;
    dstZ += HALFSIZE;


	/*
	//cud_cone
    srcZ = 0.0;
    srcX = 0.0;
    srcY = -SOD;
    float theta = (float)alpha/NPROJ*2*PI;
    cud_rotate_2d(srcX, srcY, theta);

    dstX = detectorY - HALFDET + 0.5;
    dstX *= PIXELSIZE;
    dstY = (SDD - SOD);
    dstZ = HALFDET - detectorX - 0.5;
    dstZ *= PIXELSIZE;
    cud_rotate_2d(dstX, dstY, theta);

    srcX /= SAMPLESIZE;
    srcZ /= SAMPLESIZE;
    srcY /= SAMPLESIZE;
    dstX /= SAMPLESIZE;
    dstY /= SAMPLESIZE;
    dstZ /= SAMPLESIZE;

    srcX += NX/2.0;
    srcY += NY/2.0;
    srcZ += NZ/2.0;
    dstX += NX/2.0;
    dstY += NY/2.0;
    dstZ += NZ/2.0;
	*/

    float Af = -g[alpha*NDX*NDY+detectorX*NDY+detectorY];

	int64_t line[384];
	float weight[384];
	float d[384];
    int numb;

	for (int i=0; i<MAX_RAYLEN; ++i){
		d[i]=0.f;
		line[i]=0;
		weight[i]=0.f;
		numb = 0;
	}

    cud_forward_proj(NX, NY, NZ,
                 srcX, srcY, srcZ,
                 dstX, dstY, dstZ,
                 line, weight, numb);
   
    if (d_max_numb < numb) d_max_numb = numb;
    d_ave_numb += numb;


	length[detectorX*NDY+detectorY]=numb;
	return;

    for (int i = 0; i<numb; ++i) {
        Af += f[line[i]] * weight[i];
    }
    for (int i = 0; i<numb; ++i) {
        int64_t ind = line[i];
        int64_t plain = ind%(NX*NY);
        int x = plain/NY, y = plain%NY;

        float tmp = 0.;
        float lap = 0.;
        
        if (x+1<NX) tmp += cud_sqr(v[ind])*(v[ind+NY]-v[ind]);
        else             tmp += cud_sqr(v[ind])*(       0        -v[ind]);
        
        if (y+1<NY) tmp += cud_sqr(v[ind])*(v[ind+1]-v[ind]);
        else             tmp += cud_sqr(v[ind])*(     0    -v[ind]);
        
        if (x-1>=0)      tmp -= cud_sqr(v[ind-NY])*(v[ind]-v[ind-NY]);
        else             tmp -=                        (v[ind]-0        );
        
        if (y-1>=0)      tmp -= cud_sqr(v[ind-1])*(v[ind]-v[ind-1]);
        else             tmp -=                  (v[ind]-0       );
        
        if (x+1<NX) lap += v[ind+NY];
        if (y+1<NY) lap += v[ind+1];
        if (x-1>=0)      lap += v[ind-NY];
        if (y-1>=0)      lap += v[ind-1];
        lap -= 4*v[ind];

        d[i] = -Af*weight[i]+ALPHA*(tmp+cud_sqr(EPSILON)*lap);
    }
    for (int i = 0; i<numb; i++) {
        int64_t ind = line[i];
        float tmp = f[ind] + LAMBDA_IMG * d[i];
        if (tmp<0) tmp = 0;
        f[ind] = tmp;
    }

	for (int i=0; i<MAX_RAYLEN; ++i)
		d[i]=0.f;
    for (int i = 0; i<numb; ++i) {
        int64_t ind = line[i];
        int64_t plain = ind%(NX*NY);
        int x = plain/NY, y = plain%NY;

        float a = 0.;
        float b = 0.;
        float c = 0.;

        if (x-1>=0)      a += cud_sqr(v[ind]-v[ind-NY]);
        else             a += cud_sqr(v[ind]-0        );
        
        if (y-1>=0)      a += cud_sqr(v[ind]-v[ind-1]);
        else             a += cud_sqr(v[ind]-0       );
        
        a *= v[ind];
        
        b = v[ind]-1;

        if (x+1<NX) c += v[ind+NY];
        if (y+1<NY) c += v[ind+1];
        if (x-1>=0)      c += v[ind-NY];
        if (y-1>=0)      c += v[ind-1];
        c -= 4*v[ind];
        
        d[i] = -ALPHA*a-BETA/(4*EPSILON)*b+BETA*EPSILON*c;
    }
    for (int i = 0; i<numb; i++) {
        int64_t ind = line[i];
        float tmp = v[ind] + LAMBDA_EDGE * d[i];
        if (tmp<0) tmp = 0;
        if (tmp>1) tmp = 1;
        v[ind] = tmp;
    }
}


void ct3d(Parameter &args,const CTInput &in,CTOutput &out) {

    max_numb = 0;
	ave_numb = 0;

    out.allocate();

    int64_t global_start = timer_s();


	float *d_g=NULL;
	float *d_f=NULL;
	float *d_v=NULL;

    cudaError_t err = cudaSuccess;
	err=cudaMalloc((void**)&d_f,sizeof(float)*args.NX*args.NY*args.NZ);
    if (err != cudaSuccess){
		fprintf(stderr, "Failed to allocate device vector f (error code %s)!\n", cudaGetErrorString(err));
    	exit(EXIT_FAILURE);
    }

	err=cudaMalloc((void**)&d_v,sizeof(float)*args.NX*args.NY*args.NZ);
    if (err != cudaSuccess){
		fprintf(stderr, "Failed to allocate device vector v (error code %s)!\n", cudaGetErrorString(err));
    	exit(EXIT_FAILURE);
    }

	err=cudaMalloc((void**)&d_g,sizeof(float)*args.NPROJ*args.NDX*args.NDY);
    if (err != cudaSuccess){
		fprintf(stderr, "Failed to allocate device vector g (error code %s)!\n", cudaGetErrorString(err));
    	exit(EXIT_FAILURE);
    }
	
	err=cudaMemcpy(d_f, out.img, sizeof(float)*args.NX*args.NY*args.NZ, cudaMemcpyHostToDevice);
	if (err != cudaSuccess){
		fprintf(stderr, "Failed to copy vector f from host to device (error code %s)!\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	err=cudaMemcpy(d_v, out.edge, sizeof(float)*args.NX*args.NY*args.NZ, cudaMemcpyHostToDevice);
	if (err != cudaSuccess){
		fprintf(stderr, "Failed to copy vector v from host to device (error code %s)!\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	err=cudaMemcpy(d_g, in.sino, sizeof(float)*args.NPROJ*args.NDX*args.NDY, cudaMemcpyHostToDevice);
	if (err != cudaSuccess){
		fprintf(stderr, "Failed to copy vector g from host to device (error code %s)!\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	cudaMemcpyToSymbol(d_max_numb, &max_numb, sizeof(int));
	cudaMemcpyToSymbol(d_ave_numb, &ave_numb, sizeof(float));

    float init_lambda_img = 1/args.LAMBDA_IMG;
    float init_lambda_edge = 1/args.LAMBDA_EDGE;
	Debug debug;

    for (int iters = 0; iters < args.ITERATIONS; ++iters) {

        args.LAMBDA_IMG = 1.0/(init_lambda_img+(float)iters/args.ITERATIONS*(args.AMPLIFIER-1)*init_lambda_img);
        args.LAMBDA_EDGE = 1.0/(init_lambda_edge+(float)iters/args.ITERATIONS*(args.AMPLIFIER-1)*init_lambda_edge);
        cout << format("iter = %1%, lambda_img = %2%, lambda_edge = %3%") % iters % args.LAMBDA_IMG % args.LAMBDA_EDGE <<endl;

		dim3 blocks(128, 1, 1);
		int threads=128;


		kernel<<<blocks, threads>>>(args.LAMBDA_IMG, args.LAMBDA_EDGE, d_f, d_v, d_g, args.NX, args.NY, args.NZ, args.NPROJ, args.NDX, args.NDY, args.HALFDET, args.PIXELSIZE, args.SOD, args.SDD, args.MAX_RAYLEN, args.SAMPLESIZE, args.HALFSIZE, args.ALPHA, args.BETA, args.EPSILON, debug.d_length);

		err = cudaGetLastError();

		if (err != cudaSuccess)
		{
			fprintf(stderr, "Failed to launch MumfordShah kernel (error code %s)!\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE);
		}
    }

	err=cudaMemcpy(out.img, d_f, sizeof(float)*args.NX*args.NY*args.NZ, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess){
		fprintf(stderr, "Failed to copy vector f from device to host (error code %s)!\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	err=cudaMemcpy(out.edge, d_v, sizeof(float)*args.NX*args.NY*args.NZ, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess){
		fprintf(stderr, "Failed to copy vector v from device to host (error code %s)!\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	cudaMemcpyFromSymbol(&max_numb, d_max_numb, sizeof(int));
	cudaMemcpyFromSymbol(&ave_numb, d_ave_numb, sizeof(float));
	
    int64_t global_end = timer_s();
	ave_numb /= args.NPROJ*args.NDX*args.NDY;
    cout << "===========================================" <<endl;
    cout << boost::format("TOTAL time used = %1% seconds") % (global_end-global_start) <<endl;
    cout<< boost::format("actual max raylen = %1%, average raylen = %2%") % max_numb % ave_numb << endl;

	debug.print();
	err=cudaFree(d_f);
    if (err != cudaSuccess){
       	fprintf(stderr, "Failed to free vector d_f (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
	}

	err=cudaFree(d_v);
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to free vector d_v (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
	}

	err=cudaFree(d_g);
    if (err != cudaSuccess){
       	fprintf(stderr, "Failed to free vector d_f (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
	}
    

}

int main(int argc, char** argv) {
    Parameter args;

    args.parse_config(argc, argv);
    args.print_options();

    CTInput in = CTInput(args);
    CTOutput out = CTOutput(args);

    in.read_sino(args.RAW_DATA_FILE);
    ct3d(args,in,out);
    out.write_img(args.OUTPUT_DIR);
    out.write_edge(args.OUTPUT_DIR);

    return 0;
}
