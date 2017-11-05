#include "config.h"
#include "utility.h"
#include "main.h"
#include "ct3d.h"
#include "tracing.h"
#include <stdio.h>
#include <boost/lexical_cast.hpp>
using namespace std;
using namespace boost;

using ushort = unsigned short;

double *f, *v, *g, *ed;
double SOD, SDD;//src to obj, src to board
double AT_sum;
double lambda;
int max_numb;

void load_data(double* data, double* ed, int w, int h, int p,
	const char* data_name)
{
	FILE* input = fopen(data_name, "rb");
	ushort ustmp[1];
	float ftmp[1];
	if(!input)
	{
		fprintf(stderr, "cannot open data file\n");
		exit(1);
	}
	puts("loading");
	fseek(input, 1024, SEEK_SET);//ignore header
	for(int i=0; i<w*h*p; i++)
	{
		if(!fread(ustmp, sizeof(ushort), 1, input))
		{
			fprintf(stderr, "unexpect ending\n");
			break;
		}
		data[i] = ustmp[0];
	}
	if(ed)
	{
		for(int i=0; i<w*h; i++)
		{
			if(!fread(ftmp, sizeof(float), 1, input))
			{
				fprintf(stderr, "unexpect ending\n");
				break;
			}
			ed[i] = ftmp[0];
		}
	}
	fclose(input);
	puts("load done");
}

void minIMAGE(double Af, int *line, double *weight, int numb)
{
	int i;
	for(i = 0; i<numb; i++)
	{
		Af += f[line[i]] * weight[i];
	}
	AT_sum += sqr(Af);
	for(i = 0; i<numb; i++)
	{
		int ind = line[i];
		f[ind] += lambda * (-Af * weight[i]);
		if(f[ind]<0) f[ind] = 0;
	}
}

double *cos_table;
double *sin_table;

void compute(int alpha, int detectorX, int detectorZ)
/*alpha is the rotate angle of light source
(or the object, then result will be a mirror image),
from the view of +Z to -Z, counterclockwise(+X->+Y->-X->-Y)
is positive. Initial angle is +Y
the parameters should be int, as it represents the relative
coordinate of grid in the detector
*/
{
	const double sina = sin_table[alpha], cosa = cos_table[alpha];
	double srcX = -SOD * sina,
		   srcY =  SOD * cosa,
		   oridstX = detectorX - NDETECTORX*0.5 + 0.5,
		   oridstY = SOD - SDD,
		   dstX = oridstX * cosa - oridstY * sina,
		   dstY = oridstX * sina + oridstY * cosa,
		   dstZ = detectorZ - NDETECTORZ*0.5 + 0.5;
	int *ind = new int[MAX_ELE_RAY];
	double *wgt = new double[MAX_ELE_RAY];

	double Af = -array_3d_sino(g, detectorX, detectorZ, alpha);
	//I don't know which is correct? X,Z,a? Z,X,a?
	int numb;
	{
	forward_proj(srcX, srcY, 0.0,
				 dstX, dstY, dstZ,
				 ind, wgt, numb);
	}
	{
	minIMAGE(Af, ind, wgt, numb);
	}
	delete[] ind;
	delete[] wgt;
}

const int GROUP = 16;

double *hist_val;
double *new_hist_val;

int *group_index;

double sum_val;
int count_val;
double *sub_sum;
int *sub_count;

void wrapper()
{
	for(int k = 0; k < NANGLE; ++k)
	{
		for(int i = 0; i < NDETECTORX; i++)
		{
			for(int j = 0; j < NDETECTORZ; j++)
			{
				compute(k, i, j);
			}
		}
	}
}

void ct3d(double *_image_data,/* double *_edge_data,*/ double *_sino_data)
{
//deprecated: edge data
	f = _image_data;
//	v = _edge_data;
	g = _sino_data;

	ed = new double[NX*NY];
	hist_val = new double[GROUP];
	new_hist_val = new double[GROUP];
	group_index = new int[NZ*NX*NY];
	sub_sum = new double[THREAD_NUMB];
	sub_count = new int[THREAD_NUMB];

	load_data(g, ed, NX, NY, NZ, sino_filename.c_str());
	memset(f, 0, NZ*NX*NY*sizeof(double));
//	memset(v, 0, NZ*NX*NY*sizeof(double));

	cos_table = new double[NANGLE];
	sin_table = new double[NANGLE];

	for(int i = 0; i < NANGLE; i++)
	{
		const double alpha = (double)i/NANGLE*2*PI;
		sin_table[i] = sin(alpha);
		cos_table[i] = cos(alpha);
	}

	for(int iters = 0; iters < CT_ITERATIONS; ++iters)
	{
		lambda = 1.0/(1000.0+iters*50.0);

		printf("iter = %d, lambda = %lf\n", iters, lambda);
		AT_sum = 0;
		fflush(stdout);

		wrapper();

		printf("AT_sum = %E\n", AT_sum);

		write_data_3d(f, NZ, NX, NY, OUTPUT_DIR+"/"+image_filename+lexical_cast<string>(iters));
	}
	puts("iter end");

	delete[] sin_table;
	delete[] cos_table;
	delete[] ed;
    delete [] hist_val;
    delete [] new_hist_val;
    delete [] group_index;
    delete [] sub_sum;
    delete [] sub_count;
}