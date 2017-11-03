#include "config.h"
#include "utility.h"
#include <stdio.h>

using ushort = unsigned short;

double *f, *v, *g;

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
		if(eof(input))
		{
			fprintf(stderr, "unexpect ending\n");
			break;
		}
		fread(ustmp, sizeof(ushort), 1, input);
		data[i] = ustmp[0];
	}
	if(ed)
	{
		for(int i=0; i<w*h; i++)
		{
			if(eof(input))
			{
				fprintf(stderr, "unexpect ending\n");
				break;
			}
			fread(ftmp, sizeof(float), 1, input);
			ed[i] = ftmp[0];
		}
	}
	fclose(input);
}

void ct3d(double *_image_data, double *_edge_data, double *_sino_data) {
	f = _image_data;
	v = _edge_data;
	g = _sino_data;

	ed = new float[NX*NY];
	hist_val = new double[GROUP];
	new_hist_val = new double[GROUP];
	group_index = new int[NZ*NX*NY];
	sub_sum = new double[THREAD_NUMB];
	sub_count = new int[THREAD_NUMB];

	load_data(g, ed, NX, NY, , sino_filename);
}