#include "config.h"
#include <stdio.h>

using ushort = unsigned short;

double *f, *v, *g;

void load_data(ushort* data, float* ed, int w, int h, int p,
	const char* data_name)
{
	FILE* input = fopen(data_name, "r");
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
		fscanf("%hd", data+i);
	}
	if(ed)
	{
		for(int i=0; i<w*h; i++)
		{
			fscanf("%f", ed+i);
		}
	}
	fclose(input);
}

void ct3d(double *_image_data, double *_edge_data, double *_sino_data) {
	f = _image_data;
	v = _edge_data;
	g = _sino_data;

	hist_val = new double[GROUP];
	new_hist_val = new double[GROUP];
	group_index = new int[NZ*NX*NY];
	sub_sum = new double[THREAD_NUMB];
	sub_count = new int[THREAD_NUMB];

	load_data
}