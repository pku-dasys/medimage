#include "config.h"
#include <stdio.h>

using ushort = unsigned short;

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