#include <GL/glut.h>
#include <fstream>
#include <cstdlib>
#include <cstdio>
using namespace std;
const char* dir;
int x,y,z;
volatile int curz;
float *img;
void mydisp()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels(x*2, y, GL_LUMINANCE, GL_FLOAT, img+curz*x*y*2);
	glFlush();
}
void keyhandler(int ch, int x, int y)
{
	switch(ch)
	{
		case GLUT_KEY_UP:
			if(curz >= z-1) break;
			curz++;
			printf("current z = %d\n", curz);
			break;

		case GLUT_KEY_DOWN:
			if(curz <= 0) break;
			curz--;
			printf("current z = %d\n", curz);
			break;

		default:
			break;
	}
}
void normalkeyhandler(unsigned char ch, int x, int y)
{
	switch(ch)
	{
		case 'q':
		case 27://esc
			exit(0);
	}
}
static float inf()
{
	return 1.0/0.0;
}
void readdata()
{
	char buf[1024];
	float imini = inf(), imaxi = -inf(), emini = inf(), emaxi = -inf();
	for(int i = 0; i < z; i++)
	{
		snprintf(buf, 1024, "%si_%d", dir, i);
		ifstream fimg(buf, ios::in);
		snprintf(buf, 1024, "%se_%d", dir, i);
		ifstream fedg(buf, ios::in);

		for(int j = y - 1; j >= 0; j--)
		{
			for(int k = 0; k < x; k++)
			{
				const int idx = i*x*y*2 + j*x*2 + k;
				fimg >> img[idx];
				if(img[idx] < imini)
					imini = img[idx];
				if(img[idx] > imaxi)
					imaxi = img[idx];
			}
			for(int k = x; k < 2*x; k++)
			{
				const int idx = i*x*y*2 + j*x*2 + k;
				fedg >> img[idx];
				if(img[idx] < emini)
					emini = img[idx];
				if(img[idx] > emaxi)
					emaxi = img[idx];
			}
		}
		fimg.close();
		fedg.close();
	}
	for(int i = 0; i < z; i++)
	{
		for(int j = y - 1; j >= 0; j--)
		{
			for(int k = 0; k < x; k++)
			{
				const int idx = i*x*y*2 + j*x*2 + k;
				img[idx] = (img[idx]-imini)/(imaxi-imini);
			}
			for(int k = x; k < 2*x; k++)
			{
				const int idx = i*x*y*2 + j*x*2 + k;
				img[idx] = (img[idx]-emini)/(emaxi-emini);
			}
		}
	}
}
int main(int argc, char **argv)
{
	if(argc < 5)
	{
		printf("usage: viewer dir x y z [z0]\n");
		exit(0);
	}
	dir = argv[1];
	x = atoi(argv[2]);
	y = atoi(argv[3]);
	z = atoi(argv[4]);
	if(argc > 5)
	{
		curz = atoi(argv[5]);
		argv[5] = argv[0];
		argc -= 5;
		argv += 5;
	}
	else
	{
		curz = 0;
		argv[4] = argv[0];
		argc -= 4;
		argv += 4;
	}
	img = new float[x*y*z*2];
	readdata();
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
	glutInitWindowSize(x*2, y);
	glutInitWindowPosition(200, 200);
	glutCreateWindow("simple viewer");
	glutSpecialFunc(keyhandler);
	glutKeyboardFunc(normalkeyhandler);
	glutDisplayFunc(mydisp);
	glutIdleFunc(mydisp);
	glutMainLoop();
	delete[] img;
	return 0;
}