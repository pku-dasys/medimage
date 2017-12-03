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
	glDrawPixels(x, y, GL_LUMINANCE, GL_FLOAT, img+curz*x*y);
	glFlush();
}
void keyhandler(int ch, int x, int y)
{
	switch(ch)
	{
		case GLUT_KEY_UP:
			curz++;
			break;

		case GLUT_KEY_DOWN:
			curz--;
			break;

		default:
			break;
	}
	if(curz < 0) curz = 0;
	if(curz >= z) curz = z-1;
	printf("current z = %d\n", curz);
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
	float mini = inf(), maxi = -inf();
	for(int i = 0; i < z; i++)
	{
		snprintf(buf, 1024, "%s/e_%d", dir, i);
		ifstream fin(buf, ios::in);
		for(int j = y - 1; j >= 0; j--)
		{
			for(int k = 0; k < x; k++)
			{
				const int idx = i*x*y + j*x + k;
				fin >> img[idx];
				if(img[idx] < mini)
					mini = img[idx];
				if(img[idx] > maxi)
					maxi = img[idx];
			}
		}
		fin.close();
	}
	for(int i = 0; i < z; i++)
	{
		for(int j = y - 1; j >= 0; j--)
		{
			for(int k = 0; k < x; k++)
			{
				const int idx = i*x*y + j*x + k;
				img[idx] = (img[idx]-mini)/(maxi-mini);
			}
		}
	}
}
int main(int argc, char **argv)
{
	if(argc < 5)
	{
		printf("usage: viewer dir x y z\n");
		exit(0);
	}
	dir = argv[1];
	x = atoi(argv[2]);
	y = atoi(argv[3]);
	z = atoi(argv[4]);
	curz = 0;
	img = new float[x*y*z];
	argv[4] = argv[0];
	argc -= 4;
	readdata();
	glutInit(&argc, argv+4);
	glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
	glutInitWindowSize(x, y);
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