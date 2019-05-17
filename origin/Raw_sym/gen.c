#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main() {
	srand(time(0));
	int i,j;
	for (i = 0; i<512; ++i) {
		for (j = 0; j<512; ++j) {
			int c = 16;
			int b = (i/c)%2+(j/c)%2;
			if (b==1) {
				double r = rand()%500000;
				printf("%f ",r/1000000);
			}
			else {
				double r = rand()%500000+500000;
				printf("%f ",r/1000000);
			}
			//double r = rand()%1000000;
			//printf("%f ",r/1000000);
		}
		printf("\n");
	}
	return 0;
}