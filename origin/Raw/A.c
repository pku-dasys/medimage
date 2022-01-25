#include "A.h"
#include "c_imp.h"

#include <string.h>
#include <math.h>

void A(double g[nAlpha][imageHight],double f[imageHight][imageWidth]) {
    int i;
    
    for (i = 0; i<nAlpha; ++i)
        alpha[i] = i;
    for (i = 0; i<imageHight; ++i)
        xy[i] = (double)i/(imageHight-1) - 0.5;

    for (i = 0; i<nAlpha; ++i) {
        double a = alpha[i] * M_PI / 180.0;

        RotateImage(f, rot, a, xy);
        int nx,ny;
        for (nx = 0; nx<imageHight; nx++) {
            g[i][nx] = 0;
            for (ny = 0; ny<imageWidth; ny++)
                g[i][nx] += rot[nx][ny];
        }
    }
}
