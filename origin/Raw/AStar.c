#include "AStar.h"
#include "c_imp.h"

#include <string.h>
#include <math.h>

void AStar(double f[imageHight][imageWidth],double g[nAlpha][imageWidth]) {
    int i;
    
    for (i = 0; i<nAlpha; ++i)
        alpha[i] = i;
    for (i = 0; i<imageHight; ++i)
        xy[i] = (double)i/(imageHight-1) - 0.5;

    memset(f, 0, imageHight*imageWidth*sizeof(double));

    for (i = 0; i<nAlpha; ++i) {
        double a = alpha[i] * M_PI / 180.0;

        int nx,ny;
        for (nx = 0; nx<imageHight; ++nx)
            for (ny = 0; ny<imageWidth; ++ny)
                lines[nx][ny] = g[i][nx];

        RotateImage(lines, rot, -a, xy);
        for (nx = 0; nx<imageHight; ++nx)
            for (ny = 0; ny<imageWidth; ++ny)
                f[nx][ny] += rot[nx][ny];
    }
}
