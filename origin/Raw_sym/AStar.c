#include "AStar.h"
#include "c_imp.h"

#include <string.h>
#include <math.h>

#define midray 383
#define tilex_l 128
#define tiley_l 128
#define tilex_n 4
#define tiley_n 4
#define ntiles 16

void AStar(double back[imageHight][imageWidth],double proj[nAlpha][projWidth]) {
  memset (back, 0, sizeof (double)*imageHight*imageWidth);
    //For nAlpha = 0
  int RayNum0,RayNum90,alpha;
    for(RayNum0=0; RayNum0<projWidth; RayNum0++) {
     int ind0 = RayNum0 - (projWidth / 2 );
     if((ind0>(imageHight/2-1)) || (ind0<(1-imageHight/2)))
         continue;
     int iy0 = imageHight/2 - 1 + ind0;
     int tmp0;
     for(tmp0=0; tmp0<imageWidth; tmp0++){
         back[iy0][tmp0] += proj[0][RayNum0];
     }
 }

    //For nAlpha = 90
 for(RayNum90=0; RayNum90<projWidth; RayNum90++){
    int ind90 = RayNum90 - (projWidth/2);
    if((ind90>(imageWidth/2 -1 )) || (ind90<(1-imageWidth/2)))
        continue;
    int ix90 = imageWidth/2 - 1 + ind90;
    int tmp90;
    for(tmp90=0; tmp90<imageHight; tmp90++){
        back[tmp90][ix90] += proj[90][RayNum90];
    }
}

    //For nAlpha in [1,45] & [89,45] & [91,135] & [179,135]
for (alpha = 1; alpha < 46; alpha++)
{
    double cf;
    double sf;
    double hx;
    double hy;
    double tanphi;
    double rtnphi;
    double c;
    int ix;
    int iy;
    double fp;
    int incx;
    double w;
    double yc;
    double xc;

    int next;

      //jk  5/21/2008
      //if the values for _sorx and _sory are provided, do not call posit()
    sf = sin (3.141593 / 180 * alpha);
    cf = cos (3.141593 / 180 * alpha);

    hx = ((double) 1.0) / sf;
    hy = ((double) 1.0) / cf;
    tanphi = sf / cf;
    rtnphi = cf / sf;

    int TileNum,RayNum;
      // CALCULATE THE INTERSECTION OF REFERANCE LINE AND THE RAY IN QUESTI
      //printf("hx=%f, hy=%f, alpha = %d\n", hx, hy, alpha);
    for (TileNum = 0; TileNum < ntiles; TileNum ++){
      for (RayNum = 0; RayNum < (midray + 1); RayNum++)
      {
      //int ind = RayNum - 383;
          int ind = RayNum - (projWidth / 2 );
      //For tiling
          int dx = tilex_l * (tilex_n - 1 - (TileNum%tilex_n));
          int dy = tiley_l * (TileNum/tilex_n);

          c = cf * ind + sf * ind * tanphi;
     // c = imageWidth / 2 - tanphi * imageWidth / 2 - c;
          c = (imageWidth/2 - dy) - tanphi * (imageWidth/2 - tilex_l * (TileNum%tilex_n)) - c;
      //printf("hx=%f, hy=%f, alpha = %d, c=%f, RayNum = %d\n", hx, hy, alpha, c, RayNum);

      // EQUATION OF RAY IS NOW Y=X*TANPHI+C
      // C IS THE Y INTERCEPT
      // IF INTERCEPT IS ABOVE (N) THEN THE RAY HAS MISSED ALL PIXELS

          if (c < 0)
          {
          // IF Y INTERCEPT IS NEGATIVE THEN TAKE X INTERCEPT (WHICH MUST
          // BE POSITIVE SINCE SLOPE IS ALWAYS POSITIVE

          //  RAY ENTERS FROM BOTTOM

              c = -c * rtnphi;

          // RECALCULATE C SO THAT EQUATION OF LINE IS NOW
          // X=RTANPHI*Y+C
          // IF X INTERCEPT IS GREATER TNAN N THE RAY HAS MISSED THE RECONSTRUC
          // AREA

              ix = (int) c;
              iy = 0;

          // FP IS ACTUAL DISPLACEMENT OF INTERCEPT FROM THE NEXT GREATER
          // GRID LINE

              fp = ((double) (ix + 1)) - c;

          // START AT LOW X VALUES AND INCREMENT ACROSS PICTURE UNLESS
          // SLOPE IS NEGATIVE

              ix = tilex_l - ix - 1;
              incx = -1;

          // SET INITIAL VALUES OF THE LENGTH OF THE RAY IN CURRENT CELL

              w = hy * fp;
              yc = fp * tanphi;
              xc = -fp;

          // INITIAL CONDITIONS SET FOR ENTRANCE TO LOOPS WE HAVE ENTERED
          // THE PICTURE AREA FROM THE BOTTOM

              next = 2;
          }
          else
          {

          // ENTER PICTURE FROM SIDE (LEFT)
          // CALCULATE INITIAL GRID COORDINATES (IX,IY)
              ix = tilex_l - 1;
              incx = -1;

              iy = (int) c;

          // FP IS DISTANCE FROM Y INTERCEPT TO NEXT GREATER GRID LINE

              fp = ((double) (iy + 1)) - c;

          // CALCULATE ARRAY SUBSCRIPT IN PICTURE AREA OF INITIAL CELL
          // INDEX INTO FORTRAN TYPE ARRAY

          // CALCULATE INITIAL DESCRIPTION OF RAY IN FIRST CELL FOR LOOP

              yc = -fp;
              w = hx * fp;
              xc = fp * rtnphi;
          // RAY ENTERS FROM LEFT SIDE
              next = 1;
          }

          while (ix >= 0 && ix < tilex_l && iy >= 0 && iy < tiley_l)
          {
              switch (next)
              {
                case 1:
          //  MAKE ROOM FOR NEW PICTURE ELEMENT
          // STEP ONE IN X DIRECTION AND CORRESPONDING STEP IN Y DIRECTION
                xc -= (double) 1.0;
                yc += tanphi;

          // IF XC IS NEGATIVE THE Y GRID LINE HAS BEEN CROSSED THE RAY
          // HAS LEFT THROUGHTHE TOP

                if (xc > 0)
                {
                  back[iy+dy][ix+dx] +=
                  proj[alpha][RayNum] * hy;

                  back[imageWidth-1-ix-dx][iy+dy] +=
                  proj[90+alpha][RayNum] * hy;

                  if(RayNum != midray){

                      back[imageWidth-1-iy-dy][imageWidth-1-ix-dx] +=
                      proj[alpha][projWidth-1-RayNum] * hy;

                      back[ix+dx][imageWidth-1-iy-dy] +=
                      proj[90+alpha][projWidth-1-RayNum] * hy;
                  }

                  if(alpha != 45){

                      back[ix+dx][iy+dy] +=
                      proj[90-alpha][RayNum] * hy;

                      back[imageWidth-1-iy-dy][ix+dx] +=
                      proj[180-alpha][RayNum] * hy;

                      if(RayNum != midray){

                          back[imageWidth-1-ix-dx][imageWidth-1-iy-dy] +=
                          proj[90-alpha][projWidth-1-RayNum] * hy;

                          back[iy+dy][imageWidth-1-ix-dx] +=
                          proj[180-alpha][projWidth-1-RayNum] * hy;
                      }
                  }
                  w -= hy;
                  ix += incx;
              // next = 1;
              }
              else
              {
                  back[iy+dy][ix+dx] +=
                  proj[alpha][RayNum] * w;

                  back[imageWidth-1-ix-dx][iy+dy] +=
                  proj[90+alpha][RayNum] * w;

                  if(RayNum != midray){

                      back[imageWidth-1-iy-dy][imageWidth-1-ix-dx] +=
                      proj[alpha][projWidth-1-RayNum] * w;

                      back[ix+dx][imageWidth-1-iy-dy] +=
                      proj[90+alpha][projWidth-1-RayNum] * w;
                  }

                  if(alpha != 45){

                      back[ix+dx][iy+dy] +=
                      proj[90-alpha][RayNum] * w;

                      back[imageWidth-1-iy-dy][ix+dx] +=
                      proj[180-alpha][RayNum] * w;

                      if(RayNum != midray){

                          back[imageWidth-1-ix-dx][imageWidth-1-iy-dy] +=
                          proj[90-alpha][projWidth-1-RayNum] * w;

                          back[iy+dy][imageWidth-1-ix-dx] +=
                          proj[180-alpha][projWidth-1-RayNum] * w;
                      }
                  }

                  w = hy - w;
                  iy++;
                  next = 2;
              }
              break;
              case 2:
          // STEP ONE IN Y DIRECTION AND CORRESPONDING X INCREMENT
              yc -= 1.0;
              xc += rtnphi;
          // IF YC IS NEGATIVE WH HAVE CROSSED A X-GRID LINE AND THAT
          // WE LEAVE THROUGE THE SIDE
              if (yc > 0)
              {
                  w -= hx;
              // ENTERED BOTTOM LEFT TOP ASSIGN MAXIMUM X-TYPE WEIGHT
                  back[iy+dy][ix+dx] +=
                  proj[alpha][RayNum] * hx;

                  back[imageWidth-1-ix-dx][iy+dy] +=
                  proj[90+alpha][RayNum] * hx;

                  if(RayNum != midray){

                      back[imageWidth-1-iy-dy][imageWidth-1-ix-dx] +=
                      proj[alpha][projWidth-1-RayNum] * hx;

                      back[ix+dx][imageWidth-1-iy-dy] +=
                      proj[90+alpha][projWidth-1-RayNum] * hx;
                  }

                  if(alpha != 45){

                      back[ix+dx][iy+dy] +=
                      proj[90-alpha][RayNum] * hx;

                      back[imageWidth-1-iy-dy][ix+dx] +=
                      proj[180-alpha][RayNum] * hx;

                      if(RayNum != midray){

                          back[imageWidth-1-ix-dx][imageWidth-1-iy-dy] +=
                          proj[90-alpha][projWidth-1-RayNum] * hx;

                          back[iy+dy][imageWidth-1-ix-dx] +=
                          proj[180-alpha][projWidth-1-RayNum] * hx;
                      }
                  }

                  iy++;
              //next = 2;
              }
              else
              {
                  back[iy+dy][ix+dx] +=
                  proj[alpha][RayNum] * w;

                  back[imageWidth-1-ix-dx][iy+dy] +=
                  proj[90+alpha][RayNum] * w;

                  if(RayNum != midray){

                      back[imageWidth-1-iy-dy][imageWidth-1-ix-dx] +=
                      proj[alpha][projWidth-1-RayNum] * w;

                      back[ix+dx][imageWidth-1-iy-dy] +=
                      proj[90+alpha][projWidth-1-RayNum] * w;
                  }

                  if(alpha != 45){

                      back[ix+dx][iy+dy] +=
                      proj[90-alpha][RayNum] * w;

                      back[imageWidth-1-iy-dy][ix+dx] +=
                      proj[180-alpha][RayNum] * w;

                      if(RayNum != midray){

                          back[imageWidth-1-ix-dx][imageWidth-1-iy-dy] +=
                          proj[90-alpha][projWidth-1-RayNum] * w;

                          back[iy+dy][imageWidth-1-ix-dx] +=
                          proj[180-alpha][projWidth-1-RayNum] * w;
                      }
                  }
                  w = hx - w;
                  ix += incx;
                  next = 1;
              }
              break;
          }
      }
  }
}
}
}