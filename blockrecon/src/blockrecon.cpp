#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>

#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "config.h"
#include "wray.h"
#include "utility.h"

using namespace std;

// global variables
float _f[IMGSIZE][IMGSIZE];
float _v[IMGSIZE][IMGSIZE];
float _g[NPROJ][NRAY];

// true blocking g
float global_g[NBLK][NPROJ][NRAY];

// iterative variables
float f[NBLK][BLKSIZE][BLKSIZE];
float v[NBLK][BLKSIZE][BLKSIZE];
float g[NBLK][NPROJ][NRAY];


float pad_f[NBLK][NDIR][BLKSIZE];
float pad_v[NBLK][NDIR][BLKSIZE];
float g_dot[NBLK][NPROJ][NRAY];

float g_sum[NPROJ][NRAY];

// fake constants
int cnt_ray[NPROJ][NRAY];
float w_blk[NPROJ][NRAY][NBLK];
int adj[NBLK][NDIR];
float lambda, g_relax;

void init_vars();

void minIMAGE(int _bid,int np,int nr,int *line,float *weight,int numb,float snorm) {
    int i;
    float d[IMGSIZE*2] = {0};

    float Af = g[_bid][np][nr];
    for (i = 0; i<numb; ++i) {
        int ind = line[i],x,y,bid,bx,by;
        ind2xy(ind,x,y);
        bid = xy2bid(x,y,bx,by);
        if (_bid==bid)
            Af -= f[bid][bx][by]*weight[i];
    }

    for (i = 0; i<numb; ++i) {
        int ind = line[i],x,y,bid,bx,by;
        ind2xy(ind, x,y);
        bid = xy2bid(x,y,bx,by);
        if (_bid==bid) {
            float tmp = 0.;
            float lap = 0.;

            if (bx-1>=0)      tmp -= sqr(  v[bid][bx-1][by] )*( f[bid][bx][by] - f[bid][bx-1][by] );
            else              tmp -= sqr( pad_v[bid][0][by] )*( f[bid][bx][by] - pad_f[bid][0][by] );

            if (by+1<BLKSIZE) tmp += sqr( v[bid][bx][by] )*( f[bid][bx][by+1] -  f[bid][bx][by] );
            else              tmp += sqr( v[bid][bx][by] )*( pad_f[bid][1][bx] - f[bid][bx][by] );

            if (bx+1<BLKSIZE) tmp += sqr( v[bid][bx][by] )*( f[bid][bx+1][by] -  f[bid][bx][by] );
            else              tmp += sqr( v[bid][bx][by] )*( pad_f[bid][2][by] - f[bid][bx][by] );

            if (by-1>=0)      tmp -= sqr( v[bid][bx][by-1] )*(  f[bid][bx][by] - f[bid][bx][by-1] );
            else              tmp -= sqr( pad_v[bid][3][bx] )*( f[bid][bx][by] - pad_f[bid][3][bx] );

            if (bx-1>=0)      lap += f[bid][bx-1][by];
            else              lap += pad_f[bid][0][by];
            if (by+1<BLKSIZE) lap += f[bid][bx][by+1];
            else              lap += pad_f[bid][1][bx];
            if (bx+1<BLKSIZE) lap += f[bid][bx+1][by];
            else              lap += pad_f[bid][2][by];
            if (by-1>=0)      lap += f[bid][bx][by-1];
            else              lap += pad_f[bid][3][bx];
            lap -= 4*f[bid][bx][by];

            d[i] = Af*weight[i]+ALPHA*(tmp+sqr(EPSILON)*lap);
        }
    }
    for (i = 0; i<numb; ++i) {
        int ind = line[i],x,y,bid,bx,by;
        ind2xy(ind, x,y);
        bid = xy2bid(x,y,bx,by);
        if (_bid==bid) {
            f[bid][bx][by] += lambda*d[i];
            window(f[bid][bx][by],0,255);
        }
    }
}

void minEDGE(int _bid,int np,int nr,int *line,float *weight,int numb) {
    int i;

    float d[IMGSIZE*2] = {0};

    for (i = 0; i<numb; ++i) {
        int ind = line[i],x,y,bid,bx,by;
        ind2xy(ind, x,y);
        bid = xy2bid(x,y,bx,by);

        if (_bid==bid) {
            float a = 0.;
            float b = 0.;
            float c = 0.;

            if (bx-1>=0)      a += sqr(f[bid][bx][by] - f[bid][bx-1][by]);
            else              a += sqr(f[bid][bx][by] - pad_f[bid][0][by]);

            if (by-1>=0)      a += sqr(f[bid][bx][by] - f[bid][bx][by-1]);
            else              a += sqr(f[bid][bx][by] - pad_f[bid][3][bx]);

            a *= v[bid][bx][by];

            b = v[bid][bx][by]-1;

            if (bx-1>=0)      c += v[bid][bx-1][by];
            else              c += pad_v[bid][0][by];
            if (by+1<BLKSIZE) c += v[bid][bx][by+1];
            else              c += pad_v[bid][1][bx];
            if (bx+1<BLKSIZE) c += v[bid][bx+1][by];
            else              c += pad_v[bid][2][by];
            if (by-1>=0)      c += v[bid][bx][by-1];
            else              c += pad_v[bid][3][bx];
            c -= 4*v[bid][bx][by];

            d[i] = -ALPHA*a-BETA/(4*EPSILON)*b+BETA*EPSILON*c;
        }
    }

    for (i = 0; i<numb; ++i) {
        int ind = line[i],x,y,bid,bx,by;
        ind2xy(ind, x,y);
        bid = xy2bid(x,y,bx,by);

        if (_bid==bid) {
            v[bid][bx][by] += lambda*d[i];
            window(v[bid][bx][by],0,1);
        }
    }
}

void sumG(int bid) {
    for (int np = 0; np<NPROJ; ++np)
        for (int nr = 0; nr<NRAY; ++nr) {
            g_sum[np][nr] += g_dot[bid][np][nr];
        }
}

void updateG(int bid) {
    for (int np = 0; np<NPROJ; ++np)
        for (int nr = 0; nr<NRAY; ++nr) {
            if (fabs(g_dot[bid][np][nr])<F_TOL && fabs(g_sum[np][nr])<F_TOL) g[bid][np][nr] = _g[np][nr];
            else g[bid][np][nr] = (g_dot[bid][np][nr]/(g_sum[np][nr]+0.00005))*_g[np][nr];
        }
}

void padding(int bid) {
    for (int k = 0; k<NDIR; ++k) {
        for (int i = 0; i<BLKSIZE; ++i) {
            pad_f[bid][k][i] = 0;
            pad_v[bid][k][i] = 1;
        }
    }
    for (int i = 0; i<BLKSIZE; ++i) {
        if (adj[bid][0]!=-1) {
            pad_f[bid][0][i] = f[ adj[bid][0] ][BLKSIZE-1][i];
            pad_v[bid][0][i] = v[ adj[bid][0] ][BLKSIZE-1][i];
        }
        if (adj[bid][1]!=-1) {
            pad_f[bid][1][i] = f[ adj[bid][1] ][i][    0    ];
            pad_v[bid][1][i] = v[ adj[bid][1] ][i][    0    ];
        }
        if (adj[bid][2]!=-1) {
            pad_f[bid][2][i] = f[ adj[bid][2] ][    0    ][i];
            pad_v[bid][2][i] = v[ adj[bid][2] ][    0    ][i];
        }
        if (adj[bid][3]!=-1) {
            pad_f[bid][3][i] = f[ adj[bid][3] ][i][BLKSIZE-1];
            pad_v[bid][3][i] = v[ adj[bid][3] ][i][BLKSIZE-1];
        }
    }
}

bool cmp_t(int x,int y) {
    int ax = x/NRAY, ay = x%NRAY;
    int bx = y/NRAY, by = y%NRAY;
    
    float sa = 0, sb = 0;
    for (int k = 0; k<NBLK; ++k) {
        sa += w_blk[ax][ay][k];
        sb += w_blk[bx][by][k];
    }
    sa /= NBLK; sb /= NBLK;
    float da = 0, db = 0;
    for (int k = 0; k<NBLK; ++k) {
        da += sqr(w_blk[ax][ay][k]-sa);
        db += sqr(w_blk[bx][by][k]-sb);
    }
    
    return da < db;
}

void BlockRecon() {
    int i;

    // initialization
    memset(f, 0, sizeof(float)*NBLK*BLKSIZE*BLKSIZE);
    memset(g, 0, sizeof(float)*NBLK*NPROJ*NRAY);

    for (int k = 0; k<NBLK; ++k)
        for (int i = 0; i<BLKSIZE; ++i)
            for (int j = 0; j<BLKSIZE; ++j)
                v[k][i][j] = 1.0;

    lambda = 0.001;
    
    ofstream fou("g_MSE.log");
    ofstream dou("g_delta.log");

    int t[NPROJ*NRAY],total = NPROJ*NRAY;
    for (int k = 0; k<total; ++k) {
        t[k] = k;
    }
    
    //sort(t,t+total, cmp_t);
    
    for (int k = 1; k<total; ++k) {
        int tmp = rand()%k;
        swap(t[tmp],t[k]);
    }
    

    printf("Begin MSbeam minimization ...\n");
    for (i = 1; i<=ALL_ITER; ++i) {
        printf("Iteration %d ...\n", i);
        printf("lambda = %f\n",lambda);
        fflush(stdout);
        
        g_relax = (0.618)*(i-1);

        // if (i>5) {
        // for (int extra = 0; extra<2; extra++) {
        // for (int k = 0; k<NBLK; ++k)
        //     padding(k);
        
        // for (int loop = 0; loop<total; ++loop) {
        //     int np = t[loop]/NRAY, nr = t[loop]%NRAY;
        //     int line[IMGSIZE*2];
        //     float weight[IMGSIZE*2];
        //     int numb;
        //     float snorm;

        //     wray(np, nr, line, weight, &numb, &snorm);

        //     for (int k = 0; k<NBLK; ++k) {
        //         minIMAGE(k,np,nr,line,weight,numb,snorm);
        //         minEDGE(k,np,nr,line,weight,numb);
        //     }
        // }
        // }
        // i++;
        // }
        // else {
        for (int k = 0; k<NBLK; ++k)
            padding(k);
        
        for (int loop = 0; loop<total; ++loop) {
            int np = t[loop]/NRAY, nr = t[loop]%NRAY;
            int line[IMGSIZE*2];
            float weight[IMGSIZE*2];
            int numb;
            float snorm;

            wray(np, nr, line, weight, &numb, &snorm);

            for (int k = 0; k<NBLK; ++k) {
                minIMAGE(k,np,nr,line,weight,numb,snorm);
                minEDGE(k,np,nr,line,weight,numb);
            }
        }
        // }

        for (int k = 0; k<NBLK; ++k)
            A_blk(g_dot[k], k, f[k]);

        memset( g_sum, 0, sizeof(g_sum) );
       
        for (int k = 0; k<NBLK; ++k)
            sumG(k);
        for (int k = 0; k<NBLK; ++k)
            updateG(k);

        for (int k = 0; k<NBLK; ++k) {
            float delta = 0.;
            for (int np = 0; np<NPROJ; ++np)
                for (int nr = 0; nr<NRAY; ++nr) {
                    delta += sqr(global_g[k][np][nr]-g[k][np][nr]);
                }
            delta /= NPROJ*NRAY;
            printf("Block %d's g has a MSE of %f\n",k+1,delta);
            fou<<delta<<'\t';
            fflush(stdout);
        }
        fou<<endl;
        for (int p = 90; p<=90; ++p) {
            for (int r = 384-1; r<=384+1; ++r) {
                printf("(p=%d, r=%d)\n",p,r);
                cout<<"g - global_g: ";
                for (int k = 0; k<NBLK; ++k) {
                    cout<<g[k][p][r]-global_g[k][p][r]<<"\t\t";
                    dou<<g[k][p][r]-global_g[k][p][r]<<'\t';
                }
                cout<<endl;
            }
        }
        dou<<endl;
        
        //compute_AT(i);
        
        //lambda = lambda/(1 + 500 * lambda);

        stringstream sid_g;
        sid_g<<"g_iter_"<<i<<".dat";
        cout<<"Writing G to "<<sid_g.str()<<endl;
        write_data(g_sum, sid_g.str().c_str());

        stringstream sid_f;
        sid_f<<"f_iter_"<<i<<".dat";
        cout<<"Writing F to "<<sid_f.str()<<endl;
        unblock(f, _f);
        write_file(_f, sid_f.str().c_str());

        stringstream sid_v;
        sid_v<<"v_iter_"<<i<<".dat";
        cout<<"Writing V to "<<sid_v.str()<<endl;
        unblock(v, _v);
        write_file(_v, sid_v.str().c_str());
    }
    
    fou.close();
    dou.close();
    
    printf("MSbeam minimization done.\n");
}

int main(int argc,char **argv) {
    srand(time(0));
    
    init_vars();

    if (argc>1) read_file(_f,argv[1]);
    else        read_file(_f,PHANTOM_FILE);
    /*
    for (int i = 0; i<IMGSIZE; ++i)
        for (int j = 0; j<IMGSIZE; ++j) {
            int bx,by;
            int bid = xy2bid(i,j,bx,by);
            _f[i][j] = bid;
        }
    */

    normalize(_f);

    write_file(_f, "std_f.dat");

    A(_g, _f);
    write_data(_g,"std_g.dat");

    block(_f, f);
    for (int k = 0; k<NBLK; ++k)
        A_blk(global_g[k], k, f[k]);

    BlockRecon();

    unblock(f, _f);
    unblock(v, _v);

    write_file(_f, "f.dat");
    write_file(_v, "v.dat");

    return 0;
}

void init_vars() {
    memset(cnt_ray,0,sizeof(cnt_ray));
    memset(w_blk, 0, sizeof(w_blk));
    
    for (int np = 0; np<NPROJ; ++np) {
        for (int nr = 0; nr<NRAY; ++nr) {
            int line[IMGSIZE*2];
            float weight[IMGSIZE*2];
            int numb;
            float snorm;
            wray(np, nr, line, weight, &numb, &snorm);
            
            int cc[NBLK] = {0};
            
            for (int i = 0; i<numb; ++i) {
                int ind = line[i],x,y,bid,bx,by;
                ind2xy(ind, x,y);
                bid = xy2bid(x,y,bx,by);
                ++cc[bid];
            }
            
            float sum = 0;
            for (int i = 0; i<NBLK; ++i) {
                if (cc[i]>0) ++cnt_ray[np][nr];
                sum += cc[i];
            }
            for (int i = 0; i<NBLK; ++i) {
                w_blk[np][nr][i] = cc[i]/(sum+0.00005);
            }
        }
    }

    for (int i = 0; i<DIMSIZE; ++i)
        for (int j = 0; j<DIMSIZE; ++j) {
            for (int k = 0; k<NDIR; ++k) {
                adj[i*DIMSIZE+j][k] = -1;
                int x = i+DX[k], y = j+DY[k];
                if (0<=x && x<DIMSIZE && 0<=y && y<DIMSIZE) {
                    adj[i*DIMSIZE+j][k] = x*DIMSIZE+y;
                }
            }
        }
}
