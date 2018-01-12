#include "testcase.h"

#include "utility.h"
#include "tracing.h"
#include "ct3d.h"

#include <cmath>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/optional.hpp>
#include <boost/foreach.hpp>

using namespace std;
using namespace boost::property_tree;

void wrapper_hip(const Parameter &args, float *img, ushort *prj, int a,int x,int y, int *bf,int *br) {
    float srcX,srcY,srcZ;
    float dstX,dstY,dstZ;

    if(args.BEAM == "Parallel")
        parallel(args,a,x,y,
                srcX,srcY,srcZ,dstX,dstY,dstZ);
    else if(args.BEAM == "Cone")
        cone(args,a,x,y,
                srcX,srcY,srcZ,dstX,dstY,dstZ);
    else
        assert(false);

    int64_t *ind = new int64_t[args.MAX_RAYLEN];
    float *wgt = new float[args.MAX_RAYLEN];
    int numb;
    
    forward_proj(args,
                 srcX, srcY, srcZ,
                 dstX, dstY, dstZ,
                 ind, wgt, numb);

    float sum = 0;

    for (int i = 0; i<numb; ++i) {
        int64_t idx = ind[i];
        sum += img[idx]*wgt[i];
        idx /= args.NX*args.NY;
        if (bf[idx]>x) bf[idx] = x;
        if (br[idx]<x) br[idx] = x;
    }

    prj[x*args.NDY+y] = ushort(sum);
    delete[] ind;
    delete[] wgt;
}

static float put_img(ptree& contents, float* img, const Parameter &args)
{
    BOOST_FOREACH(ptree::value_type& pairs, contents)
    {
        int sx=0, dx=args.NX,
            sy=0, dy=args.NY,
            sz=0, dz=args.NZ;

        ptree& content = pairs.second;

        string type;
        float cx, cy, cz, rho;
        boost::optional<ptree&> constraint_op=content.get_child_optional("constraint");
        if(constraint_op)
        {
            ptree& constraint = constraint_op.get();
            boost::optional<ptree&> res_op;
            res_op=constraint.get_child_optional("x");
            if(res_op)
            {
                ptree& res = res_op.get();
                if(res.get<string>("rel") == "GT") sx = max(sx, (int)(res.get<float>("num")/args.SAMPLESIZE));
                if(res.get<string>("rel") == "LT") dx = min(dx, (int)(res.get<float>("num")/args.SAMPLESIZE));
            }
            res_op=constraint.get_child_optional("y");
            if(res_op)
            {
                ptree& res = res_op.get();
                if(res.get<string>("rel") == "GT") sy = max(sy, (int)(res.get<float>("num")/args.SAMPLESIZE));
                if(res.get<string>("rel") == "LT") dy = min(dy, (int)(res.get<float>("num")/args.SAMPLESIZE));
            }
            res_op=constraint.get_child_optional("z");
            if(res_op)
            {
                ptree& res = res_op.get();
                if(res.get<string>("rel") == "GT") sz = max(sz, (int)(res.get<float>("num")/args.SAMPLESIZE));
                if(res.get<string>("rel") == "LT") dz = min(dz, (int)(res.get<float>("num")/args.SAMPLESIZE));
            }
        }
        try
        {
            type=content.get<string>("Type");
            cx=content.get<float>("x");
            cy=content.get<float>("y");
            cz=content.get<float>("z");
            rho=content.get<float>("rho");
        }
        catch (boost::property_tree::ptree_error &e) {
            printf("JSON file corrupted.\n");
            exit(1);
        }
        if(type == "Sphere")
        {
            float r=content.get<float>("r");
            sx = max(sx, (int)((cx-r)/args.SAMPLESIZE));
            dx = min(dx, (int)((cx+r)/args.SAMPLESIZE));
            sy = max(sy, (int)((cy-r)/args.SAMPLESIZE));
            dy = min(dy, (int)((cy+r)/args.SAMPLESIZE));
            sz = max(sz, (int)((cz-r)/args.SAMPLESIZE));
            dz = min(dz, (int)((cz+r)/args.SAMPLESIZE));
            for(int i = sx; i < dx; i++)
                for(int j = sy; j < dy; j++)
                    for(int k = sz; k < dz; k++)
                        if(sqr(args.SAMPLESIZE*i-cx)+sqr(args.SAMPLESIZE*j-cy)+sqr(args.SAMPLESIZE*k-cz) <= sqr(r))
                            img[(k*args.NX+i)*args.NY+j] = rho;
        }
        else if(type == "Box")
        {
            float ddx=content.get<float>("dx")/2,
                  ddy=content.get<float>("dy")/2,
                  ddz=content.get<float>("dz")/2;
            sx = max(sx, (int)((cx-ddx)/args.SAMPLESIZE));
            dx = min(dx, (int)((cx+ddx)/args.SAMPLESIZE));
            sy = max(sy, (int)((cy-ddy)/args.SAMPLESIZE));
            dy = min(dy, (int)((cy+ddy)/args.SAMPLESIZE));
            sz = max(sz, (int)((cz-ddz)/args.SAMPLESIZE));
            dz = min(dz, (int)((cz+ddz)/args.SAMPLESIZE));
            for(int i = sx; i < dx; i++)
                for(int j = sy; j < dy; j++)
                    for(int k = sz; k < dz; k++)
                        img[(k*args.NX+i)*args.NY+j] = rho;
        }
        else if(type == "Cylinder_x")
        {
            float r=content.get<float>("r"),
                  l=content.get<float>("l")/2;
            sx = max(sx, (int)((cx-l)/args.SAMPLESIZE));
            dx = min(dx, (int)((cx+l)/args.SAMPLESIZE));
            sy = max(sy, (int)((cy-r)/args.SAMPLESIZE));
            dy = min(dy, (int)((cy+r)/args.SAMPLESIZE));
            sz = max(sz, (int)((cz-r)/args.SAMPLESIZE));
            dz = min(dz, (int)((cz+r)/args.SAMPLESIZE));
            for(int i = sx; i < dx; i++)
                for(int j = sy; j < dy; j++)
                    for(int k = sz; k < dz; k++)
                        if(sqr(args.SAMPLESIZE*j-cy)+sqr(args.SAMPLESIZE*k-cz) <= sqr(r))
                            img[(k*args.NX+i)*args.NY+j] = rho;
        }
        else if(type == "Cylinder_y")
        {
            float r=content.get<float>("r"),
                  l=content.get<float>("l");
            sx = max(sx, (int)((cx-r)/args.SAMPLESIZE));
            dx = min(dx, (int)((cx+r)/args.SAMPLESIZE));
            sy = max(sy, (int)((cy-l)/args.SAMPLESIZE));
            dy = min(dy, (int)((cy+l)/args.SAMPLESIZE));
            sz = max(sz, (int)((cz-r)/args.SAMPLESIZE));
            dz = min(dz, (int)((cz+r)/args.SAMPLESIZE));
            for(int i = sx; i < dx; i++)
                for(int j = sy; j < dy; j++)
                    for(int k = sz; k < dz; k++)
                        if(sqr(args.SAMPLESIZE*i-cx)+sqr(args.SAMPLESIZE*k-cz) <= sqr(r))
                            img[(k*args.NX+i)*args.NY+j] = rho;
        }
        else if(type == "Cylinder_z")
        {
            float r=content.get<float>("r"),
                  l=content.get<float>("l");
            sx = max(sx, (int)((cx-r)/args.SAMPLESIZE));
            dx = min(dx, (int)((cx+r)/args.SAMPLESIZE));
            sy = max(sy, (int)((cy-r)/args.SAMPLESIZE));
            dy = min(dy, (int)((cy+r)/args.SAMPLESIZE));
            sz = max(sz, (int)((cz-l)/args.SAMPLESIZE));
            dz = min(dz, (int)((cz+l)/args.SAMPLESIZE));
            for(int i = sx; i < dx; i++)
                for(int j = sy; j < dy; j++)
                    for(int k = sz; k < dz; k++)
                        if(sqr(args.SAMPLESIZE*j-cy)+sqr(args.SAMPLESIZE*i-cx) <= sqr(r))
                            img[(k*args.NX+i)*args.NY+j] = rho;
        }
        else if(type == "Ellipt_Cyl_x")
        {
            float l  =content.get<float>("l")/2,
                  ddy=content.get<float>("dy"),
                  ddz=content.get<float>("dz");
            sx = max(sx, (int)((cx-l)/args.SAMPLESIZE));
            dx = min(dx, (int)((cx+l)/args.SAMPLESIZE));
            sy = max(sy, (int)((cy-ddy)/args.SAMPLESIZE));
            dy = min(dy, (int)((cy+ddy)/args.SAMPLESIZE));
            sz = max(sz, (int)((cz-ddz)/args.SAMPLESIZE));
            dz = min(dz, (int)((cz+ddz)/args.SAMPLESIZE));
            for(int i = sx; i < dx; i++)
                for(int j = sy; j < dy; j++)
                    for(int k = sz; k < dz; k++)
                        if(sqr((args.SAMPLESIZE*j-cy)/ddy)+sqr((args.SAMPLESIZE*k-cz)/ddz) <= 1.0)
                            img[(k*args.NX+i)*args.NY+j] = rho;
        }
        else if(type == "Ellipt_Cyl_y")
        {
            float l  =content.get<float>("l")/2,
                  ddx=content.get<float>("dx"),
                  ddz=content.get<float>("dz");
            sx = max(sx, (int)((cx-ddx)/args.SAMPLESIZE));
            dx = min(dx, (int)((cx+ddx)/args.SAMPLESIZE));
            sy = max(sy, (int)((cy-l)/args.SAMPLESIZE));
            dy = min(dy, (int)((cy+l)/args.SAMPLESIZE));
            sz = max(sz, (int)((cz-ddz)/args.SAMPLESIZE));
            dz = min(dz, (int)((cz+ddz)/args.SAMPLESIZE));
            for(int i = sx; i < dx; i++)
                for(int j = sy; j < dy; j++)
                    for(int k = sz; k < dz; k++)
                        if(sqr((args.SAMPLESIZE*i-cx)/ddx)+sqr((args.SAMPLESIZE*k-cz)/ddz) <= 1.0)
                            img[(k*args.NX+i)*args.NY+j] = rho;
        }
        else if(type == "Ellipt_Cyl_z")
        {
            float l  =content.get<float>("l")/2,
                  ddx=content.get<float>("dx"),
                  ddy=content.get<float>("dy");
            sx = max(sx, (int)((cx-ddx)/args.SAMPLESIZE));
            dx = min(dx, (int)((cx+ddx)/args.SAMPLESIZE));
            sy = max(sy, (int)((cy-ddy)/args.SAMPLESIZE));
            dy = min(dy, (int)((cy+ddy)/args.SAMPLESIZE));
            sz = max(sz, (int)((cz-l)/args.SAMPLESIZE));
            dz = min(dz, (int)((cz+l)/args.SAMPLESIZE));
            for(int i = sx; i < dx; i++)
                for(int j = sy; j < dy; j++)
                    for(int k = sz; k < dz; k++)
                        if(sqr((args.SAMPLESIZE*i-cx)/ddx)+sqr((args.SAMPLESIZE*j-cy)/ddy) <= 1.0)
                            img[(k*args.NX+i)*args.NY+j] = rho;
        }
        else if(type == "Ellipsoid")
        {
            float ddx=content.get<float>("dx"),
                  ddy=content.get<float>("dy"),
                  ddz=content.get<float>("dz");
            sx = max(sx, (int)((cx-ddx)/args.SAMPLESIZE));
            dx = min(dx, (int)((cx+ddx)/args.SAMPLESIZE));
            sy = max(sy, (int)((cy-ddy)/args.SAMPLESIZE));
            dy = min(dy, (int)((cy+ddy)/args.SAMPLESIZE));
            sz = max(sz, (int)((cz-ddz)/args.SAMPLESIZE));
            dz = min(dz, (int)((cz+ddz)/args.SAMPLESIZE));
            for(int i = sx; i < dx; i++)
                for(int j = sy; j < dy; j++)
                    for(int k = sz; k < dz; k++)
                        if(sqr((args.SAMPLESIZE*i-cx)/ddx)
                          +sqr((args.SAMPLESIZE*j-cy)/ddy)
                          +sqr((args.SAMPLESIZE*k-cz)/ddz) <= 1.0)
                            img[(k*args.NX+i)*args.NY+j] = rho;
        }
        else if(type == "Ellipsoid_free")
        {
            float ddx=content.get<float>("dx"),
                  ddy=content.get<float>("dy"),
                  ddz=content.get<float>("dz");
            float axis[3][3] = {{1,0,0},{0,1,0},{0,0,1}},
                  inv[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
            float ma = max(max(ddx,ddy),ddz)*1.05;
            sx = max(sx, (int)((cx-ma)/args.SAMPLESIZE));
            dx = min(dx, (int)((cx+ma)/args.SAMPLESIZE));
            sy = max(sy, (int)((cy-ma)/args.SAMPLESIZE));
            dy = min(dy, (int)((cy+ma)/args.SAMPLESIZE));
            sz = max(sz, (int)((cz-ma)/args.SAMPLESIZE));
            dz = min(dz, (int)((cz+ma)/args.SAMPLESIZE));
            boost::optional<ptree&> res_op;
            res_op = content.get_child_optional("a_x");
            if(res_op)
            {
                int i = 0;
                float norm = 0;
                ptree& res = res_op.get();
                BOOST_FOREACH(ptree::value_type& pairs, res)
                {
                    axis[i][0] = pairs.second.get_value<float>();
                    norm += sqr(axis[i++][0]);
                }
                norm = sqrt(norm);
                for(i = 0; i < 3; i++)
                    axis[i][0]/=norm;
            }
            res_op = content.get_child_optional("a_y");
            if(res_op)
            {
                int i = 0;
                float norm = 0;
                ptree& res = res_op.get();
                BOOST_FOREACH(ptree::value_type& pairs, res)
                {
                    axis[i][1] = pairs.second.get_value<float>();
                    norm += sqr(axis[i++][1]);
                }
                norm = sqrt(norm);
                for(i = 0; i < 3; i++)
                    axis[i][1]/=norm;
            }
            res_op = content.get_child_optional("a_z");
            if(res_op)
            {
                int i = 0;
                float norm = 0;
                ptree& res = res_op.get();
                BOOST_FOREACH(ptree::value_type& pairs, res)
                {
                    axis[i][2] = pairs.second.get_value<float>();
                    norm += sqr(axis[i++][2]);
                }
                norm = sqrt(norm);
                for(i = 0; i < 3; i++)
                    axis[i][2]/=norm;
            }
            for(int i = 0; i < 3; i++)
            {
                if(fabs(axis[i][i])<1e-5)
                {
                    int tar = i;
                    for(int j = i+1; j < 3; j++)
                    {
                        if(fabs(axis[j][i])>1e-5)
                            tar = j;
                    }
                    if(tar==i)
                    {
                        fprintf(stderr, "matrix singular\n");
                        exit(1);
                    }
                    for(int j = 0; j < 3; j++)
                    {
                        swap(axis[i][j], axis[tar][j]);
                        swap(inv[i][j], inv[tar][j]);
                    }
                }
                for(int j = 0; j < 3; j++)
                {
                    if(j == i) continue;
                    axis[i][j] /= axis[i][i];
                    inv[i][j] /= axis[i][i];
                }
                inv[i][i] /= axis[i][i];
                axis[i][i] = 1;
                for(int j = 0; j < 3; j++)
                {
                    if(j == i) continue;
                    for(int k = 0; k < 3; k++)
                    {
                        if(k == i) continue;
                        axis[j][k] -= axis[i][k] * axis[j][i];
                        inv[j][k] -= inv[i][k] * axis[j][i];
                    }
                    inv[j][i] -= inv[i][i] * axis[j][i];
                    axis[j][i] = 0;
                }
            }
            for(int i = sx; i < dx; i++)
                for(int j = sy; j < dy; j++)
                    for(int k = sz; k < dz; k++)
                    {
                        float x = (args.SAMPLESIZE*i-cx),
                              y = (args.SAMPLESIZE*j-cy),
                              z = (args.SAMPLESIZE*k-cz);
                        float xx = inv[0][0]*x+inv[0][1]*y+inv[0][2]*z,
                              yy = inv[1][0]*x+inv[1][1]*y+inv[1][2]*z,
                              zz = inv[2][0]*x+inv[2][1]*y+inv[2][2]*z;
                        if(sqr(xx/ddx)+sqr(yy/ddy)+sqr(zz/ddz) <= 1.0)
                            img[(k*args.NX+i)*args.NY+j] = rho;
                    }
        }
        else if(type == "Ellipt_Cyl")
        {
            float ddx=content.get<float>("dx"),
                  ddy=content.get<float>("dy"),
                  l  =content.get<float>("l")/2;
            float axis[3][3] = {{1,0,0},{0,1,0},{0,0,1}},
                  inv[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
            float ma = sqrt(sqr(max(ddx,ddy))+sqr(l));
            sx = max(sx, (int)((cx-ma)/args.SAMPLESIZE));
            dx = min(dx, (int)((cx+ma)/args.SAMPLESIZE));
            sy = max(sy, (int)((cy-ma)/args.SAMPLESIZE));
            dy = min(dy, (int)((cy+ma)/args.SAMPLESIZE));
            sz = max(sz, (int)((cz-ma)/args.SAMPLESIZE));
            dz = min(dz, (int)((cz+ma)/args.SAMPLESIZE));
            boost::optional<ptree&> res_op;
            res_op = content.get_child_optional("a_x");
            if(res_op)
            {
                int i = 0;
                float norm = 0;
                ptree& res = res_op.get();
                BOOST_FOREACH(ptree::value_type& pairs, res)
                {
                    axis[i][0] = pairs.second.get_value<float>();
                    norm += sqr(axis[i++][0]);
                }
                norm = sqrt(norm);
                for(i = 0; i < 3; i++)
                    axis[i][0]/=norm;
            }
            res_op = content.get_child_optional("a_y");
            if(res_op)
            {
                int i = 0;
                float norm = 0;
                ptree& res = res_op.get();
                BOOST_FOREACH(ptree::value_type& pairs, res)
                {
                    axis[i][1] = pairs.second.get_value<float>();
                    norm += sqr(axis[i++][1]);
                }
                norm = sqrt(norm);
                for(i = 0; i < 3; i++)
                    axis[i][1]/=norm;
            }
            {
                float norm = 0;
                for(int i = 0; i < 3; i++)
                {
                    axis[i][2] = axis[(i+1)%3][0]*axis[(i+2)%3][1] - axis[(i+2)%3][0]*axis[(i+1)%3][1];
                    norm += sqr(axis[i][2]);
                }
                norm = sqrt(norm);
                for(int i = 0; i < 3; i++)
                    axis[i][2]/=norm;
            }
            for(int i = 0; i < 3; i++)
            {
                if(fabs(axis[i][i])<1e-5)
                {
                    int tar = i;
                    for(int j = i+1; j < 3; j++)
                    {
                        if(fabs(axis[j][i])>1e-5)
                            tar = j;
                    }
                    if(tar==i)
                    {
                        fprintf(stderr, "matrix singular\n");
                        exit(1);
                    }
                    for(int j = 0; j < 3; j++)
                    {
                        swap(axis[i][j], axis[tar][j]);
                        swap(inv[i][j], inv[tar][j]);
                    }
                }
                for(int j = 0; j < 3; j++)
                {
                    if(j == i) continue;
                    axis[i][j] /= axis[i][i];
                    inv[i][j] /= axis[i][i];
                }
                inv[i][i] /= axis[i][i];
                axis[i][i] = 1;
                for(int j = 0; j < 3; j++)
                {
                    if(j == i) continue;
                    for(int k = 0; k < 3; k++)
                    {
                        if(k == i) continue;
                        axis[j][k] -= axis[i][k] * axis[j][i];
                        inv[j][k] -= inv[i][k] * axis[j][i];
                    }
                    inv[j][i] -= inv[i][i] * axis[j][i];
                    axis[j][i] = 0;
                }
            }
            for(int i = sx; i < dx; i++)
                for(int j = sy; j < dy; j++)
                    for(int k = sz; k < dz; k++)
                    {
                        float x = (args.SAMPLESIZE*i-cx),
                              y = (args.SAMPLESIZE*j-cy),
                              z = (args.SAMPLESIZE*k-cz);
                        float xx = inv[0][0]*x+inv[0][1]*y+inv[0][2]*z,
                              yy = inv[1][0]*x+inv[1][1]*y+inv[1][2]*z,
                              zz = inv[2][0]*x+inv[2][1]*y+inv[2][2]*z;
                        if(sqr(xx/ddx)+sqr(yy/ddy) <= 1.0 && fabs(zz) <= l)
                            img[(k*args.NX+i)*args.NY+j] = rho;
                    }
        }
        else
        {
            printf("unknown type\n");
            exit(1);
        }
    }
}



void hip(const Parameter &args, const char* hip_file_path) {

    int *back_front = new int[args.NZ];
    int *back_rear = new int[args.NZ];

    for (int i = 0; i<args.NZ; ++i)
        back_front[i] = args.NDX+1;
    for (int i = 0; i<args.NZ; ++i)
        back_rear[i] = -1;

    boost::filesystem::path dr_file(args.RAW_DATA_FILE);

    if (boost::filesystem::exists(dr_file))
        return;

    float *img = new float[args.NX*args.NY*args.NZ];
    ushort *prj = new ushort[args.NPROJ*args.NDX*args.NDY];

    ptree root;
    try {
        read_json(hip_file_path, root);
    }
    catch (boost::property_tree::ptree_error &e) {
        printf("Could not read from the JSON file.\n");
        cout << e.what() << endl;
        exit(1);
    }

    memset(img, 0, args.NX*args.NY*args.NZ*sizeof(float));
    put_img(root, img, args);

    for(int a = 0; a < args.NPROJ; a++) {
        for (int ndx = 0; ndx<args.NDX; ++ndx) {
            for (int ndy = 0; ndy<args.NDY; ++ndy) {
                wrapper_hip(args,img,prj+a*args.NDX*args.NDY,a,ndx,ndy, back_front, back_rear);
            }
        }
    }

    ofstream fou;
    fou.open(args.RAW_DATA_FILE, ios::binary);
    char empty[1024] = {};
    fou.write(empty,1024);
    int64_t size = args.NDX * args.NDY;
    for (int k = 0; k<args.NPROJ; ++k)
        fou.write((char*)(prj+k*args.NDX*args.NDY), sizeof(ushort) * size);
    fou.close();

    //fou.open(args.PRETRACING_FILE);
    //for (int i = 0; i<args.NZ; ++i)
    //    fou<<back_front[i]<<' '<<back_rear[i]<<endl;
    //fou.close();

    for(int i=0; i<args.NZ; i++)
    {
        char buf[32];
        sprintf(buf, "pretracing/i_%d", i);
        fou.open(buf);
        for(int j = 0; j<args.NX; j++)
        {
            for(int k = 0; k < args.NY; k++)
                fou << img[(i*args.NX+j)*args.NY+k] << ' ';
            fou << endl;
        }
        fou.close();
    }

    delete [] back_front;
    delete [] back_rear;

    delete [] img;
    delete [] prj;
}

