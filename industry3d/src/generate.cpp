#include <cmath>
#include <fstream>
#include <iostream>
#include <cstring>

#include "tracing.h"
#include "ct3d.h"
#include "utility.h"

#include "testcase.h"

using namespace std;



void generate_parallel_circle() {
    char **json = new char*[2];
    json[0] = new char[256]{};
    json[1] = new char[256]{};

    char parallel_circle_json[] = "test_360x128x128_128_parallel_circle.json";

    strcpy(json[1], parallel_circle_json);
    Parameter args;
    args.parse_config(2,json);
    args.derive();

    parallel_circle(args);
}
void generate_parallel_sphere()
{
    char **json = new char*[2];
    json[0] = new char[256]{};
    json[1] = new char[256]{};

    char parallel_sphere_json[] = "test_360x128x128_128_parallel_sphere.json";

    strcpy(json[1], parallel_sphere_json);
    Parameter args;
    args.parse_config(2,json);
    args.derive();

    parallel_sphere(args);
}
void generate_parallel_ellipse()
{
    char **json = new char*[2];
    json[0] = new char[256]{};
    json[1] = new char[256]{};

    char parallel_ellipse_json[] = "test_360x128x128_128_parallel_ellipse.json";

    strcpy(json[1], parallel_ellipse_json);
    Parameter args;
    args.parse_config(2,json);
    args.derive();

    parallel_ellipse(args);
}
void generate_parallel_taiji()
{
    char **json = new char*[2];
    json[0] = new char[256]{};
    json[1] = new char[256]{};

    char parallel_taiji_json[] = "test_360x128x128_128_parallel_taiji.json";

    strcpy(json[1], parallel_taiji_json);
    Parameter args;
    args.parse_config(2,json);
    args.derive();

    taiji(args);
}
void generate_cone_taiji()
{
    char **json = new char*[2];
    json[0] = new char[256]{};
    json[1] = new char[256]{};

    char cone_taiji_json[] = "test_360x256x256_128_cone_taiji.json";

    strcpy(json[1], cone_taiji_json);
    Parameter args;
    args.parse_config(2,json);
    args.derive();

    taiji(args);
}
void generate_cone_90_taiji()
{
    char **json = new char*[2];
    json[0] = new char[256]{};
    json[1] = new char[256]{};

    char cone_taiji_json[] = "test_90x256x256_128_cone_taiji.json";

    strcpy(json[1], cone_taiji_json);
    Parameter args;
    args.parse_config(2,json);
    args.derive();

    taiji(args);
}
void generate_cone_180_taiji()
{
    char **json = new char*[2];
    json[0] = new char[256]{};
    json[1] = new char[256]{};

    char cone_taiji_json[] = "test_180x256x256_128_cone_taiji.json";

    strcpy(json[1], cone_taiji_json);
    Parameter args;
    args.parse_config(2,json);
    args.derive();

    taiji(args);
}
void generate_cone_shepp_logan()
{
    char **json = new char*[2];
    json[0] = new char[256]{};
    json[1] = new char[256]{};

    char cone_shepp_logan_json[] = "test_360x256x256_128_cone_shepp_logan.json";

    strcpy(json[1], cone_shepp_logan_json);
    Parameter args;
    args.parse_config(2,json);
    args.derive();

    shepp_logan(args);
}
void generate_cone_90_shepp_logan()
{
    char **json = new char*[2];
    json[0] = new char[256]{};
    json[1] = new char[256]{};

    char cone_shepp_logan_json[] = "test_90x256x256_128_cone_shepp_logan.json";

    strcpy(json[1], cone_shepp_logan_json);
    Parameter args;
    args.parse_config(2,json);
    args.derive();

    shepp_logan(args);
}

void generate_cone_90_shepp_logan_straight()
{
    char **json = new char*[2];
    json[0] = new char[256]{};
    json[1] = new char[256]{};

    char cone_shepp_logan_json[] = "test_90x256x256_128_cone_shepp_logan_str.json";

    strcpy(json[1], cone_shepp_logan_json);
    Parameter args;
    args.parse_config(2,json);
    args.derive();

    shepp_logan_straight(args);
}

void generate_cone_180_shepp_logan_straight()
{
    char **json = new char*[2];
    json[0] = new char[256]{};
    json[1] = new char[256]{};

    char cone_shepp_logan_json[] = "test_180x256x256_128_cone_shepp_logan_str.json";

    strcpy(json[1], cone_shepp_logan_json);
    Parameter args;
    args.parse_config(2,json);
    args.derive();

    shepp_logan_straight(args);
}
void generate_parallel_jaw()
{
    char **json = new char*[2];
    json[0] = new char[256]{};
    json[1] = new char[256]{};

    char parallel_jaw_json[] = "test_360x180x180_120x180x100_parallel_jaw.json";

    strcpy(json[1], parallel_jaw_json);
    Parameter args;
    args.parse_config(2,json);
    args.derive();

    jaw(args);
}
void generate_parallel_hip()
{
    char **json = new char*[2];
    json[0] = new char[256]{};
    json[1] = new char[256]{};

    char parallel_hip_json[] = "test_360x150x350_350x150x150_parallel_hip.json";

    strcpy(json[1], parallel_hip_json);
    Parameter args;
    args.parse_config(2,json);
    args.derive();

    hip(args);
}
void generate_parallel_head()
{
    char **json = new char*[2];
    json[0] = new char[256]{};
    json[1] = new char[256]{};

    char parallel_head_json[] = "test_360x256x256_256_parallel_head.json";

    strcpy(json[1], parallel_head_json);
    Parameter args;
    args.parse_config(2,json);
    args.derive();

    head(args);
}
void generate_parallel_abdomen()
{
    char **json = new char*[2];
    json[0] = new char[256]{};
    json[1] = new char[256]{};

    char parallel_abdomen_json[] = "test_360x256x256_256_parallel_abdomen.json";

    strcpy(json[1], parallel_abdomen_json);
    Parameter args;
    args.parse_config(2,json);
    args.derive();

    abdomen(args);
}
int main(int argc, char** argv) {

    //generate_parallel_circle();
    //generate_parallel_sphere();
    //generate_parallel_ellipse();
    //generate_parallel_taiji();
    //generate_cone_taiji();
    //generate_cone_90_taiji();
    //generate_cone_180_taiji();
    //generate_cone_shepp_logan();
    //generate_cone_90_shepp_logan();
    //generate_cone_180_shepp_logan_straight();
    //generate_cone_90_shepp_logan_straight();
    //generate_parallel_jaw();
    //generate_parallel_hip();
    //generate_parallel_head();
    generate_parallel_abdomen();
    return 0;
}
