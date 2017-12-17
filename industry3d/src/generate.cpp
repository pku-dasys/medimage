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

int main(int argc, char** argv) {

    //generate_parallel_circle();
    //generate_parallel_sphere();
    //generate_parallel_ellipse();
    //generate_parallel_taiji();
    generate_cone_taiji();
    return 0;
}
