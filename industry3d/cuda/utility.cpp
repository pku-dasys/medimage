
bool equals_draw(float x,int y) {
    return fabs(x-y)<1.25;
}

float sqr(float x) { return x*x; }

void rotate_axis_3d(float &x,float &y,float &z,float ux,float uy,float uz,float theta) {
    float c = 1/sqrt(sqr(ux)+sqr(uy)+sqr(uz));
    ux *= c;
    uy *= c;
    uz *= c;

    float X[3][1] = {{x},{y},{z}};
    float costheta = cos(theta), sintheta = sin(theta);
    float R[3][3] = {
        {costheta+sqr(ux)*(1-costheta), ux*uy*(1-costheta)-uz*sintheta, ux*uz*(1-costheta)+uy*sintheta},
        {uy*ux*(1-costheta)+uz*sintheta, costheta+sqr(uy)*(1-costheta), uy*uz*(1-costheta)-ux*sintheta},
        {uz*ux*(1-costheta)-uy*sintheta, uz*uy*(1-costheta)+ux*sintheta, costheta+sqr(uz)*(1-costheta)}
    };

    float d[3][1] = {};
    for (int k = 0; k<3; ++k)
        for (int i = 0; i<3; ++i)
            for (int j = 0; j<1; ++j)
                d[i][j] += R[i][k]*X[k][j];
    
    x = d[0][0];
    y = d[1][0];
    z = d[2][0];
}

void rotate_2d(float &x,float &y,float theta) {
    float X[2][1] = {{x},{y}};
    float costheta = cos(theta), sintheta = sin(theta);
    float R[2][2] = {
        {costheta, -sintheta},
        {sintheta, costheta}
    };

    float d[2][1] = {};
    for (int k = 0; k<2; ++k)
        for (int i = 0; i<2; ++i)
            for (int j = 0; j<1; ++j)
                d[i][j] += R[i][k]*X[k][j];
    
    x = d[0][0];
    y = d[1][0];
}
