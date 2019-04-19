#ifndef _RAY_H_
#define _RAY_H_

#include "config.h"

#pragma SDS data access_pattern(wgt:SEQUENTIAL, \
                                f0:SEQUENTIAL, f1:SEQUENTIAL, f2:SEQUENTIAL, f3:SEQUENTIAL, f4:SEQUENTIAL,\
                                v0:SEQUENTIAL, v1:SEQUENTIAL, v2:SEQUENTIAL, v3:SEQUENTIAL, v4:SEQUENTIAL,\
                                commit_f:SEQUENTIAL, commit_v:SEQUENTIAL)
#pragma SDS data mem_attribute(wgt:PHYSICAL_CONTIGUOUS, \
                               f0:PHYSICAL_CONTIGUOUS, f1:PHYSICAL_CONTIGUOUS, f2:PHYSICAL_CONTIGUOUS, f3:PHYSICAL_CONTIGUOUS, f4:PHYSICAL_CONTIGUOUS,\
                               v0:PHYSICAL_CONTIGUOUS, v1:PHYSICAL_CONTIGUOUS, v2:PHYSICAL_CONTIGUOUS, v3:PHYSICAL_CONTIGUOUS, v4:PHYSICAL_CONTIGUOUS,\
                               commit_f:PHYSICAL_CONTIGUOUS, commit_v:PHYSICAL_CONTIGUOUS)
void ray(Data lambda_f, Data lambda_v, Data ALPHA, Data BETA, Data EPSILON,
        Data minus_g, Number numb, Weight wgt[MAX_ELE_RAY],
        Data f0[MAX_ELE_RAY], Data f1[MAX_ELE_RAY], Data f2[MAX_ELE_RAY], Data f3[MAX_ELE_RAY], Data f4[MAX_ELE_RAY],
        Data v0[MAX_ELE_RAY], Data v1[MAX_ELE_RAY], Data v2[MAX_ELE_RAY], Data v3[MAX_ELE_RAY], Data v4[MAX_ELE_RAY],
        Data commit_f[MAX_ELE_RAY], Data commit_v[MAX_ELE_RAY]);

#endif
