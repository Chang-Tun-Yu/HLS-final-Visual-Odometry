#ifndef MYMATCH_H
#define MYMATCH_H

#include <iostream>
#include <stdint.h>

#include "size.h"
#include "delaunator.hpp"
#include "matcher.h"

using namespace std;

#define ABS(a, b) ((a > b)? a-b:b-a);

void myCreateIndexVector (int32_t m[MAX_FEATURE_ARRAY_SZIE],int32_t& n,int32_t k[BIN_NUM][MAX_FP_IN_BIN], int32_t k_num[BIN_NUM],const int32_t &u_bin_num,const int32_t &v_bin_num, int32_t m_num[BIN_NUM]);

inline void myFindMatch (int32_t m1[MAX_FEATURE_ARRAY_SZIE],const int32_t &i1,int32_t m2[MAX_FEATURE_ARRAY_SZIE],const int32_t &step_size,int32_t k2[BIN_NUM][MAX_FP_IN_BIN], int32_t k2_num[BIN_NUM], 
                                const int32_t &u_bin_num,const int32_t &v_bin_num,const int32_t &stat_bin,
                                int32_t& min_ind,int32_t stage);

void myMatching (int32_t m1p[MAX_FEATURE_ARRAY_SZIE],int32_t m1c[MAX_FEATURE_ARRAY_SZIE], int32_t n1p[BIN_NUM],int32_t n1c[BIN_NUM],Matcher::p_match p_matched[POINT_L], int32_t& p_matched_num);

#endif