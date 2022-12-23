#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <emmintrin.h>
#include <algorithm>
#include <vector>

#include "matrix.h"
#include "size.h"

using namespace std;

#define ABS(a, b) ((a > b)? a-b:b-a);

void Matcher::myCreateIndexVector (int32_t* m,int32_t n,int32_t k[BIN_NUM][MAX_FP_IN_BIN], int32_t k_num[BIN_NUM],const int32_t &u_bin_num,const int32_t &v_bin_num) {

  // descriptor step size
  int32_t step_size = sizeof(Matcher::maximum)/sizeof(int32_t);
  int32_t bin_idx;
  // for all points do
  for (int32_t i=0; i<n; i++) {
    // cout << i << endl;
    // extract coordinates and class
    int32_t u = *(m+step_size*i+0); // u-coordinate
    int32_t v = *(m+step_size*i+1); // v-coordinate
    int32_t c = *(m+step_size*i+3); // class
    
    // compute row and column of bin to which this observation belongs
    int32_t u_bin = min((int32_t)floor((float)u/(float)param.match_binsize),u_bin_num-1);
    int32_t v_bin = min((int32_t)floor((float)v/(float)param.match_binsize),v_bin_num-1);
    
    // save index
    // k[(c*v_bin_num+v_bin)*u_bin_num+u_bin].push_back(i);
    bin_idx = (c*v_bin_num+v_bin)*u_bin_num+u_bin;
    k[bin_idx][k_num[bin_idx]] = i;
    k_num[bin_idx] += 1;
  }
}

inline void Matcher::myFindMatch (int32_t* m1,const int32_t &i1,int32_t* m2,const int32_t &step_size,int32_t k2[BIN_NUM][MAX_FP_IN_BIN], int32_t k2_num[BIN_NUM], 
                                const int32_t &u_bin_num,const int32_t &v_bin_num,const int32_t &stat_bin,
                                int32_t& min_ind,int32_t stage) {
  
  // init and load image coordinates + feature
  min_ind          = 0;
  double  min_cost = 10000000;
  int32_t u1       = *(m1+step_size*i1+0);
  int32_t v1       = *(m1+step_size*i1+1);
  int32_t c        = *(m1+step_size*i1+3);
  uint8_t d1[32];
  uint8_t d2[32];
  int16_t  psum;
  memcpy(d1, m1+step_size*i1+4, 32);
  
//   __m128i xmm1     = _mm_load_si128((__m128i*)(m1+step_size*i1+4));
//   __m128i xmm2     = _mm_load_si128((__m128i*)(m1+step_size*i1+8));
  
  float u_min,u_max,v_min,v_max;
  
  u_min = u1-param.match_radius;
  u_max = u1+param.match_radius;
  v_min = v1-param.match_radius;
  v_max = v1+param.match_radius;

  // bins of interest
  int32_t u_bin_min = min(max((int32_t)floor(u_min/(float)param.match_binsize),0),u_bin_num-1);
  int32_t u_bin_max = min(max((int32_t)floor(u_max/(float)param.match_binsize),0),u_bin_num-1);
  int32_t v_bin_min = min(max((int32_t)floor(v_min/(float)param.match_binsize),0),v_bin_num-1);
  int32_t v_bin_max = min(max((int32_t)floor(v_max/(float)param.match_binsize),0),v_bin_num-1);
  
  // for all bins of interest do
  for (int32_t u_bin=u_bin_min; u_bin<=u_bin_max; u_bin++) {
    for (int32_t v_bin=v_bin_min; v_bin<=v_bin_max; v_bin++) {
      int32_t k2_ind = (c*v_bin_num+v_bin)*u_bin_num+u_bin;
    //   for (vector<int32_t>::const_iterator i2_it=k2[k2_ind].begin(); i2_it!=k2[k2_ind].end(); i2_it++) {
      for (int i2=0; i2 < k2_num[k2_ind]; i2++) {
        // int32_t u2   = *(m2+step_size*(*i2_it)+0);
        // int32_t v2   = *(m2+step_size*(*i2_it)+1);
        int32_t u2   = *(m2+step_size*(k2[k2_ind][i2])+0);
        int32_t v2   = *(m2+step_size*(k2[k2_ind][i2])+1);
        
        if (u2>=u_min && u2<=u_max && v2>=v_min && v2<=v_max) {
        //   __m128i xmm3 = _mm_load_si128((__m128i*)(m2+step_size*(*i2_it)+4));
        //   __m128i xmm4 = _mm_load_si128((__m128i*)(m2+step_size*(*i2_it)+8));              
            memcpy(d2, m2+step_size*(k2[k2_ind][i2])+4, 32);
        //   xmm3 = _mm_sad_epu8 (xmm1,xmm3);
        //   xmm4 = _mm_sad_epu8 (xmm2,xmm4);
        //   xmm4 = _mm_add_epi16(xmm3,xmm4);
            psum = 0;
            for (int i = 0; i < 32; i++) {
                psum += (int16_t) ABS(d1[i], d2[i]);
            }
            double cost = (double) psum;
        //   double cost = (double)(_mm_extract_epi16(xmm4,0)+_mm_extract_epi16(xmm4,4));
          
        //   if (u_>=0 && v_>=0) {
        //     double du = (double)u2-u_;
        //     double dv = (double)v2-v_;
        //     double dist = sqrt(du*du+dv*dv);
        //     cost += 4*dist;
        //   }
          
          if (cost<min_cost) {
            min_ind  = k2[k2_ind][i2];
            min_cost = cost;
          }
        }
      }
    }
  }
}

void Matcher::myMatching (int32_t *m1p,int32_t *m1c, int32_t n1p,int32_t n1c, vector<Matcher::p_match> &p_matched) {
// cout << "my" << endl;
  // descriptor step size (number of int32_t elements in struct)
  int32_t step_size = sizeof(Matcher::maximum)/sizeof(int32_t);
  
  // compute number of bins
  int32_t u_bin_num = (int32_t)ceil((float)dims_c[0]/(float)param.match_binsize);
  int32_t v_bin_num = (int32_t)ceil((float)dims_c[1]/(float)param.match_binsize);
//   cout << "u_bin_num " <<u_bin_num << " v_bin_num "  <<  v_bin_num << endl;
//   u_bin_num 21 v_bin_num 6


  int32_t bin_num   = 4*v_bin_num*u_bin_num; // 4 classes
  
  // allocate memory for index vectors (needed for efficient search)
//   vector<int32_t> *k1p = new vector<int32_t>[bin_num];
//   vector<int32_t> *k1c = new vector<int32_t>[bin_num];
  int32_t k1p[BIN_NUM][MAX_FP_IN_BIN];
  int32_t k1c[BIN_NUM][MAX_FP_IN_BIN];
  int32_t num1p[BIN_NUM] = {0};
  int32_t num1c[BIN_NUM] = {0};
//   cout << "init" << endl;
  // loop variables
//   int32_t* M = (int32_t*)calloc(dims_c[0]*dims_c[1],sizeof(int32_t));
  int32_t M[IMG_SIZE] = {0};
  int32_t i1p,i1c,i1c2;
  int32_t u1p,v1p,u1c,v1c;
  

  /////////////////////////////////////////////////////
  // method: flow

    
  // create position/class bin index vectors
  myCreateIndexVector(m1p,n1p,k1p, num1p,u_bin_num,v_bin_num);
  myCreateIndexVector(m1c,n1c,k1c, num1c,u_bin_num,v_bin_num);
//   cout << "create" << endl;
  // for all points do
  for (i1c=0; i1c<n1c; i1c++) {

    // coordinates in previous left image
    u1c = *(m1c+step_size*i1c+0);
    v1c = *(m1c+step_size*i1c+1);

    // compute row and column of statistics bin to which this observation belongs
    int32_t u_bin = min((int32_t)floor((float)u1c/(float)param.match_binsize),u_bin_num-1);
    int32_t v_bin = min((int32_t)floor((float)v1c/(float)param.match_binsize),v_bin_num-1);
    int32_t stat_bin = v_bin*u_bin_num+u_bin;

    // match forward/backward
    myFindMatch(m1c,i1c,m1p,step_size,k1p, num1p,u_bin_num,v_bin_num,stat_bin,i1p, 0);
    myFindMatch(m1p,i1p,m1c,step_size,k1c, num1c,u_bin_num,v_bin_num,stat_bin,i1c2,1);

    // circle closure success?
    if (i1c2==i1c) {

      // extract coordinates
      u1p = *(m1p+step_size*i1p+0);
      v1p = *(m1p+step_size*i1p+1);

      // add match if this pixel isn't matched yet
      if (M[getAddressOffsetImage(u1c,v1c,dims_c[0])]==0) {
        p_matched.push_back(Matcher::p_match(u1p,v1p,i1p,-1,-1,-1,u1c,v1c,i1c,-1,-1,-1));
        M[getAddressOffsetImage(u1c,v1c,dims_c[0])] = 1;
      }
    }
  }
  
}