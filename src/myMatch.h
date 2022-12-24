#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <emmintrin.h>
#include <algorithm>
#include <vector>
#include <cassert>

#include "matrix.h"
#include "size.h"

using namespace std;

#define ABS(a, b) ((a > b)? a-b:b-a);

void Matcher::myCreateIndexVector (int32_t m[MAX_FEATURE_ARRAY_SZIE],int32_t& n,int32_t k[BIN_NUM][MAX_FP_IN_BIN], int32_t k_num[BIN_NUM],const int32_t &u_bin_num,const int32_t &v_bin_num, int32_t m_num[BIN_NUM]) {
  // 
  int32_t fp_cnt = 0;
  int32_t max_tmp[MAX_FEATURE_ARRAY_SZIE];
  n = 0;
  // for each bin
  int32_t bin_idx;
  
  for (int _bin = 0; _bin < BIN_NUM; _bin++) {
    int c = _bin % 4;
    int spatial_idx = _bin / 4;

    // copy to tmp
    // construct k
    for (int e = 0; e < m_num[_bin]; e++) {
        for (int s=0; s < 12; s++) {
            max_tmp[12*fp_cnt+s] = m[MAX2_BIN_OFFSET*_bin+12*e+s];
        }
        // check
        // assert(m[MAX2_BIN_OFFSET*_bin+12*e+3]==c);
        // assert((m[MAX2_BIN_OFFSET*_bin+12*e+0]/BIN_W + U_BIN_NUM*m[MAX2_BIN_OFFSET*_bin+12*e+1]/BIN_H)==spatial_idx);
        // cout << "u " << m[MAX2_BIN_OFFSET*_bin+12*e+0] << " v " << m[MAX2_BIN_OFFSET*_bin+12*e+1] << " c " << m[MAX2_BIN_OFFSET*_bin+12*e+3] 
        //     << " ubin " << m[MAX2_BIN_OFFSET*_bin+12*e+0]/BIN_W << " vbin " << m[MAX2_BIN_OFFSET*_bin+12*e+1]/BIN_H << endl;
        bin_idx = c*U_BIN_NUM*V_BIN_NUM + (spatial_idx);
        k[bin_idx][k_num[bin_idx]] = fp_cnt;
        k_num[bin_idx] += 1;
        fp_cnt++;
    }
  }


  n = fp_cnt;

  // copy tmp back to m
  for (int i=0; i < n; i++) {
    for (int s=0; s <12; s++) {
        m[i*12+s] = max_tmp[i*12+s];
    }
  }




//   // descriptor step size
//   int32_t step_size = sizeof(Matcher::maximum)/sizeof(int32_t);
//   // for all points do
//   for (int32_t i=0; i<n; i++) {
//     // cout << i << endl;
//     // extract coordinates and class
//     int32_t u = *(m+step_size*i+0); // u-coordinate
//     int32_t v = *(m+step_size*i+1); // v-coordinate
//     int32_t c = *(m+step_size*i+3); // class
    
//     // compute row and column of bin to which this observation belongs
//     // int32_t u_bin = min((int32_t)floor((float)u/(float)BIN_W),u_bin_num-1);
//     // int32_t v_bin = min((int32_t)floor((float)v/(float)BIN_H),v_bin_num-1);
//     int32_t u_bin = u/BIN_W;
//     int32_t v_bin = v/BIN_H;
   
//     // save index
//     // k[(c*v_bin_num+v_bin)*u_bin_num+u_bin].push_back(i);
//     bin_idx = (c*v_bin_num+v_bin)*u_bin_num+u_binc;
    
//     k[bin_idx][k_num[bin_idx]] = i;
//     k_num[bin_idx] += 1;
//   }
}

inline void Matcher::myFindMatch (int32_t m1[MAX_FEATURE_ARRAY_SZIE],const int32_t &i1,int32_t m2[MAX_FEATURE_ARRAY_SZIE],const int32_t &step_size,int32_t k2[BIN_NUM][MAX_FP_IN_BIN], int32_t k2_num[BIN_NUM], 
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
  
  u_min = u1-SEARCH_RAD_U;
  u_max = u1+SEARCH_RAD_U;
  v_min = v1-SEARCH_RAD_V;
  v_max = v1+SEARCH_RAD_V;

  // bins of interest
//   int32_t u_bin_min = min(max((int32_t)floor(u_min/(float)BIN_W),0),u_bin_num-1);
//   int32_t u_bin_max = min(max((int32_t)floor(u_max/(float)BIN_W),0),u_bin_num-1);
//   int32_t v_bin_min = min(max((int32_t)floor(v_min/(float)BIN_H),0),v_bin_num-1);
//   int32_t v_bin_max = min(max((int32_t)floor(v_max/(float)BIN_H),0),v_bin_num-1);
  int32_t u_bin_min = min(max( (int32_t)u_min/BIN_W, 0), U_BIN_NUM-1);
  int32_t u_bin_max = min(max( (int32_t)u_max/BIN_W, 0), U_BIN_NUM-1);
  int32_t v_bin_min = min(max( (int32_t)v_min/BIN_H, 0), V_BIN_NUM-1);
  int32_t v_bin_max = min(max( (int32_t)v_max/BIN_H, 0), V_BIN_NUM-1);
  
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

void Matcher::myMatching (int32_t m1p[MAX_FEATURE_ARRAY_SZIE],int32_t m1c[MAX_FEATURE_ARRAY_SZIE], int32_t n1p[BIN_NUM],int32_t n1c[BIN_NUM],Matcher::p_match p_matched[POINT_L], int32_t& p_matched_num) {
// cout << "my" << endl;
// static int bin_max;
  // descriptor step size (number of int32_t elements in struct)
  int32_t step_size = sizeof(Matcher::maximum)/sizeof(int32_t);
  
  // compute number of bins
//   int32_t u_bin_num = (int32_t)ceil((float)dims_c[0]/(float)BIN_W);
//   int32_t v_bin_num = (int32_t)ceil((float)dims_c[1]/(float)BIN_H);
  int32_t u_bin_num = U_BIN_NUM;
  int32_t v_bin_num = V_BIN_NUM;
//   cout << "u_bin_num " <<u_bin_num << " v_bin_num "  <<  v_bin_num << endl;
//   u_bin_num 21 v_bin_num 6


//   int32_t bin_num   = 4*v_bin_num*u_bin_num; // 4 classes
  int32_t bin_num   = BIN_NUM; // 4 classes
  
  // allocate memory for index vectors (needed for efficient search)
//   vector<int32_t> *k1p = new vector<int32_t>[bin_num];
//   vector<int32_t> *k1c = new vector<int32_t>[bin_num];
  int32_t k1p[BIN_NUM][MAX_FP_IN_BIN];
  int32_t k1c[BIN_NUM][MAX_FP_IN_BIN];
  int32_t num1p[BIN_NUM] = {0};
  int32_t num1c[BIN_NUM] = {0};
  // loop variables
//   int32_t* M = (int32_t*)calloc(dims_c[0]*dims_c[1],sizeof(int32_t));
  int32_t M[IMG_SIZE] = {0};
  int32_t i1p,i1c,i1c2;
  int32_t u1p,v1p,u1c,v1c;
  

  /////////////////////////////////////////////////////
  // method: flow
    int nump = 0;
    int numc = 0;
    
  // create position/class bin index vectors
//   myCreateIndexVector(m1p,n1p,k1p, num1p,u_bin_num,v_bin_num);
//   myCreateIndexVector(m1c,n1c,k1c, num1c,u_bin_num,v_bin_num);  
  // TODO: matching change to bucket style

  myCreateIndexVector(m1p,nump, k1p, num1p,u_bin_num,v_bin_num,n1p);
  myCreateIndexVector(m1c,numc, k1c, num1c,u_bin_num,v_bin_num,n1c);
  // FIND BIN MAX
//   int max = 0;
//   for (int i=0; i < BIN_NUM; i++) {
//     if (num1c[i] > max) {
//         max = num1c[i];
//         // cout << max << endl;
//     }
//   }
//   if (max > bin_max) {
//     bin_max = max;
//     cout << "bin_max " << bin_max << endl;
//   }
//   cout << "create" << endl;
  // for all points do
//   for (i1c=0; i1c<n1c; i1c++) {

  for (i1c=0; i1c<numc; i1c++) {
    // cout << "i1c"<<i1c << endl;

    // coordinates in previous left image
    u1c = *(m1c+step_size*i1c+0);
    v1c = *(m1c+step_size*i1c+1);

    // compute row and column of statistics bin to which this observation belongs
    int32_t u_bin = u1c/BIN_W;
    int32_t v_bin = v1c/BIN_H;
    // cout << u_bin << " " << v_bin <<endl;
    int32_t stat_bin = v_bin*u_bin_num+u_bin;

    // match forward/backward
    myFindMatch(m1c,i1c,m1p,step_size,k1p, num1p,u_bin_num,v_bin_num,stat_bin,i1p, 0);
    myFindMatch(m1p,i1p,m1c,step_size,k1c, num1c,u_bin_num,v_bin_num,stat_bin,i1c2,1);
    // cout <<" find finish" << endl;
    // circle closure success?
    if (i1c2==i1c) {

      // extract coordinates
      u1p = *(m1p+step_size*i1p+0);
      v1p = *(m1p+step_size*i1p+1);

      // add match if this pixel isn't matched yet
      if (M[getAddressOffsetImage(u1c,v1c,dims_c[0])]==0) {
        // p_matched.push_back(Matcher::p_match(u1p,v1p,i1p,-1,-1,-1,u1c,v1c,i1c,-1,-1,-1));
        p_matched[p_matched_num] = Matcher::p_match(u1p,v1p,i1p,-1,-1,-1,u1c,v1c,i1c,-1,-1,-1);
        p_matched_num += 1;
        M[getAddressOffsetImage(u1c,v1c,dims_c[0])] = 1;
      }
    }
  }
  
}