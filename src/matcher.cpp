/*
Copyright 2012. All rights reserved.
Institute of Measurement and Control Systems
Karlsruhe Institute of Technology, Germany

This file is part of libviso2.
Authors: Andreas Geiger

libviso2 is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or any later version.

libviso2 is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
libviso2; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301, USA 
*/

#include "matcher.h"
#include "triangle.h"
#include "filter.h"

using namespace std;

#include "myComputeFeature.h"
//////////////////////
// PUBLIC FUNCTIONS //
//////////////////////

// constructor (with default parameters)
Matcher::Matcher(parameters param) : param(param) {

  // init match ring buffer to zero
  n1p2 = 0;
  m2p2 = 0; n2p2 = 0;
  n1c2 = 0;
  m2c2 = 0; n2c2 = 0;
  // margin needed to compute descriptor + sobel responses
  margin = 8+1;
  
}

// deconstructor
Matcher::~Matcher() {
  if (m2p2)         _mm_free(m2p2);
  if (m2c2)         _mm_free(m2c2);
}

void Matcher::pushBack (uint8_t *I1,uint8_t* I2,int32_t* dims,const bool replace) {
  cout << endl << "push back" << endl;
  // image dimensions
  int32_t width  = dims[0];
  int32_t height = dims[1];
  int32_t bpl    = dims[2];

  // sanity check
  if (width<=0 || height<=0 || bpl<width || I1==0) {
    cerr << "ERROR: Image dimension mismatch!" << endl;
    return;
  }

  if (replace) {
    if (m2c2)        _mm_free(m2c2);
  } else {
    if (m2p2)        _mm_free(m2p2);
    for (int i =0; i < 288000; i++) {
      m1p2[i] = m1c2[i];
    }
    n1p2 = n1c2;
    m2p2 = m2c2; n2p2 = n2c2;

    dims_p[0]   = dims_c[0];
    dims_p[1]   = dims_c[1];
    dims_p[2]   = dims_c[2];
  }

  // set new dims (bytes per line must be multiple of 16)
  dims_c[0] = width;
  dims_c[1] = height;
  dims_c[2] = width + 15-(width-1)%16;


  // compute new features for current frame
  // computeFeatures(I1c,dims_c,m1c1,n1c1,m1c2,n1c2,I1c_du,I1c_dv,I1c_du_full,I1c_dv_full);
  myComputeFeatures(I1,dims_c, m1c2,n1c2);
}

void Matcher::matchFeatures(int32_t method, Matrix *Tr_delta) {
  
  //////////////////
  // sanity check //
  //////////////////

  cout << "match" << endl;

  // clear old matches
  p_matched_1.clear();
  p_matched_2.clear();

  // matching(m1p2,m2p2,m1c2,m2c2,n1p2,n2p2,n1c2,n2c2,p_matched_2,method,false,Tr_delta);
  matching(m1p2,0,m1c2,0,n1p2,0,n1c2,0,p_matched_2,method,false,Tr_delta);
  removeOutliers(p_matched_2,method);
  cout << "match done" << endl;
}

void Matcher::bucketFeatures(int32_t max_features,float bucket_width,float bucket_height) {

  // find max values
  float u_max = 0;
  float v_max = 0;
  for (vector<p_match>::iterator it = p_matched_2.begin(); it!=p_matched_2.end(); it++) {
    if (it->u1c>u_max) u_max=it->u1c;
    if (it->v1c>v_max) v_max=it->v1c;
  }

  // allocate number of buckets needed
  int32_t bucket_cols = (int32_t)floor(u_max/bucket_width)+1;
  int32_t bucket_rows = (int32_t)floor(v_max/bucket_height)+1;
  vector<p_match> *buckets = new vector<p_match>[bucket_cols*bucket_rows];

  // assign matches to their buckets
  for (vector<p_match>::iterator it=p_matched_2.begin(); it!=p_matched_2.end(); it++) {
    int32_t u = (int32_t)floor(it->u1c/bucket_width);
    int32_t v = (int32_t)floor(it->v1c/bucket_height);
    buckets[v*bucket_cols+u].push_back(*it);
  }
  
  // refill p_matched from buckets
  p_matched_2.clear();
  for (int32_t i=0; i<bucket_cols*bucket_rows; i++) {
    
    // shuffle bucket indices randomly
    std::random_shuffle(buckets[i].begin(),buckets[i].end());
    
    // add up to max_features features from this bucket to p_matched
    int32_t k=0;
    for (vector<p_match>::iterator it=buckets[i].begin(); it!=buckets[i].end(); it++) {
      p_matched_2.push_back(*it);
      k++;
      if (k>=max_features)
        break;
    }
  }

  // free buckets
  delete []buckets;
}


///////////////////////
// PRIVATE FUNCTIONS //
///////////////////////

void Matcher::createIndexVector (int32_t* m,int32_t n,vector<int32_t> *k,const int32_t &u_bin_num,const int32_t &v_bin_num) {

  // descriptor step size
  int32_t step_size = sizeof(Matcher::maximum)/sizeof(int32_t);
  
  // for all points do
  for (int32_t i=0; i<n; i++) {
    
    // extract coordinates and class
    int32_t u = *(m+step_size*i+0); // u-coordinate
    int32_t v = *(m+step_size*i+1); // v-coordinate
    int32_t c = *(m+step_size*i+3); // class
    
    // compute row and column of bin to which this observation belongs
    int32_t u_bin = min((int32_t)floor((float)u/(float)param.match_binsize),u_bin_num-1);
    int32_t v_bin = min((int32_t)floor((float)v/(float)param.match_binsize),v_bin_num-1);
    
    // save index
    k[(c*v_bin_num+v_bin)*u_bin_num+u_bin].push_back(i);
  }
}

inline void Matcher::findMatch (int32_t* m1,const int32_t &i1,int32_t* m2,const int32_t &step_size,vector<int32_t> *k2,
                                const int32_t &u_bin_num,const int32_t &v_bin_num,const int32_t &stat_bin,
                                int32_t& min_ind,int32_t stage,bool flow,bool use_prior,double u_,double v_) {
  
  // init and load image coordinates + feature
  min_ind          = 0;
  double  min_cost = 10000000;
  int32_t u1       = *(m1+step_size*i1+0);
  int32_t v1       = *(m1+step_size*i1+1);
  int32_t c        = *(m1+step_size*i1+3);
  __m128i xmm1     = _mm_load_si128((__m128i*)(m1+step_size*i1+4));
  __m128i xmm2     = _mm_load_si128((__m128i*)(m1+step_size*i1+8));
  
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
      for (vector<int32_t>::const_iterator i2_it=k2[k2_ind].begin(); i2_it!=k2[k2_ind].end(); i2_it++) {
        int32_t u2   = *(m2+step_size*(*i2_it)+0);
        int32_t v2   = *(m2+step_size*(*i2_it)+1);
        if (u2>=u_min && u2<=u_max && v2>=v_min && v2<=v_max) {
          __m128i xmm3 = _mm_load_si128((__m128i*)(m2+step_size*(*i2_it)+4));
          __m128i xmm4 = _mm_load_si128((__m128i*)(m2+step_size*(*i2_it)+8));                    
          xmm3 = _mm_sad_epu8 (xmm1,xmm3);
          xmm4 = _mm_sad_epu8 (xmm2,xmm4);
          xmm4 = _mm_add_epi16(xmm3,xmm4);
          double cost = (double)(_mm_extract_epi16(xmm4,0)+_mm_extract_epi16(xmm4,4));
          
          if (u_>=0 && v_>=0) {
            double du = (double)u2-u_;
            double dv = (double)v2-v_;
            double dist = sqrt(du*du+dv*dv);
            cost += 4*dist;
          }
          
          if (cost<min_cost) {
            min_ind  = *i2_it;
            min_cost = cost;
          }
        }
      }
    }
  }
}

void Matcher::matching (int32_t *m1p,int32_t *m2p,int32_t *m1c,int32_t *m2c,
                        int32_t n1p,int32_t n2p,int32_t n1c,int32_t n2c,
                        vector<Matcher::p_match> &p_matched,int32_t method,bool use_prior,Matrix *Tr_delta) {

  // descriptor step size (number of int32_t elements in struct)
  int32_t step_size = sizeof(Matcher::maximum)/sizeof(int32_t);
  
  // compute number of bins
  int32_t u_bin_num = (int32_t)ceil((float)dims_c[0]/(float)param.match_binsize);
  int32_t v_bin_num = (int32_t)ceil((float)dims_c[1]/(float)param.match_binsize);
  int32_t bin_num   = 4*v_bin_num*u_bin_num; // 4 classes
  
  // allocate memory for index vectors (needed for efficient search)
  vector<int32_t> *k1p = new vector<int32_t>[bin_num];
  vector<int32_t> *k2p = new vector<int32_t>[bin_num];
  vector<int32_t> *k1c = new vector<int32_t>[bin_num];
  vector<int32_t> *k2c = new vector<int32_t>[bin_num];
  
  // loop variables
  int32_t* M = (int32_t*)calloc(dims_c[0]*dims_c[1],sizeof(int32_t));
  int32_t i1p,i2p,i1c,i2c,i1c2,i1p2;
  int32_t u1p,v1p,u2p,v2p,u1c,v1c,u2c,v2c;
  
  double t00,t01,t02,t03,t10,t11,t12,t13,t20,t21,t22,t23;

  /////////////////////////////////////////////////////
  // method: flow

    
  // create position/class bin index vectors
  createIndexVector(m1p,n1p,k1p,u_bin_num,v_bin_num);
  createIndexVector(m1c,n1c,k1c,u_bin_num,v_bin_num);
  
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
    findMatch(m1c,i1c,m1p,step_size,k1p,u_bin_num,v_bin_num,stat_bin,i1p, 0,true,use_prior);
    findMatch(m1p,i1p,m1c,step_size,k1c,u_bin_num,v_bin_num,stat_bin,i1c2,1,true,use_prior);

    // circle closure success?
    if (i1c2==i1c) {

      // extract coordinates
      u1p = *(m1p+step_size*i1p+0);
      v1p = *(m1p+step_size*i1p+1);

      // add match if this pixel isn't matched yet
      if (*(M+getAddressOffsetImage(u1c,v1c,dims_c[0]))==0) {
        p_matched.push_back(Matcher::p_match(u1p,v1p,i1p,-1,-1,-1,u1c,v1c,i1c,-1,-1,-1));
        *(M+getAddressOffsetImage(u1c,v1c,dims_c[0])) = 1;
      }
    }
  }
  
  // free memory
  free(M);
  delete []k1p;
  delete []k2p;
  delete []k1c;
  delete []k2c;
}

void Matcher::removeOutliers (vector<Matcher::p_match> &p_matched,int32_t method) {

  // do we have enough points for outlier removal?
  if (p_matched.size()<=3)
    return;

  // input/output structure for triangulation
  struct triangulateio in, out;

  // inputs
  in.numberofpoints = p_matched.size();
  in.pointlist = (float*)malloc(in.numberofpoints*2*sizeof(float));
  int32_t k=0;
  
  // create copy of p_matched, init vector with number of support points
  // and fill triangle point vector for delaunay triangulation
  vector<Matcher::p_match> p_matched_copy;  
  vector<int32_t> num_support;
  for (vector<Matcher::p_match>::iterator it=p_matched.begin(); it!=p_matched.end(); it++) {
    p_matched_copy.push_back(*it);
    num_support.push_back(0);
    in.pointlist[k++] = it->u1c;
    in.pointlist[k++] = it->v1c;
  }
  
  // input parameters
  in.numberofpointattributes = 0;
  in.pointattributelist      = NULL;
  in.pointmarkerlist         = NULL;
  in.numberofsegments        = 0;
  in.numberofholes           = 0;
  in.numberofregions         = 0;
  in.regionlist              = NULL;
  
  // outputs
  out.pointlist              = NULL;
  out.pointattributelist     = NULL;
  out.pointmarkerlist        = NULL;
  out.trianglelist           = NULL;
  out.triangleattributelist  = NULL;
  out.neighborlist           = NULL;
  out.segmentlist            = NULL;
  out.segmentmarkerlist      = NULL;
  out.edgelist               = NULL;
  out.edgemarkerlist         = NULL;

  // do triangulation (z=zero-based, n=neighbors, Q=quiet, B=no boundary markers)
  // attention: not using the B switch or using the n switch creates a memory leak (=> use valgrind!)
  char parameters[] = "zQB";
  triangulate(parameters, &in, &out, NULL);
  
  // for all triangles do
  for (int32_t i=0; i<out.numberoftriangles; i++) {
    
    // extract triangle corner points
    int32_t p1 = out.trianglelist[i*3+0];
    int32_t p2 = out.trianglelist[i*3+1];
    int32_t p3 = out.trianglelist[i*3+2];
    
    // method: flow
    if (method==0) {
      
      // 1. corner disparity and flow
      float p1_flow_u = p_matched_copy[p1].u1c-p_matched_copy[p1].u1p;
      float p1_flow_v = p_matched_copy[p1].v1c-p_matched_copy[p1].v1p;

      // 2. corner disparity and flow
      float p2_flow_u = p_matched_copy[p2].u1c-p_matched_copy[p2].u1p;
      float p2_flow_v = p_matched_copy[p2].v1c-p_matched_copy[p2].v1p;

      // 3. corner disparity and flow
      float p3_flow_u = p_matched_copy[p3].u1c-p_matched_copy[p3].u1p;
      float p3_flow_v = p_matched_copy[p3].v1c-p_matched_copy[p3].v1p;

      // consistency of 1. edge
      if (fabs(p1_flow_u-p2_flow_u)+fabs(p1_flow_v-p2_flow_v)<param.outlier_flow_tolerance) {
        num_support[p1]++;
        num_support[p2]++;
      }

      // consistency of 2. edge
      if (fabs(p2_flow_u-p3_flow_u)+fabs(p2_flow_v-p3_flow_v)<param.outlier_flow_tolerance) {
        num_support[p2]++;
        num_support[p3]++;
      }

      // consistency of 3. edge
      if (fabs(p1_flow_u-p3_flow_u)+fabs(p1_flow_v-p3_flow_v)<param.outlier_flow_tolerance) {
        num_support[p1]++;
        num_support[p3]++;
      }
      
    // method: stereo
    } else if (method==1) {
      
      // 1. corner disparity and flow
      float p1_disp   = p_matched_copy[p1].u1c-p_matched_copy[p1].u2c;

      // 2. corner disparity and flow
      float p2_disp   = p_matched_copy[p2].u1c-p_matched_copy[p2].u2c;

      // 3. corner disparity and flow
      float p3_disp   = p_matched_copy[p3].u1c-p_matched_copy[p3].u2c;

      // consistency of 1. edge
      if (fabs(p1_disp-p2_disp)<param.outlier_disp_tolerance) {
        num_support[p1]++;
        num_support[p2]++;
      }

      // consistency of 2. edge
      if (fabs(p2_disp-p3_disp)<param.outlier_disp_tolerance) {
        num_support[p2]++;
        num_support[p3]++;
      }

      // consistency of 3. edge
      if (fabs(p1_disp-p3_disp)<param.outlier_disp_tolerance) {
        num_support[p1]++;
        num_support[p3]++;
      }
      
    // method: quad matching
    } else {
      
      // 1. corner disparity and flow
      float p1_flow_u = p_matched_copy[p1].u1c-p_matched_copy[p1].u1p;
      float p1_flow_v = p_matched_copy[p1].v1c-p_matched_copy[p1].v1p;
      float p1_disp   = p_matched_copy[p1].u1p-p_matched_copy[p1].u2p;

      // 2. corner disparity and flow
      float p2_flow_u = p_matched_copy[p2].u1c-p_matched_copy[p2].u1p;
      float p2_flow_v = p_matched_copy[p2].v1c-p_matched_copy[p2].v1p;
      float p2_disp   = p_matched_copy[p2].u1p-p_matched_copy[p2].u2p;

      // 3. corner disparity and flow
      float p3_flow_u = p_matched_copy[p3].u1c-p_matched_copy[p3].u1p;
      float p3_flow_v = p_matched_copy[p3].v1c-p_matched_copy[p3].v1p;
      float p3_disp   = p_matched_copy[p3].u1p-p_matched_copy[p3].u2p;

      // consistency of 1. edge
      if (fabs(p1_disp-p2_disp)<param.outlier_disp_tolerance && fabs(p1_flow_u-p2_flow_u)+fabs(p1_flow_v-p2_flow_v)<param.outlier_flow_tolerance) {
        num_support[p1]++;
        num_support[p2]++;
      }

      // consistency of 2. edge
      if (fabs(p2_disp-p3_disp)<param.outlier_disp_tolerance && fabs(p2_flow_u-p3_flow_u)+fabs(p2_flow_v-p3_flow_v)<param.outlier_flow_tolerance) {
        num_support[p2]++;
        num_support[p3]++;
      }

      // consistency of 3. edge
      if (fabs(p1_disp-p3_disp)<param.outlier_disp_tolerance && fabs(p1_flow_u-p3_flow_u)+fabs(p1_flow_v-p3_flow_v)<param.outlier_flow_tolerance) {
        num_support[p1]++;
        num_support[p3]++;
      }
    }
  }
  
  // refill p_matched
  p_matched.clear();
  for (int i=0; i<in.numberofpoints; i++)
    if (num_support[i]>=4)
      p_matched.push_back(p_matched_copy[i]);
  
  // free memory used for triangulation
  free(in.pointlist);
  free(out.pointlist);
  free(out.trianglelist);
}

float Matcher::mean(const uint8_t* I,const int32_t &bpl,const int32_t &u_min,const int32_t &u_max,const int32_t &v_min,const int32_t &v_max) {
  float mean = 0;
  for (int32_t v=v_min; v<=v_max; v++)
    for (int32_t u=u_min; u<=u_max; u++)
      mean += (float)*(I+getAddressOffsetImage(u,v,bpl));
  return
    mean /= (float)((u_max-u_min+1)*(v_max-v_min+1));
}


// ============================================================================================ // 
//
//
// Compute Feature Function
//
//
// ============================================================================================ // 


void Matcher::nonMaximumSuppression (int16_t* I_f1,int16_t* I_f2,const int32_t* dims,vector<Matcher::maximum> &maxima,int32_t nms_n) {
  
  // extract parameters
  int32_t width  = dims[0];
  int32_t height = dims[1];
  int32_t bpl    = dims[2];
  int32_t n      = nms_n;
  int32_t tau    = param.nms_tau;
  
  // loop variables
  int32_t f1mini,f1minj,f1maxi,f1maxj,f2mini,f2minj,f2maxi,f2maxj;
  int32_t f1minval,f1maxval,f2minval,f2maxval,currval;
  int32_t addr;
  // cout << "margin: " << margin << endl;
  // cout << "width: " << width << "height: " <<height<<  "bpl: " << bpl << endl;
  for (int32_t j=n+margin; j<height-n-margin;j+=n+1) {
    for (int32_t i=n+margin; i<width-n-margin;i+=n+1) { 
      // cout << i<< " " << j << endl;
      f1mini = i; f1minj = j; f1maxi = i; f1maxj = j;
      f2mini = i; f2minj = j; f2maxi = i; f2maxj = j;
      
      addr     = getAddressOffsetImage(i,j,bpl);
      f1minval = *(I_f1+addr);
      f1maxval = f1minval;
      f2minval = *(I_f2+addr);
      f2maxval = f2minval;

      for (int32_t j2=j; j2<=j+n; j2++) {
        for (int32_t i2=i; i2<=i+n; i2++) {
          addr    = getAddressOffsetImage(i2,j2,bpl);
          currval = *(I_f1+addr);
          if (currval<f1minval) {
            f1mini   = i2;
            f1minj   = j2;
            f1minval = currval;
          } else if (currval>f1maxval) {
            f1maxi   = i2;
            f1maxj   = j2;
            f1maxval = currval;
          }
          currval = *(I_f2+addr);
          if (currval<f2minval) {
            f2mini   = i2;
            f2minj   = j2;
            f2minval = currval;
          } else if (currval>f2maxval) {
            f2maxi   = i2;
            f2maxj   = j2;
            f2maxval = currval;
          }
        }
      }
      
      // f1 minimum
      for (int32_t j2=f1minj-n; j2<=min(f1minj+n,height-1-margin); j2++) {
        for (int32_t i2=f1mini-n; i2<=min(f1mini+n,width-1-margin); i2++) {
          currval = *(I_f1+getAddressOffsetImage(i2,j2,bpl));
          if (currval<f1minval && (i2<i || i2>i+n || j2<j || j2>j+n))
            goto failed_f1min;            
        }
      }
      if (f1minval<=-tau)
        maxima.push_back(Matcher::maximum(f1mini,f1minj,f1minval,0));
      failed_f1min:;

      // f1 maximum
      for (int32_t j2=f1maxj-n; j2<=min(f1maxj+n,height-1-margin); j2++) {
        for (int32_t i2=f1maxi-n; i2<=min(f1maxi+n,width-1-margin); i2++) {
          currval = *(I_f1+getAddressOffsetImage(i2,j2,bpl));
          if (currval>f1maxval && (i2<i || i2>i+n || j2<j || j2>j+n))
            goto failed_f1max;
        }
      }
      if (f1maxval>=tau)
        maxima.push_back(Matcher::maximum(f1maxi,f1maxj,f1maxval,1));
      failed_f1max:;
      
      // f2 minimum
      for (int32_t j2=f2minj-n; j2<=min(f2minj+n,height-1-margin); j2++) {
        for (int32_t i2=f2mini-n; i2<=min(f2mini+n,width-1-margin); i2++) {
          currval = *(I_f2+getAddressOffsetImage(i2,j2,bpl));
          if (currval<f2minval && (i2<i || i2>i+n || j2<j || j2>j+n))
            goto failed_f2min;
        }
      }
      if (f2minval<=-tau)
        maxima.push_back(Matcher::maximum(f2mini,f2minj,f2minval,2));
      failed_f2min:;

      // f2 maximum
      for (int32_t j2=f2maxj-n; j2<=min(f2maxj+n,height-1-margin); j2++) {
        for (int32_t i2=f2maxi-n; i2<=min(f2maxi+n,width-1-margin); i2++) {
          currval = *(I_f2+getAddressOffsetImage(i2,j2,bpl));
          if (currval>f2maxval && (i2<i || i2>i+n || j2<j || j2>j+n))
            goto failed_f2max;
        }
      }
      if (f2maxval>=tau)
        maxima.push_back(Matcher::maximum(f2maxi,f2maxj,f2maxval,3));
      failed_f2max:;
    }
  }
}

inline void Matcher::computeDescriptor (const uint8_t* I_du,const uint8_t* I_dv,const int32_t &bpl,const int32_t &u,const int32_t &v,uint8_t *desc_addr) {
  
    // get address indices
  int32_t addr_m1 = getAddressOffsetImage(u,v-1,bpl);
  int32_t addr_m3 = addr_m1-2*bpl;
  int32_t addr_m5 = addr_m3-2*bpl;
  int32_t addr_p1 = addr_m1+2*bpl;
  int32_t addr_p3 = addr_p1+2*bpl;
  int32_t addr_p5 = addr_p3+2*bpl;
  
  // compute descriptor
  uint32_t k = 0;
  desc_addr[k++] = I_du[addr_m1-3];
  desc_addr[k++] = I_dv[addr_m1-3];
  desc_addr[k++] = I_du[addr_p1-3];
  desc_addr[k++] = I_dv[addr_p1-3];
  desc_addr[k++] = I_du[addr_m1-1];
  desc_addr[k++] = I_dv[addr_m1-1];
  desc_addr[k++] = I_du[addr_p1-1];
  desc_addr[k++] = I_dv[addr_p1-1];
  desc_addr[k++] = I_du[addr_m1+3];
  desc_addr[k++] = I_dv[addr_m1+3];
  desc_addr[k++] = I_du[addr_p1+3];
  desc_addr[k++] = I_dv[addr_p1+3];
  desc_addr[k++] = I_du[addr_m1+1];
  desc_addr[k++] = I_dv[addr_m1+1];
  desc_addr[k++] = I_du[addr_p1+1];
  desc_addr[k++] = I_dv[addr_p1+1];
  desc_addr[k++] = I_du[addr_m5-1];
  desc_addr[k++] = I_dv[addr_m5-1];
  desc_addr[k++] = I_du[addr_p5-1];
  desc_addr[k++] = I_dv[addr_p5-1];
  desc_addr[k++] = I_du[addr_m5+1];
  desc_addr[k++] = I_dv[addr_m5+1];
  desc_addr[k++] = I_du[addr_p5+1];
  desc_addr[k++] = I_dv[addr_p5+1];
  desc_addr[k++] = I_du[addr_m3-5];
  desc_addr[k++] = I_dv[addr_m3-5];
  desc_addr[k++] = I_du[addr_p3-5];
  desc_addr[k++] = I_dv[addr_p3-5];
  desc_addr[k++] = I_du[addr_m3+5];
  desc_addr[k++] = I_dv[addr_m3+5];
  desc_addr[k++] = I_du[addr_p3+5];
  desc_addr[k++] = I_dv[addr_p3+5];
}

inline void Matcher::computeSmallDescriptor (const uint8_t* I_du,const uint8_t* I_dv,const int32_t &bpl,const int32_t &u,const int32_t &v,uint8_t *desc_addr) {
  
  // get address indices
  int32_t addr2 = getAddressOffsetImage(u,v,bpl);
  int32_t addr1 = addr2-bpl;
  int32_t addr0 = addr1-bpl;
  int32_t addr3 = addr2+bpl;
  int32_t addr4 = addr3+bpl;
  
  // compute ELAS-descriptor
  uint32_t k = 0;
  desc_addr[k++] = I_du[addr0];
  desc_addr[k++] = I_du[addr1-2];
  desc_addr[k++] = I_du[addr1];
  desc_addr[k++] = I_du[addr1+2];
  desc_addr[k++] = I_du[addr2-1];
  desc_addr[k++] = I_du[addr2];
  desc_addr[k++] = I_du[addr2];
  desc_addr[k++] = I_du[addr2+1];
  desc_addr[k++] = I_du[addr3-2];
  desc_addr[k++] = I_du[addr3];
  desc_addr[k++] = I_du[addr3+2];
  desc_addr[k++] = I_du[addr4];
  desc_addr[k++] = I_dv[addr1];
  desc_addr[k++] = I_dv[addr2-1];
  desc_addr[k++] = I_dv[addr2+1];
  desc_addr[k++] = I_dv[addr3];
}

void Matcher::computeDescriptors (uint8_t* I_du,uint8_t* I_dv,const int32_t bpl,std::vector<Matcher::maximum> &maxima) {
  
  // loop variables
  int32_t u,v;
  uint8_t *desc_addr;
  
  // for all maxima do
  for (vector<Matcher::maximum>::iterator it=maxima.begin(); it!=maxima.end(); it++) {
    u = (*it).u;
    v = (*it).v;
    desc_addr = (uint8_t*)(&((*it).d1));
    computeDescriptor(I_du,I_dv,bpl,u,v,desc_addr);    
  }
}

inline uint8_t Matcher::saturate (int16_t in) {
  if (in<0)   return 0;
  if (in>255) return 255;
  return in;
}

void Matcher::getHalfResolutionDimensions(const int32_t *dims,int32_t *dims_half) {
  dims_half[0] = dims[0]/2;
  dims_half[1] = dims[1]/2;
  dims_half[2] = dims_half[0]+15-(dims_half[0]-1)%16;
}

uint8_t* Matcher::createHalfResolutionImage(uint8_t *I,const int32_t* dims) {
  int32_t dims_half[3];
  getHalfResolutionDimensions(dims,dims_half);
  uint8_t* I_half = (uint8_t*)_mm_malloc(dims_half[2]*dims_half[1]*sizeof(uint8_t),16);
  for (int32_t v=0; v<dims_half[1]; v++)
    for (int32_t u=0; u<dims_half[0]; u++)
      I_half[v*dims_half[2]+u] =  (uint8_t)(((int32_t)I[(v*2+0)*dims[2]+u*2+0]+
                                             (int32_t)I[(v*2+0)*dims[2]+u*2+1]+
                                             (int32_t)I[(v*2+1)*dims[2]+u*2+0]+
                                             (int32_t)I[(v*2+1)*dims[2]+u*2+1])/4);
  return I_half;
}

void Matcher::computeFeatures (uint8_t *I,const int32_t* dims,int32_t* &max1,int32_t &num1,int32_t* &max2,int32_t &num2,uint8_t* &I_du,uint8_t* &I_dv,uint8_t* &I_du_full,uint8_t* &I_dv_full) {
  
  int16_t *I_f1;
  int16_t *I_f2;
  
  int32_t dims_matching[3];
  memcpy(dims_matching,dims,3*sizeof(int32_t));
  
  // allocate memory for sobel images and filter images
  if (!param.half_resolution) {
    // cout << "check0" << endl;
    I_du = (uint8_t*)_mm_malloc(dims[2]*dims[1]*sizeof(uint8_t*),16);
    I_dv = (uint8_t*)_mm_malloc(dims[2]*dims[1]*sizeof(uint8_t*),16);
    I_f1 = (int16_t*)_mm_malloc(dims[2]*dims[1]*sizeof(int16_t),16);
    I_f2 = (int16_t*)_mm_malloc(dims[2]*dims[1]*sizeof(int16_t),16);
    filter::sobel5x5(I,I_du,I_dv,dims[2],dims[1]);
    filter::blob5x5(I,I_f1,dims[2],dims[1]);
    filter::checkerboard5x5(I,I_f2,dims[2],dims[1]);
  } else {
    uint8_t* I_matching = createHalfResolutionImage(I,dims);
    getHalfResolutionDimensions(dims,dims_matching);
    I_du      = (uint8_t*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(uint8_t*),16);
    I_dv      = (uint8_t*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(uint8_t*),16);
    I_f1      = (int16_t*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(int16_t),16);
    I_f2      = (int16_t*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(int16_t),16);
    I_du_full = (uint8_t*)_mm_malloc(dims[2]*dims[1]*sizeof(uint8_t*),16);
    I_dv_full = (uint8_t*)_mm_malloc(dims[2]*dims[1]*sizeof(uint8_t*),16);
    filter::sobel5x5(I_matching,I_du,I_dv,dims_matching[2],dims_matching[1]);
    filter::sobel5x5(I,I_du_full,I_dv_full,dims[2],dims[1]);
    filter::blob5x5(I_matching,I_f1,dims_matching[2],dims_matching[1]);
    filter::checkerboard5x5(I_matching,I_f2,dims_matching[2],dims_matching[1]);
    _mm_free(I_matching);
  }
  
  // extract sparse maxima (1st pass) via non-maximum suppression
  vector<Matcher::maximum> maxima1;
  if (param.multi_stage) {
    int32_t nms_n_sparse = param.nms_n*4;  // param.nms_n*3
    if (nms_n_sparse>10)
      nms_n_sparse = max(param.nms_n,10);
      // cout << "1st pass n: " << nms_n_sparse << endl;
    nonMaximumSuppression(I_f1,I_f2,dims_matching,maxima1,nms_n_sparse);
    computeDescriptors(I_du,I_dv,dims_matching[2],maxima1);
  }
  
  // extract dense maxima (2nd pass) via non-maximum suppression
  vector<Matcher::maximum> maxima2;
  // cout << "2nd pass n: " << param.nms_n << endl;

  nonMaximumSuppression(I_f1,I_f2,dims_matching,maxima2,param.nms_n);
  computeDescriptors(I_du,I_dv,dims_matching[2],maxima2);

  // release filter images
  _mm_free(I_f1);
  _mm_free(I_f2);  
  
  // get number of interest points and init maxima pointer to NULL
  num1 = maxima1.size();
  num2 = maxima2.size();
  max1 = 0;
  max2 = 0;
  
  int32_t s = 1;
  if (param.half_resolution)
    s = 2;

  // return sparse maxima as 16-bytes aligned memory
  if (num1!=0) {
    max1 = (int32_t*)_mm_malloc(sizeof(Matcher::maximum)*num1,16);
    int32_t k=0;
    for (vector<Matcher::maximum>::iterator it=maxima1.begin(); it!=maxima1.end(); it++) {
      *(max1+k++) = it->u*s;  *(max1+k++) = it->v*s;  *(max1+k++) = 0;        *(max1+k++) = it->c;
      *(max1+k++) = it->d1;   *(max1+k++) = it->d2;   *(max1+k++) = it->d3;   *(max1+k++) = it->d4;
      *(max1+k++) = it->d5;   *(max1+k++) = it->d6;   *(max1+k++) = it->d7;   *(max1+k++) = it->d8;
    }
  }
  
  // return dense maxima as 16-bytes aligned memory
  if (num2!=0) {
    max2 = (int32_t*)_mm_malloc(sizeof(Matcher::maximum)*num2,16);
    int32_t k=0;
    for (vector<Matcher::maximum>::iterator it=maxima2.begin(); it!=maxima2.end(); it++) {
      *(max2+k++) = it->u*s;  *(max2+k++) = it->v*s;  *(max2+k++) = 0;        *(max2+k++) = it->c;
      *(max2+k++) = it->d1;   *(max2+k++) = it->d2;   *(max2+k++) = it->d3;   *(max2+k++) = it->d4;
      *(max2+k++) = it->d5;   *(max2+k++) = it->d6;   *(max2+k++) = it->d7;   *(max2+k++) = it->d8;
    }
  }
}

