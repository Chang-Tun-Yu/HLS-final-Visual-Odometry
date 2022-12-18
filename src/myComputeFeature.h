
void Matcher::myGetHalfResolutionDimensions(const int32_t *dims,int32_t *dims_half) {
  dims_half[0] = dims[0]/2;
  dims_half[1] = dims[1]/2;
  dims_half[2] = dims_half[0]+15-(dims_half[0]-1)%16;
}

uint8_t* Matcher::myCreateHalfResolutionImage(uint8_t *I,const int32_t* dims) {
  int32_t dims_half[3];
  myGetHalfResolutionDimensions(dims,dims_half);
  uint8_t* I_half = (uint8_t*)_mm_malloc(dims_half[2]*dims_half[1]*sizeof(uint8_t),16);
  for (int32_t v=0; v<dims_half[1]; v++)
    for (int32_t u=0; u<dims_half[0]; u++)
      I_half[v*dims_half[2]+u] =  (uint8_t)(((int32_t)I[(v*2+0)*dims[2]+u*2+0]+
                                             (int32_t)I[(v*2+0)*dims[2]+u*2+1]+
                                             (int32_t)I[(v*2+1)*dims[2]+u*2+0]+
                                             (int32_t)I[(v*2+1)*dims[2]+u*2+1])/4);
  return I_half;
}

void Matcher::myComputeDescriptor (const uint8_t* I_du,const uint8_t* I_dv,const int32_t &bpl,const int32_t &u,const int32_t &v,uint8_t *desc_addr) {
  
    // get address indices
  int32_t addr_m1 = myGetAddressOffsetImage(u,v-1,bpl);
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

void Matcher::myComputeDescriptors (uint8_t* I_du,uint8_t* I_dv,const int32_t bpl,std::vector<Matcher::maximum> &maxima) {
      
  // loop variables
  int32_t u,v;
  uint8_t *desc_addr;
  
  // for all maxima do
  for (vector<Matcher::maximum>::iterator it=maxima.begin(); it!=maxima.end(); it++) {
    u = (*it).u;
    v = (*it).v;
    desc_addr = (uint8_t*)(&((*it).d1));
    myComputeDescriptor(I_du,I_dv,bpl,u,v,desc_addr);    
  }
}

void  Matcher::mySobel5x5 ( const uint8_t* in, uint8_t* out_v, uint8_t* out_h, int w, int h ) {
    const int16_t v_weight[5][5] = {
        {1, 2, 0, -2, -1},
        {4, 8, 0, -8, -4},
        {6, 12, 0, -12, -6},
        {4, 8, 0, -8, -4},
        {1, 2, 0, -2, -1}
    };
    const int16_t h_weight[5][5] = {
        {1, 4, 6, 4, 1},
        {2, 8, 12, 8, 2}, 
        {0, 0, 0, 0, 0},
        {-2, -8, -12, -8, -2}, 
        {-1, -4, -6, -4, -1}
    };
    for (int y = 2; y < h-2; y++) {
        for (int x = 2; x < w-2; x++) {
          // cout << x << "    " << y << endl;
            int16_t v_tmp = 0;
            int16_t h_tmp = 0;
            for (int ky = -2; ky <= 2; ky++) {
                for (int kx = -2; kx <=2 ; kx++) {
                    v_tmp += in[myGetAddressOffsetImage(x+kx, y+ky, w)] * v_weight[ky+2][kx+2];
                    h_tmp += in[myGetAddressOffsetImage(x+kx, y+ky, w)] * h_weight[ky+2][kx+2];
                }
            }
            // todo: saturate to uint8
            v_tmp /= 128;
            h_tmp /= 128;
            v_tmp += 128;
            h_tmp += 128;
            if (v_tmp < 0) {
                out_v[myGetAddressOffsetImage(x, y, w)] = 0;
            }
            else if (v_tmp > 255) {
                out_v[myGetAddressOffsetImage(x, y, w)] = 255;
            }
            else {
                out_v[myGetAddressOffsetImage(x, y, w)] = (uint8_t) v_tmp;
            }

            if (h_tmp < 0) {
                out_h[myGetAddressOffsetImage(x, y, w)] = 0;
            }
            else if (h_tmp > 255) {
                out_h[myGetAddressOffsetImage(x, y, w)] = 255;
            }
            else {
                out_h[myGetAddressOffsetImage(x, y, w)] = (uint8_t) h_tmp;
            }
        }
    }
}

void  Matcher::myCheckerboard5x5 ( const uint8_t* in, int16_t* out, int w, int h ) {
    const int16_t weight[5][5] = {
        {-1, -1, 0, 1, 1},
        {-1, -1, 0, 1, 1},
        {0, 0, 0, 0, 0},
        {1, 1, 0, -1, -1},
        {1, 1, 0, -1, -1}
    };
    for (int y = 2; y < h-2; y++) {
        for (int x = 2; x < w-2; x++) {
            int16_t tmp = 0;
            for (int ky = -2; ky <= 2; ky++) {
                for (int kx = -2; kx <=2 ; kx++) {
                    tmp += in[myGetAddressOffsetImage(x+kx, y+ky, w)] * weight[ky+2][kx+2];
                }
            }
            out[myGetAddressOffsetImage(x, y, w)] = tmp;
        }
    }
}

void Matcher::myBlob5x5( const uint8_t* in, int16_t* out, int w, int h ) {
    const int16_t weight[5][5] = {
        {-1, -1, -1, -1, -1},
        {-1, 1, 1, 1, -1},
        {-1, 1, 8, 1, -1},
        {-1, 1, 1, 1, -1},
        {-1, -1, -1, -1, -1}
    };
    for (int y = 2; y < h-2; y++) {
        for (int x = 2; x < w-2; x++) {
            int16_t tmp = 0;
            for (int ky = -2; ky <= 2; ky++) {
                for (int kx = -2; kx <=2 ; kx++) {
                    tmp += in[myGetAddressOffsetImage(x+kx, y+ky, w)] * weight[ky+2][kx+2];
                }
            }
            out[myGetAddressOffsetImage(x, y, w)] = tmp;
        }
    }
}

void Matcher::myNonMaximumSuppression (int16_t* I_f1,int16_t* I_f2,const int32_t* dims,vector<Matcher::maximum> &maxima1,vector<Matcher::maximum> &maxima2) {
    int32_t width  = dims[0];
    int32_t height = dims[1];
    int32_t bpl    = dims[2];
    int32_t tau    = 50;
    
    int16_t f1map[27][27];
    int16_t f2map[27][27];
    int32_t f1mini[10],f1minj[10],f1maxi[10],f1maxj[10],f2mini[10],f2minj[10],f2maxi[10],f2maxj[10];
    int32_t f1minval[10],f1maxval[10],f2minval[10],f2maxval[10];
    // cout << "[Non] INIT" << endl;
    for (int j=8+9; j < height-8-9; j+=8+1) {
        for (int i=8+9; i < width-8-9; i+=8+1) {
            // get map
            for (int jj = 0; jj < 27; jj++) {
                for (int ii = 0; ii < 27; ii++) {
                    f1map[jj][ii] = I_f1[myGetAddressOffsetImage(i + ii - 9, j + jj - 9, bpl)];
                    f2map[jj][ii] = I_f2[myGetAddressOffsetImage(i + ii - 9, j + jj - 9, bpl)];
                }
            }
    // cout << "[Non] COPY" << endl;
            
            // 3*3 filter
            int cnt = 0;
            for (int jj = 9; jj < 18; jj += 3) {
                for (int ii = 9; ii < 18; ii += 3) {
                    // init
                    f1minval[cnt] = 32767;
                    f2minval[cnt] = 32767;
                    f1maxval[cnt] = -32767;
                    f2maxval[cnt] = -32767;
                    for (int ky=0; ky < 3; ky++) {
                        for (int kx=0; kx < 3; kx++) {
                            // f1 min
                            if (f1map[jj+ky][ii+kx] < f1minval[cnt]) {
                                f1mini[cnt] = ii+kx;
                                f1minj[cnt] = jj+ky;
                                f1minval[cnt] = f1map[jj+ky][ii+kx];
                            }

                            // f1 max
                            if (f1map[jj+ky][ii+kx] > f1maxval[cnt]) {
                                f1maxi[cnt] = ii+kx;
                                f1maxj[cnt] = jj+ky;
                                f1maxval[cnt] = f1map[jj+ky][ii+kx];
                            }

                            // f2 min
                            if (f2map[jj+ky][ii+kx] < f2minval[cnt]) {
                                f2mini[cnt] = ii+kx;
                                f2minj[cnt] = jj+ky;
                                f2minval[cnt] = f2map[jj+ky][ii+kx];
                            }

                            // f2 max
                            if (f2map[jj+ky][ii+kx] > f2maxval[cnt]) {
                                f2maxi[cnt] = ii+kx;
                                f2maxj[cnt] = jj+ky;
                                f2maxval[cnt] = f2map[jj+ky][ii+kx];
                            }
                        }
                    }
                    cnt++;
                }
            }
    // cout << "[Non] 3*3 F" << endl;

            // 9*9 filter
            // init
            f1minval[9] = 32767;
            f2minval[9] = 32767;
            f1maxval[9] = -32767;
            f2maxval[9] = -32767;

            for (int k = 0; k < 9; k++) {
                // f1 min
                if (f1minval[k] < f1minval[9]) {
                    f1mini[9] = f1mini[k];
                    f1minj[9] = f1minj[k];
                    f1minval[9] = f1minval[k];
                }

                // f1 max
                if (f1maxval[k] > f1maxval[9]) {
                    f1maxi[9] = f1maxi[k];
                    f1maxj[9] = f1maxj[k];
                    f1maxval[9] = f1maxval[k];
                }

                // f2 min
                if (f2minval[k] < f2minval[9]) {
                    f2mini[9] = f2mini[k];
                    f2minj[9] = f2minj[k];
                    f2minval[9] = f2minval[k];
                }

                // f2 max
                if (f2maxval[k] > f2maxval[9]) {
                    f2maxi[9] = f2maxi[k];
                    f2maxj[9] = f2maxj[k];
                    f2maxval[9] = f2maxval[k];
                }
            }
    // cout << "[Non] 9*9 F" << endl;

            // check 3*3
            bool f1min_valid = true;
            bool f1max_valid = true;
            bool f2min_valid = true;
            bool f2max_valid = true;
            for (int k = 0; k < 9; k++) {
                f1min_valid = true;
                f1max_valid = true;
                f2min_valid = true;
                f2max_valid = true;
                // f1min
                // cout << f1minj[k] << "   " <<  f1mini[k] << endl;
                // cout << f1maxj[k] << "   " <<  f1maxi[k] << endl;
                // cout << f2minj[k] << "   " <<  f2mini[k] << endl;
                // cout << f2maxj[k] << "   " <<  f2maxi[k] << endl;
                for (int ky = -2; ky <= 2; ky++) {
                    for (int kx = -2; kx <= 2; kx++) {
                        if (f1map[f1minj[k]+ky][f1mini[k]+kx] < f1minval[k])
                            f1min_valid = false; 
                    }
                }

                // f1max
                for (int ky = -2; ky <= 2; ky++) {
                    for (int kx = -2; kx <= 2; kx++) {
                        if (f1map[f1maxj[k]+ky][f1maxi[k]+kx] > f1maxval[k])
                            f1max_valid = false; 
                    }
                }

                // f2min
                for (int ky = -2; ky <= 2; ky++) {
                    for (int kx = -2; kx <= 2; kx++) {
                        if (f2map[f2minj[k]+ky][f2mini[k]+kx] < f2minval[k])
                            f2min_valid = false; 
                    }
                }

                // f2max
                for (int ky = -2; ky <= 2; ky++) {
                    for (int kx = -2; kx <= 2; kx++) {
                        if (f2map[f2maxj[k]+ky][f2maxi[k]+kx] > f2maxval[k])
                            f2max_valid = false; 
                    }
                }
    // cout << "[Non] 9*9 bEF PUSH" << endl;

                // push back
                // f1min
                if (f1min_valid && (f1minval[k] <= -tau)) {
    // cout << "[Non] 9*9 bEF PUSH 0 " << endl;

                    maxima2.push_back(Matcher::maximum(i-9+f1mini[k],j-9+f1minj[k],f1minval[k],0));

                }
    // cout << "[Non] 9*9 bEF PUSH 0 done" << endl;
    

                // f1max
                if (f1max_valid && (f1maxval[k] >= tau)) {
    // cout << "[Non] 9*9 bEF PUSH 1" << endl;

                    maxima2.push_back(Matcher::maximum(i-9+f1maxi[k],j-9+f1maxj[k],f1maxval[k],1));
                }
    // cout << "[Non] 9*9 bEF PUSH 1 done" << endl;

                // f2min
                if (f2min_valid && (f2minval[k] <= -tau)) {
    // cout << "[Non] 9*9 bEF PUSH 2" << endl;

                    maxima2.push_back(Matcher::maximum(i-9+f2mini[k],j-9+f2minj[k],f2minval[k],2));
                }
    // cout << "[Non] 9*9 bEF PUSH 2 done" << endl;

                // f2max
                if (f2max_valid && (f2maxval[k] >= tau)) {
    // cout << "[Non] 9*9 bEF PUSH 3" << endl;

                    maxima2.push_back(Matcher::maximum(i-9+f2maxi[k],j-9+f2maxj[k],f2maxval[k],3));
                }
    // cout << "[Non] 9*9 bEF PUSH 3 done" << endl;

            }

    // cout << "[Non] 3*3 C" << endl;

            // check 9*9
            f1min_valid = true;
            f1max_valid = true;
            f2min_valid = true;
            f2max_valid = true;
            // f1min
            for (int ky = -8; ky <= 8; ky++) {
                for (int kx = -8; kx <= 8; kx++) {
                    if (f1map[f1minj[9]+ky][f1mini[9]+kx] < f1minval[9])
                        f1min_valid = false; 
                }
            }

            // f1max
            for (int ky = -8; ky <= 8; ky++) {
                for (int kx = -8; kx <= 8; kx++) {
                    if (f1map[f1maxj[9]+ky][f1maxi[9]+kx] > f1maxval[9])
                        f1max_valid = false; 
                }
            }

            // f2min
            for (int ky = -8; ky <= 8; ky++) {
                for (int kx = -8; kx <= 8; kx++) {
                    if (f2map[f2minj[9]+ky][f2mini[9]+kx] < f2minval[9])
                        f2min_valid = false; 
                }
            }

            // f2max
            for (int ky = -8; ky <= 8; ky++) {
                for (int kx = -8; kx <= 8; kx++) {
                    if (f2map[f2maxj[9]+ky][f2maxi[9]+kx] > f2maxval[9])
                        f2max_valid = false; 
                }
            }

            // push back
            // f1min
            if (f1min_valid && (f1minval[9] <= -tau)) {
                maxima1.push_back(Matcher::maximum(i-9+f1mini[9],j-9+f1minj[9],f1minval[9],0));
            }
            // f1max
            if (f1max_valid && (f1maxval[9] >= tau)) {
                maxima1.push_back(Matcher::maximum(i-9+f1maxi[9],j-9+f1maxj[9],f1maxval[9],1));
            }

            // f2min
            if (f2min_valid && (f2minval[9] <= -tau)) {
                maxima1.push_back(Matcher::maximum(i-9+f2mini[9],j-9+f2minj[9],f2minval[9],2));
            }
            // f2max
            if (f2max_valid && (f2maxval[9] >= tau)) {
                maxima1.push_back(Matcher::maximum(i-9+f2maxi[9],j-9+f2maxj[9],f2maxval[9],3));
            }
    // cout << "[Non] 9*9 C" << endl;

        }
    }
    // cout << "[Non] END" << endl;

}

void Matcher::myComputeFeatures (uint8_t *I,const int32_t* dims,int32_t* &max1,int32_t &num1,int32_t* &max2,int32_t &num2,uint8_t* &I_du,uint8_t* &I_dv,uint8_t* &I_du_full,uint8_t* &I_dv_full) {
  cout << endl << " MY VERSION" << endl;
  int16_t *I_f1;
  int16_t *I_f2;
  
  int32_t dims_matching[3];
  memcpy(dims_matching,dims,3*sizeof(int32_t));
  
  // allocate memory for sobel images and filter images

    uint8_t* I_matching = myCreateHalfResolutionImage(I,dims);
    myGetHalfResolutionDimensions(dims,dims_matching);
    I_du      = (uint8_t*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(uint8_t*),16);
    I_dv      = (uint8_t*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(uint8_t*),16);
    I_f1      = (int16_t*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(int16_t),16);
    I_f2      = (int16_t*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(int16_t),16);
    I_du_full = (uint8_t*)_mm_malloc(dims[2]*dims[1]*sizeof(uint8_t*),16);
    I_dv_full = (uint8_t*)_mm_malloc(dims[2]*dims[1]*sizeof(uint8_t*),16);
    // cout << endl<<"start Sobel!" << endl;
   mySobel5x5(I_matching,I_du,I_dv,dims_matching[2],dims_matching[1]);
    mySobel5x5(I,I_du_full,I_dv_full,dims[2],dims[1]);
    // cout << endl<<"start myBlob5x5!" << endl;

    myBlob5x5(I_matching,I_f1,dims_matching[2],dims_matching[1]);
    // cout << endl<<"start myCheckerboard5x5!" << endl;

    myCheckerboard5x5(I_matching,I_f2,dims_matching[2],dims_matching[1]);
    _mm_free(I_matching);
  
  
  // extract sparse maxima (1st pass) via non-maximum suppression
  vector<Matcher::maximum> maxima1;
  vector<Matcher::maximum> maxima2;
  // cout << endl<<"start myNonMaximumSuppression" << endl;
    myNonMaximumSuppression(I_f1,I_f2,dims_matching,maxima1, maxima2);
    myComputeDescriptors(I_du,I_dv,dims_matching[2],maxima1);

  
  // extract dense maxima (2nd pass) via non-maximum suppression


  myComputeDescriptors(I_du,I_dv,dims_matching[2],maxima2);

  // release filter images
  _mm_free(I_f1);
  _mm_free(I_f2);  
  
  // get number of interest points and init maxima pointer to NULL
  num1 = maxima1.size();
  num2 = maxima2.size();
  max1 = 0;
  max2 = 0;
  
  int32_t s = 1;
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
  //  cout << endl<<"end compute" << endl;
}