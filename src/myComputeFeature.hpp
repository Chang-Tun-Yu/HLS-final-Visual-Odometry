#include <iostream>
#include <stdint.h>


inline int32_t myGetAddressOffsetImage (const int32_t u,const int32_t v,const int32_t width);

void myComputeDescriptor (const uint8_t I_du[290816],const uint8_t I_dv[290816],const int32_t bpl,const int32_t u,const int32_t v,uint8_t desc_addr[32]);


void  mySobel5x5 ( const uint8_t in[290816], uint8_t out_v[290816], uint8_t out_h[290816], int w, int h );

void  myCheckerboard5x5 ( const uint8_t in[290816], int16_t out[290816], int w, int h );

void myBlob5x5 ( const uint8_t in[290816], int16_t out[290816], int w, int h );

void myNonMaximumSuppression_and_ComputeDescriptors (int16_t I_f1[290816],int16_t I_f2[290816],const int32_t dims[290816], uint8_t I_du[290816],uint8_t I_dv[290816], int32_t max2[288000], int32_t &num2);

void myComputeFeatures (uint8_t I[290816], int32_t max2[288000],int32_t &num2);