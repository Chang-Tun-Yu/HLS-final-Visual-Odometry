#include <iostream>
#include <stdint.h>


inline int32_t myGetAddressOffsetImage (const int32_t& u,const int32_t& v,const int32_t& width);

void myComputeDescriptor (const uint8_t* I_du,const uint8_t* I_dv,const int32_t &bpl,const int32_t &u,const int32_t &v,uint8_t *desc_addr);


void  mySobel5x5 ( const uint8_t* in, uint8_t* out_v, uint8_t* out_h, int w, int h );

void  myCheckerboard5x5 ( const uint8_t* in, int16_t* out, int w, int h );

void myBlob5x5 ( const uint8_t* in, int16_t* out, int w, int h );

void myNonMaximumSuppression_and_ComputeDescriptors (int16_t* I_f1,int16_t* I_f2,const int32_t* dims, uint8_t* I_du,uint8_t* I_dv, int32_t* max2, int32_t &num2);

void myComputeFeatures (uint8_t *I,const int32_t* dims, int32_t* max2,int32_t &num2);