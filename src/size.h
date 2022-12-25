#ifndef __SIZE_H__
#define __SIZE_H__

#define IMG_W                   1024
#define IMG_H                   284
#define IMG_SIZE                290816
#define MAX_FEATURE_POINT_SZIE  24000


#define BIN_W                   60
#define BIN_H                   45
#define U_BIN_NUM               17
#define V_BIN_NUM               7
#define SEARCH_RAD_U            BIN_W*3
#define SEARCH_RAD_V            BIN_H*3

#define BIN_NUM                 U_BIN_NUM*V_BIN_NUM*4
#define MAX_FP_IN_BIN           160 // 158
#define MAX2_BIN_OFFSET         MAX_FP_IN_BIN*12 
#define MAX_FEATURE_ARRAY_SZIE  BIN_NUM*MAX_FP_IN_BIN*12



#endif