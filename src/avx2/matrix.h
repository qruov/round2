#pragma once
#include "qruov_misc.h"
#include "Fql.h"

typedef uint8_t QRUOV_SEED  [QRUOV_SEED_LEN] ;
typedef uint8_t QRUOV_SALT  [QRUOV_SALT_LEN] ;

TYPEDEF_STRUCT(QRUOV_SIGNATURE,
  QRUOV_SALT r           ;
  Fql        s [QRUOV_N] ;
) ;

/*
#define ALIGN_8BIT   __attribute((__aligned__( 1))) ; //   8 bit // SHIFT = 0 ref
#define ALIGN_16BIT  __attribute((__aligned__( 2))) ; //  16 bit // SHIFT = 1
#define ALIGN_32BIT  __attribute((__aligned__( 4))) ; //  32 bit // SHIFT = 2
#define ALIGN_64BIT  __attribute((__aligned__( 8))) ; //  64 bit // SHIFT = 3
#define ALIGN_128BIT __attribute((__aligned__(16))) ; // 128 bit // SHIFT = 4
#define ALIGN_256BIT __attribute((__aligned__(32))) ; // 256 bit // SHIFT = 5 avx2
#define ALIGN_512BIT __attribute((__aligned__(64))) ; // 512 bit // SHIFT = 6 avx512
*/

#define SHIFT 5
#define BYTES2BLOCKS(BYTES)          (((BYTES)>>(SHIFT))+(((BYTES)&((1<<SHIFT)-1))?1:0))
#define BLOCK_ALIGNED_BYTES(BYTES)   (BYTES2BLOCKS(BYTES)<<(SHIFT))

typedef Fq         vector_V        [BLOCK_ALIGNED_BYTES(QRUOV_V)]  ALIGN_256BIT ;
typedef Fq         vector_M        [BLOCK_ALIGNED_BYTES(QRUOV_M)]  ALIGN_256BIT ;

typedef vector_V   VECTOR_V        [QRUOV_L]          ;
typedef vector_M   VECTOR_M        [QRUOV_L]          ;

typedef VECTOR_V   MATRIX_VxV      [QRUOV_V]          ;
typedef VECTOR_V   MATRIX_MxV      [QRUOV_M]          ;
typedef VECTOR_M   MATRIX_VxM      [QRUOV_V]          ;
typedef VECTOR_M   MATRIX_MxM      [QRUOV_M]          ;

typedef MATRIX_VxM QRUOV_Sd                           ;
typedef MATRIX_MxV QRUOV_SdT                          ;
typedef MATRIX_VxV QRUOV_P1        [QRUOV_m]          ;
typedef MATRIX_VxM QRUOV_P2        [QRUOV_m]          ;
typedef MATRIX_MxV QRUOV_P2T       [QRUOV_m]          ;
typedef MATRIX_MxM QRUOV_P3        [QRUOV_m]          ;

extern void VECTOR_M_ADD(const VECTOR_M A, const VECTOR_M B, VECTOR_M C) ;
extern void VECTOR_M_SUB(const VECTOR_M A, const VECTOR_M B, VECTOR_M C) ;
extern void VECTOR_V_SUB(const VECTOR_V A, const VECTOR_V B, VECTOR_V C) ;
extern void VECTOR_V_dot_VECTOR_V (const VECTOR_V A, const VECTOR_V B, Fq C[QRUOV_L]) ;
extern void VECTOR_M_dot_VECTOR_M (const VECTOR_M A, const VECTOR_M B, Fq C[QRUOV_L]) ;
extern void VECTOR_V_MUL_SYMMETRIC_MATRIX_VxV(const VECTOR_V A, const MATRIX_VxV B, VECTOR_V C) ;
extern void MATRIX_TRANSPOSE_VxM(const MATRIX_VxM A, MATRIX_MxV C) ;
extern void VECTOR_V_MUL_MATRIX_VxM(const VECTOR_V A, const MATRIX_VxM B, VECTOR_M C) ;
extern void MATRIX_MxV_MUL_SYMMETRIC_MATRIX_VxV(const MATRIX_MxV A, const MATRIX_VxV B, MATRIX_MxV C) ;
extern void MATRIX_MUL_MxV_VxM(const MATRIX_MxV A, const MATRIX_VxM B, MATRIX_MxM C) ;
extern void MATRIX_ADD_MxM(const MATRIX_MxM A, const MATRIX_MxM B, MATRIX_MxM C) ;
extern void MATRIX_MUL_ADD_MxV_VxM(const MATRIX_MxV A, const MATRIX_VxM B, MATRIX_MxM C) ;
extern void MATRIX_SUB_MxV(const MATRIX_MxV A, const MATRIX_MxV B, MATRIX_MxV C) ;
extern void SIG_GEN(const VECTOR_M oil, const MATRIX_VxM Sd, const VECTOR_V vineger, QRUOV_SIGNATURE sig) ;
extern uint8_t VERIFY_i(const MATRIX_VxV Pi1, const MATRIX_VxM Pi2, const MATRIX_MxM Pi3, const VECTOR_M oil, const VECTOR_V vineger, const Fq msg_i) ;
