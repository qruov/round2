#include "qruov.h"
#include <x86intrin.h>

/*
BYTES2BLOCKS(QRUOV_M)
(1<<SHIFT)

inline static Fq_vec Fq_vec_add(Fq_vec a, Fq_vec b){
    __m256i c = _mm256_add_epi32(a, b) ;
    __m256i d = _mm256_set1_epi32(QRUOV_q-1) ;
    __m256i e = _mm256_cmpgt_epi32(c,d) ;
    __m256i f = _mm256_srli_epi32(e, 32-QRUOV_ceil_log_2_q) ;
    __m256i g = _mm256_sub_epi32(c, f) ;
    return g ;
}
*/

inline static void vector_M_add(const vector_M A, const vector_M B, vector_M C){
  for(int j = 0 ; j < QRUOV_M ; j++) C[j] = Fq_add(A[j], B[j]) ;
}

inline static void vector_M_sub(const vector_M A, const vector_M B, vector_M C){
  for(int j = 0 ; j < QRUOV_M ; j++) C[j] = Fq_sub(A[j], B[j]) ;
}

inline static void vector_V_sub(const vector_V A, const vector_V B, vector_V C){
  for(int j = 0 ; j < QRUOV_V ; j++) C[j] = Fq_sub(A[j], B[j]) ;
}

inline static int64_t vector_V_dot_vector_V (const vector_V A, const vector_V B){
  int64_t C = 0 ;
  for(int i = 0 ; i < QRUOV_V ; i++) C += (int64_t)A[i] * (int64_t)B[i] ;
  return C ;
}

inline static int64_t vector_M_dot_vector_M (const vector_M A, const vector_M B){
  int64_t C = 0 ;
  for(int i = 0 ; i < QRUOV_M ; i++) C += (int64_t)A[i] * (int64_t)B[i] ;
  return C ;
}

// =============================================================================

void VECTOR_M_ADD(const VECTOR_M A, const VECTOR_M B, VECTOR_M C){
  for(int i = 0 ; i < QRUOV_L ; i++) vector_M_add(A[i], B[i], C[i]) ;
}

void VECTOR_M_SUB(const VECTOR_M A, const VECTOR_M B, VECTOR_M C){
  for(int i = 0 ; i < QRUOV_L ; i++) vector_M_sub(A[i], B[i], C[i]) ;
}

void VECTOR_V_SUB(const VECTOR_V A, const VECTOR_V B, VECTOR_V C){
  for(int i = 0 ; i < QRUOV_L ; i++) vector_V_sub(A[i], B[i], C[i]) ;
}

void VECTOR_V_dot_VECTOR_V (
  const VECTOR_V A,         // input
  const VECTOR_V B,         // input
        Fq       C[QRUOV_L] // output
){
  int64_t T[2*QRUOV_L - 1] ;
  memset(T, 0, sizeof(T)) ;

  for(int i = 0 ; i < QRUOV_L ; i++)
    for(int j = 0 ; j < QRUOV_L ; j++)
      T[i+j] += vector_V_dot_vector_V (A[i], B[j]) ;

  for(int i = 2*QRUOV_L-2; i >= QRUOV_L; i--){
    T[i-QRUOV_L]          += QRUOV_fc0 * T[i] ;
    T[i-QRUOV_L+QRUOV_fe] += QRUOV_fc  * T[i] ;
  }

  for(int i = 0; i < QRUOV_L; i++) C[i] = (Fq)(T[i] % QRUOV_q) ;
}


void VECTOR_M_dot_VECTOR_M (
  const VECTOR_M A,         // input
  const VECTOR_M B,         // input
        Fq       C[QRUOV_L] // output
){
  int64_t T[2*QRUOV_L - 1] ;
  memset(T, 0, sizeof(T)) ;

  for(int i = 0 ; i < QRUOV_L ; i++)
    for(int j = 0 ; j < QRUOV_L ; j++)
      T[i+j] += vector_M_dot_vector_M (A[i], B[j]) ;

  for(int i = 2*QRUOV_L-2; i >= QRUOV_L; i--){
    T[i-QRUOV_L]          += QRUOV_fc0 * T[i] ;
    T[i-QRUOV_L+QRUOV_fe] += QRUOV_fc  * T[i] ;
  }

  for(int i = 0; i < QRUOV_L; i++) C[i] = (Fq)(T[i] % QRUOV_q) ;
}

void VECTOR_V_MUL_SYMMETRIC_MATRIX_VxV(const VECTOR_V A, const MATRIX_VxV B, VECTOR_V C){
  Fq tmp [QRUOV_L] ;
  for(int i = 0 ; i < QRUOV_V ; i++){
    VECTOR_V_dot_VECTOR_V (A, B[i], tmp) ;
    for(int j = 0; j < QRUOV_L; j++) C[j][i] = tmp[j] ;
  }
}

void MATRIX_TRANSPOSE_VxM(const MATRIX_VxM A, MATRIX_MxV C){
  for(int i=0;i<QRUOV_V;i++){
    for(int k=0;k<QRUOV_L;k++){
      for(int j=0;j<QRUOV_M;j++){
        C[j][k][i] = A[i][k][j] ;
      }
    }
  }
}

void VECTOR_V_MUL_MATRIX_VxM(const VECTOR_V A, const MATRIX_VxM B, VECTOR_M C){
  MATRIX_MxV BT ;
  MATRIX_TRANSPOSE_VxM(B, BT) ;
  Fq tmp [QRUOV_L] ;
  for(int i = 0 ; i < QRUOV_M ; i++){
    VECTOR_V_dot_VECTOR_V (A, BT[i], tmp) ;
    for(int j = 0; j < QRUOV_L; j++) C[j][i] = tmp[j] ;
  }
}

void MATRIX_MxV_MUL_SYMMETRIC_MATRIX_VxV(const MATRIX_MxV A, const MATRIX_VxV B, MATRIX_MxV C){
  for(int i = 0 ; i < QRUOV_M ; i++){
    VECTOR_V_MUL_SYMMETRIC_MATRIX_VxV(A[i], B, C[i]) ;
  }
}

void MATRIX_MUL_MxV_VxM(const MATRIX_MxV A, const MATRIX_VxM B, MATRIX_MxM C){
  for(int i = 0 ; i < QRUOV_M ; i++){
    VECTOR_V_MUL_MATRIX_VxM(A[i], B, C[i]) ;
  }
}

void MATRIX_ADD_MxM(const MATRIX_MxM A, const MATRIX_MxM B, MATRIX_MxM C){
  for(int i = 0 ; i < QRUOV_M ; i++){
    VECTOR_M_ADD(A[i], B[i], C[i]) ;
  }
}

void MATRIX_MUL_ADD_MxV_VxM(const MATRIX_MxV A, const MATRIX_VxM B, MATRIX_MxM C){
  MATRIX_MxM T ;
  MATRIX_MUL_MxV_VxM(A, B, T) ;
  MATRIX_ADD_MxM(C, T, C) ;
}

void MATRIX_SUB_MxV(const MATRIX_MxV A, const MATRIX_MxV B, MATRIX_MxV C){
  for(int i = 0 ; i < QRUOV_M ; i++){
    VECTOR_V_SUB(A[i], B[i], C[i]) ;
  }
}

void SIG_GEN(const VECTOR_M oil, const MATRIX_VxM Sd, const VECTOR_V vineger, QRUOV_SIGNATURE sig){
  for(int i=0;i<QRUOV_V;i++){
    Fq u [QRUOV_L] ;
    for(int j=0;j<QRUOV_L;j++) u[j] = vineger[j][i] ;
    Fq t [QRUOV_L] ;
    VECTOR_M_dot_VECTOR_M(oil, Sd[i], t) ;
    sig->s[i] = Fql_sub(
      Fq2Fql_uint8_t_star(u),
      Fq2Fql_uint8_t_star(t)
    ) ;
  }
  for(int i=QRUOV_V;i<QRUOV_N;i++){
    Fq u [QRUOV_L] ;
    for(int j=0;j<QRUOV_L;j++) u[j] = oil[j][i-QRUOV_V] ;
    sig->s[i] = Fq2Fql_uint8_t_star(u) ;
  }
}

uint8_t VERIFY_i(const MATRIX_VxV Pi1, const MATRIX_VxM Pi2, const MATRIX_MxM Pi3, const VECTOR_M oil, const VECTOR_V vineger, const Fq msg_i){
  int j,k ;

  VECTOR_V tmp_v ;
  VECTOR_M tmp_o ;

  Fq t [QRUOV_L] ;
  Fq u [QRUOV_L] ;

  for(j=0;j<QRUOV_V;j++){
      VECTOR_M_dot_VECTOR_M(oil, Pi2[j], t) ;
      VECTOR_V_dot_VECTOR_V(vineger, Pi1[j], u) ;
      for(k=0;k<QRUOV_L;k++) tmp_v[k][j] = Fq_add(Fq_add(t[k],t[k]),u[k]) ;
  }

  for(j=0;j<QRUOV_M;j++){
      VECTOR_M_dot_VECTOR_M(oil, Pi3[j], t) ;
      for(k=0;k<QRUOV_L;k++) tmp_o[k][j] = t[k] ;
  }

  VECTOR_V_dot_VECTOR_V(vineger, tmp_v, t) ;
  VECTOR_M_dot_VECTOR_M(oil    , tmp_o, u) ;

  return msg_i == Fq_add(t[QRUOV_perm(0)],u[QRUOV_perm(0)]) ;
}
