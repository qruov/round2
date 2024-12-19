#include "qruov.h"
#include <x86intrin.h>

/*
BYTES2BLOCKS(QRUOV_M)
(1<<SHIFT)
*/

typedef __m512i Fq_vec ;

inline static Fq_vec Fq_vec_add(Fq_vec a, Fq_vec b){
  __m512i   q = _mm512_set1_epi8(QRUOV_q)      ; // Fq_vec_PRINT(q) ;
  __m512i   c = _mm512_adds_epu8(a, b)         ; // Fq_vec_PRINT(c) ;
  __mmask64 m = _mm512_cmple_epu8_mask(q, c)   ; // printf("m = %016llx\n",m);
  __m512i   d = _mm512_mask_subs_epu8(c,m,c,q) ; // Fq_vec_PRINT(d) ;
  return d ;
}

inline static Fq_vec Fq_vec_sub(Fq_vec a, Fq_vec b){
  __m512i   q = _mm512_set1_epi8(QRUOV_q)      ; // Fq_vec_PRINT(q) ;
  __mmask64 m = _mm512_cmplt_epu8_mask(a, b)   ; // printf("m = %016llx\n",m);
  __m512i   c = _mm512_mask_adds_epu8(a,m,a,q) ; // Fq_vec_PRINT(c) ;
  __m512i   d = _mm512_subs_epu8(c, b)         ; // Fq_vec_PRINT(d) ;
  return d ;
}

inline static Fq_vec Fq_vec_dot_accumulate(Fq_vec a, Fq_vec b, Fq_vec c){
  return _mm512_dpbusd_epi32(c,a,b) ;
}

inline static int64_t Fq_vec_dot_accumulate_final(Fq_vec c){
  __m256i lo0 = _mm512_extracti64x4_epi64(c, 0) ;
  __m256i hi0 = _mm512_extracti64x4_epi64(c, 1) ;
  __m256i f   = _mm256_add_epi32(lo0, hi0)     ;
// --
  __m128i lo = _mm256_castsi256_si128(f) ;
  __m128i hi = _mm256_extracti128_si256(f, 1) ;
  __m128i s1 = _mm_hadd_epi32(lo, hi) ;
  __m128i s2 = _mm_hadd_epi32(s1, s1) ;
  __m128i s3 = _mm_hadd_epi32(s2, s2) ;
  return ((uint64_t)_mm_extract_epi32(s3, 0)) ;
}

#if 0


inline static void Fq_vec_print(char * header, Fq_vec a_){
  uint8_t * a = (uint8_t *) & a_ ;
  printf("%s",header);
  for(int i=(1<<SHIFT)-1; i>=0; i--) printf("%02x ", a[i]);
  printf("\n");
}

#define   Fq_vec_PRINT(a)    Fq_vec_print( #a " = ", a )

#endif

inline static void vector_M_add_avx2(const vector_M A, const vector_M B, vector_M C){
  Fq_vec * a = (Fq_vec *) A ;
  Fq_vec * b = (Fq_vec *) B ;
  Fq_vec * c = (Fq_vec *) C ;
  for(int j = 0 ; j < BYTES2BLOCKS(QRUOV_M) ; j++) *c++ = Fq_vec_add(*a++, *b++) ; // , c++) ;
}

#if 1
#  define vector_M_add vector_M_add_avx2
#else
inline static void vector_M_add(const vector_M A, const vector_M B, vector_M C){
  vector_M A2 ; memcpy(A2, A, sizeof(vector_M)) ;
  vector_M B2 ; memcpy(B2, B, sizeof(vector_M)) ;
  vector_M C2 ;
  for(int j = 0 ; j < QRUOV_M ; j++) C[j] = Fq_add(A[j], B[j]) ;
#if 0
  { // debug
    vector_M_add_avx2(A2, B2, C2) ;
    for(int j = 0 ; j < QRUOV_M ; j++) if(C[j] != C2[j]){
      printf("A = ") ; for(int k = 0 ; k < QRUOV_M ; k++) printf("%02x ", A2[k]) ; printf("\n") ;
      printf("B = ") ; for(int k = 0 ; k < QRUOV_M ; k++) printf("%02x ", B2[k]) ; printf("\n") ;
      printf("C = ") ; for(int k = 0 ; k < QRUOV_M ; k++) printf("%02x ", C [k]) ; printf("\n") ;
      printf("C2= ") ; for(int k = 0 ; k < QRUOV_M ; k++) printf("%02x ", C2[k]) ; printf("\n") ;
      return ;
    }
  }
#endif
}
#endif

#if 1
inline static void vector_M_sub(const vector_M A, const vector_M B, vector_M C){
  Fq_vec * a = (Fq_vec *) A ;
  Fq_vec * b = (Fq_vec *) B ;
  Fq_vec * c = (Fq_vec *) C ;
  for(int j = 0 ; j < BYTES2BLOCKS(QRUOV_M) ; j++) *c++ = Fq_vec_sub(*a++, *b++) ; // , c++) ;
}
#else
inline static void vector_M_sub(const vector_M A, const vector_M B, vector_M C){
  for(int j = 0 ; j < QRUOV_M ; j++) C[j] = Fq_sub(A[j], B[j]) ;
}
#endif


#if 1
inline static void vector_V_sub(const vector_V A, const vector_V B, vector_V C){
  Fq_vec * a = (Fq_vec *) A ;
  Fq_vec * b = (Fq_vec *) B ;
  Fq_vec * c = (Fq_vec *) C ;
  for(int j = 0 ; j < BYTES2BLOCKS(QRUOV_V) ; j++) *c++ = Fq_vec_sub(*a++, *b++) ; // , c++) ;
}
#else
inline static void vector_V_sub(const vector_V A, const vector_V B, vector_V C){
  for(int j = 0 ; j < QRUOV_V ; j++) C[j] = Fq_sub(A[j], B[j]) ;
}
#endif

#if 1
inline static int64_t vector_V_dot_vector_V (const vector_V A, const vector_V B){
  Fq_vec * a = (Fq_vec *) A ;
  Fq_vec * b = (Fq_vec *) B ;
  Fq_vec   c ; c = _mm512_xor_epi32(c,c) ;
  for(int j = 0 ; j < BYTES2BLOCKS(QRUOV_V) ; j++) c = Fq_vec_dot_accumulate(*a++, *b++, c) ;
  int64_t C = Fq_vec_dot_accumulate_final(c) ;
  return C ;
}
#else
inline static int64_t vector_V_dot_vector_V (const vector_V A, const vector_V B){
  int64_t C = 0 ;
  for(int i = 0 ; i < QRUOV_V ; i++) C += (int64_t)A[i] * (int64_t)B[i] ;
  return C ;
}
#endif

#if 1
inline static int64_t vector_M_dot_vector_M(const vector_V A, const vector_V B){
  Fq_vec * a = (Fq_vec *) A ;
  Fq_vec * b = (Fq_vec *) B ;
  Fq_vec   c ; c = _mm512_xor_epi32(c,c) ;
  for(int j = 0 ; j < BYTES2BLOCKS(QRUOV_M) ; j++) c = Fq_vec_dot_accumulate(*a++, *b++, c) ;
  int64_t C = Fq_vec_dot_accumulate_final(c) ;
  return C ;
}
#else
inline static int64_t vector_M_dot_vector_M (const vector_M A, const vector_M B){
  int64_t C = 0 ;
  for(int i = 0 ; i < QRUOV_M ; i++) C += (int64_t)A[i] * (int64_t)B[i] ;
  return C ;
}
#endif

inline static void print_vector_M (char * header, const vector_M A){
  printf("%s",header) ;
  for(int k = 0 ; k < QRUOV_M ; k++) printf("%02x ", A[k]) ;
  printf(" | ");
  for(int k = QRUOV_M ; k < BLOCK_ALIGNED_BYTES(QRUOV_M) ; k++) printf("%02x ", A[k]) ;
  printf("\n");
}

#define PRINT_vector_M(a) print_vector_M( #a " = ", a )

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
    for(int k = 0; k < QRUOV_L; k++) C[k][i] = tmp[k] ;
  }
  for(int k=0;k<QRUOV_L;k++) memset(C[k]+QRUOV_V, 0, BLOCK_ALIGNED_BYTES(QRUOV_V)-QRUOV_V) ;
}

void MATRIX_TRANSPOSE_VxM(const MATRIX_VxM A, MATRIX_MxV C){
  for(int i=0;i<QRUOV_V;i++){
    for(int k=0;k<QRUOV_L;k++){
      for(int j=0;j<QRUOV_M;j++){
        C[j][k][i] = A[i][k][j] ;
      }
    }
  }
  for(int j=0;j<QRUOV_M;j++){
    for(int k=0;k<QRUOV_L;k++) memset(C[j][k]+QRUOV_V, 0, BLOCK_ALIGNED_BYTES(QRUOV_V)-QRUOV_V) ;
  }
}

void VECTOR_V_MUL_MATRIX_VxM(const VECTOR_V A, const MATRIX_VxM B, VECTOR_M C){
  MATRIX_MxV BT ;
  MATRIX_TRANSPOSE_VxM(B, BT) ;
  Fq tmp [QRUOV_L] ;
  for(int i = 0 ; i < QRUOV_M ; i++){
    VECTOR_V_dot_VECTOR_V (A, BT[i], tmp) ;
    for(int k = 0; k < QRUOV_L; k++) C[k][i] = tmp[k] ;
  }
  for(int k=0;k<QRUOV_L;k++) memset(C[k]+QRUOV_M, 0, BLOCK_ALIGNED_BYTES(QRUOV_M)-QRUOV_M) ;
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
