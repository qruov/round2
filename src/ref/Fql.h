#pragma once
#include "qruov_misc.h"
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
typedef          __int128  INT128_T ;
typedef unsigned __int128 UINT128_T ;

#define Fq_reduction(X)    Fq_reduction_1(X)
#define Fql_reduction(X)   Fql_reduction_1(X)
#define Fql_acc_refresh(X) Fql_acc_refresh_1(X)
#define Fql_acc_reduce(X)  Fql_acc_reduce_1(X)
#define Fql_mul(X,Y)       Fql_mul_1(X,Y)

// ============================================================================
// F_q  (q = 2^c - 1)
// ============================================================================

typedef uint8_t Fq ;

// ============================================================================
// Fq_reduction
// ============================================================================

inline static int Fq_reduction_0(int Z){ return Z % QRUOV_q ; }

inline static int Fq_reduction_1(int Z){
      Z = (Z & QRUOV_q) + ((Z & ~QRUOV_q) >> QRUOV_ceil_log_2_q) ;
  int C = ((Z+1) & ~QRUOV_q) ;
      Z += (C>>QRUOV_ceil_log_2_q) ;
      Z -= C ;
  return Z ;
}

inline static int Fq_reduction_debug(int Z);

// ============================================================================
// Fq add/sub ...
// ============================================================================

inline static Fq Fq_add(Fq X, Fq Y){ return (Fq)Fq_reduction((int)X+(int)Y) ; }
inline static Fq Fq_sub(Fq X, Fq Y){ return (Fq)Fq_reduction((int)X-(int)Y+QRUOV_q) ; }
inline static Fq Fq_mul(Fq X, Fq Y){ return (Fq)Fq_reduction((int)X*(int)Y) ; }
inline static Fq Fq_inv(Fq X){ extern Fq Fq_inv_table[QRUOV_q] ; return Fq_inv_table[X] ; }

// ============================================================================
// hardware
// ============================================================================
// ws: word size
// wm: word mask

#if QRUOV_L == 3 
#  define Fql_ws               (24)
#  define Fql_wm               ((1ULL<<Fql_ws)-1)
#  define Fql_mask_(n)         (((uint64_t)QRUOV_q)<<(Fql_ws*(n)))
#  define Fql_mask             (Fql_mask_(0)|Fql_mask_(1)|Fql_mask_(2))
#  define Fql_mask_one_(n)     (((uint64_t)1ULL)<<(Fql_ws*(n)))
#  define Fql_mask_one         (Fql_mask_one_(0)|Fql_mask_one_(1)|Fql_mask_one_(2))
#  define Fql_2_mask_(n)       (((uint64_t)((1<<(QRUOV_ceil_log_2_q*2))-1))<<(Fql_ws*(n)))
#  define Fql_2_mask           (Fql_2_mask_(0)|Fql_2_mask_(1)|Fql_2_mask_(2))
#  define Fql_acc_mask_(n)     (((UINT128_T)QRUOV_q)<<(Fql_ws*(n)))
#  define Fql_acc_mask         (Fql_acc_mask_(0)|Fql_acc_mask_(1)|Fql_acc_mask_(2)|Fql_acc_mask_(3)|Fql_acc_mask_(4))
#  define Fql_acc_mask_one_(n) (((UINT128_T)1ULL)<<(Fql_ws*(n)))
#  define Fql_acc_mask_one     (Fql_acc_mask_one_(0)|Fql_acc_mask_one_(1)|Fql_acc_mask_one_(2)|Fql_acc_mask_one_(3)|Fql_acc_mask_one_(4))
#  define Fql_2_acc_mask_(n)   (((UINT128_T)((1<<(QRUOV_ceil_log_2_q*2))-1))<<(Fql_ws*(n)))
#  define Fql_2_acc_mask       (Fql_2_acc_mask_(0)|Fql_2_acc_mask_(1)|Fql_2_acc_mask_(2)|Fql_2_acc_mask_(3)|Fql_2_acc_mask_(4))
#  define Fql_U_SIZE 1
#elif QRUOV_L == 10 
#  define Fql_ws               (16)
#  define Fql_wm               ((1ULL<<Fql_ws)-1)
#  define Fql_mask_(n)         (((uint64_t)QRUOV_q)<<(Fql_ws*(n)))
#  define Fql_mask             (Fql_mask_(0)|Fql_mask_(1)|Fql_mask_(2)|Fql_mask_(3))
#  define Fql_mask_one_(n)     (((uint64_t)1ULL)<<(Fql_ws*(n)))
#  define Fql_mask_one         (Fql_mask_one_(0)|Fql_mask_one_(1)|Fql_mask_one_(2)|Fql_mask_one_(3))
#  define Fql_2_mask_(n)       (((uint64_t)((1<<(QRUOV_ceil_log_2_q*2))-1))<<(Fql_ws*(n)))
#  define Fql_2_mask           (Fql_2_mask_(0)|Fql_2_mask_(1)|Fql_2_mask_(2)|Fql_2_mask_(3))
#  define Fql_U_SIZE 3
#endif
#define Fql_AU_SIZE  (2*(Fql_U_SIZE))

typedef union Fql_union_t {
  uint64_t c64[Fql_U_SIZE*1] ;
  uint32_t c32[Fql_U_SIZE*2] ;
  uint16_t c16[Fql_U_SIZE*4] ;
  uint8_t  c8 [Fql_U_SIZE*8] ;
} Fql_union ;

typedef union Fql_acc_union_t {
  uint64_t  c64[Fql_AU_SIZE*1] ;
  uint32_t  c32[Fql_AU_SIZE*2] ;
  uint16_t  c16[Fql_AU_SIZE*4] ;
  uint8_t   c8 [Fql_AU_SIZE*8] ;
  Fql_union c                  ;
} Fql_acc_union ;

#if QRUOV_L == 3
typedef   uint64_t     Fql ;
typedef   UINT128_T    Fql_acc ;
#  define Fql_zero     ((Fql)0)
#  define Fql_acc_zero ((Fql_acc)0)
#elif QRUOV_L == 10
typedef Fql_union      Fql ;
typedef Fql_acc_union  Fql_acc ;
extern  Fql            Fql_zero ;
extern  Fql_acc        Fql_acc_zero ;
#endif

// ============================================================================
// F_q^L House keeping
// ============================================================================

inline static void Fql_fprint_n(FILE *stream, int n, char * header, void * A_){
  Fql_acc_union * A = (Fql_acc_union *) A_ ;
  fprintf(stream, "%s",header) ;
  for(int i=n-1;i>=0;i--)fprintf(stream, "%016lx", A->c64[i]) ;
  fprintf(stream, "\n") ;
}

inline static void Fql_print_n  (int n, char * header, void * A_){ Fql_fprint_n(stderr, n,   header, A_) ; }
inline static void Fql_print    (       char * header, Fql A    ){ Fql_print_n (Fql_U_SIZE,  header, &A) ; }
inline static void Fql_acc_print(       char * header, Fql_acc A){ Fql_print_n (Fql_AU_SIZE, header, &A) ; }

#define Fql_PRINT(a)      Fql_print(#a " = ", a)
#define Fql_acc_PRINT(a)  Fql_acc_print(#a " = ", a)

inline static int Fql_eq(Fql a, Fql b){ return memcmp(&a, &b, sizeof(Fql)) == 0 ; }
inline static int Fql_ne(Fql a, Fql b){ return ! Fql_eq(a, b) ; }
inline static int Fql_acc_eq(Fql_acc a, Fql_acc b){ return memcmp(&a, &b, sizeof(Fql_acc)) == 0 ; }
inline static int Fql_acc_ne(Fql_acc a, Fql_acc b){ return ! Fql_acc_eq(a, b) ; }

#if QRUOV_L ==  3
inline static Fq  Fql2Fq(Fql Z, int i){ return ((Z >> ((Fql_ws)*i)) & QRUOV_q) ; }
inline static Fql Fq2Fql(Fql z0, Fql z1, Fql z2){ return z0|(z1<<Fql_ws)|(z2<<(Fql_ws*2)) ; }
inline static Fql Fq2Fql_uint8_t_star(uint8_t c[QRUOV_L]){ return Fq2Fql(c[0],c[1],c[2]) ; }
inline static Fql_acc Fq2Fql_acc(Fql z0, Fql z1, Fql z2, Fql z3, Fql z4){
  return ((UINT128_T)z0<<(Fql_ws*0))|
         ((UINT128_T)z1<<(Fql_ws*1))|
	 ((UINT128_T)z2<<(Fql_ws*2))|
	 ((UINT128_T)z3<<(Fql_ws*3))|
	 ((UINT128_T)z4<<(Fql_ws*4));
}

#elif QRUOV_L == 10
#  if   __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#    define WORD_ORDER(i)   (i)
#  elif   __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
#    define WORD_ORDER(i)   (((i>>2)<<2)+(3-(i&3)))
#  else 
#    error "unsupported WORD_ORDER()"
#  endif
inline static Fq  Fql2Fq(Fql Z, int i){ return Z.c16[WORD_ORDER(i)] ; }
inline static Fql Fq2Fql(uint16_t c[QRUOV_L]){
  Fql Z ;
  Z.c64[Fql_U_SIZE-1] = 0 ;
  for(int i=0; i<QRUOV_L; i++) Z.c16[WORD_ORDER(i)] = c[i] ;
  return Z ;
}
inline static Fql Fq2Fql_uint8_t_star(uint8_t c[QRUOV_L]){
  Fql Z ;
  Z.c64[Fql_U_SIZE-1] = 0 ;
  for(int i=0; i<QRUOV_L; i++) Z.c16[WORD_ORDER(i)] = c[i] ;
  return Z ;
}
inline static Fql_acc Fq2Fql_acc(uint16_t c[2*QRUOV_L-1]){
  Fql_acc Z ;
  Z.c64[Fql_AU_SIZE-1] = 0 ;
  for(int i=0; i<2*QRUOV_L-1; i++) Z.c16[WORD_ORDER(i)] = c[i] ;
  return Z ;
}

#endif

#if   (QRUOV_q == 127) && (QRUOV_L == 3)
#  include "Fql_L3.h"
#elif (QRUOV_q ==  31) && (QRUOV_L == 3)
#  include "Fql_L3.h"
#elif (QRUOV_q ==  31) && (QRUOV_L == 10)
#  include "Fql_L10.h"
#  include "Fql_q31L10.h"
#elif (QRUOV_q ==   7) && (QRUOV_L == 10)
#  include "Fql_L10.h"
#  include "Fql_q7L10.h"
#else
#  error "unknown (QRUOV_q, QRUOV_L)"
#endif
