//--------------------------------------------------------------------
// Reed-Solomon C library - rs.h
//
// Wojciech Kaczmarski, SP5WWP
// 6 September 2026
//--------------------------------------------------------------------
#pragma once

#ifdef __cplusplus
extern "C"
{
#endif
#include <stdint.h>
#include <math.h>

#define MAX_ITEMS 255

  typedef struct
  {
    uint8_t m;            // RS code over GF(2**m)
    uint8_t n;            // codeword length (in symbols)
    uint8_t t;            // erroneous symbols that can be corrected
    uint8_t k;            // message length, k=n-2*t
    uint8_t *p;           // primitive polynomial
    int alpha[MAX_ITEMS]; // powers of alpha array
    int index[MAX_ITEMS]; //?
    int g[MAX_ITEMS];     // generator polynomial
  } rs_t;

  typedef enum
  {
    RS_INIT_OK,
    RS_INIT_CW_DT_MISMATCH,
    RS_INIT_T_ODD,
    RS_INIT_K_ZERO,
    RS_INIT_CW_LEN_POW_2,
  } rs_init_t;

  typedef enum
  {
    RS_NO_ERROR,
    RS_CORRECTED,
    RS_UNCORRECTABLE,
  } rs_status_t;

  void gen_GF(rs_t *rs);
  void gen_poly(rs_t *rs);
  rs_init_t init_RS(rs_t *rs, uint8_t cw_len, uint8_t dt_len, uint8_t *poly);
  void encode_RS(rs_t *rs, uint8_t *out, uint8_t *inp);
  rs_status_t decode_RS(rs_t *rs, uint8_t *inp);

#ifdef __cplusplus
}
#endif