//--------------------------------------------------------------------
// Reed-Solomon C library - rs.h
//
// Wojciech Kaczmarski, SP5WWP
// 29 September 2024
//--------------------------------------------------------------------
#pragma once

#ifdef __cplusplus
extern "C" {
#endif
#include <stdint.h>
#include <math.h>

#define MAX_ITEMS   255

typedef struct
{
    uint8_t m;              //RS code over GF(2**m)
    uint8_t n;              //codeword length (in symbols)
    uint8_t t;              //erroneous symbols that can be corrected
    uint8_t k;              //message length, k=n-2*t
    uint8_t *p;             //primitive polynomial
    int alpha[MAX_ITEMS];   //powers of alpha array
    int index[MAX_ITEMS];   //?
    int g[MAX_ITEMS];       //generator polynomial
} rs_t;

void gen_GF(rs_t *rs);
void gen_poly(rs_t *rs);
void init_RS(rs_t *rs, uint8_t cw_len, uint8_t dt_len, uint8_t *poly);
void encode_RS(rs_t *rs, uint8_t *out, uint8_t *inp);
void decode_RS(rs_t *rs, int8_t *inp);

#ifdef __cplusplus
}
#endif