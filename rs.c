//--------------------------------------------------------------------
// Reed-Solomon C library - rs.c
//
// Based on Simon Rockliff's code, 26th June 1991
//
// Wojciech Kaczmarski, SP5WWP
// 29 September 2024
//--------------------------------------------------------------------
#include "rs.h"

/*
    Generate GF(2**mm) from the irreducible polynomial p(X) in p[0]..p[mm]
    lookup tables:  index->polynomial form   alpha[] contains j=alpha**i;
                    polynomial form -> index form  index[j=alpha**i] = i
    alpha=2 is the primitive element of GF(2**mm)
*/
void gen_GF(rs_t *rs)
{
    uint8_t mask = 1;

    rs->alpha[rs->m] = 0;
    for(uint8_t i=0; i<rs->m; i++)
    {
        rs->alpha[i] = mask;
        rs->index[rs->alpha[i]] = i;
        if(rs->p[i])
            rs->alpha[rs->m] ^= mask;
        mask <<= 1;
    }

    rs->index[rs->alpha[rs->m]] = rs->m;
    mask >>= 1;
    for(uint8_t i=rs->m+1; i<rs->n; i++)
    {
        if (rs->alpha[i-1] >= mask)
            rs->alpha[i] = rs->alpha[rs->m] ^ ((rs->alpha[i-1]^mask)<<1) ;
        else rs->alpha[i] = rs->alpha[i-1]<<1 ;
            rs->index[rs->alpha[i]] = i;
    }

    rs->index[0] = -1;
}

/*
    Obtain the generator polynomial of the tt-error correcting, length
    nn=(2**m -1) Reed Solomon code  from the product of (X+alpha**i), i=1..2*t
*/
void gen_poly(rs_t *rs)
{
    rs->g[0] = 2;    //primitive element alpha = 2  for GF(2**mm)
    rs->g[1] = 1;    //g(x) = (X+alpha) initially

    for(uint8_t i=2; i<=(rs->n)-(rs->k); i++)
    {
        rs->g[i] = 1;

        for(uint8_t j=i-1; j>0; j--)
        {
            if(rs->g[j])
                rs->g[j] = rs->g[j-1]^ rs->alpha[(rs->index[rs->g[j]]+i)%rs->n];
            else
                rs->g[j] = rs->g[j-1];
        }

        rs->g[0] = rs->alpha[(rs->index[rs->g[0]]+i)%(rs->n)];     //g[0] can never be zero
    }

    //convert g[] to index form for quicker encoding
    for(uint8_t i=0; i<=(rs->n)-(rs->k); i++)
        rs->g[i] = rs->index[rs->g[i]];
}

/*
    Init the RS de/encoder
*/
void init_RS(rs_t *rs, uint8_t cw_len, uint8_t dt_len, uint8_t *poly)
{
    //TODO: add required assertions here

    rs->m=log2(cw_len+1);
    rs->n=cw_len;
    rs->t=(cw_len-dt_len)/2;
    rs->k=cw_len-2*rs->t;
    rs->p=poly;

    gen_GF(rs);
    gen_poly(rs);
}

/*
    Take the string of symbols in inp[i], i=0..(k-1) and encode systematically
    to produce 2*tt parity symbols in out[0]..out[2*tt-1]
    inp[] is input and out[] is output in polynomial form.
    Encoding is done by using a feedback shift register with appropriate
    connections specified by the elements of g[], which was generated above.
    Codeword is   c(X) = inp(X)*X**(nn-kk)+ b(X)
*/
void encode_RS(rs_t *rs, uint8_t *out, uint8_t *inp)
{
    int feedback;

    for(uint8_t i=0; i<(rs->n)-(rs->k); i++)
        out[i] = 0;

    for(int16_t i=(rs->k)-1; i>=0; i--)
    {
        feedback = rs->index[inp[i]^out[(rs->n)-(rs->k)-1]] ;
        if(feedback != -1)
        {
            for(int16_t j=(rs->n)-(rs->k)-1; j>0; j--)
            {
                if (rs->g[j] != -1)
                    out[j] = out[j-1]^rs->alpha[(rs->g[j]+feedback)%(rs->n)] ;
                else
                    out[j] = out[j-1] ;
            }

            out[0] = rs->alpha[(rs->g[0]+feedback)%(rs->n)] ;
        }
        else
        {
            for(int16_t j=(rs->n)-(rs->k)-1; j>0; j--)
                out[j] = out[j-1] ;
            out[0] = 0 ;
        }
    }
}

