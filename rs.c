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
  for (uint8_t i = 0; i < rs->m; i++)
  {
    rs->alpha[i] = mask;
    rs->index[rs->alpha[i]] = i;
    if (rs->p[i])
      rs->alpha[rs->m] ^= mask;
    mask <<= 1;
  }

  rs->index[rs->alpha[rs->m]] = rs->m;
  mask >>= 1;
  for (uint8_t i = rs->m + 1; i < rs->n; i++)
  {
    if (rs->alpha[i - 1] >= mask)
      rs->alpha[i] = rs->alpha[rs->m] ^ ((rs->alpha[i - 1] ^ mask) << 1);
    else
      rs->alpha[i] = rs->alpha[i - 1] << 1;
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
  rs->g[0] = 2; // primitive element alpha = 2  for GF(2**mm)
  rs->g[1] = 1; // g(x) = (X+alpha) initially

  for (uint8_t i = 2; i <= (rs->n) - (rs->k); i++)
  {
    rs->g[i] = 1;

    for (uint8_t j = i - 1; j > 0; j--)
    {
      if (rs->g[j])
        rs->g[j] = rs->g[j - 1] ^ rs->alpha[(rs->index[rs->g[j]] + i) % rs->n];
      else
        rs->g[j] = rs->g[j - 1];
    }

    rs->g[0] = rs->alpha[(rs->index[rs->g[0]] + i) % (rs->n)]; // g[0] can never be zero
  }

  // convert g[] to index form for quicker encoding
  for (uint8_t i = 0; i <= (rs->n) - (rs->k); i++)
    rs->g[i] = rs->index[rs->g[i]];
}

/*
    Init the RS de/encoder
*/
void init_RS(rs_t *rs, uint8_t cw_len, uint8_t dt_len, uint8_t *poly)
{
  // TODO: add required assertions here

  rs->m = log2(cw_len + 1);
  rs->n = cw_len;
  rs->t = (cw_len - dt_len) / 2;
  rs->k = cw_len - 2 * rs->t;
  rs->p = poly;

  gen_GF(rs);
  gen_poly(rs);
}

/*
    Take the string of symbols in inp[i], i=0..(k-1) and encode systematically
    to produce 2*t parity symbols in out[0]..out[2*t-1]
    inp[] is input and out[] is output in polynomial form.
    Encoding is done by using a feedback shift register with appropriate
    connections specified by the elements of g[], which was generated above.
    Codeword is   c(X) = inp(X)*X**(n-k)+ b(X)
*/
void encode_RS(rs_t *rs, uint8_t *out, uint8_t *inp)
{
  int feedback;

  for (uint8_t i = 0; i < (rs->n) - (rs->k); i++)
    out[i] = 0;

  for (int16_t i = (rs->k) - 1; i >= 0; i--)
  {
    feedback = rs->index[inp[i] ^ out[(rs->n) - (rs->k) - 1]];
    if (feedback != -1)
    {
      for (int16_t j = (rs->n) - (rs->k) - 1; j > 0; j--)
      {
        if (rs->g[j] != -1)
          out[j] = out[j - 1] ^ rs->alpha[(rs->g[j] + feedback) % (rs->n)];
        else
          out[j] = out[j - 1];
      }

      out[0] = rs->alpha[(rs->g[0] + feedback) % (rs->n)];
    }
    else
    {
      for (int16_t j = (rs->n) - (rs->k) - 1; j > 0; j--)
        out[j] = out[j - 1];
      out[0] = 0;
    }
  }
}

/*
    Assume we have received bits grouped into mm-bit symbols in ino[i],
    i=0..(n-1),  and inp[i] is index form (ie as powers of alpha).
    We first compute the 2*t syndromes by substituting alpha**i into rec(X) and
    evaluating, storing the syndromes in s[i], i=1..2t (leave s[0] zero) .
    Then we use the Berlekamp iteration to find the error location polynomial
    elp[i].   If the degree of the elp is >t, we cannot correct all the errors
    and hence just put out the information symbols uncorrected. If the degree of
    elp is <=t, we substitute alpha**i , i=1..n into the elp to get the roots,
    hence the inverse roots, the error location numbers. If the number of errors
    located does not equal the degree of the elp, we have more than t errors
    and cannot correct them.  Otherwise, we then solve for the error value at
    the error location and correct the error.  The procedure is that found in
    Lin and Costello. For the cases where the number of errors is known to be too
    large to correct, the information symbols as received are output (the
    advantage of systematic encoding is that hopefully some of the information
    symbols will be okay and that if we are in luck, the errors are in the
    parity part of the transmitted codeword).  Of course, these insoluble cases
    can be returned as error flags to the calling routine if desired.

    Edit: this function accepts data in the polynomial form and overwrites the input buffer
*/
void decode_RS(rs_t *rs, int8_t *inp)
{
  int elp[(rs->n) - (rs->k) + 2][(rs->n) - (rs->k)], d[(rs->n) - (rs->k) + 2], l[(rs->n) - (rs->k) + 2], u_lu[(rs->n) - (rs->k) + 2], s[(rs->n) - (rs->k) + 1];
  int count = 0, syn_error = 0, root[rs->t], loc[rs->t], z[rs->t + 1], err[rs->n], reg[rs->t + 1];

  // convert to polynomial form
  for(uint8_t i = 0; i < (rs->n); i++)
    inp[i] = rs->index[inp[i]];

  /* first form the syndromes */
  for(uint8_t i = 1; i <= (rs->n) - (rs->k); i++)
  {
    s[i] = 0;
    for(uint8_t j = 0; j < rs->n; j++)
      if(inp[j] != -1)
        s[i] ^= rs->alpha[(inp[j] + i * j) % rs->n]; /* inp[j] in index form */
                                                     /* convert syndrome from polynomial form to index form  */
    if(s[i] != 0)
      syn_error = 1; /* set flag if non-zero syndrome => error */
    s[i] = rs->index[s[i]];
  };

  if(syn_error) /* if errors, try and correct */
  {
    /* compute the error location polynomial via the Berlekamp iterative algorithm,
       following the terminology of Lin and Costello :   d[u] is the 'mu'th
       discrepancy, where u='mu'+1 and 'mu' (the Greek lers->ter!) is the step number
       ranging from -1 to 2*rs->t (see L&C),  l[u] is the
       degree of the elp at that step, and u_l[u] is the difference between the
       step number and the degree of the elp.
    */
    /* initialise table entries */
    d[0] = 0;      /* index form */
    d[1] = s[1];   /* index form */
    elp[0][0] = 0; /* index form */
    elp[1][0] = 1; /* polynomial form */
    for(uint8_t i = 1; i < (rs->n) - (rs->k); i++)
    {
      elp[0][i] = -1; /* index form */
      elp[1][i] = 0;  /* polynomial form */
    }
    l[0] = 0;
    l[1] = 0;
    u_lu[0] = -1;
    u_lu[1] = 0;
    uint8_t u = 0;

    do
    {
      u++;
      if(d[u] == -1)
      {
        l[u + 1] = l[u];
        for(uint8_t i = 0; i <= l[u]; i++)
        {
          elp[u + 1][i] = elp[u][i];
          elp[u][i] = rs->index[elp[u][i]];
        }
      }
      else
      /* search for words with greatest u_lu[q] for which d[q]!=0 */
      {
        uint8_t q = u - 1;
        while((d[q] == -1) && (q > 0))
          q--;
        /* have found first non-zero d[q]  */
        if(q > 0)
        {
          uint8_t j = q;
          do
          {
            j--;
            if ((d[j] != -1) && (u_lu[q] < u_lu[j]))
              q = j;
          }while (j > 0);
        };

        /* have now found q such that d[u]!=0 and u_lu[q] is maximum */
        /* store degree of new elp polynomial */
        if(l[u] > l[q] + u - q)
          l[u + 1] = l[u];
        else
          l[u + 1] = l[q] + u - q;

        /* form new elp(x) */
        for(uint8_t i = 0; i < (rs->n) - (rs->k); i++)
        {
          elp[u + 1][i] = 0;
        }
        for(uint8_t i = 0; i <= l[q]; i++)
        {
          if(elp[q][i] != -1)
            elp[u + 1][i + u - q] = rs->alpha[(d[u] + rs->n - d[q] + elp[q][i]) % rs->n];
        }
        for(uint8_t i = 0; i <= l[u]; i++)
        {
          elp[u + 1][i] ^= elp[u][i];
          elp[u][i] = rs->index[elp[u][i]]; /*convert old elp value to index*/
        }
      }
      u_lu[u + 1] = u - l[u + 1];

      /* form (u+1)th discrepancy */
      if(u < (rs->n) - (rs->k)) /* no discrepancy computed on last iteration */
      {
        if(s[u + 1] != -1)
          d[u + 1] = rs->alpha[s[u + 1]];
        else
          d[u + 1] = 0;
        for(uint8_t i = 1; i <= l[u + 1]; i++)
        {
          if((s[u + 1 - i] != -1) && (elp[u + 1][i] != 0))
            d[u + 1] ^= rs->alpha[(s[u + 1 - i] + rs->index[elp[u + 1][i]]) % rs->n];
        }
        d[u + 1] = rs->index[d[u + 1]]; /* put d[u+1] into index form */
      }
    }while((u < (rs->n) - (rs->k)) && (l[u + 1] <= rs->t));

    u++;
    if(l[u] <= rs->t) /* can correct error */
    {
      /* put elp into index form */
      for(uint8_t i = 0; i <= l[u]; i++)
        elp[u][i] = rs->index[elp[u][i]];

      /* find roots of the error location polynomial */
      for(uint8_t i = 1; i <= l[u]; i++)
      {
        reg[i] = elp[u][i];
      }
      count = 0;
      for(uint8_t i = 1; i <= rs->n; i++)
      {
        uint8_t q = 1;
        for(uint8_t j = 1; j <= l[u]; j++)
          if(reg[j] != -1)
          {
            reg[j] = (reg[j] + j) % rs->n;
            q ^= rs->alpha[reg[j]];
          };
        if(!q) /* store root and error location number indices */
        {
          root[count] = i;
          loc[count] = rs->n - i;
          count++;
        }
      }
      if(count == l[u]) /* no. roots = degree of elp hence <= rs->t errors */
      {
        /* form polynomial z(x) */
        for(uint8_t i = 1; i <= l[u]; i++) /* Z[0] = 1 always - do not need */
        {
          if((s[i] != -1) && (elp[u][i] != -1))
            z[i] = rs->alpha[s[i]] ^ rs->alpha[elp[u][i]];
          else if ((s[i] != -1) && (elp[u][i] == -1))
            z[i] = rs->alpha[s[i]];
          else if ((s[i] == -1) && (elp[u][i] != -1))
            z[i] = rs->alpha[elp[u][i]];
          else
            z[i] = 0;
          for(uint8_t j = 1; j < i; j++)
          {
            if ((s[j] != -1) && (elp[u][i - j] != -1))
              z[i] ^= rs->alpha[(elp[u][i - j] + s[j]) % rs->n];
          }
          z[i] = rs->index[z[i]]; /* put into index form */
        };

        /* evaluate errors at locations given by error location numbers loc[i] */
        for(uint8_t i = 0; i < rs->n; i++)
        {
          err[i] = 0;
          if (inp[i] != -1) /* convert inp[] to polynomial form */
            inp[i] = rs->alpha[inp[i]];
          else
            inp[i] = 0;
        }
        for(uint8_t i = 0; i < l[u]; i++) /* compute numerator of error term first */
        {
          err[loc[i]] = 1; /* accounts for z[0] */
          for(uint8_t j = 1; j <= l[u]; j++)
          {
            if (z[j] != -1)
              err[loc[i]] ^= rs->alpha[(z[j] + j * root[i]) % rs->n];
          }
          if(err[loc[i]] != 0)
          {
            err[loc[i]] = rs->index[err[loc[i]]];
            uint8_t q = 0; /* form denominator of error term */
            for(uint8_t j = 0; j < l[u]; j++)
            {
              if(j != i)
                q += rs->index[1 ^ rs->alpha[(loc[j] + root[i]) % rs->n]];
            }
            q = q % rs->n;
            err[loc[i]] = rs->alpha[(err[loc[i]] - q + rs->n) % rs->n];
            inp[loc[i]] ^= err[loc[i]]; /*inp[i] must be in polynomial form */
          }
        }
      }
      else                          /* no. roots != degree of elp => >rs->t errors and cars->not solve */
      {
        for(uint8_t i = 0; i < rs->n; i++) /* could return error flag if desired */
        {
          if(inp[i] != -1)         /* convert inp[] to polynomial form */
            inp[i] = rs->alpha[inp[i]];
          else
            inp[i] = 0; /* just output received codeword as is */
        }
      }
    }
    else                          /* elp has degree has degree >rs->t hence cars->not solve */
    {
      for(uint8_t i = 0; i < rs->n; i++) /* could return error flag if desired */
      {
        if(inp[i] != -1)         /* convert inp[] to polynomial form */
          inp[i] = rs->alpha[inp[i]];
        else
          inp[i] = 0; /* just output received codeword as is */
      }
    }
  }
  else /* no non-zero syndromes => no errors: output received codeword */
  {
    for(uint8_t i = 0; i < rs->n; i++)
    {
      if(inp[i] != -1) /* convert inp[] to polynomial form */
        inp[i] = rs->alpha[inp[i]];
      else
        inp[i] = 0;
    }
  }
}
