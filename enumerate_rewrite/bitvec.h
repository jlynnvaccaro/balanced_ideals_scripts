/*********************************************************************
 * Filename:      bitvec.h
 *
 * Description:   Bit vectors implemented as uint64_t arrays, size
 *                fixed at compile time (in weylconfig.h).  Supports
 *                efficient set operations: Union, difference, count.
 *                Uses SSE4 64-bit popcount instruction.
 *
 * Author:        David Dumas <david@dumas.io>
 *
 * This program is free software distributed under the MIT license.
 * See the file LICENSE for details.
 ********************************************************************/

#ifndef __BITVEC_H__
#define __BITVEC_H__ 1

#include <string.h>

#include <inttypes.h>

#include <stdio.h>
#include <stdlib.h>

// FIRSTBITS(n) only yields useful result when 0 <= n < 64
#define BV_RANK 10
#define BV_BLOCKSIZE 64

#define FIRSTBITS(n) (((uint64_t)1 << (n)) - 1l)
#define BIT(n) (((uint64_t)1 << (n)))
#define ALLBITS ((uint64_t)-1)
#define BLOCK(n) ((n)/64)
#define INDEX(n) ((n)%64)

typedef struct {
  uint64_t v[BV_RANK];
} bitvec_t;

static inline void bv_clear_bit(bitvec_t *x, int k)
{
  x->v[BLOCK(k)] &= ~BIT(INDEX(k));
}

static inline void bv_set_bit(bitvec_t *x, int k)
{
  x->v[BLOCK(k)] |= BIT(INDEX(k));
}

static inline int bv_get_bit(const bitvec_t *x, int k)
{
  return (x->v[BLOCK(k)] >> INDEX(k)) & 0x1;
}

static inline void bv_clear(bitvec_t *x)
{
  int i;
  for (i=0;i<BV_RANK;i++)
    x->v[i] = 0;
}

static inline int bv_is_zero(const bitvec_t *x)
{
  int i;
  for (i=0;i<BV_RANK;i++)
    if (x->v[i])
      return 0;

  return 1;
}

static inline void bv_print(FILE *f, const bitvec_t *x, int len)
{
  for(int i = 0; i < len; i++) {
    fputc(bv_get_bit(x, i) ? '1' : '0', f);
    //    if(i % BLOCKSIZE == BLOCKSIZE - 1)
    //    fputc('-',f);
  }
}

static inline void bv_print_nice(FILE *f, const bitvec_t *pos, const bitvec_t *neg, int special, int len)
{
  for(int i = 0; i < len; i++) {
    if(i == special)
      fputc('X', f);
    else if(bv_get_bit(pos, i) && !bv_get_bit(neg, i))
      fputc('1', f);
    else if(!bv_get_bit(pos, i) && bv_get_bit(neg, i))
      fputc('0', f);
    else if(!bv_get_bit(pos, i) && !bv_get_bit(neg, i))
      fputc(' ', f);
    else
      fputc('-', f);
  }
}

static inline void bv_union(const bitvec_t *x, const bitvec_t *y, bitvec_t *result)
{
  int i;
  for (i=0; i < BV_RANK; i++) {
    result->v[i] = x->v[i] | y->v[i];
  }
}

static inline void bv_intersection(const bitvec_t *x, const bitvec_t *y, bitvec_t *result)
{
  int i;
  for (i=0; i < BV_RANK; i++) {
    result->v[i] = x->v[i] & y->v[i];
  }
}

static inline void bv_difference(const bitvec_t *x, const bitvec_t *y, bitvec_t *result)
{
  int i;
  for (i=0; i < BV_RANK; i++) {
    result->v[i] = x->v[i] & ~y->v[i];
  }
}

static inline int bv_disjoint(const bitvec_t *x, const bitvec_t *y)
{
  for(int i = 0; i < BV_RANK; i++)
    if(x->v[i] & y->v[i])
      return 0;

  return 1;
}

static inline int bv_full(const bitvec_t *x, int len)
{
  int i;
  for(i = 0; i < BLOCK(len); i++)
    if(x->v[i] != ALLBITS)
      return 0;

  return (x->v[i] & FIRSTBITS(INDEX(len))) == FIRSTBITS(INDEX(len));
}

// set bits in range start...end (including start and excluding end)
static inline void bv_set_range(bitvec_t *x, int start, int end)
{
  if(BLOCK(start) == BLOCK(end))
    x->v[BLOCK(start)] |= ~FIRSTBITS(INDEX(start)) & FIRSTBITS(INDEX(end));
  else {
    x->v[BLOCK(start)] |= ~FIRSTBITS(INDEX(start));
    for(int i = BLOCK(start) + 1; i < BLOCK(end); i++)
      x->v[i] = ALLBITS;
    x->v[BLOCK(end)] |= FIRSTBITS(INDEX(end));
  }
}

// set bits in range start...end (including start and excluding end), except if they are set in mask
static inline void bv_set_range_except(bitvec_t *x, const bitvec_t *mask, int start, int end)
{
  if(BLOCK(start) == BLOCK(end))
    x->v[BLOCK(start)] |= ~FIRSTBITS(INDEX(start)) & FIRSTBITS(INDEX(end)) & ~mask->v[BLOCK(start)];
  else {
    x->v[BLOCK(start)] |= ~FIRSTBITS(INDEX(start)) & ~mask->v[BLOCK(start)];
    for(int i = BLOCK(start) + 1; i < BLOCK(end); i++)
      x->v[i] |= ~mask->v[i];
    x->v[BLOCK(end)] |= FIRSTBITS(INDEX(end)) & ~mask->v[BLOCK(end)];
  }
}

// find least significant 0 bit starting from position start (included)
static inline int bv_next_zero(const bitvec_t *x, int start)
{
  int position;

  position = ffsll(~(x->v[BLOCK(start)] | FIRSTBITS(INDEX(start))));

  if(position)
    return BLOCK(start)*BV_BLOCKSIZE + position - 1; // found zero in same chunk

  for(int i = BLOCK(start) + 1; i < BV_RANK; i++) {
    position = ffsll(~x->v[i]);
    if(position) // found a 0
      return i*BV_BLOCKSIZE + position - 1;
  }

  return BV_RANK*BV_BLOCKSIZE; // found nothing
}

static inline void bv_copy(const bitvec_t *from, bitvec_t *to)
{
  for(int i = 0; i < BV_RANK; i++)
    to->v[i] = from->v[i];
}

static inline void bv_negate(const bitvec_t *from, bitvec_t *to)
{
  for(int i = 0; i < BV_RANK; i++)
    to->v[i] = ~from->v[i];
}

static inline int bv_count_bits(const bitvec_t *vec,int len)
{
  int count = 0;
  for(int i=0; i<len;i++)
  {
    count += bv_get_bit((vec), i);
  }
  return count;
}

#endif /* __BITVEC_H__ */
