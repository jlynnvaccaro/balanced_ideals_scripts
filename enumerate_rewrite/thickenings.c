#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>

#include "thickenings.h"
#include "weyl.h"
#include "queue.h"

/*
  This function enumerates all balanced ideals satisfying certain constraints, given by its arguments pos, neg and next_neg

  - info: constant information which just gets passed on to recursive calls, mainly contains the principal ideals
  - pos: a set of elements which have to be positive (that is, in the ideal)
  - neg: a set of elements which have to be negative (not in the ideal)
  - next_neg: this element has to be the first negative one not already contained in neg; if next_neg is info.size/2, then everything not in neg has to be positive
  - already_known: not a constraint, but just a hint to speed things up; tells the function that the first already_known elements are set either in neg or in pos; must be less or equal to next_neg
  - returns number of balanced ideals found

  uses the bitvector functions bv_union, bv_copy, bv_set_range_except, bv_disjoint, bv_next_zero

 */

static long enumerate_tree(const enumeration_info_t *info, const bitvec_t *pos, const bitvec_t *neg, int next_neg, int already_known, int level)
{
  static long totcount = 0;
  bitvec_t newpos, newneg, known;
  int next_next_neg;
  long count = 0;

  // the omission of next_neg means inclusion of info->size - 1 - next_neg
  // add its principal ideal to pos and the opposite to neg
  if(next_neg != info->size/2) {
    bv_union(&info->principal_pos[info->size - 1 - next_neg], pos, &newpos);
    bv_union(&info->principal_neg[info->size - 1 - next_neg], neg, &newneg);
  } else { // or, if there is no next_neg, just copy
    bv_copy(pos, &newpos);
    bv_copy(neg, &newneg);
  }

  // everything before next_neg which was unknown should be set to positive; to speed this up, we can start with already_known
  bv_set_range_except(&newpos, neg, already_known, next_neg);

  #ifdef _DEBUG
  bv_print_nice(stderr, &newpos, &newneg, -1, info->size/2);
  fprintf(stderr, "\n");
  #endif

  // check if this leads to any conflicts (equivalently, violates slimness)
  if(!bv_disjoint(&newpos, &newneg))
    return 0;

  // what do we know so far?
  bv_union(&newpos, &newneg, &known);

  next_next_neg = bv_next_zero(&known, next_neg + 1);

  if(next_next_neg >= info->size/2) {
    // there is no unknown left, so we found a balanced ideal
    if(info->callback)
      info->callback(&newpos, info->size, info);
    return 1;
  }

  do {
    count += enumerate_tree(info, &newpos, &newneg, next_next_neg, next_neg + 1, level + 1);
    next_next_neg = bv_next_zero(&known, next_next_neg + 1);
  } while(next_next_neg <= info->size/2);

  return count;
}

static void generate_principal_ideals(doublequotient_t *dq, bitvec_t *pos, bitvec_t *neg, int *is_slim)
{
  queue_t queue;
  int current;
  doublecoset_list_t *edge;
  int size = dq->count;

  // generate principal ideals
  int *principal = (int*)malloc(size*sizeof(int));
  for(int i = 0; i < size; i++) {
    memset(principal, 0, size*sizeof(int));
    principal[i] = 1;
    queue_init(&queue);
    queue_put(&queue, i);
    while((current = queue_get(&queue)) != -1)
      for(edge = dq->cosets[current].bruhat_lower; edge; edge = edge->next)
	if(!principal[edge->to->index]) {
	  principal[edge->to->index] = 1;
	  queue_put(&queue, edge->to->index);
	}

    // copy the first half into bitvectors
    bv_clear(&pos[i]);
    bv_clear(&neg[i]);
    is_slim[i] = 1;
    for(int j = 0; j < size/2; j++)
      if(principal[j])
	bv_set_bit(&pos[i], j);
    for(int j = 0; j < size/2; j++)
      if(principal[size - 1 - j]) {
	bv_set_bit(&neg[i], j);
	if(bv_get_bit(&pos[i], j))
	  is_slim[i] = 0;
      }

#ifdef _DEBUG
    if(is_slim[i]) {
      fprintf(stderr, " ids: [0");
      for(int j = 1; j < size; j++)
	if(principal[j])
	  fprintf(stderr, ", %d", dq->cosets[j].min->id);
      fprintf(stderr, "]\n");
    }
#endif

  }
  free(principal);

  // output principal ideals
#ifdef _DEBUG
  for(int i = 0; i < size; i++) {
    fprintf(stderr, "%2d: ", i);
    bv_print_nice(stderr, &pos[i], &neg[i], -1, size/2);
    fprintf(stderr, "\n");
  }
  fprintf(stderr,"\n");
#endif

}

/*
   enumerates all balanced ideals

   - graph: hasse diagram of the bruhat order (of double cosets) with opposition pairing
   - size: number of nodes in graph
   - callback to call when a balanced ideal was found
   - arbitrary data for callback function

   returns the number of balanced ideals
*/

long enumerate_balanced_thickenings(doublequotient_t *dq, enumeration_callback callback, void *callback_data)
{
  long count = 0;
  enumeration_info_t info;

  info.size = dq->count;
  info.callback = callback;
  info.callback_data = callback_data;
  info.principal_pos = (bitvec_t*)malloc(info.size*sizeof(bitvec_t));
  info.principal_neg = (bitvec_t*)malloc(info.size*sizeof(bitvec_t));
  info.principal_is_slim = (int*)malloc(info.size*sizeof(int));

  // the algorithm only works if the opposition pairing does not stabilize any element
  // if this happens, there can be no balanced thickenings
  for(int i = 0; i < dq->count; i++)
    if(dq->cosets[i].opposite->min->id == dq->cosets[i].min->id)
      return 0;

  // we can only handle bitvectors up to BV_BLOCKSIZE*BV_RANK bits, but we only store half of the weyl group
  ERROR(info.size > 2*BV_BLOCKSIZE*BV_RANK, "We can handle at most %d cosets. Increase BV_RANK if more is needed.\n", 2*BV_BLOCKSIZE*BV_RANK);

  generate_principal_ideals(dq, info.principal_pos, info.principal_neg, info.principal_is_slim);

  // enumerate balanced ideals
  bitvec_t pos, neg;
  bv_clear(&pos);
  bv_clear(&neg);
  for(int i = 0; i <= info.size/2; i++)
    count += enumerate_tree(&info, &pos, &neg, i, 0, 0);

  free(info.principal_is_slim);
  free(info.principal_pos);
  free(info.principal_neg);

  return count;
}
