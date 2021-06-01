#include "weyl.h"
#include "queue.h"

#include <stdio.h>
#include <memory.h>
#include <stdlib.h>

#define BIT(n) ((uint64_t)1 << (n))

typedef struct {
  weylid_t id;
  int position;
} weylid_lookup_t;

static void generate_left_and_ids(semisimple_type_t type, weylgroup_element_t *group);
static int search(const void *key, const void *base, size_t nmem, size_t size, int (*compar) (const void *, const void *, void *), void *arg);
static int compare_root_vectors(int rank, const int *x, const int *y);
static int compare_root_vectors_qsort(const void *x, const void *y, void *arg);
static int compare_weylid(const void *x, const void *y);
static int compare_weylid_lookup(const void *x, const void *y);
static int lookup_id(weylid_t id, weylid_lookup_t *list, int len);
static weylid_t multiply_generator(int s, weylid_t w, const int *simple, const int *mapping, int rank, int positive);
static void reflect_root_vector(const int *cartan, int rank, int i, int *old, int *new);
static weylgroup_element_t* apply_word(int *word, int len, weylgroup_element_t *current);
static weylgroup_element_t* apply_word_reverse(int *word, int len, weylgroup_element_t *current);

/******** generate_left_and_ids and a pile of helper functions **************/

static void generate_left_and_ids(semisimple_type_t type, weylgroup_element_t *group)
{
  int rank = weyl_rank(type);
  int order = weyl_order(type);
  int positive = weyl_positive(type);

  queue_t queue;
  int current;
  int roots_known, elements, length_elements, nextids_count;
  int *cartan_matrix;
  int *root_vectors;
  int *vector;
  int *simple_roots;
  int *root_mapping;
  weylid_t *ids, *edges, *nextids;
  weylid_lookup_t *lookup;

  // allocate temporary stuff

  cartan_matrix =      (int*)malloc(rank*rank      *sizeof(int));
  root_vectors =       (int*)malloc(2*positive*rank*sizeof(int));
  vector =             (int*)malloc(rank           *sizeof(int));
  root_mapping =       (int*)malloc(positive*rank  *sizeof(int));
  simple_roots =       (int*)malloc(rank           *sizeof(int));
  ids =           (weylid_t*)malloc(order          *sizeof(weylid_t));
  edges =         (weylid_t*)malloc(rank*order     *sizeof(weylid_t));
  nextids =       (weylid_t*)malloc(rank*order     *sizeof(weylid_t));
  lookup = (weylid_lookup_t*)malloc(order          *sizeof(weylid_lookup_t));

  // get all information on the cartan type
  LOG("Get Cartan matrix.\n");

  weyl_cartan_matrix(type, cartan_matrix);

  // enumerate roots, first the simple ones, then all others by reflecting
  LOG("Enumerate roots.\n");

  memset(root_vectors, 0, 2*positive*rank*sizeof(int));
  roots_known = 0;

  queue_init(&queue);
  for(int i = 0; i < rank; i++) {
    root_vectors[rank*i + i] = 1; // (r_i)_j = delta_ij
    queue_put(&queue, i);
    roots_known++;
  }

  while((current = queue_get(&queue)) != -1) {
    for(int i = 0; i < rank; i++) {
      reflect_root_vector(cartan_matrix, rank, i, &root_vectors[rank*current], vector);
      int j;
      for(j = 0; j < roots_known; j++)
	if(compare_root_vectors(rank, &root_vectors[rank*j], vector) == 0)
	  break;
      if(j == roots_known) {
	memcpy(&root_vectors[rank*roots_known], vector, rank*sizeof(int));
	queue_put(&queue, roots_known);
	roots_known++;
      }
    }
  }

  ERROR(roots_known != 2*positive, "Number of roots does not match!\n");

  // sort roots and restrict to positives
  LOG("Sort roots.\n");

  qsort_r(root_vectors, 2*positive, rank*sizeof(int), compare_root_vectors_qsort, &rank);
  memcpy(root_vectors, &root_vectors[positive*rank], positive*rank*sizeof(int)); // this just copies the second part of the list onto the first; source and destination are disjoint!

  // generate root_mapping, which gives the action of root reflections on positive roots (-1 if result is not a positive root)
  LOG("Compute root reflections.\n");

  for(int i = 0; i < positive; i++) {
    for(int j = 0; j < rank; j++) {
      reflect_root_vector(cartan_matrix, rank, j, &root_vectors[rank*i], vector);
      root_mapping[i*rank+j] =
	search(vector, root_vectors, positive, rank*sizeof(int), compare_root_vectors_qsort, &rank);
    }
  }

  // find simple roots in the list
  LOG("Find simple roots.\n");

  for(int i = 0; i < rank; i++) {
    memset(vector, 0, rank*sizeof(int));
    vector[i] = 1;
    simple_roots[i] = search(vector, root_vectors, positive, rank*sizeof(int), compare_root_vectors_qsort, &rank);
  }

  // enumerate weyl group elements using difference sets
  LOG("Enumerate Weyl group elements.\n");

  nextids[0] = 0;
  nextids_count = 1;
  elements = 0;
  for(int len = 0; len <= positive; len++) {
    length_elements = 0;

    // find unique ids in edges added in the last iteration
    qsort(nextids, nextids_count, sizeof(weylid_t), compare_weylid);
    for(int i = 0; i < nextids_count; i++)
      if(i == 0 || nextids[i] != nextids[i-1])
	ids[elements + length_elements++] = nextids[i];

    // add new edges
    nextids_count = 0;
    for(int i = elements; i < elements + length_elements; i++)
      for(int j = 0; j < rank; j++) {
	edges[i*rank+j] = multiply_generator(j, ids[i], simple_roots, root_mapping, rank, positive);
	if(!(ids[i] & BIT(simple_roots[j]))) // the new element is longer then the old one
	  nextids[nextids_count++] = edges[i*rank+j];
      }

    elements += length_elements;
  }

  // translate the ids to list positions (i.e. local continuous ids)
  LOG("Reorder Weyl group elements.\n");

  for(int i = 0; i < order; i++) {
    lookup[i].id = ids[i];
    lookup[i].position = i;
  }
  qsort(lookup, order, sizeof(weylid_lookup_t), compare_weylid_lookup);

  // fill in results
  LOG("Compute left multiplication.\n");

  for(int i = 0; i < order; i++) {
    group[i].id = ids[i];
    for(int j = 0; j < rank; j++)
      group[i].left[j] = group + lookup_id(edges[i*rank+j], lookup, order);
  }

  // free temporary stuff

  free(cartan_matrix);
  free(root_vectors);
  free(vector);
  free(root_mapping);
  free(simple_roots);
  free(ids);
  free(edges);
  free(nextids);
  free(lookup);
}

// glibc search function, but with user pointer and returning index (or -1 if not found)
static int search(const void *key, const void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *, void *), void *arg)
{
  size_t l, u, idx;
  const void *p;
  int comparison;

  l = 0;
  u = nmemb;
  while (l < u) {
    idx = (l + u) / 2;
    p = (void *) (((const char *) base) + (idx * size));
    comparison = (*compar) (key, p, arg);
    if (comparison < 0)
      u = idx;
    else if (comparison > 0)
      l = idx + 1;
    else
      return idx;
  }

  return -1;
}

// maybe we want a different ordering here?
static int compare_root_vectors(int rank, const int *x, const int *y)
{
  for(int i = 0; i < rank; i++)
    if(x[i] != y[i])
      return x[i] - y[i];

  return 0;
}

static int compare_root_vectors_qsort(const void *x, const void *y, void *arg)
{
  return compare_root_vectors(*((int*)arg), x, y);
}

static int compare_weylid(const void *x, const void *y)
{
  weylid_t u = *((weylid_t*)x);
  weylid_t v = *((weylid_t*)y);

  return u > v ? 1 : u < v ? -1 : 0;
}

static int compare_weylid_lookup(const void *x, const void *y)
{
  weylid_t u = ((weylid_lookup_t*)x)->id;
  weylid_t v = ((weylid_lookup_t*)y)->id;

  return u > v ? 1 : u < v ? -1 : 0;
}

static int lookup_id(weylid_t id, weylid_lookup_t *list, int len)
{
  weylid_lookup_t key;
  key.id = id;
  weylid_lookup_t *p = (weylid_lookup_t*)bsearch(&key, list, len, sizeof(weylid_lookup_t), compare_weylid_lookup);
  return p->position;
}

static weylid_t multiply_generator(int s, weylid_t w, const int* simple, const int* mapping, int rank, int positive)
{
  weylid_t sw = 0;

  for(int i = 0; i < positive; i++) {
    if(w & BIT(i))
      if(mapping[i*rank+s] != -1)
	sw |= BIT(mapping[i*rank+s]);
  }

  if(w & BIT(simple[s]))
    return sw;
  else
    return sw | BIT(simple[s]);
}

static void reflect_root_vector(const int *cartan, int rank, int i, int *old, int *new)
{
  memcpy(new, old, rank*sizeof(int));
  for(int j = 0; j < rank; j++)
    new[i] -= cartan[i*rank + j]*old[j];
}

/************* Weyl group infos ************************/

static int weyl_exists(simple_type_t type)
{
  if(type.series < 'A' || type.series > 'G' || type.rank < 1 ||
     type.series == 'B' && type.rank < 2 ||
     type.series == 'C' && type.rank < 2 ||
     type.series == 'D' && type.rank < 3 ||
     type.series == 'E' && type.rank != 6 && type.rank != 7 && type.rank != 8 ||
     type.series == 'F' && type.rank != 4 ||
     type.series == 'G' && type.rank != 2)
    return 0;
  else
    return 1;
}

int weyl_rank(semisimple_type_t type)
{
  // Total rank of the Weyl group, which is the sum of the factor ranks.
  int rank = 0;
  for(int i = 0; i < type.n; i++)
    rank += type.factors[i].rank;
  return rank;
}

int weyl_order(semisimple_type_t type)
{
  // Total order of the Weyl group, which is the product of the orders of the factors. E.g. order(An) = (n+1)!, order(Bn) = (n+1)! * 2^(n+1)
  int order = 1;
  for(int i = 0; i < type.n; i++) {
    ERROR(!weyl_exists(type.factors[i]), "A Weyl group of type %c%d does not exist!\n", type.factors[i].series, type.factors[i].rank);

    switch(type.factors[i].series) {
    case 'A':
      for(int j = 1; j <= type.factors[i].rank + 1; j++)
	order *= j;
      break;

    case 'B': case 'C':
      for(int j = 1; j <= type.factors[i].rank; j++)
	order *= 2*j;
      break;

    case 'D':
      for(int j = 2; j <= type.factors[i].rank; j++)
	order *= 2*j;
      break;

    case 'E':
      if(type.factors[i].rank == 6)
	order *= 51840;
      else if(type.factors[i].rank == 7)
	order *= 2903040;
      else if(type.factors[i].rank == 8)
	order *= 696729600;
      break;

    case 'F':
      order *= 1152;
      break;

    case 'G':
      order *= 12;
      break;
    }
  }

  return order;
}

int weyl_positive(semisimple_type_t type)
{
  // Maximum length, which is the sum of the maximum length of the summands. This is the length for w_0 in W.
  int positive = 0;

  for(int i = 0; i < type.n; i++) {
    ERROR(!weyl_exists(type.factors[i]), "A Weyl group of type %c%d does not exist!\n", type.factors[i].series, type.factors[i].rank);

    switch(type.factors[i].series) {
    case 'A':
      positive += (type.factors[i].rank * (type.factors[i].rank + 1)) / 2;
      break;

    case 'B': case 'C':
      positive += type.factors[i].rank * type.factors[i].rank;
      break;

    case 'D':
      positive += type.factors[i].rank * (type.factors[i].rank - 1);
      break;

    case 'E':
      if(type.factors[i].rank == 6)
	positive += 36;
      else if(type.factors[i].rank == 7)
	positive += 63;
      else if(type.factors[i].rank == 8)
	positive += 120;
      break;

    case 'F':
      positive += 24;
      break;

    case 'G':
      positive += 6;
      break;
    }
  }

  return positive;
}

int weyl_opposition(semisimple_type_t type, int simple_root)
{
  int offset = 0;
  int factor = 0;
  int r, iota_r;

  for(factor = 0; factor < type.n; factor++)
    if(simple_root < offset + type.factors[factor].rank)
      break;
    else
      offset += type.factors[factor].rank;
  r = simple_root - offset;

  ERROR(!weyl_exists(type.factors[factor]), "A Weyl group of type %c%d does not exist!\n", type.factors[factor].series, type.factors[factor].rank);

  switch(type.factors[factor].series) {
  case 'B': case 'C': case 'F': case 'G':
    iota_r = r;
    break;

  case 'A':
    iota_r = type.factors[factor].rank - 1 - r;
    break;

  case 'D':
    if(type.factors[factor].rank % 2 == 0)
      iota_r = r;
    else
      iota_r = r == 0 ? 1 : r == 1 ? 0 : r;
    break;

  case 'E':
    if(type.factors[factor].rank != 6)
      iota_r = r;
    else
      iota_r = r == 2 || r == 3 ? r : 5 - r;
    break;
  }

  return iota_r + offset;
}

void weyl_cartan_matrix(semisimple_type_t type, int *m)
{
  int offset = 0;
  int rank = weyl_rank(type);

  int **A = (int**)malloc(rank*sizeof(int*));

  memset(m, 0, rank*rank*sizeof(int));
  for(int i = 0; i < rank; i++)
    m[i*rank+i] = 2;

  for(int k = 0; k < type.n; k++) {
    ERROR(!weyl_exists(type.factors[k]), "A Weyl group of type %c%d does not exist!\n", type.factors[k].series, type.factors[k].rank);

    for(int i = 0; i < type.factors[k].rank; i++)  // A is the submatrix corresponding to the current simple factor
      A[i] = &m[(i+offset)*rank + offset];

    for(int i = 1; i < type.factors[k].rank; i++) {
      A[i][i-1] = -1;
      A[i-1][i] = -1;
    }

    switch(type.factors[k].series) {
    case 'A':
      break;
    case 'B':
      A[0][1] = -2;
      break;
    case 'C':
      A[1][0] = -2;
      break;
    case 'D':
      A[0][1] = A[1][0] = 0;
      A[0][2] = A[2][0] = -1;
      break;
    case 'E':
      A[1][2] = A[2][1] = 0;
      A[1][3] = A[3][1] = -1;
      break;
    case 'F':
      A[2][1] = -2;
      break;
    case 'G':
      A[1][0] = -3;
      break;
    }

    offset += type.factors[k].rank;
  }

  free(A);
}

/************ weyl_generate etc. ********************/

static weylgroup_element_t* apply_word(int *word, int len, weylgroup_element_t *current)
{
  for(int k = len - 1; k >= 0; k--) // apply group element from right to left
    current = current->left[word[k]];

  return current;
}

static weylgroup_element_t* apply_word_reverse(int *word, int len, weylgroup_element_t *current)
{
  for(int k = 0; k < len; k++) // apply group element from left to right (i.e. apply inverse)
    current = current->left[word[k]];

  return current;
}

weylgroup_t *weyl_generate(semisimple_type_t type)
{
  int rank = weyl_rank(type);
  int order = weyl_order(type);
  int positive = weyl_positive(type);

  ERROR(positive > 64, "We can't handle root systems with more than 64 positive roots!\n");

  // allocate result

  weylgroup_element_t *group = (weylgroup_element_t*)malloc(order*sizeof(weylgroup_element_t));
  weylgroup_t *result = malloc(sizeof(weylgroup_t));
  result->type = type;
  result->elements = group;
  result->lists = (weylgroup_element_t**)malloc(2*order*rank*sizeof(weylgroup_element_t*));

  for(int i = 0; i < order; i++) {
    group[i].left = result->lists + 2*i*rank;
    group[i].right = result->lists + (2*i+1)*rank;
    group[i].coset = (doublecoset_t*)0;
    group[i].index = i;
  }

  // the main part
  LOG("Start generating Weyl group.\n");

  generate_left_and_ids(type, group);

  // word length is just the number of 1s in the binary id
  LOG("Find word lengths.\n");

  for(int i = 0; i < order; i++) {
    group[i].wordlength = 0;
    for(int j = 0; j < positive; j++)
      if(group[i].id & BIT(j))
	group[i].wordlength++;
  }

  // allocate letters

  int total_wordlength = 0;
  for(int i = 0; i < order; i++)
    total_wordlength += group[i].wordlength;
  result->letters = (int*)malloc(total_wordlength*sizeof(int));
  total_wordlength = 0;
  for(int i = 0; i < order; i++) {
    group[i].word = result->letters + total_wordlength;
    total_wordlength += group[i].wordlength;
  }

  // find shortest words (using that the elements are already ordered by word length)
  LOG("Find shortest words.\n");

  memset(result->letters, -1, total_wordlength*sizeof(int));
  for(int i = 0; i < order - 1; i++) {
    weylgroup_element_t *this = &group[i];
    for(int j = 0; j < rank; j++) {
      weylgroup_element_t *that = group[i].left[j];
      if(that->wordlength > this->wordlength && that->word[0] == -1) {
	memcpy(that->word + 1, this->word, this->wordlength*sizeof(int));
	that->word[0] = j;
      }
    }
  }

  // generate right edges
  LOG("Compute right multiplication.\n");

  for(int i = 0; i < order; i++)
    for(int j = 0; j < rank; j++)
      group[i].right[j] = apply_word(group[i].word, group[i].wordlength, group[0].left[j]);

  // find opposites
  LOG("Find opposites.\n");

  weylgroup_element_t *longest = &group[order-1];
  for(int i = 0; i < order; i++)
    group[i].opposite = apply_word(longest->word, longest->wordlength, &group[i]);

  // check for root reflections
  LOG("Find root reflections.\n");

  for(int i = 0; i < order; i++)
    group[i].is_root_reflection = 0;
  for(int i = 0; i < order; i++)
    for(int j = 0; j < rank; j++) // we want to calculate word^{-1} * j * word; this is a root reflection
      apply_word_reverse(group[i].word, group[i].wordlength, group[i].left[j]) -> is_root_reflection = 1; // TODO: What does this code do?

  return result;
}

void weyl_destroy(weylgroup_t *group)
{
  free(group->elements);
  free(group->lists);
  free(group->letters);
  free(group);
}

doublequotient_t *weyl_generate_bruhat(semisimple_type_t type, int left_invariance, int right_invariance)
{
  int rank = weyl_rank(type);
  int order = weyl_order(type);
  int positive = weyl_positive(type);
  int count;

  int is_minimum, is_maximum;

  weylgroup_t *wgroup = weyl_generate(type);
  weylgroup_element_t *group = wgroup->elements;
  doublecoset_t *cosets;

  for(int i = 0; i < rank; i++) {
    int oppi = weyl_opposition(type, i);
    if(left_invariance & BIT(i) && !(left_invariance & BIT(oppi)) ||
       left_invariance & BIT(oppi) && !(left_invariance & BIT(i)))
      ERROR(1, "The specified left invariance is not invariant under the opposition involution!\n");
  }

  doublequotient_t *result = (doublequotient_t*)malloc(sizeof(doublequotient_t));
  result->type = type;
  result->left_invariance = left_invariance;
  result->right_invariance = right_invariance;
  result->group = wgroup->elements;
  result->grouplists = wgroup->lists;
  result->groupletters = wgroup->letters;

  free(wgroup); // dissolved in result and not needed anymore

  LOG("Count cosets.\n"); // count cosets by finding the minimum length element in every coset

  count = 0;
  for(int i = 0; i < order; i++) {
    is_minimum = 1;
    for(int j = 0; j < rank; j++)
      if(left_invariance  & BIT(j) && group[i].left[j]->wordlength  < group[i].wordlength ||
	 right_invariance & BIT(j) && group[i].right[j]->wordlength < group[i].wordlength)
	is_minimum = 0;
    if(is_minimum)
      count++;
  }
  result->count = count;

  // alloc more stuff

  cosets = result->cosets = (doublecoset_t*)malloc(count*sizeof(doublecoset_t));
  for(int i = 0; i < count; i++) {
    cosets[i].bruhat_lower = cosets[i].bruhat_higher = (doublecoset_list_t*)0;
  }
  result->lists = (doublecoset_list_t*)malloc(2*count*positive*sizeof(doublecoset_list_t)); // 2 times, for bruhat lower and higher

  LOG("Find minimal length elements in cosets.\n"); // basically same code as above

  count = 0;
  for(int i = 0; i < order; i++) {
    is_minimum = 1;
    for(int j = 0; j < rank; j++)
      if(left_invariance  & BIT(j) && group[i].left[j]->wordlength  < group[i].wordlength ||
	 right_invariance & BIT(j) && group[i].right[j]->wordlength < group[i].wordlength)
	is_minimum = 0;
    if(is_minimum) {
      cosets[count].min = &group[i];
      group[i].coset = &cosets[count];
      count++;
    }
  }

  LOG("Generate quotient map.\n");

  for(int i = 0; i < order; i++) {
    for(int j = 0; j < rank; j++) {
      if(left_invariance & BIT(j) && group[i].left[j]->wordlength > group[i].wordlength)
	group[i].left[j]->coset = group[i].coset;
      if(right_invariance & BIT(j) && group[i].right[j]->wordlength > group[i].wordlength)
	group[i].right[j]->coset = group[i].coset;
    }
  }

  LOG("Find maximal length elements.\n");

  for(int i = 0; i < order; i++) {
    is_maximum = 1;
    for(int j = 0; j < rank; j++)
      if(left_invariance  & BIT(j) && group[i].left[j]->wordlength  > group[i].wordlength ||
	 right_invariance & BIT(j) && group[i].right[j]->wordlength > group[i].wordlength)
	is_maximum = 0;
    if(is_maximum) {
      group[i].coset->max = &group[i];
    }
  }

  LOG("Find opposites.\n");

  for(int i = 0; i < count; i++)
    cosets[i].opposite = cosets[i].min->opposite->coset;

  LOG("Sort opposites.\n");

  int *old2newindices = (int*)malloc(count*sizeof(int));
  int *new2oldindices = (int*)malloc(count*sizeof(int));

  // give the cosets some temporary indices
  for(int i = 0; i < count; i++)
    cosets[i].index = i;

  // generate a nice ordering, where element j is opposite to n-j, except the self-opposite ones, which are in the middle
  int j = 0;
  for(int i = 0; i < count; i++)
    if(i < cosets[i].opposite->index) {
      old2newindices[i] = j;
      old2newindices[cosets[i].opposite->index] = count-1-j;
      j++;
    }
  for(int i = 0; i < count; i++)
    if(i == cosets[i].opposite->index)
      old2newindices[i] = j++;

  for(int i = 0; i < count; i++)
    new2oldindices[old2newindices[i]] = i;

  // rewrite everything in the new ordering
  doublecoset_t *oldcosets = (doublecoset_t*)malloc(count*sizeof(doublecoset_t));
  memcpy(oldcosets, cosets, count*sizeof(doublecoset_t));
  for(int i = 0; i < count; i++) {
    cosets[i].min = oldcosets[new2oldindices[i]].min;
    cosets[i].max = oldcosets[new2oldindices[i]].max;
    cosets[i].opposite = cosets + old2newindices[oldcosets[new2oldindices[i]].opposite->index];
    //    cosets[i].bruhat_lower = oldcosets[new2oldindices[i]].bruhat_lower;
    //    cosets[i].bruhat_higher = oldcosets[new2oldindices[i]].bruhat_higher;
    //    for(doublecoset_list_t *current = cosets[i].bruhat_lower; current; current = current -> next)
    //      current->to = &cosets[old2newindices[current->to->index]];
    //    for(doublecoset_list_t *current = cosets[i].bruhat_higher; current; current = current -> next)
    //      current->to = &cosets[old2newindices[current->to->index]];
  }
  for(int i = 0; i < order; i++)
    group[i].coset = old2newindices[group[i].coset->index] + cosets;
  for(int i = 0; i < count; i++) // do this in the end, so we can use the "index" attribute before to translate pointers to indices
    cosets[i].index = i;

  free(old2newindices);
  free(new2oldindices);
  free(oldcosets);

  LOG("Find bruhat order.\n");

  int edgecount = 0;
  for(int i = 0; i < order; i++) {
    if(group[i].is_root_reflection) {
      for(int j = 0; j < count; j++) {
	weylgroup_element_t *this = cosets[j].min;
	weylgroup_element_t *that = apply_word(group[i].word, group[i].wordlength, cosets[j].min);
	if(this->wordlength > that->wordlength) { // this is higher in bruhat order than that
	  doublecoset_list_t *new = &result->lists[edgecount++];
	  new->next = this->coset->bruhat_lower;
	  this->coset->bruhat_lower = new;
	  new->to = that->coset;
	}
      }
    }
  }

  LOG("Perform transitive reduction.\n"); // eliminating redudant order relations

  doublecoset_t *origin;
  doublecoset_list_t *current;
  doublecoset_list_t *prev;
  queue_t queue;
  int idx;
  int *seen = malloc(count*sizeof(int));
  for(int i = 0; i < count; i++) {
    memset(seen, 0, count*sizeof(int));
    queue_init(&queue);

    for(int len = 1; len <= cosets[i].min->wordlength; len++) {

      // remove all edges originating from i of length len which connect to something already seen using shorter edges
      origin = &cosets[i];
      prev = (doublecoset_list_t*)0;

      for(current = origin->bruhat_lower; current; current = current->next) {
	if(origin->min->wordlength - current->to->min->wordlength != len) {
	  prev = current;
	} else if(seen[current->to->index]) {
	  if(prev)
	    prev->next = current->next;
	  else
	    origin->bruhat_lower = current->next;
	} else {
	  prev = current;
	  seen[current->to->index] = 1;
	  queue_put(&queue, current->to->index);
	}
      }

      // see which nodes we can reach using only edges up to length len, mark them as seen
      while((idx = queue_get(&queue)) != -1) {
	current = cosets[idx].bruhat_lower;
	for(current = cosets[idx].bruhat_lower; current; current = current->next) {
	  if(!seen[current->to->index]) {
	    seen[current->to->index] = 1;
	    queue_put(&queue, current->to->index);
	  }
	}
      }
    }
  }

  free(seen);

  LOG("Revert bruhat order.\n");

  for(int i = 0; i < count; i++) {
    for(current = cosets[i].bruhat_lower; current; current = current->next) {
      doublecoset_list_t *new = &result->lists[edgecount++];
      new->to = &cosets[i];
      new->next = current->to->bruhat_higher;
      current->to->bruhat_higher = new;
    }
  }

  return result;
}

void weyl_destroy_bruhat(doublequotient_t *dq)
{
  free(dq->group);
  free(dq->grouplists);
  free(dq->groupletters);
  free(dq->cosets);
  free(dq->lists);
  free(dq);
}
