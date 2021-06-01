#include "thickenings.h"
#include "weyl.h"
#include "queue.h"

#include <strings.h>
#include <stdio.h>

char stringbuffer[100];
char stringbuffer2[100];

typedef struct {
  doublequotient_t *dq;
  int rank;
  int order;
  int positive;
  int *buffer;
  int level;
} info_t;

static char* alphabetize(weylgroup_element_t *e, char *str)
{
  if(e->wordlength == 0)
    sprintf(str, "1");
  else {
    for(int j = 0; j < e->wordlength; j++)
      str[j] = e->word[j] + 'a';
    str[e->wordlength] = 0;
  }

  return str;
}

void balanced_thickening_callback(const bitvec_t *pos, int size, const enumeration_info_t *ei)
{
  static long totcount = 0;

  if(ei->callback_data) {
    info_t *info = (info_t*)ei->callback_data;

    unsigned long right_invariance = FIRSTBITS(info->rank);
    unsigned long left_invariance = FIRSTBITS(info->rank);

    int bit1, bit2left, bit2right, left, right;

    for(int i = 0; i < size; i++) {
      bit1 = i < size/2 ? bv_get_bit(pos, i) : !bv_get_bit(pos, size - 1 - i);
      for(int j = 0; j < info->rank; j++) {
	left = info->dq->cosets[i].min->left[j]->coset->index;
	right = info->dq->cosets[i].min->right[j]->coset->index;
	bit2left = left < size/2 ? bv_get_bit(pos, left) : !bv_get_bit(pos, size - 1 - left);
	bit2right = right < size/2 ? bv_get_bit(pos, right) : !bv_get_bit(pos, size - 1 - right);
	if(bit1 != bit2left)
	  left_invariance &= ~BIT(j);
	if(bit1 != bit2right)
	  right_invariance &= ~BIT(j);
      }
    }

    printf("%4ld left: ", totcount++);
    for(int j = 0; j < info->rank; j++)
      printf("%c", left_invariance & (1 << j) ? j + 'a' : ' ');
    printf(" right: ");
    for(int j = 0; j < info->rank; j++)
      printf("%c", right_invariance & (1 << j) ? j + 'a' : ' ');

    if(info->buffer) {
      bitvec_t low, high;
      bv_copy(pos, &low);
      bv_negate(pos, &high);

      printf(" gen: ");


      for(int i = 0; i < size/2; i++) {
	if(!bv_get_bit(&high, i))
	  continue;

	printf("%s ", alphabetize(info->dq->cosets[size-1-i].min, stringbuffer));

	bv_difference(&high, &ei->principal_neg[size-1-i], &high);
	bv_difference(&low,  &ei->principal_pos[size-1-i], &low);
      }

      for(int i = size/2 - 1; i >= 0; i--) {
	if(!bv_get_bit(&low, i))
	  continue;

	printf("%s ", alphabetize(info->dq->cosets[i].min, stringbuffer));

	bv_difference(&low, &ei->principal_pos[i], &low);
      }
    }

    int max_length = 0;
    for(int i = 0; i < size/2; i++) {
      if(bv_get_bit(pos, i)) {
	if(info->dq->cosets[i].max->wordlength > max_length)
	  max_length = info->dq->cosets[i].max->wordlength;
      } else {
	if(info->dq->cosets[size-i-1].max->wordlength > max_length)
	  max_length = info->dq->cosets[size-i-1].max->wordlength;
      }
    }

    printf("\n");
  }
}

void balanced_thickening_simple_callback(const bitvec_t *pos, int size, const enumeration_info_t *ei)
{
  long *count = (long*)ei->callback_data;

  if((++(*count)) % 100000000 == 0) {
    bv_print(stderr, pos, size/2);
    fprintf(stderr, "\n");
  }
}

int main(int argc, const char *argv[])
{
  semisimple_type_t type;
  unsigned long right_invariance, left_invariance;
  int rank, order, positive;
  int fixpoints;

  doublequotient_t *dq;

  const char *alphabet = "abcdefghijklmnopqrstuvwxyz";

  // read arguments

  ERROR(argc < 2, "Too few arguments!\n\nUsage is '%s A2 A3' or '%s A2 A3 -abc -abc' with\nA2,A3 simple Weyl factors and abc,abc left/right invariance.\n\nTo adjust output detail, set environment variable OUTPUT_LEVEL (1-4).\n",argv[0],argv[0]);

  // Count the number of simple factors in the semisimple Weyl group
  type.n = 0;
  for(int i = 0; i < argc - 1; i++) {
    // Skip any arguments that don't start with a letter A-G.
    if(argv[i+1][0] < 'A' || argv[i+1][0] > 'G')
      break;
    type.n++;
  }

  // Allocate memory, then read in the actual simple factors by letter/number, e.g. A5 is series 'A' and rank '5'. Series is A-G and the max rank is 9.
  type.factors = (simple_type_t*)malloc(type.n*sizeof(simple_type_t));
  for(int i = 0; i < type.n; i++) {
    type.factors[i].series = argv[i+1][0];
    type.factors[i].rank = argv[i+1][1] - '0';
    ERROR(argv[i+1][0] < 'A' || argv[i+1][0] > 'G' || argv[i+1][1] < '1' || argv[i+1][1] > '9', "Arguments must be Xn with X out of A-G and n out of 1-9\n");
  }

  left_invariance = right_invariance = 0;
  // Additional command line arguments that were not factors are the left/right invariance.
  if(argc - type.n >= 3) {
    if(strcmp(argv[type.n + 1], "-") != 0){
      printf("%s\n",argv[type.n+1]);
      for(int i = 0; i < strlen(argv[type.n + 1]); i++)
	      left_invariance |= (1 << (argv[type.n + 1][i] - 'a'));
    }
    if(strcmp(argv[type.n + 2], "-") != 0){
      for(int i = 0; i < strlen(argv[type.n + 2]); i++)
	      right_invariance |= (1 << (argv[type.n + 2][i] - 'a'));
    }
  }

  // generate graph
  // dq is the Weyl graph

  dq = weyl_generate_bruhat(type, left_invariance, right_invariance);

  // print stuff

  // The system output_level sets the level of detail. This would be if you ran it with a debugger?
  int output_level = 2;
  if(getenv("OUTPUT_LEVEL"))
    output_level = atoi(getenv("OUTPUT_LEVEL"));

  rank = weyl_rank(type);                // number of simple roots
  order = weyl_order(type);              // number of Weyl group elements
  positive = weyl_positive(type);        // number of positive roots

  if(output_level >= 1) {
    if(left_invariance) {
      printf("<");
      for(int j = 0; j < rank; j++)
	if(left_invariance & BIT(j))
	  fputc(alphabet[j], stdout);
      printf("> \\ ");
    }

    for(int i = 0; i < type.n; i++)
      printf("%s%c%d", i == 0 ? "" : " x ", type.factors[i].series, type.factors[i].rank);

    if(right_invariance) {
      printf(" / <");
      for(int j = 0; j < rank; j++)
	if(right_invariance & BIT(j))
	  fputc(alphabet[j], stdout);
      printf(">");
    }
    fprintf(stdout, "\n");

    // Top message printout
    fprintf(stdout, "Rank: %d\tOrder: %d\tPositive Roots: %d\tCosets: %d\n\n", rank, order, positive, dq->count);
  }

  if(output_level >= 3) {
    fprintf(stdout, "Shortest coset representatives: \n");
    for(int i = 0, wl = 0; i < dq->count; i++) {
      if(dq->cosets[i].min->wordlength > wl) {
	printf("\n");
	wl = dq->cosets[i].min->wordlength;
      }
      //      fprintf(stdout, "%s(%d) ", alphabetize(dq->cosets[i].min, stringbuffer), dq->cosets[i].max->wordlength);
      fprintf(stdout, "%s ", alphabetize(dq->cosets[i].min, stringbuffer));
    }
    fprintf(stdout, "\n\n");
  }

  if(output_level >= 4) {
    fprintf(stdout, "Bruhat order in graphviz format:\n");
    fprintf(stdout, "digraph test123 {\n");
    for(int i = 0; i < dq->count; i++)
      for(doublecoset_list_t *current = dq->cosets[i].bruhat_lower; current; current = current->next)
	fprintf(stdout, "%s -> %s;\n",
		alphabetize(dq->cosets[i].min, stringbuffer),
		alphabetize(current->to->min, stringbuffer2));
    fprintf(stdout, "}\n\n");
  }

  if(output_level >= 4) {
    fprintf(stdout, "Opposites:\n");
    for(int i = 0; i < dq->count; i++)
      fprintf(stdout, "%s <-> %s\n",
	      alphabetize(dq->cosets[i].min, stringbuffer),
	      alphabetize(dq->cosets[i].opposite->min, stringbuffer2));
    fprintf(stdout, "\n");
  }

  // Check if there were no balanced ideals
  fixpoints = 0;
  for(int i = 0; i < dq->count; i++)
    if(dq->cosets[i].opposite == &dq->cosets[i]) {
      if(output_level >= 1) {
	if(fixpoints == 0)
	  fprintf(stdout, "No balanced ideals since the longest element fixes the following cosets:");
	fprintf(stdout, " %s", alphabetize(dq->cosets[i].min, stringbuffer));
      }
      fixpoints++;
    }
  if(output_level >= 1 && fixpoints)
    fprintf(stdout, "\n\n");

  // If there were balanced ideals then print a message
  if(!fixpoints) {
    int *buffer = (int*)malloc(dq->count*sizeof(int));

    info_t info;
    info.dq = dq;
    info.rank = weyl_rank(type);
    info.order = weyl_order(type);
    info.positive = weyl_positive(type);
    info.buffer = buffer;
    info.level = output_level;

    ERROR(dq->count > 2*BV_BLOCKSIZE*BV_RANK, "We can handle at most %d cosets. Increase BV_RANK if more is needed.\n", 2*BV_BLOCKSIZE*BV_RANK);

    long count;
    if(output_level >= 2) {
      fprintf(stdout, "Balanced ideals:\n");
      count = enumerate_balanced_thickenings(dq, balanced_thickening_callback, &info);
      fprintf(stdout, "\n");
    } else {
      long outputcount = 0;
      count = enumerate_balanced_thickenings(dq, balanced_thickening_simple_callback, &outputcount);
    }

    if(output_level >= 1)
      fprintf(stdout, "Found %ld balanced ideal%s\n", count, count == 1 ? "" : "s");
  }
  // Deconstruct the dq
  weyl_destroy_bruhat(dq);
  free(type.factors);

  return 0;
}
