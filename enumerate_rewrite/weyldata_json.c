#include "thickenings.h"
#include "weyl.h"
#include "queue.h"

#include <strings.h>
#include <stdio.h>
#include <time.h>

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

void json_balanced_thickening_callback(const bitvec_t *pos, int size, const enumeration_info_t *ei)
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
    if (totcount==0){
      printf("[");
    } else {
      printf(",\n");
    }

    printf("{\"id\":%ld, \"left\": [", totcount++);
    int second = 0;
    for(int j = 0; j < info->rank; j++){
      if (left_invariance & (1 << j)){
        if (second==1){
          printf(",");
        }
        printf("\"%c\"",j + 'a');
        second = 1;
      }
    }
    printf("], \"right\": [");
    second = 0;
    for(int j = 0; j < info->rank; j++){
      if (right_invariance & (1 << j) ){
        if (second==1){
          printf(",");
        }
        printf("\"%c\"", j + 'a');
        second = 1;
      }
    }
    second = 0;
    if(info->buffer) {
      bitvec_t low, high;
      bv_copy(pos, &low);
      bv_negate(pos, &high);
      printf("], \"gen\": [");
      
      for(int i = 0; i < size/2; i++) {
	      if(!bv_get_bit(&high, i))
	        continue;
        if (second==1)
          printf(",");
	      printf("\"%s\"", alphabetize(info->dq->cosets[size-1-i].min, stringbuffer));
        second = 1;

	      bv_difference(&high, &ei->principal_neg[size-1-i], &high);
	      bv_difference(&low,  &ei->principal_pos[size-1-i], &low);
      }

      for(int i = size/2 - 1; i >= 0; i--) {
	      if(!bv_get_bit(&low, i))
	        continue;
        if (second==1)
          printf(",");
        printf("\"%s\"", alphabetize(info->dq->cosets[i].min, stringbuffer));
        second = 1;

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

    printf("]}");
  }
}

void json_balanced_thickening_simple_callback(const bitvec_t *pos, int size, const enumeration_info_t *ei)
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
  int rank, order, positive;
  int fixpoints;
  const char* commands[3];
  commands[0] = "all";
  commands[1] = "elts";
  commands[2] = "ideals";

  doublequotient_t *dq;

  const char *alphabet = "abcdefghijklmnopqrstuvwxyz";

  // read arguments

  ERROR(argc < 3, "Too few arguments!\n\nUsage is \"%s A2A3\" with\nA2,A3 simple Weyl factors.\n\n",argv[0]);

  
  // Count the number of simple factors in the semisimple Weyl group
  type.n = strlen(argv[2])/2;
  // fprintf(stdout, "type.n=%d\n\n",type.n);

  // Allocate memory, then read in the actual simple factors by letter/number, e.g. A5 is series 'A' and rank '5'. Series is A-G and the max rank is 9.
  type.factors = (simple_type_t*)malloc(type.n*sizeof(simple_type_t));
  int index_shift = 0;
  int new_rank = 0;
  int new_n = 0;
  for(int i = 0; i < type.n; i++) {
    if (2*i+index_shift+1>2*type.n){
      break;
    }
    type.factors[i].series = argv[2][2*i+index_shift];
    new_rank = argv[2][2*i+1+index_shift] - '0';
    while(2*i+1+index_shift<2*type.n && argv[2][2*i+1+index_shift+1]>='0' && argv[2][2*i+1+index_shift+1] <= '9'){
      new_rank = new_rank*10 + (argv[2][2*i+1+index_shift+1] - '0');
      index_shift = index_shift+1;
    }
    type.factors[i].rank = new_rank;
    new_n = i+1;
    // fprintf(stdout, "type.factors[%d].series=%c, type.factors[%d].rank=%d\n\n",i,type.factors[i].series,i,type.factors[i].rank);
    // ERROR(argv[1][2*i] < 'A' || argv[1][2*i] > 'G' || argv[1][2*i+1] < '1' || argv[1][2*i+1] > '9', "Arguments must be Xn with X out of A-G and n out of 1-9\n");
  }
  type.n = new_n;

  rank = weyl_rank(type);                // number of simple roots
  order = weyl_order(type);              // number of Weyl group elements
  positive = weyl_positive(type);        // number of positive roots

  
  // If command is graphviz, then print only the graphviz
  if(strcmp(argv[1],"graphviz")==0) {
    dq = weyl_generate_bruhat(type, 0, 0);
    fprintf(stdout, "digraph %s {\n",argv[2]);
    for(int i = 0; i < dq->count; i++){
      for(doublecoset_list_t *current = dq->cosets[i].bruhat_lower; current; current = current->next){
	      fprintf(stdout, "%s -> %s;\n",
		    alphabetize(dq->cosets[i].min, stringbuffer),
		    alphabetize(current->to->min, stringbuffer2));
      }
    }
    fprintf(stdout, "}\n");

    // Deconstruct and return
    weyl_destroy_bruhat(dq);
    free(type.factors);
    return 0;
  }

  // If command is anything else, start with the general JSON output.

  // Create the cartan matrix
  int *cartan_matrix;
  cartan_matrix = (int*)malloc(rank*rank*sizeof(int));
  weyl_cartan_matrix(type, cartan_matrix); //cartan matrix
  
  // Create the JSON-formatted timestamp
  time_t now;
  struct tm * local;
  char buffer [80];
  time(&now);
  local = localtime(&now);
  strftime(buffer,80,"%FT%X.000Z",local);

  // Output the general JSON stuff
  fprintf(stdout,"{");
  fprintf(stdout,"\"timestamp\": \"%s\",\n",buffer); // TODO: Make it more like JSON
  fprintf(stdout,"\"creator\": \"%s\",\n",argv[0]);
  fprintf(stdout,"\"version\": \"0.0.1\",\n");
  fprintf(stdout, "\"cartan_type\": \"%s\",\n",argv[2]);
  fprintf(stdout, "\"summands\": [");
  for (int i=0; i<type.n; i++){
    fprintf(stdout, "\"%c%d\"",type.factors[i].series,type.factors[i].rank);
    if (i<type.n-1) {
      fprintf(stdout, ",");
    }
  }
  fprintf(stdout, "],\n");
  fprintf(stdout, "\"rank\": %d,\n\"weyl_order\": %d,\n\"max_len\": %d,\n", rank, order, positive);
  
  // Print out the cartan matrix, formatted nicely
  fprintf(stdout, "\"cartan_matrix\":\n[");
  for (int i=0; i<rank; i++) {
    fprintf(stdout, "[ ");
    for (int j=0; j<rank; j++){
        // Make the spacing nice
        if (cartan_matrix[rank*i+j]>=0) {
            fprintf(stdout, " ");
        }
        fprintf(stdout, "%d",cartan_matrix[rank*i+j]);
        if (j<rank-1) {
          fprintf(stdout, ", ");
        } else {
          fprintf(stdout, " ");
        }
    }
    if (i<rank-1){
      fprintf(stdout, "],\n ");
    } else {
      fprintf(stdout, "]]");
    }
  }

  // print out weylgroup elements
  if (strcmp(argv[1],"elts")==0) {
    weylgroup_t *wgroup = weyl_generate(type); // TODO: This makes the code take much longer to run
    fprintf(stdout, ",\n\"elements\": [");
    for (int i=0; i<order; i++){
      if (i!= 0){
        fprintf(stdout, ", ");
      }
      fprintf(stdout, "\"%s\"", alphabetize(&wgroup->elements[i], stringbuffer));
    }
    fprintf(stdout,"]");
    weyl_destroy(wgroup);
  }
  
  // Print out the balanced ideals
  if (strcmp(argv[1],"ideals")==0) {
    dq = weyl_generate_bruhat(type, 0, 0);
    // Check if there were no balanced ideals
    fixpoints = 0;
    for(int i = 0; i < dq->count; i++)
      if(dq->cosets[i].opposite == &dq->cosets[i]) {
        if(fixpoints == 0)
          fprintf(stdout, "No balanced ideals since the longest element fixes the following cosets:");
        fprintf(stdout, " %s", alphabetize(dq->cosets[i].min, stringbuffer));
        fixpoints++;
      }
    if(fixpoints)
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

      ERROR(dq->count > 2*BV_BLOCKSIZE*BV_RANK, "We can handle at most %d cosets. Increase BV_RANK if more is needed.\n", 2*BV_BLOCKSIZE*BV_RANK);

      long count;
      fprintf(stdout, ",\n\"balanced_ideals\":\n");
      count = enumerate_balanced_thickenings(dq, json_balanced_thickening_callback, &info);
      fprintf(stdout, "],");

      fprintf(stdout, "\n\"num_balanced_ideals\":%ld", count);
    }
    // Deconstruct the dq
    weyl_destroy_bruhat(dq);
  }
  fprintf(stdout, "\n}\n");

  // free memory back
  free(type.factors);

  return 0;
}
