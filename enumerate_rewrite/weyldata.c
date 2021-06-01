#include "thickenings.h"
#include "weyl.h"
#include "queue.h"

#include <strings.h>
#include <stdio.h>

int main(int argc, const char *argv[])
{
  semisimple_type_t type;
  unsigned long right_invariance, left_invariance;
  int rank, order, positive;
  int fixpoints;

  doublequotient_t *dq;

  const char *alphabet = "abcdefghijklmnopqrstuvwxyz";

  // read arguments

  ERROR(argc < 2, "Too few arguments!\n\nUsage is '%s A2A3' or '%s A2A3 -abc -abc' with\nA2,A3 simple Weyl factors and abc,abc left/right invariance.\n\n",argv[0],argv[0]);
  
  // Count the number of simple factors in the semisimple Weyl group
  type.n = strlen(argv[1])/2;
  fprintf(stdout, "type.n=%d\n\n",type.n);

  // Allocate memory, then read in the actual simple factors by letter/number, e.g. A5 is series 'A' and rank '5'. Series is A-G and the max rank is 9.
  type.factors = (simple_type_t*)malloc(type.n*sizeof(simple_type_t));
  for(int i = 0; i < type.n; i++) {
    type.factors[i].series = argv[1][2*i];
    type.factors[i].rank = argv[1][2*i+1] - '0';
    fprintf(stdout, "type.factors[%d].series=%c, type.factors[%d].rank=%d\n\n",i,type.factors[i].series,i,type.factors[i].rank);
    ERROR(argv[1][2*i] < 'A' || argv[1][2*i] > 'G' || argv[1][2*i+1] < '1' || argv[1][2*i+1] > '9', "Arguments must be Xn with X out of A-G and n out of 1-9\n");
  }

  left_invariance = right_invariance = 0;
  // Additional command line arguments that were not factors are the left/right invariance.
  if(argc >= 3) {
    if(strcmp(argv[2], "-") != 0){
      printf("%s\n",argv[type.n+1]);
      for(int i = 0; i < strlen(argv[type.n + 1]); i++)
	      left_invariance |= (1 << (argv[type.n + 1][i] - 'a'));
    }
    if(strcmp(argv[3], "-") != 0){
      for(int i = 0; i < strlen(argv[type.n + 2]); i++)
	      right_invariance |= (1 << (argv[type.n + 2][i] - 'a'));
    }
  }

  // generate graph
  // dq is the Weyl graph

  dq = weyl_generate_bruhat(type, left_invariance, right_invariance);

  // print stuff

  rank = weyl_rank(type);                // number of simple roots
  order = weyl_order(type);              // number of Weyl group elements
  positive = weyl_positive(type);        // number of positive roots

  int *cartan_matrix;
  cartan_matrix = (int*)malloc(rank*rank*sizeof(int));
  weyl_cartan_matrix(type, cartan_matrix); //cartan matrix

  fprintf(stdout, "Rank: %d\nOrder: %d\nPositive Roots: %d\nCosets: %d\n\n", rank, order, positive, dq->count);
  
  fprintf(stdout, "Cartan matrix for %s:\n",argv[1]);
  for (int i=0; i<rank; i++) {
      fprintf(stdout, "[ ");
      for (int j=0; j<rank; j++){
          // Make the spacing nice
          // TODO: Could also be rank*k+j, which will make the transpose.
          if (cartan_matrix[rank*i+j]>=0) {
              fprintf(stdout, " ");
          }
          fprintf(stdout, "%d ",cartan_matrix[rank*i+j]);
      }
      fprintf(stdout, "]\n");
  }
  fprintf(stdout, "\n");


  // Deconstruct the dq
  weyl_destroy_bruhat(dq);
  free(type.factors);

  return 0;
}
