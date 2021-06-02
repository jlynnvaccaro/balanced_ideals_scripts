#include "thickenings.h"
#include "weyl.h"
#include "queue.h"

#include <strings.h>
#include <stdio.h>
#include <time.h>

int main(int argc, const char *argv[])
{
  semisimple_type_t type;
  int rank, order, positive;

  doublequotient_t *dq;

  const char *alphabet = "abcdefghijklmnopqrstuvwxyz";

  // read arguments

  ERROR(argc < 2, "Too few arguments!\n\nUsage is \"%s A2A3\" with\nA2,A3 simple Weyl factors.\n\n",argv[0]);
  
  // Count the number of simple factors in the semisimple Weyl group
  type.n = strlen(argv[1])/2;
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
    type.factors[i].series = argv[1][2*i+index_shift];
    new_rank = argv[1][2*i+1+index_shift] - '0';
    while(2*i+1+index_shift<2*type.n && argv[1][2*i+1+index_shift+1]>='0' && argv[1][2*i+1+index_shift+1] <= '9'){
      new_rank = new_rank*10 + (argv[1][2*i+1+index_shift+1] - '0');
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

  int *cartan_matrix;
  cartan_matrix = (int*)malloc(rank*rank*sizeof(int));
  weyl_cartan_matrix(type, cartan_matrix); //cartan matrix

  // JSON OUTPUT BEGINS HERE
  time_t now;
  struct tm * local;
  char buffer [80];
  time(&now);
  local = localtime(&now);
  strftime(buffer,80,"%FT%X.000Z",local);

  fprintf(stdout,"{");
  fprintf(stdout,"\"timestamp\":\"%s\",\n",buffer); // TODO: Make it more like JSON
  fprintf(stdout,"\"creator\":\"%s\",\n",argv[0]);
  fprintf(stdout,"\"version\":\"0.0.1\",\n");
  fprintf(stdout, "\"cartan_type\":\"%s\",\n",argv[1]);
  fprintf(stdout, "\"summands\":[");
  for (int i=0; i<type.n; i++){
    fprintf(stdout, "\"%c%d\"",type.factors[i].series,type.factors[i].rank);
    if (i<type.n-1) {
      fprintf(stdout, ",");
    }
  }
  fprintf(stdout, "],\n");
  fprintf(stdout, "\"rank\": %d,\n\"weyl_order\": %d,\n\"max_len\": %d,\n", rank, order, positive);
  
  fprintf(stdout, "\"cartan_matrix\":\n[");
  for (int i=0; i<rank; i++) {
      fprintf(stdout, "[ ");
      for (int j=0; j<rank; j++){
          // Make the spacing nice
          // TODO: Could also be rank*j+i, which will make the transpose.
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
        fprintf(stdout, "]]\n");
      }
  }
  fprintf(stdout, "}\n");

  // free memory back
  free(type.factors);

  return 0;
}
