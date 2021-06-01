#include "weyl.h"
#include "queue.h"

#include <strings.h>
#include <stdio.h>
#include <memory.h>

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

int main(int argc, const char *argv[])
{
	semisimple_type_t type;
	unsigned long right_invariance, left_invariance;
	doublequotient_t *dq;
	const char *alphabet = "abcdefghijklmnopqrstuvwxyz";
	char stringbuffer[100];
	char stringbuffer2[100];

	ERROR(argc < 2, "Too few arguments!\n");

	type.n = 0;
	for(int i = 0; i < argc - 1; i++) {
		if(argv[i+1][0] < 'A' || argv[i+1][0] > 'G')
			break;
		type.n++;
	}

	type.factors = (simple_type_t*)malloc(type.n*sizeof(simple_type_t));
	for(int i = 0; i < type.n; i++) {
		type.factors[i].series = argv[i+1][0];
		type.factors[i].rank = argv[i+1][1] - '0';
		ERROR(argv[i+1][0] < 'A' || argv[i+1][0] > 'G' || argv[i+1][1] < '1' || argv[i+1][1] > '9', "Arguments must be Xn with X out of A-G and n out of 1-9\n");
	}

	left_invariance = right_invariance = 0;

	if(argc - type.n >= 3) {
		if(strcmp(argv[type.n + 1], "-") != 0)
			for(int i = 0; i < strlen(argv[type.n + 1]); i++)
				left_invariance |= (1 << (argv[type.n + 1][i] - 'a'));
		if(strcmp(argv[type.n + 2], "-") != 0)
			for(int i = 0; i < strlen(argv[type.n + 2]); i++)
				right_invariance |= (1 << (argv[type.n + 2][i] - 'a'));
	}

	// generate graph

	dq = weyl_generate_bruhat(type, left_invariance, right_invariance);

    fprintf(stdout, "digraph test123 {\n");
    for(int i = 0; i < dq->count; i++)
      for(doublecoset_list_t *current = dq->cosets[i].bruhat_lower; current; current = current->next)
        fprintf(stdout, "%s -> %s;\n",
                alphabetize(dq->cosets[i].min, stringbuffer),
                alphabetize(current->to->min, stringbuffer2));
    fprintf(stdout, "}\n\n");

	// clean up
	weyl_destroy_bruhat(dq);
	free(type.factors);
}
