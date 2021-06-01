#ifndef THICKENINGS_H
#define THICKENINGS_H

#include "bitvec.h"
#include "weyl.h"

#define DEBUG(msg, ...) do{fprintf(stderr, msg, ##__VA_ARGS__); }while(0)

struct enumeration_info {
  int size;        // the size of the weyl group. We store however only the first size/2 elements
  bitvec_t *principal_pos;
  bitvec_t *principal_neg;
  int *principal_is_slim;
  void (*callback)(const bitvec_t *, int, const struct enumeration_info *);
  void *callback_data;
};

typedef void (*enumeration_callback)(const bitvec_t *, int, const struct enumeration_info *);
typedef struct enumeration_info enumeration_info_t;

// enumerating balanced thickenings
long enumerate_balanced_thickenings(doublequotient_t *dq, enumeration_callback callback, void *callback_data);

#endif
