#ifndef QUEUE_H
#define QUEUE_H

#include <stdio.h>
#include <stdlib.h>

#define QUEUE_SIZE 5000

#define ERROR(condition, msg, ...) if(condition){fprintf(stderr, msg, ##__VA_ARGS__); exit(1);}
#ifdef _DEBUG
#define LOG(msg, ...) fprintf(stderr, msg, ##__VA_ARGS__)
#else
#define LOG(msg, ...)
#endif

typedef struct {
  unsigned int start;
  unsigned int end;
  int data[QUEUE_SIZE];
} queue_t;

static void queue_init(queue_t *q)
{
  q->start = q->end = 0;
}

static void queue_put(queue_t *q, int x)
{
  q->data[q->end++] = x;
  q->end %= QUEUE_SIZE;

  ERROR(q->start == q->end, "The queue is full! Increase QUEUE_SIZE\n");
}

static int queue_get(queue_t *q)
{
  if(q->start == q->end)
    return -1;

  int result = q->data[q->start++];
  q->start %= QUEUE_SIZE;

  return result;
}

#endif
