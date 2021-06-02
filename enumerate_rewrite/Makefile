HEADERS=weyl.h thickenings.h queue.h bitvec.h

#SPECIAL_OPTIONS=-O0 -g -D_DEBUG
#SPECIAL_OPTIONS=-O3 -pg -funroll-loops -fno-inline
SPECIAL_OPTIONS=-O3 -flto -funroll-loops -Winline

OPTIONS=-m64 -march=native -mtune=native -std=gnu99 -D_GNU_SOURCE $(SPECIAL_OPTIONS)
NAME=enumerate-rewrite

all: enumerate graph weyldata weyldata_json

$(NAME).tar.bz2: $(NAME) $(HEADERS) enumerate.c weyl.c thickenings.c
	tar cjhf $(NAME).tar.bz2 $(NAME)/weyldata.c $(NAME)/weyldata_json.c $(NAME)/enumerate.c $(NAME)/weyl.c $(NAME)/thickenings.c $(NAME)/weyl.h $(NAME)/thickenings.h $(NAME)/queue.h $(NAME)/bitvec.h $(NAME)/Makefile $(NAME)/graph.c

$(NAME):
	ln -s . $(NAME)

weyldata_json: weyldata_json.o weyl.o thickenings.o
	gcc $(OPTIONS) -o weyldata_json weyldata_json.o thickenings.o weyl.o

weyldata: weyldata.o weyl.o thickenings.o
	gcc $(OPTIONS) -o weyldata weyldata.o thickenings.o weyl.o

enumerate: enumerate.o weyl.o thickenings.o
	gcc $(OPTIONS) -o enumerate enumerate.o thickenings.o weyl.o

graph: graph.o weyl.o
	gcc $(OPTIONS) -o graph graph.o weyl.o

weyldata_json.o: weyldata_json.c $(HEADERS)
	gcc $(OPTIONS) -c weyldata_json.c

weyldata.o: weyldata.c $(HEADERS)
	gcc $(OPTIONS) -c weyldata.c

enumerate.o: enumerate.c $(HEADERS)
	gcc $(OPTIONS) -c enumerate.c

thickenings.o: thickenings.c $(HEADERS)
	gcc $(OPTIONS) -c thickenings.c

weyl.o: weyl.c $(HEADERS)
	gcc $(OPTIONS) -c weyl.c

graph.o: graph.c $(HEADERS)
	gcc $(OPTIONS) -c graph.c

clean:
	rm -f enumerate graph weyldata weyldata_json thickenings.o weyldata.o weyldata_json.o weyl.o enumerate.o graph.o $(NAME) $(NAME).tar.bz2
