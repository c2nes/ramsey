CC=gcc
CFLAGS= --std=c99 -Wall -pedantic -O2 -funroll-loops
#CFLAGS= --std=c99 -Wall -pedantic -pg -g

PRGMS=find_cliques extend_graph

all: $(PRGMS)

clean:
	rm $(PRGMS)

..c:
	$(CC) $(CFLAGS) -o $@ $<

.PHONY: all clean
