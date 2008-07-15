# Definition de compilation
all: Liveness
clean:
	rm -f *.o
install:
	mkdir -p ../bin/
	cp liveness ../bin/

CXX		= g++
CC		= gcc
CFLAGS	= -O3 -Wall

# Compilation
CommonTools.o: CommonTools.h CommonTools.c
	$(CC) $(CFLAGS) -c CommonTools.c
Liveness.o: Liveness.c
	$(CC) $(CFLAGS) -c Liveness.c
Liveness: CommonTools.o Liveness.o
	$(CC) $(CFLAGS) -o liveness CommonTools.o Liveness.o -lm
