# Definition de compilation
all: liveness
clear: rm *.o

CXX		= g++
CC		= gcc
CFLAGS	= -O3 -Wall

# Compilation
commonTools.o: commonTools.h commonTools.c
	$(CC) $(CFLAGS) -c commonTools.c
liveness.o: liveness.c
	$(CC) $(CFLAGS) -c liveness.c
liveness: commonTools.o liveness.o
	$(CC) $(CFLAGS) -o liveness commonTools.o liveness.o -lm
