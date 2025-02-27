 # the compiler: gcc for C program, define as g++ for C++
CC = gcc

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
#CFLAGS = -g -Wall
CFLAGS = -Wall -std=c99 -O3
TARGET = GEMF
all: $(TARGET)

$(TARGET): gemfc_nrm.c nrm.o para.o common.o pcg_basic.o
	rm -rf $(TARGET)
	$(CC) $(CFLAGS) -o $(TARGET) gemfc_nrm.c nrm.o para.o common.o pcg_basic.o -lm

nrm.o:  nrm.c nrm.h
	$(CC) $(CFLAGS) -c nrm.c
para.o:  para.c para.h
	$(CC) $(CFLAGS) -c para.c
common.o:  common.c common.h
	$(CC) $(CFLAGS) -c common.c
pcg_basic.o: pcg_basic.c pcg_basic.h
	$(CC) $(CFLAGS) -c pcg_basic.c

clean:
	rm -rf $(TARGET)
	rm -rf nrm.o
	rm -rf common.o
	rm -rf para.o
	rm -rf pcg_basic.o


