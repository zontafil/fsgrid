CC=mpiCC
CFLAGS=-O2 -std=c++11 -march=native -g -Wall

test: test.cpp fsgrid.hpp
	$(CC) $(CFLAGS) -o $@ $<
