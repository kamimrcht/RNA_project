.PHONY: all clean

CC=/usr/bin/g++

CFLAGS=  -Wall  -Ofast -std=c++11 -march=native -pthread -fmax-errors=3 -flto -pipe
LDFLAGS=-pthread

ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O3  -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O3 -g
LDFLAGS=-g
endif

SRCS = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp,%,$(SRCS))

all: $(OBJS)

clean:
	rm -fv $(OBJS)

rebuild: clean $(OBJS)

%: %.cpp
	$(CC) $(CFLAGS) $^ -o $@  $(LDFLAGS)
