# Makefile for testparticle.c

CC = gcc
CFLAGS = -I/Users/richardanderson/rebound -I/opt/homebrew/opt/gsl/include
LDFLAGS = -L/Users/richardanderson/rebound -lrebound -L/opt/homebrew/opt/gsl/lib -lgsl -lgslcblas -lm

TARGET = out
SRC = testparticle.c
DYLD_PATH = /Users/richardanderson/rebound

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

run: $(TARGET)
	DYLD_LIBRARY_PATH=$(DYLD_PATH) ./$(TARGET) 183.303 100000

clean:
	rm -f $(TARGET)

