CC = gcc
SRC = test.c
TARGET = out

REBOUND_PATH = /Users/richardanderson/rebound
GSL_PATH = /opt/homebrew/opt/gsl

CFLAGS = -I$(REBOUND_PATH)/src -I$(GSL_PATH)/include
LDFLAGS = -L$(REBOUND_PATH) -lrebound -L$(GSL_PATH)/lib -lgsl -lgslcblas -lm

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

run: $(TARGET)
	DYLD_LIBRARY_PATH=$(REBOUND_PATH) ./$(TARGET)

clean:
	rm -f $(TARGET)

