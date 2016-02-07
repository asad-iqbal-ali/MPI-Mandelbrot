OUTPUT=btes
CFLAGS=-Wall -g

all: $(OUTPUT)

$(OUTPUT): ex2.c
	mpicc -o $(OUTPUT) $(CFLAGS) ex2.c

clean:
	-rm $(OUTPUT)
