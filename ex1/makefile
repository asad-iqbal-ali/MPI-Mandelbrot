OUTPUT=mandelbrot
CFLAGS=-Wall -g

all: $(OUTPUT)

$(OUTPUT): ex1.c
	mpicc -o $(OUTPUT) $(CFLAGS) ex1.c

clean:
	-rm $(OUTPUT)
