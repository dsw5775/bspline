CC       = nvcc
CCFLAGS  = -I./
REGFLAGS = -O2  
SOURCES = testbspline.cu bspline.cu
BIN = bspline 

all: $(BIN)

bspline: $(SOURCES)
	$(CC) $(REGFLAGS) -o $@ $^

clean:
	$(RM) $(BIN)

