#CC       = pgc++
CC       = g++
CCFLAGS  = -I./
#REGFLAGS = -fast 
#ACCFLAGS = -acc -ta=nvidia -Minfo=accel 
#ACCFLAGS = -acc -ta=tesla:managed -Minfo=accel
#OMPFLAGS = -fast -mp 
OMPFLAGS = -fopenmp
SOURCES = testbspline.cpp bspline.cpp
#BIN =  bspline bspline_omp bspline_acc
BIN =  bspline bspline_omp

all: $(BIN)

bspline: $(SOURCES)
	$(CC) $(REGFLAGS) -o $@ $^

#bspline_acc: $(SOURCES)
#	$(CC) $(CCFLAGS) $(ACCFLAGS) -o $@ $^

bspline_omp: $(SOURCES)
	$(CC) $(CCFLAGS) $(OMPFLAGS) -o $@ $^

clean:
	$(RM) $(BIN) *.o

