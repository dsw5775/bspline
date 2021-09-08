#CC       = pgc++
FC       = gfortran
FCFLAGS  = -Ofast
#REGFLAGS = -fast 
#ACCFLAGS = -acc -ta=nvidia -Minfo=accel 
#ACCFLAGS = -acc -ta=tesla:managed -Minfo=accel
#OMPFLAGS = -fast -mp 
OMPFLAGS = -fopenmp
SOURCES = bspline.f90 testbspline.f90
#BIN =  bspline bspline_omp bspline_acc
BIN =  bspline bspline_omp

all: $(BIN)

bspline: $(SOURCES)
	$(FC) $(FCFLAGS) $(REGFLAGS) -o $@ $^

bspline_omp: $(SOURCES)
	$(FC) $(FCFLAGS) $(OMPFLAGS) -o $@ $^

clean:
	$(RM) $(BIN) *.o *.mod

