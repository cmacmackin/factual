# libraries needed for linking, unused in the examples
SOURCES := $(shell find -name '*.f90')
SOURCESPP := $(shell find -name '*.F90')
OBJS = $(SOURCES:%.f90=%.o) $(SOURCESPP:%.F90=%.o)

# List of executables to be built within the package
INC = ../mod
FCFLAGS += $(COVFLAGS)

# "make" builds all
all: lib

lib: init $(ARCHIVE)

$(ARCHIVE): $(OBJS)
	ar rcs $@ $^

cheb1d_fields.o: cheb1d_fields.F90 abstract_fields.o chebyshev.o
factual.o: factual.f90 cheb1d_fields.o abstract_fields.o

gcov:
	gcov *.gcno -pb

%: $(ODIR)/%.o
	$(F90) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(F90) $(FCFLAGS) -o $@ -c $<

%.o: %.F90
	$(F90) $(FCFLAGS) -o $@ -c $<

clean:
	/bin/rm  -rf *.o ./tmp $(INC)/*.mod *.gcno *.gcda *.gcov *.s *~

init:
	mkdir -p $(INC)

echo:
	echo "$(FCFLAGS)"
	echo $(INC)
