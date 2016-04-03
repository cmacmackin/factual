VENDOR ?= GNU
VENDOR_ = $(shell echo $(VENDOR) | tr A-Z a-z)

ifeq ($(VENDOR_),gnu)
F90 = gfortran-5
PFUNIT = /opt/pfunit-gfortran-5
FCFLAGS = -Og -g -I$(INC) -J$(INC)
LDFLAGS = -Og -g 
COVFLAGS = -fprofile-arcs -ftest-coverage
else ifeq ($(VENDOR_),intel)
F90 = ifort
PFUNIT = /opt/pfunit-ifort
FCFLAGS = -O0 -g -I$(INC) -module $(INC) -traceback -assume realloc_lhs
LDFLAGS = -O0 -g -traceback
COVFLAGS = 
endif

#PFUNIT = /opt/pfunit-gfortran-5
export PFUNIT
export F90
export FCFLAGS
export LDFLAGS
export COVFLAGS

# flags for debugging or for maximum performance, comment as necessary
PWD = $(shell pwd)
INC = $(PWD)/mod

ARCHIVE = libfactual.a
SRC = ./src
TEST = ./tests
EXE = tests.x
#INC = ./mod
TESTINC = $(TEST)/mod

export ARCHIVE

# "make" builds all
all: lib

$(ARCHIVE): lib

lib:
	make -C $(SRC) lib
	cp $(SRC)/$(ARCHIVE) .

tests: $(EXE)
	./$(EXE) 

SUT: lib
	make -C $(TEST) tests

gcov: tests
	make -C $(SRC) gcov

%: $(ODIR)/%.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f
	$(FC) $(FCFLAGS) -o $@ -c $<

%.o: %.F
	$(FC) $(FCFLAGS) -o $@ -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -o $@ -c $<

%.o: %.f95
	$(FC) $(FCFLAGS) -o $@ -c $<

%.o: %.F95
	$(FC) $(FCFLAGS) -o $@ -c $<

clean:
	make -C $(SRC) clean
	make -C $(TEST) clean
	rm -f $(ARCHIVE)
	rm -f *.gcov

init:
	mkdir -p $(ODIR)
	mkdir -p $(INC)

echo:
	make -C $(SRC) echo


$(EXE): SUT
	$(F90) -o $@ -I$(PFUNIT)/mod -I$(PFUNIT)/include \
		-I$(TEST) -I$(TESTINC) \
		$(PFUNIT)/include/driver.F90 $(TEST)/*.o $(FCFLAGS) \
		 -L. -lfactual -L$(PFUNIT)/lib -lpfunit $(COVFLAGS)

