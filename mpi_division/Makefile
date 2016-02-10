FC=mpif90 -O0 -g -traceback -check bounds  
#FOPT=-Wall -Wextra #-fopenmp
OUT=exe_simu

all: mod_random.o mod_constants.o mod_variables.o main.o
	$(FC) $(FOPT) $^ -o $(OUT)

#=========Executable=========#
main.o: main.f90  
	$(FC) $(FOPT) -c $^

mod_constants.o: mod_constants.f90
	$(FC) $(FOPT) -c $^

mod_variables.o: mod_variables.f90
	$(FC) $(FOPT) -c $^

mod_random.o: mod_random.f90
	 $(FC) $(FOPT) -c $^
#==========Autre cible==========#
clean:
	\rm -rf *.o $(OUT) *.mod

run: all
	time ./$(OUT)

rerun: clean all run


