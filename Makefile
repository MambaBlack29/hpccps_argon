FC = gfortran 
FFLAGS = -O3 -fopenmp
TARGET = mdrun 

OBJECTS = main.o init.o force.o integrate.o temp.o 

$(TARGET) : $(OBJECTS) 
	$(FC) $(FFLAGS) -o $@ $^  
%.o : %.f90 
	$(FC) $(FFLAGS) -c $< 
	
.PHONY :  clean 

clean : 
	rm -f *.o *.mod $(TARGET)

