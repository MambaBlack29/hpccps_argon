FC = gfortran 
FFLAGS = -O3 
TARGET = mdrun.run

OBJECTS = main.o init.o force.o integrate.o temp.o verlet.o

$(TARGET) : $(OBJECTS) 
	$(FC) $(FFLAGS) -o $@ $^  
%.o : %.f90 
	$(FC) $(FFLAGS) -c $< 
	
.PHONY :  clean 

clean : 
	rm -f *.o *.mod $(TARGET)

