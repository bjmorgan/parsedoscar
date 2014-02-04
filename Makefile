OBJS := parsedoscar.o class_atom.o class_orbital.o class_species.o class_spin.o class_dos_file.o
EXECUTABLE := parsedoscar
COMPILER := ifort

all: app

debug: FFLAGS = -check all
debug: app

opt: FFLAGS = -fast -O2
opt: app

app: $(OBJS)
	$(COMPILER) $(FFLAGS) $(OBJS) -o $(EXECUTABLE)

parsedoscar.o: parsedoscar.f90 class_atom.mod class_orbital.mod class_species.mod class_spin.mod
	$(COMPILER) $(FFLAGS) -c parsedoscar.f90 

class_atom.o: class_atom.f90 class_orbital.mod class_dos_file.mod
	$(COMPILER) $(FFLAGS) -c class_atom.f90 

class_atom.mod: class_atom.f90 class_orbital.mod class_dos_file.mod
	$(COMPILER) $(FFLAGS) -c class_atom.f90 
	
class_orbital.o: class_orbital.f90 class_spin.mod
	$(COMPILER) $(FFLAGS) -c class_orbital.f90
	
class_orbital.mod: class_orbital.f90 class_spin.mod
	$(COMPILER) $(FFLAGS) -c class_orbital.f90
	
class_species.o: class_species.f90 class_atom.mod class_dos_file.mod
	$(COMPILER) $(FFLAGS) -c class_species.f90
	
class_species.mod: class_species.f90 class_atom.mod class_dos_file.mod
	$(COMPILER) $(FFLAGS) -c class_species.f90

class_spin.o: class_spin.f90
	$(COMPILER) $(FFLAGS) -c class_spin.f90
	
class_spin.mod: class_spin.f90
	$(COMPILER) $(FFLAGS) -c class_spin.f90

class_dos_file.o: class_dos_file.f90
	$(COMPILER) $(FFLAGS) -c class_dos_file.f90
	
class_dos_file.mod: class_dos_file.f90
	$(COMPILER) $(FFLAGS) -c class_dos_file.f90

clean:
	rm -f *.o *.mod $(EXECUTABLE)
