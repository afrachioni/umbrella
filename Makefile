.PHONY: help clean clean-all clean-%
EXE = driver_$@
SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)
ifeq (,$(findstring gitversion.cpp,$(SRC)))
	SRC += gitversion.cpp
endif

help:
	@echo ''
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'
	@echo ''
	@echo 'make machine             build umbrellaese where machine specifies'
	@echo '                         one of these makefiles from src/MAKE:'
	@echo ''
	@files="`ls MAKE/Makefile.*`"; \
		for file in $$files; do head -1 $$file; done
	@echo ''

.DEFAULT:
	@test -f MAKE/Makefile.$@
	@$(MAKE) -s $(MFLAGS) gitversion.cpp
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@cp MAKE/Makefile.$@ Obj_$@
	@$(MAKE) $(MFLAGS) -f Makefile.$@ ../driver_$@ \
		"OBJ = $(OBJ)" "SRC = $(SRC)" "TARGET = $@" -C Obj_$@

gitversion.cpp: .git/HEAD .git/index
	echo "// This file is generated automatically by the make process," > $@
	echo "// to provide information about the most recent Git commit." >> $@
	echo "// Editing or deleting it could confuse the build system." >> $@
	echo "const char *gitversion = \"$(shell git rev-parse HEAD)\";" >> $@
	echo "const char *gitmessage = \"$(shell git log --format=%s%b -n 1)\";" >> $@
	echo "const char *gitbranch = \"$(shell git rev-parse --abbrev-ref HEAD)\";" >> $@

%.d: ../%.cpp
	$(CXX) -M $(CXXFLAGS) $(LAMMPS_INC) $< -MF $@

%.o : ../%.cpp
	$(CXX) $(CXXFLAGS) $(LAMMPS_INC) -c $< -o $@
-include $(SRC:%.cpp=%.d)

../driver_$(TARGET): $(OBJ)
	$(CXX) $(LAMMPS_LIB) $(OBJ) $(LIBS) -o $@

clean:
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'

clean-all:
	rm -rf Obj_*

clean-%:
	rm -rf Obj_$(@:clean-%=%)
