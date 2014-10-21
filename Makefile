EXE = driver_$@
#vpath %.cpp ..
#vpath %.h ..
#SRC = $(wildcard *.cpp)
#INC = $(wildcard *.h)
#OBJ = $(SRC:.cpp=.o)
#"OBJ = $(OBJ)" "SRC = $(SRC)"

.DEFAULT:
	@test -f MAKE/Makefile.$@
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@cp MAKE/Makefile.$@ Obj_$@
	@cd Obj_$@; $(MAKE) $(MFLAGS) -f Makefile.$@ $@

clean:
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'

clean-all:
	rm -rf Obj_*

clean-%:
	rm -rf Obj_$(@:clean-%=%)
