EXE = driver_$@
SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

.DEFAULT:
	@test -f MAKE/Makefile.$@
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@cp MAKE/Makefile.$@ Obj_$@
	@cd Obj_$@; $(MAKE) $(MFLAGS) -f Makefile.$@ ../driver_$@ \
		"OBJ = $(OBJ)" "SRC = $(SRC)" "TARGET = $@"

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

clean:
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'

clean-all:
	rm -rf Obj_*

clean-%:
	rm -rf Obj_$(@:clean-%=%)
