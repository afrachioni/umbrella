CXX = mpic++
CPPFLAGS = 
INC = -I ~/myLAMMPS_2/src -I ~/myLAMMPS_2/lib/meam
LIB_PATHS = -L ~/myLAMMPS_2/src -L ~/myLAMMPS_2/lib/meam
LIBS = -llammps_vulcan -lmeam -lifcore -lsvml -limf
SRC = $(wildcard *.cpp)
OBJ = $(SRC:%.cpp=%.o)

%.d: %.cpp
	$(CXX) -M $(CPPFLAGS) $(INC) $< -MF $@

%.o : %.cpp 
	$(CXX) $(CPPFLAGS) $(INC) -c $< -o $@
-include $(SRC:%.cpp=%.d)

driver_vulcan: $(OBJ)
	$(CXX) $(LIB_PATHS) $(OBJ) $(LIBS) -o $@

vulcan: driver_vulcan
