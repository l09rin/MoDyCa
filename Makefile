CXX = g++
CXXFLAGS = -Wall -std=c++11
MATHLIB = -lm
THREADLIB = -lpthread

DEBUG_FLAGS = -g -fno-inline
RELEASE_FLAGS = -O3
# Check if DEBUG is set to YES
ifeq ($(DEBUG), YES)
    CXXFLAGS += $(DEBUG_FLAGS)
else
    CXXFLAGS += $(RELEASE_FLAGS)
endif

SRC_MDMC = $(wildcard src/numerical_simulation/*.cpp src/numerical_simulation/*.h src/numerical_simulation/*.tpp src/numerical_simulation/*.ipp src/numerical_simulation/Integrators/*.cpp src/numerical_simulation/Interactions/*.tpp)

_SRC_MYLIBS = matrix.cpp record.cpp
SRC_MYLIBS = $(patsubst %,src/%,$(_SRC_MYLIBS))
OBJ_MYLIBS = $(SRC_MYLIBS:.cpp=.o)

_DEPS_MYLIBS = matrix.h record.h system_comm.h system_comm.tpp
DEPS_MYLIBS = $(patsubst %,src/%,$(_DEPS_MYLIBS))

SRC_MAIN = src/calculate_GofR.cpp src/histogram.cpp src/mcmd.cpp
OBJ = $(SRC_MAIN:.cpp=.o)
_EXE = $(subst src/,,$(SRC_MAIN))
EXE = $(_EXE:.cpp=.exe)

all : $(EXE)  ## Build all objects and main executables
	@echo "Built all the targets:"
	@echo $(EXE)

$(OBJ_MYLIBS): %.o: %.cpp ${DEPS_MYLIBS}  # Compile library object files
	$(CXX) $(CXXFLAGS) -c $< -o $@


calculate_GofR.o : src/calculate_GofR.cpp $(SRC_MDMC) $(SRC_MYLIBS) $(DEPS_MYLIBS)  # Compile pair distribution function code object
	$(CXX) $(CXXFLAGS) -c src/$(@:.o=.cpp) -o $@

histogram.o : src/histogram.cpp $(SRC_MYLIBS) $(DEPS_MYLIBS)  # Compile histogram code object
	$(CXX) $(CXXFLAGS) -c src/$(@:.o=.cpp) -o $@

mcmd.o : src/mcmd.cpp $(SRC_MDMC) $(SRC_MYLIBS) $(DEPS_MYLIBS)  # Compile mcmd object
	$(CXX) $(CXXFLAGS) -c src/$(@:.o=.cpp) -o $@



calculate_GofR.exe : src/calculate_GofR.o $(OBJ_MYLIBS)  ## Link final executable of the utility to calculate pair distribution functions
	$(CXX) $(CXXFLAGS) src/$(@:.exe=.o) $(OBJ_MYLIBS) $(MATHLIB) $(THREADLIB) -o $@

histogram.exe : src/histogram.o $(OBJ_MYLIBS)  ## Link final executable of the utility to generate histograms
	$(CXX) $(CXXFLAGS) src/$(@:.exe=.o) $(OBJ_MYLIBS) $(MATHLIB) -o $@

mcmd.exe : src/mcmd.o $(OBJ_MYLIBS)  ## Link final executable of the Molecular Dynamics and MonteCarlo code
	$(CXX) $(CXXFLAGS) src/$(@:.exe=.o) $(OBJ_MYLIBS) $(MATHLIB) -o $@


clean-obj :  # Remove all object files
	rm $(OBJ) $(OBJ_MYLIBS)

clean :  ## Remove objects and executables
	rm $(EXE) $(OBJ) $(OBJ_MYLIBS)

help:  ## Show this help
	@echo "Available targets:"
	@grep -E '^[a-zA-Z._-]+ *:.*?##' $(MAKEFILE_LIST) | \
	    awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

.PHONY : all clean clean-obj help
