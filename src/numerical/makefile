CPP=tsd_multivariate.cpp main.cpp
TEST_CPP=tsd_multivariate.cpp tsd_testsuite.cpp

HPP=tsd_multivariate.hpp
EXE=tsd_multivariate.exe
TEST_EXE=tsd_testsuite.exe

CXX=g++

# flags for any compilation
CXXFLAGS=-Wall -std=c++17 -O3 -lgsl -lgslcblas
LXXFLAGS=-lgsl -lgslcblas

# flags to compile boost etc
IFLAGS=-I /usr/local/Cellar/boost/1.75.0_1/include
LFLAGS=-L /usr/local/lib/
lFlags=/usr/local/lib/libboost_unit_test_framework-mt.a

IFLAGS_GSL=-I/opt/homebrew/include
LFLAGS_GSL=-L/opt/homebrew/lib

# build everything
#all : $(EXE) $(TEST_EXE)
all : $(EXE)


# build the normal executable
$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) $(IFLAGS_GSL) $(LFLAGS_GSL) -o $(EXE) $(CPP)  $(LXXFLAGS)

# build the testsuite
$(TEST_EXE) : $(TEST_CPP) $(HPP)
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o $(TEST_EXE) $(TEST_CPP) $(LFLAGS) $(lFlags)

clean :
	rm -rf $(EXE)
	rm -rf $(TEST_EXE)
