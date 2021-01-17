CPP=tsd_multivariate.cpp main.cpp
TEST_CPP=tsd_multivariate.cpp testsuite.cpp

HPP=tsd_multivariate.hpp
EXE=tsd_multivariate.exe
TEST_EXE=tsd_testsuite.exe

CXX=g++

CXXFLAGS=-Wall -std=c++17 -O3 -lgsl -lgslcblas

IFLAGS=-I /usr/local/Cellar/boost/1.74.0/include/boost/test/unit_test.hpp
LFLAGS=-L /usr/local/lib/
lFlags=-lboost_unit_test_framework-mt

all : $(EXE) $(TEST_EXE)

$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

$(TEST_EXE) : $(TEST_CPP) $(HPP)
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o $(TEST_EXE) $(TEST_CPP) $(LFLAGS) $(lFlags)

clean :
	rm -rf $(EXE)
	rm -rf $(TEST_EXE)
