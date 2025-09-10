EXE=tsd_seasonal.exe
CPP=tsd_seasonal.cpp main.cpp individual.cpp patch.cpp
HPP=tsd_seasonal.hpp parameters.hpp individual.hpp patch.hpp
CXX=g++
CXXFLAGS=-Wall -O3 -std=c++20

$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

clean :
	rm -rf $(EXE)
