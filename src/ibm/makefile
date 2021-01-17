# makefile to compile the
# individual based simulation into
# a runnable executable file

EXE=sexsel_intragenomic.exe
CPP=sexsel_intragenomic.cpp individual.cpp patch.cpp
HPP=individual.hpp patch.hpp
LIBS=yaml-cpp
LFLAGS=-L$(HOME)/usr/local/lib64
IFLAGS=-I$(HOME)/usr/local/include
CFLAGS=-ggdb -std=c++17 -Wall -O3
all : $(EXE)

HOST=$(shell hostname | grep athena)
compileCMD=""
ifeq ($(HOST),)
compileCMD=g++ $(CFLAGS) -o $@ $(CPP) -l$(LIBS)
else
compileCMD=g++ $(CFLAGS) $(LFLAGS) $(IFLAGS) -o $@ $(CPP) -l$(LIBS)
endif

# conditional compilation 
$(EXE): $(CPP) $(HPP)
	$(compileCMD)


.PHONY :
	clean

clean :
	rm -rf $(EXE)
