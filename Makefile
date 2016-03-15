CXX=mpiCC
CXXFLAGS= -std=c++11 -march=native -g -Wall

test: test.cpp fsgrid.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
ddtest: ddtest.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	-rm test
