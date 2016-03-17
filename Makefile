CXX=CC
CXXFLAGS= -std=c++11 -march=native -g -Wall

benchmark: benchmark.cpp fsgrid.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
test: test.cpp fsgrid.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
ddtest: ddtest.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	-rm test
