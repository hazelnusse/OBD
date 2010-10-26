CXX = g++
CXXFLAGS = -Wall -O3

all: whipplesim whippleeig

whipplesim: whipplesim.o whipple.o
	$(CXX) -Wall -o $@ -lgsl -lcblas $^

whippleeig: whipple.o whippleeig.o whippleutils.o gslVecUtils.o
	$(CXX) -Wall -o $@ -lgsl -lcblas $^

whipple.o: whipple.cpp whipple.h
whippleeig.o: whippleeig.cpp
whipplesim.o: whipplesim.cpp
whippleutils.o: whippleutils.cpp whippleutils.h
gslVecUtils.o: gslVecUtils.cpp gslVecUtils.h

clean: 
	rm -rf *.o whipplesim whippleeig *.dat *.pyc *.pyo

.PHONY: all clean
