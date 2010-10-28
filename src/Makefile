CXX = g++
CXXFLAGS = -Wall -O3

all: whipplesim whippleeig whippleeig_gyro

whipplesim: whipplesim.o whippleutils.o whipple.o
	$(CXX) -Wall -o $@ -lgsl -lcblas $^

whippleeig: whipple.o whippleeig.o whippleutils.o gslVecUtils.o
	$(CXX) -Wall -o $@ -lgsl -lcblas $^

#whipplesim_gyro: whipplesim_gyro.o whipple_gyro.o whippleutils_gyro.o
#	$(CXX) -Wall -o $@ -lgsl -lcblas $^

whippleeig_gyro: whipple_gyro.o whippleeig_gyro.o whippleutils_gyro.o gslVecUtils.o
	$(CXX) -Wall -o $@ -lgsl -lcblas $^

whipple.o: whipple.cpp whipple.h
whippleeig.o: whippleeig.cpp
whipplesim.o: whipplesim.cpp
whippleutils.o: whippleutils.cpp whippleutils.h
gslVecUtils.o: gslVecUtils.cpp gslVecUtils.h

whipple_gyro.o: whipple_gyro.cpp whipple_gyro.h
whippleeig_gyro.o: whippleeig_gyro.cpp
#whipplesim_gyro.o: whipplesim_gyro.cpp
whippleutils_gyro.o: whippleutils_gyro.cpp whippleutils_gyro.h

clean: 
	rm -rf *.o whipplesim whippleeig whipplesim_gyro whippleeig_gyro

.PHONY: all clean
