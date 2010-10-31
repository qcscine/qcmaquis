BOOSTDIR = $(HOME)/dev/boost_1_23_0
#BOOSTDIR = $(HOME)/dev/boost_1_22_0
#BOOSTDIR = $(HOME)/dev/boost_1_21_1
CXX = KCC
#CXX = icc
#CXX = /usr/local/gcc-3.0/bin/g++
#CXX = g++
#OPTFLAGS = -O3 -tpp6
OPTFLAGS = +K2 -k -g
#OPTFLAGS = -w +K0 -g -k -v
#OPTFLAGS = -g -O7
#OPTFLAGS = -g -O5 -finline-functions -funroll-all-loops
CXXFLAGS = $(OPTFLAGS) --linux_strict -I $(BOOSTDIR)
#CXXFLAGS = $(OPTFLAGS) -I $(BOOSTDIR)
#CXXFLAGS = $(OPTFLAGS) -I $(BOOSTDIR) --static

all : test_lint

test_lint : test_lint.o
	$(CXX) -o $@ test_lint.o $(CXXFLAGS) 

%.o : %.cpp
	$(CXX) -c $< $(CXXFLAGS)

clean:
	rm -rf *~ *.o ti_files/ *.s
