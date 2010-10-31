CXX=KCC
CXXOPTFLAGS=+K0 -g
CXXFLAGS=$(CXXOPTFLAGS) --strict
LDFLAGS=

all : testpolynom

testpolynom : testpolynom.o polynomial.hpp
	$(CXX) -o $@ testpolynom.o $(LDFLAGS) $(CXXFLAGS)

%.o : %.C
	$(CXX) -c $< $(CXXFLAGS)

clean:
	rm -rf *~ *.o ti_files/ testpolynom
