all: zero_exp high_t_exp

zero_exp: 
	g++ -O2 -Wall -pedantic zero_expansion_testcode.cpp -I./ -I../dmrg/ -lboost_graph -lgmpxx -lgmp -o zero_exp
high_t_exp:
	g++ -O2 -Wall -pedantic high_t_expansion_testcode.cpp -I./ -I../dmrg/ -lboost_graph -lgmpxx -lgmp -o high_t_exp

.PHONY: clean
clean:
	-rm zero_exp
	-rm high_t_exp
