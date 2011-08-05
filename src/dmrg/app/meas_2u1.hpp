

template<class Matrix>
struct pre_measurements<Matrix, TwoU1> {
    
	void operator() (const Lattice& lattice, BaseParameters& model,
					 std::vector<Measurement_Term<Matrix, TwoU1> >& terms,
					 typename Measurement_Term<Matrix, TwoU1>::op_t& ident)
	{
        
	}
};

