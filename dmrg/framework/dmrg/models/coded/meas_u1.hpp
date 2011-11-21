

template<class Matrix>
struct pre_measurements<Matrix, U1> {
	typedef Measurement_Term<Matrix, U1> mterm_t;
	typedef typename mterm_t::op_t op_t;
    
	void pre_ff (const Lattice& lattice, BaseParameters& model,
				 std::vector<Measurement_Term<Matrix, U1> >& terms,
				 op_t & ident)
	{
		op_t dens, create, destroy, sign;
        
        dens.insert_block(Matrix(1, 1, 1), 1, 1);
        create.insert_block(Matrix(1, 1, 1), 0, 1);
        destroy.insert_block(Matrix(1, 1, 1), 1, 0);
        
        sign.insert_block(Matrix(1, 1, 1), 0, 0);
        sign.insert_block(Matrix(1, 1, -1), 1, 1);
        
        ident.insert_block(Matrix(1, 1, 1), 0, 0);
        ident.insert_block(Matrix(1, 1, 1), 1, 1);
        
        {
        	mterm_t term;
        	term.name = "Density";
        	term.type = mterm_t::Local;
        	term.operators.push_back( std::make_pair(dens, false) );
            
        	terms.push_back(term);
        }
        {
        	mterm_t term;
        	term.name = "DensityCorrelation";
        	term.type = mterm_t::HalfCorrelation;
        	term.fill_operator = ident;
        	term.operators.push_back( std::make_pair(dens, false) );
        	term.operators.push_back( std::make_pair(dens, false) );
            
        	terms.push_back(term);
        }
        {
        	mterm_t term;
        	term.name = "OneBodyDM";
        	term.type = mterm_t::HalfCorrelation;
        	term.fill_operator = sign;
        	term.operators.push_back( std::make_pair(create, true) );
        	term.operators.push_back( std::make_pair(destroy, true) );
            
        	terms.push_back(term);
        }
        
	}
    
    
	void operator() (const Lattice& lattice, BaseParameters& model,
					 std::vector<Measurement_Term<Matrix, U1> >& terms,
					 typename Measurement_Term<Matrix, U1>::op_t& ident)
	{
        if (model.get<std::string>("MODEL") == std::string("FreeFermions"))
            pre_ff(lattice, model, terms, ident);
	}
};

