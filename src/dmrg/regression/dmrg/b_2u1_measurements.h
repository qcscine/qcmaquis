#ifndef B_2U1_MEASUREMENTS_H
#define B_2U1_MEASUREMENTS_H

template<class Matrix>
struct measure_<Matrix, TwoU1>
{
    void measure_blbq(MPS<Matrix, TwoU1> & mps,
                      b_adj::Adjacency & adj,
                      b_mpos::Hamiltonian<Matrix, TwoU1> & H,
                      BaseParameters & model,
                      alps::hdf5::archive & ar)
    {
        TwoU1::charge A, B, C;
        B[0] = 1;
        C[1] = 1;
        
        std::vector<double> magns;
        
        std::vector<block_matrix<Matrix, TwoU1> > ops(4);
        std::vector<std::string> names;
        
        block_matrix<Matrix, TwoU1> ident = H.get_free();
        
        ops[0].insert_block(Matrix(1, 1, 1), A, A);
        ops[0].insert_block(Matrix(1, 1, 0), B, B);
        ops[0].insert_block(Matrix(1, 1, -1),C, C);
        names.push_back("Magnetization");
        
        ops[1].insert_block(Matrix(1, 1, 1), A, A);
        names.push_back("ColorDensity1");
        ops[2].insert_block(Matrix(1, 1, 1), B, B);
        names.push_back("ColorDensity2");
        ops[3].insert_block(Matrix(1, 1, 1), C, C);
        names.push_back("ColorDensity3");
        
        for (int i = 0; i < 4; ++i)
        {
            magns.clear();
            
            for (std::size_t p = 0; p < adj.size(); ++p)
            {
                b_mpos::MPOMaker<Matrix, TwoU1> mpom(adj, H);
                std::vector<std::pair<std::size_t, block_matrix<Matrix, TwoU1> > > v;
                v.push_back( std::make_pair( p, ops[i] ) );
                mpom.add_term(v);
                MPO<Matrix, TwoU1> mpo = mpom.create_mpo();
                
                double val = expval(mps, mpo, 0);
                magns.push_back(val);
            }
            
            std::string n = std::string("/spectrum/results/") + names[i] + std::string("/mean/value");
            ar << alps::make_pvp(n, magns);
        }
        
        for (int i = 0; i < 4; ++i)
        {
            std::vector<block_matrix<Matrix, TwoU1> > corr;
            corr.push_back( ops[i] );
            corr.push_back( ops[i] );
            
            std::string name = std::string("/spectrum/results/") + names[i] + std::string("Correlation");
            
            measure_2pt_correlation(mps, adj, ident, ident, ar,
                                    corr, name);
        }
        
        {
            block_matrix<Matrix, TwoU1> ident, aa, bb, cc, ab, ba, ac, ca, bc, cb;
            
#define define_op(name, I, J) name.insert_block(Matrix(1, 1, 1), I, J)
            define_op(ident, A, A);
            define_op(ident, B, B);
            define_op(ident, C, C);
            
            define_op(aa, A, A);
            define_op(bb, B, B);
            define_op(cc, C, C);
            
            define_op(ab, A, B);
            define_op(ba, B, A);
            
            define_op(ac, A, C);
            define_op(ca, C, A);
            
            define_op(bc, B, C);
            define_op(cb, C, B);
#undef define_op
            
            std::vector<double> bond_e;
            std::vector<std::string> bond_names;
            
            double sum = 0;
            
            for (int p = 0; p < adj.size(); ++p) {
                std::vector<int> fneighs = adj.forward(p);
                for (std::vector<int>::iterator it = fneighs.begin();
                     it != fneighs.end(); ++it)
                {
                    b_mpos::MPOMaker<Matrix, TwoU1> mpom(adj, H);
                    std::vector<std::pair<std::size_t, block_matrix<Matrix, TwoU1> > > v;
#define term(a,b) { v.push_back(make_pair(p, a)); v.push_back(make_pair(*it, b)); mpom.add_term(v); v.clear(); }
                    
                    term(aa, aa);
                    term(bb, bb);
                    term(cc, cc);
                    
                    term(ab, ba);
                    term(ac, ca);
                    term(bc, cb);
                    
                    term(ba, ab);
                    term(ca, ac);
                    term(cb, bc);
#undef term
                    
                    MPO<Matrix, TwoU1> mpo = mpom.create_mpo();
                    
                    double val = expval(mps, mpo, 0);
                    sum += val;
                    
                    std::ostringstream name;
                    name << "( " << p << " ) -- ( " << *it << " )";
                    bond_names.push_back(name.str());
                    bond_e.push_back(val);
                }
            }
            
            cout << "Bond energy sum: " << sum << endl;
            
            ar << alps::make_pvp("/spectrum/results/BondEnergies/mean/value", bond_e);
            ar << alps::make_pvp("/spectrum/results/BondEnergies/labels", bond_names);
        }
    }
    
    void measure_fermi_hubbard(MPS<Matrix, TwoU1> & mps,
                               b_adj::Adjacency & adj,
                               b_mpos::Hamiltonian<Matrix, TwoU1> & H,
                               BaseParameters & model,
                               alps::hdf5::archive & ar)
    {
        block_matrix<Matrix, TwoU1> create_up, create_down;
        block_matrix<Matrix, TwoU1> destroy_up, destroy_down;
        block_matrix<Matrix, TwoU1> count_up, count_down;            
        block_matrix<Matrix, TwoU1> sign_up, sign_down;            
        block_matrix<Matrix, TwoU1> ident, sign;
        
        TwoU1::charge A(0), B(0), C(0), D(1);
        B[0]=1; C[1]=1;
        
        ident.insert_block(Matrix(1, 1, 1), A, A);
        ident.insert_block(Matrix(1, 1, 1), B, B);
        ident.insert_block(Matrix(1, 1, 1), C, C);
        ident.insert_block(Matrix(1, 1, 1), D, D);
        
        // sign for single type operator (eg. up'(i)*up(j))
        sign_up.insert_block(Matrix(1, 1, 1), A, A);
        sign_up.insert_block(Matrix(1, 1, -1), B, B);
        sign_up.insert_block(Matrix(1, 1, 1), C, C);
        sign_up.insert_block(Matrix(1, 1, -1), D, D);
        
        sign_down.insert_block(Matrix(1, 1, 1), A, A);
        sign_down.insert_block(Matrix(1, 1, 1), B, B);
        sign_down.insert_block(Matrix(1, 1, -1), C, C);
        sign_down.insert_block(Matrix(1, 1, -1), D, D);
 
        // sign for mixed operator (eg. up(i)*down(i)*up'(j)*down'(j))
        sign.insert_block(Matrix(1, 1, 1), A, A);
        sign.insert_block(Matrix(1, 1, -1), B, B);
        sign.insert_block(Matrix(1, 1, -1), C, C);
        sign.insert_block(Matrix(1, 1, 1), D, D);

        create_up.insert_block(Matrix(1, 1, 1), A, B);
        create_up.insert_block(Matrix(1, 1, 1), C, D);
        create_down.insert_block(Matrix(1, 1, 1), A, C);
        create_down.insert_block(Matrix(1, 1, 1), B, D);
        
        destroy_up.insert_block(Matrix(1, 1, 1), B, A);
        destroy_up.insert_block(Matrix(1, 1, 1), D, C);
        destroy_down.insert_block(Matrix(1, 1, 1), C, A);
        destroy_down.insert_block(Matrix(1, 1, 1), D, B);
        
        count_up.insert_block(Matrix(1, 1, 1), B, B);
        count_up.insert_block(Matrix(1, 1, 1), D, D);
        count_down.insert_block(Matrix(1, 1, 1), C, C);
        count_down.insert_block(Matrix(1, 1, 1), D, D);
                
        std::vector<double> density_up;
        std::vector<double> density_down;
        for (std::size_t p = 0; p < adj.size(); ++p)
        {
            double val;
            std::vector<std::pair<std::size_t, block_matrix<Matrix, TwoU1> > > v;
            
            // density up
            b_mpos::MPOMaker<Matrix, TwoU1> mpom_up(adj, H);
            v.clear();
            v.push_back( std::make_pair( p, count_up ) );
            mpom_up.add_term(v);
            MPO<Matrix, TwoU1> mpo_up = mpom_up.create_mpo();
            val = expval(mps, mpo_up, 0);
            density_up.push_back(val);
            
            // density down
            b_mpos::MPOMaker<Matrix, TwoU1> mpom_down(adj, H);
            v.clear();
            v.push_back( std::make_pair( p, count_down ) );
            mpom_down.add_term(v);
            MPO<Matrix, TwoU1> mpo_down = mpom_down.create_mpo();
            val = expval(mps, mpo_down, 0);
            density_down.push_back(val);
            
        }
        ar << alps::make_pvp("/spectrum/results/DensityUp/mean/value", density_up);
        ar << alps::make_pvp("/spectrum/results/DensityDown/mean/value", density_down);

        double total_count = 0.;
        total_count += std::accumulate(density_up.begin(), density_up.end(), 0.);
        total_count += std::accumulate(density_down.begin(), density_down.end(), 0.);
        ar << alps:: make_pvp("/spectrum/results/Filling/mean/value", total_count);
        
        
        cout << "Calculating z-z spin correlation." << std::endl;
        std::vector<block_matrix<Matrix, TwoU1> > zz_corr;
        zz_corr.push_back( count_up );
        zz_corr.push_back( count_up );
        measure_2pt_correlation(mps, adj, ident, ident,
                                ar, zz_corr,
                                "/spectrum/results/ZZSpinCorrelation1");
        zz_corr.clear();
        zz_corr.push_back( count_up );
        zz_corr.push_back( -1.*count_down );
        measure_2pt_correlation(mps, adj, ident, ident,
                                ar, zz_corr,
                                "/spectrum/results/ZZSpinCorrelation2");
        zz_corr.clear();
        zz_corr.push_back( -1.*count_down );
        zz_corr.push_back( count_up );
        measure_2pt_correlation(mps, adj, ident, ident,
                                ar, zz_corr,
                                "/spectrum/results/ZZSpinCorrelation3");
        zz_corr.clear();
        zz_corr.push_back( count_down );
        zz_corr.push_back( count_down );
        measure_2pt_correlation(mps, adj, ident, ident,
                                ar, zz_corr,
                                "/spectrum/results/ZZSpinCorrelation4");
        
        cout << "Calculating one body density matrix." << std::endl;
        std::vector<block_matrix<Matrix, TwoU1> > onebody;
        onebody.push_back( create_up );
        onebody.push_back( destroy_up );
        measure_2pt_correlation(mps, adj, ident, sign_up,
                                ar, onebody,
                                "/spectrum/results/OneBodyDMUp");
        onebody.clear();
        onebody.push_back( create_down );
        onebody.push_back( destroy_down );
        measure_2pt_correlation(mps, adj, ident, sign_down,
                                ar, onebody,
                                "/spectrum/results/OneBodyDMDown");
       
        
//        cout << "Calculating rung-rung pair field correlation." << std::endl;
//        std::vector<block_matrix<Matrix, TwoU1> > rr_pairfield;
//        rr_pairfield.push_back( destroy_down );
//        rr_pairfield.push_back( destroy_up );
//        rr_pairfield.push_back( create_up );
//        rr_pairfield.push_back( create_down );
//        measure_4ptnn_correlation(mps, adj, ident, sign,
//                                  ar, rr_pairfield,
//                                  "/spectrum/results/RungRungPairFieldCorrelation1");
//        rr_pairfield.clear();
//        rr_pairfield.push_back( -1.*destroy_down );
//        rr_pairfield.push_back( destroy_up );
//        rr_pairfield.push_back( create_down );
//        rr_pairfield.push_back( create_up );
//        measure_4ptnn_correlation(mps, adj, ident, sign,
//                                  ar, rr_pairfield,
//                                  "/spectrum/results/RungRungPairFieldCorrelation2");
//        rr_pairfield.clear();
//        rr_pairfield.push_back( -1.*destroy_up );
//        rr_pairfield.push_back( destroy_down );
//        rr_pairfield.push_back( create_up );
//        rr_pairfield.push_back( create_down );
//        measure_4ptnn_correlation(mps, adj, ident, sign,
//                                  ar, rr_pairfield,
//                                  "/spectrum/results/RungRungPairFieldCorrelation3");
//        rr_pairfield.clear();
//        rr_pairfield.push_back( destroy_up );
//        rr_pairfield.push_back( destroy_down );
//        rr_pairfield.push_back( create_down );
//        rr_pairfield.push_back( create_up );
//        measure_4ptnn_correlation(mps, adj, ident, sign,
//                                  ar, rr_pairfield,
//                                  "/spectrum/results/RungRungPairFieldCorrelation4");
//        rr_pairfield.clear();
//        rr_pairfield.push_back( create_up );
//        rr_pairfield.push_back( create_down );
//        rr_pairfield.push_back( destroy_down );
//        rr_pairfield.push_back( destroy_up );
//        measure_4ptnn_correlation(mps, adj, ident, sign,
//                                  ar, rr_pairfield,
//                                  "/spectrum/results/RungRungPairFieldCorrelation1_inv");
//        rr_pairfield.clear();
//        rr_pairfield.push_back( create_down );
//        rr_pairfield.push_back( create_up );
//        rr_pairfield.push_back( -1.*destroy_down );
//        rr_pairfield.push_back( destroy_up );
//        measure_4ptnn_correlation(mps, adj, ident, sign,
//                                  ar, rr_pairfield,
//                                  "/spectrum/results/RungRungPairFieldCorrelation2_inv");
//        rr_pairfield.clear();
//        rr_pairfield.push_back( create_up );
//        rr_pairfield.push_back( create_down );
//        rr_pairfield.push_back( -1.*destroy_up );
//        rr_pairfield.push_back( destroy_down );
//        measure_4ptnn_correlation(mps, adj, ident, sign,
//                                  ar, rr_pairfield,
//                                  "/spectrum/results/RungRungPairFieldCorrelation3_inv");
//        rr_pairfield.clear();
//        rr_pairfield.push_back( create_down );
//        rr_pairfield.push_back( create_up );
//        rr_pairfield.push_back( destroy_up );
//        rr_pairfield.push_back( destroy_down );
//        measure_4ptnn_correlation(mps, adj, ident, sign,
//                                  ar, rr_pairfield,
//                                  "/spectrum/results/RungRungPairFieldCorrelation4_inv");

    }

    
    void operator()(MPS<Matrix, TwoU1> & mps,
                    b_adj::Adjacency & adj,
                    b_mpos::Hamiltonian<Matrix, TwoU1> & H,
                    BaseParameters & model,
                    alps::hdf5::archive & ar)
    {
        if (model.get<std::string>("model") == std::string("biquadratic"))
            measure_blbq(mps, adj, H, model, ar);
        else if (model.get<std::string>("model") == std::string("fermi_hubbard"))
            measure_fermi_hubbard(mps, adj, H, model, ar);
    }
};

#endif
