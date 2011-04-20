#ifndef B_U1_MEASUREMENTS_H
#define B_U1_MEASUREMENTS_H

template<class Matrix>
struct measure_<Matrix, U1>
{
    void measure_blbq(MPS<Matrix, U1> & mps,
                      b_adj::Adjacency & adj,
                      b_mpos::Hamiltonian<Matrix, U1> & H,
                      BaseParameters & model,
                      alps::hdf5::oarchive & ar)
    {
        std::vector<double> magns;
        
        std::vector<block_matrix<Matrix, U1> > ops(4);
        std::vector<std::string> names;
        
        block_matrix<Matrix, U1> ident = H.get_free();
        
        ops[0].insert_block(Matrix(1, 1, 1), 1, 1);
        ops[0].insert_block(Matrix(1, 1, 0), 0, 0);
        ops[0].insert_block(Matrix(1, 1, -1), -1, -1);
        names.push_back("Magnetization");
        
        ops[1].insert_block(Matrix(1, 1, 1), -1, -1);
        names.push_back("ColorDensity1");
        ops[2].insert_block(Matrix(1, 1, 1), 0, 0);
        names.push_back("ColorDensity2");
        ops[3].insert_block(Matrix(1, 1, 1), 1, 1);
        names.push_back("ColorDensity3");
        
        for (int i = 0; i < 4; ++i)
        {
            magns.clear();
            
            for (std::size_t p = 0; p < adj.size(); ++p)
            {
                b_mpos::MPOMaker<Matrix, U1> mpom(adj, H);
                std::vector<std::pair<std::size_t, block_matrix<Matrix, U1> > > v;
                v.push_back( std::make_pair( p, ops[i] ) );
                mpom.add_term(v);
                MPO<Matrix, U1> mpo = mpom.create_mpo();
                
                double val = expval(mps, mpo, 0);
                magns.push_back(val);
            }
            
            std::string n = std::string("/spectrum/results/") + names[i] + std::string("/mean/value");
            ar << alps::make_pvp(n, magns);
        }
        
        for (int i = 0; i < 4; ++i)
        {
            std::vector<block_matrix<Matrix, U1> > corr;
            corr.push_back( ops[i] );
            corr.push_back( ops[i] );
            
            std::string name = std::string("/spectrum/results/") + names[i] + std::string("Correlation");
            
            measure_2pt_correlation(mps, adj, ident, ident, ar,
                                    corr, name);
        }
        
        {
            block_matrix<Matrix, U1> ident, splus, sminus, sz, spp, smm, spm, smp, szz, szp, spz, szm, smz;
            
            ident.insert_block(Matrix(1, 1, 1), -1, -1);
            ident.insert_block(Matrix(1, 1, 1), 0, 0);
            ident.insert_block(Matrix(1, 1, 1), 1, 1);
            
            splus.insert_block(Matrix(1, 1, 1), -1, 0);
            splus.insert_block(Matrix(1, 1, 1), 0, 1);
            
            sminus.insert_block(Matrix(1, 1, 1), 1, 0);
            sminus.insert_block(Matrix(1, 1, 1), 0, -1);
            
            sz.insert_block(Matrix(1, 1, 1), 1, 1);
            sz.insert_block(Matrix(1, 1, 0), 0, 0);
            sz.insert_block(Matrix(1, 1, -1), -1, -1);
            
            gemm(splus, splus, spp);
            gemm(sminus, sminus, smm);
            gemm(splus, sminus, spm);
            gemm(sminus, splus, smp);
            gemm(sz, sz, szz);
            gemm(sz, splus, szp);
            gemm(splus, sz, spz);
            gemm(sz, sminus, szm);
            gemm(sminus, sz, smz);
            
            std::vector<double> bond_e;
            std::vector<std::string> bond_names;
            
            double Jbl = cos(M_PI * model.get<double>("theta"));
            double Jbq = sin(M_PI * model.get<double>("theta"));
            
            double sum = 0;
            
            for (int p = 0; p < adj.size(); ++p) {
                std::vector<int> fneighs = adj.forward(p);
                for (std::vector<int>::iterator it = fneighs.begin();
                     it != fneighs.end(); ++it)
                {
                    b_mpos::MPOMaker<Matrix, U1> mpom(adj, H);
                    std::vector<std::pair<std::size_t, block_matrix<Matrix, U1> > > v;
                    
#define term(a,b) { v.push_back( make_pair(p, a) ); v.push_back( make_pair(*it, b) ); mpom.add_term(v); v.clear(); }    
                    //                    term(ident, ident);
                    term(Jbl*splus, sminus);
                    term(Jbl*sminus, splus);
                    term(Jbl*sz, sz);
                    
                    term(Jbq*spp, smm);
                    term(Jbq*smm, spp);
                    term(Jbq*szz, szz);
                    
                    term(Jbq*spm, smp);
                    term(Jbq*smp, spm);
                    
                    term(Jbq*spz, smz);
                    term(Jbq*szp, szm);
                    
                    term(Jbq*smz, spz);
                    term(Jbq*szm, szp);            
#undef term
                    
                    MPO<Matrix, U1> mpo = mpom.create_mpo();
                    
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
    
    void measure_superf(MPS<Matrix, U1> & mps,
                        b_adj::Adjacency & adj,
                        b_mpos::Hamiltonian<Matrix, U1> & H,
                        BaseParameters & model,
                        alps::hdf5::oarchive & ar)
    {
        std::vector<double> magns;
        
        block_matrix<Matrix, U1> dens;
        
        dens.insert_block(Matrix(1, 1, 1), 1, 1);
        
        for (std::size_t p = 0; p < adj.size(); ++p)
        {
            b_mpos::MPOMaker<Matrix, U1> mpom(adj, H);
            std::vector<std::pair<std::size_t, block_matrix<Matrix, U1> > > v;
            v.push_back( std::make_pair( p, dens ) );
            mpom.add_term(v);
            MPO<Matrix, U1> mpo = mpom.create_mpo();
            
            double val = expval(mps, mpo, 0);
            magns.push_back(val);
        }
        
        ar << alps::make_pvp("/spectrum/results/Density/mean/value", magns);
        
        std::vector<double> corrs;
        
        for (std::size_t p = 0; p < adj.size(); ++p)
        {
            std::vector<int> neighs = adj.forward(p);
            for (std::vector<int>::iterator it = neighs.begin();
                 it != neighs.end(); ++it)
            {
                b_mpos::MPOMaker<Matrix, U1> mpom(adj, H);
                std::vector<std::pair<std::size_t, block_matrix<Matrix, U1> > > v;
                v.push_back( std::make_pair( p, dens ) );
                v.push_back( std::make_pair(*it, dens) );
                mpom.add_term(v);
                MPO<Matrix, U1> mpo = mpom.create_mpo();
                
                double val = expval(mps, mpo, 0);
                corrs.push_back(val);
            }
        }
        
        ar << alps::make_pvp("/spectrum/results/NNDensityCorrelation/mean/value", corrs);
    }
    
    void measure_ff(MPS<Matrix, U1> & mps,
                    b_adj::Adjacency & adj,
                    b_mpos::Hamiltonian<Matrix, U1> & H,
                    BaseParameters & model,
                    alps::hdf5::oarchive & ar)
    {
        block_matrix<Matrix, U1> dens, create, destroy, sign, ident;
        
        dens.insert_block(Matrix(1, 1, 1), 1, 1);
        create.insert_block(Matrix(1, 1, 1), 0, 1);
        destroy.insert_block(Matrix(1, 1, 1), 1, 0);
        
        sign.insert_block(Matrix(1, 1, 1), 0, 0);
        sign.insert_block(Matrix(1, 1, -1), 1, 1);
        
        ident.insert_block(Matrix(1, 1, 1), 0, 0);
        ident.insert_block(Matrix(1, 1, 1), 1, 1);
        
        std::vector<double> density;
        for (std::size_t p = 0; p < adj.size(); ++p)
        {
            b_mpos::MPOMaker<Matrix, U1> mpom(adj, H);
            std::vector<std::pair<std::size_t, block_matrix<Matrix, U1> > > v;
            v.push_back( std::make_pair( p, dens ) );
            mpom.add_term(v);
            MPO<Matrix, U1> mpo = mpom.create_mpo();
            
            double val = expval(mps, mpo, 0);
            density.push_back(val);
        }
        ar << alps::make_pvp("/spectrum/results/Density/mean/value", density);
        
        std::vector<block_matrix<Matrix, U1> > density_corr;
        density_corr.push_back( dens );
        density_corr.push_back( dens );
        measure_2pt_correlation(mps, adj, ident, ident,
                                ar, density_corr,
                                "/spectrum/results/DensityCorrelation");
        
        std::vector<block_matrix<Matrix, U1> > onebody;
        onebody.push_back( create );
        onebody.push_back( destroy );
        measure_2pt_correlation(mps, adj, ident, sign,
                                ar, onebody,
                                "/spectrum/results/OneBodyDM");
    }
    
    
    void operator()(MPS<Matrix, U1> & mps,
                    b_adj::Adjacency & adj,
                    b_mpos::Hamiltonian<Matrix, U1> & H,
                    BaseParameters & model,
                    alps::hdf5::oarchive & ar)
    {
        if (model.get<std::string>("model") == std::string("biquadratic"))
            measure_blbq(mps, adj, H, model, ar);
        else if (model.get<std::string>("model") == std::string("superfermion"))
            measure_superf(mps, adj, H, model, ar);
        else if (model.get<std::string>("model") == std::string("FreeFermions"))
            measure_ff(mps, adj, H, model, ar);
    }
};

#endif
