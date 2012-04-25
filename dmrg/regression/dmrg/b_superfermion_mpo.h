#ifndef SUPERFERMION_MPO_H
#define SUPERFERMION_MPO_H

#include "b_adjacency.h"
#include "b_generate_mpo.h"
#include "b_DmrgParameters.h"
#include "dmrg/mp_tensors/mpo.h"

namespace b_mpos
{
    template<class Matrix>
    class AlternateSuperfermions : public Hamiltonian<Matrix, U1>
    {
        typedef typename Hamiltonian<Matrix, U1>::op_t op_t;
        typedef std::pair<op_t, op_t> op_pair;
        typedef std::vector<op_pair> op_pairs;
        
    public:
        AlternateSuperfermions(BaseParameters & parms_, int sweep_)
        : parms(parms_)
        {
            Matrix m(2, 2);
            
            m(1, 0) = 1;
            m(0, 1) = 1;
            jump.insert_block(m, 1, 1);
            
            resize(m, 2, 1);
            
            m *= 0;
            m(0, 0) = 1;
            destroy_up.insert_block(m, 1, 0);
            
            m *= 0;
            m(1, 0) = 1;
            destroy_down.insert_block(m, 1, 0);
            
            resize(m, 1, 2);
            
            m *= 0;
            m(0, 0) = 1;
            create_up.insert_block(m, 0, 1);
            
            m *= 0;
            m(0, 1) = 1;
            create_down.insert_block(m, 0, 1);
            
            resize(m, 2, 2);
            
            m *= 0;
            m(0, 0) = 1;
            count_up.insert_block(m, 1, 1);
            
            m *= 0;
            m(1, 1) = 1;
            count_down.insert_block(m, 1, 1);
            
            m(0, 0) = 1;
            ident.insert_block(m, 1, 1);

            ident.insert_block(Matrix(1, 1, 1), 0, 0);
            
            omcount.insert_block(Matrix(1, 1, 1), 0, 0);
            
            omcount_up = ident;
            omcount_up -= count_up;
            omcount_down = ident;
            omcount_down -= count_down;
        }
        
        op_t get_identity() { return ident; }
        op_t get_free() { return ident; } // haha
        
        int num_2site_ops() { return 0; }
        op_pair get_2site_op(int) { return op_pair(); }
        int num_1site_ops() { return 0; }
        op_t get_1site_op(int) { return op_t(); }
        
        Index<U1> get_phys()
        {
            Index<U1> phys;
            phys.insert(std::make_pair(0, 1));
            phys.insert(std::make_pair(1, 2));
            return phys;
        }
        
        void push_extra_terms(MPOMaker<Matrix, U1> & mm,
                              b_adj::Adjacency & adj)
        {
            double penalty = parms.get<double>("penalty");
            
            for (int p = 0; p < adj.size(); ++p)
            {   
                vector<pair<size_t, op_t> > term;
                
                // penalty term
                vector<int> neighs = adj.forward(p);
                for (vector<int>::iterator it = neighs.begin();
                     it != neighs.end(); ++it)
                {
                    term.clear();
                    term.push_back(make_pair(p, penalty*count_up));
                    term.push_back(make_pair(*it, count_up));
                    mm.add_term(term);
                    
                    term.clear();
                    term.push_back(make_pair(p, penalty*count_down));
                    term.push_back(make_pair(*it, count_down));
                    mm.add_term(term);
                }
                
                neighs = adj.all(p);
                
                // vertical hopping terms
                term.clear();
                term.push_back(make_pair(p, jump));
                for (vector<int>::iterator n = neighs.begin();
                     n != neighs.end(); ++n)
                    term.push_back(make_pair(*n, omcount));
                mm.add_term(term);
                term.clear();
                
                // horiz hopping terms
                for (vector<int>::iterator hopto = neighs.begin();
                     hopto != neighs.end(); ++hopto)
                {
                    term.clear();
                    
                    term.push_back(make_pair(p, destroy_up));
                    term.push_back(make_pair(*hopto, create_up));
                    
                    set<size_t> projected;
                    
                    vector<int> projsites = adj.all(*hopto);
                    for (vector<int>::iterator it = neighs.begin();
                         it != neighs.end(); ++it)
                        if ((*it != p) && (*it != *hopto)) {
                            term.push_back(make_pair(*it, omcount_up));
                            projected.insert(*it);
                        }
                    for (vector<int>::iterator it = projsites.begin();
                         it != projsites.end(); ++it)
                        if ((*it != p) && (*it != *hopto)) {
                            term.push_back(make_pair(*it, omcount_up));
                            projected.insert(*it);
                        }
                    
                    mm.add_term(term);
                    term.clear();
                    projected.clear();
                    
                    term.push_back(make_pair(p, destroy_down));
                    term.push_back(make_pair(*hopto, create_down));
                    
                    for (vector<int>::iterator it = neighs.begin();
                         it != neighs.end(); ++it)
                        if ((*it != p) && (*it != *hopto)) {
                            term.push_back(make_pair(*it, omcount_down));
                            projected.insert(*it);
                        }
                    for (vector<int>::iterator it = projsites.begin();
                         it != projsites.end(); ++it)
                        if ((*it != p) && (*it != *hopto)) {
                            term.push_back(make_pair(*it, omcount_down));
                            projected.insert(*it);
                        }
                    
                    mm.add_term(term);
                }
                
                // potential
                term.clear();
                term.push_back(make_pair(p, omcount_down));
                for (vector<int>::iterator it = neighs.begin();
                     it != neighs.end(); ++it)
                    term.push_back(make_pair(*it, omcount_up));
                mm.add_term(term);
                
                term.clear();
                term.push_back(make_pair(p, omcount_up));
                for (vector<int>::iterator it = neighs.begin();
                     it != neighs.end(); ++it)
                    term.push_back(make_pair(*it, omcount_down));
                mm.add_term(term);
            }
        }
        
    private:
        block_matrix<Matrix, U1> jump, destroy_up, destroy_down, create_up, create_down;
        block_matrix<Matrix, U1> count_up, count_down, ident;
        block_matrix<Matrix, U1> omcount, omcount_up, omcount_down;
        
        BaseParameters & parms;
    };
    
    template<class Matrix>
    class Superfermions : public Hamiltonian<Matrix, U1>
    {
        typedef typename Hamiltonian<Matrix, U1>::op_t op_t;
        typedef std::pair<op_t, op_t> op_pair;
        typedef std::vector<op_pair> op_pairs;
        
    public:
        Superfermions(BaseParameters & parms_, int sweep_)
        : sweep(sweep_)
        , parms(parms_)
        {
            create.insert_block(Matrix(1, 1, 1), 0, 1);
            destroy.insert_block(Matrix(1, 1, 1), 1, 0);
            
            count.insert_block(Matrix(1, 1, 1), 1, 1);
            omcount.insert_block(Matrix(1, 1, 1), 0, 0);
            
            ident.insert_block(Matrix(1, 1, 1), 0, 0);
            ident.insert_block(Matrix(1, 1, 1), 1, 1);
            
            sign.insert_block(Matrix(1, 1, 1), 0, 0);
            sign.insert_block(Matrix(1, 1, -1), 1, 1);
        }
        
        op_t get_identity() { return ident; }
        op_t get_free() { return ident; } // haha
        
        int num_2site_ops() { return 0; }
        op_pair get_2site_op(int) { return op_pair(); }
        int num_1site_ops() { return 0; }
        op_t get_1site_op(int) { return op_t(); }
        
        void push_extra_terms(MPOMaker<Matrix, U1> & mm,
                              b_adj::Adjacency & adj)
        {
            double penalty = parms.get<double>("penalty");
            
            for (int p = 0; p < adj.size(); ++p)
            {   
                vector<pair<size_t, op_t> > term;
                
                // penalty term
                vector<int> neighs = adj.forward(p);
                for (vector<int>::iterator it = neighs.begin();
                     it != neighs.end(); ++it)
                {
                    term.clear();
                    term.push_back(make_pair(p, penalty*count));
                    term.push_back(make_pair(*it, count));
                    mm.add_term(term);
                }
                
                neighs = adj.all(p);
                
                // hopping
                for (vector<int>::iterator hopto = neighs.begin();
                     hopto != neighs.end(); ++hopto)
                {
                    double phase = 1;
                    if (parms.get<double>("twist") == 1 && adj.wraps_pbc(p, *hopto))
                        phase = -1;
                    
                    term.clear();
                    term.push_back(make_pair(p, phase*destroy));
                    term.push_back(make_pair(*hopto, create));
                    
                    set<size_t> projected;
                    
                    vector<int> projsites = adj.all(*hopto);
                    for (vector<int>::iterator it = neighs.begin();
                         it != neighs.end(); ++it)
                        if ((*it != p) && (*it != *hopto)) {
                            term.push_back(make_pair(*it, omcount));
                            projected.insert(*it);
                        }
                    for (vector<int>::iterator it = projsites.begin();
                         it != projsites.end(); ++it)
                        if ((*it != p) && (*it != *hopto)) {
                            term.push_back(make_pair(*it, omcount));
                            projected.insert(*it);
                        }
                    
                    for (size_t f0 = std::min(p, *hopto)+1;
                         f0 < std::max(p, *hopto);
                         ++f0)
                        if (projected.count(f0) == 0)
                            term.push_back(make_pair(f0, sign));
                    
                    mm.add_term(term);
                }
                
                // potential
                term.clear();
                for (vector<int>::iterator it = neighs.begin();
                     it != neighs.end(); ++it)
                    term.push_back(make_pair(*it, omcount));
                mm.add_term(term);
            }
        }
        
        Index<U1> get_phys()
        {
            Index<U1> phys;
            phys.insert(std::make_pair(0, 1));
            phys.insert(std::make_pair(1, 1));
            return phys;
        }
        
    private:
        int sweep;
        BaseParameters & parms;
        
        block_matrix<Matrix, U1> create, destroy, count, omcount, ident, sign;
    };
}
#endif
