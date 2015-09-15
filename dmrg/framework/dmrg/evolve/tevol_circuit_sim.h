/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef APP_DMRG_TEVOL_CIRCUIT_SIM_H
#define APP_DMRG_TEVOL_CIRCUIT_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "dmrg/evolve/te_utils.hpp"
#include "dmrg/utils/results_collector.h"

#include "dmrg/mp_tensors/state_mps.h"

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup>
restrict_with_op(MPSTensor<Matrix, SymmGroup> const & mps,
                 typename operator_selector<Matrix, SymmGroup>::type const & op,
                 Index<SymmGroup> const & out_phys_i)
{
    typedef typename SymmGroup::charge charge;
    
    mps.make_left_paired();
    block_matrix<Matrix, SymmGroup> const & vec = mps.data();
    
    Index<SymmGroup> const& phys_i  = mps.site_dim();
    Index<SymmGroup> const& left_i  = mps.row_dim();
    Index<SymmGroup> const& right_i = mps.col_dim();
    
    ProductBasis<SymmGroup> left_pb(phys_i, left_i);
    ProductBasis<SymmGroup> new_left_pb(out_phys_i, left_i);
    
    block_matrix<Matrix, SymmGroup> ret_vec;
    Index<SymmGroup> out_left_i;
    
    for (size_t ls = 0; ls < left_i.size(); ++ls)
        for (size_t rs = 0; rs < right_i.size(); ++rs)
            for (size_t s1 = 0; s1 < phys_i.size(); ++s1)
                for (size_t s2 = 0; s2 < out_phys_i.size(); ++s2) {
                    
                    charge lc  = left_i[ls].first, rc  = right_i[rs].first;
                    charge s1c = phys_i[s1].first, s2c = out_phys_i[s2].first;
                    
                    if (! op.has_block(s1c, s2c) )
                        continue;
                    
                    charge phys_diff = SymmGroup::fuse(s2c, -s1c);
                    
                    charge left_vec_charge  = SymmGroup::fuse(s1c, lc);
                    charge right_vec_charge = rc;
                    
                    if (! vec.has_block(left_vec_charge, right_vec_charge) )
                        continue;
                    
                    if (! out_left_i.has(lc))
                        out_left_i.insert(left_i[ls]);
                    
                    charge left_out_charge  = SymmGroup::fuse(s2c, lc);
                    charge right_out_charge = SymmGroup::fuse(rc, phys_diff);
                    
                    if (! ret_vec.has_block(left_out_charge, right_out_charge) )
                        ret_vec.insert_block(new Matrix(left_pb.size(s2c, lc), right_i[rs].second, 0), left_out_charge, right_out_charge);
                    
                    Matrix & oblock         = ret_vec(left_out_charge, right_out_charge);
                    Matrix const & iblock   = vec(left_vec_charge, right_vec_charge);
                    Matrix const & op_block = op(s1c, s2c);
                    
                    size_t i_l_offset = left_pb(s1c, lc);
                    size_t i_r_offset = 0;
                    size_t l_offset = new_left_pb(s2c, lc);
                    size_t r_offset = 0;
                    size_t i_op_offset = 0;
                    size_t op_offset = 0;
                    
#ifdef USE_AMBIENT
                    printf("UNOPTIMIZED FUNCTION (MULTIPLY WITH OP!)\n");
#endif
                    for (size_t ss1 = 0; ss1 < phys_i[s1].second; ++ss1) {
                        for (size_t ss2 = 0; ss2 < out_phys_i[s2].second; ++ss2) {
                            size_t o_l_start = l_offset + ss2*left_i[ls].second;
                            size_t o_r_start = r_offset;
                            
                            size_t i_l_start = i_l_offset + ss1*left_i[ls].second;
                            size_t i_r_start = i_r_offset;
                            
                            typename Matrix::value_type const & op_val = op_block(i_op_offset + ss1,  op_offset + ss2);
                            
                            // TODO: replace by kernel
                            for (size_t rr = 0; rr < right_i[rs].second; ++rr) {
                                for (size_t ll = 0; ll < left_i[ls].second; ++ll) {
                                    oblock(o_l_start+ll, o_r_start+rr) += iblock(i_l_start+ll, i_r_start+rr) * op_val;
                                }
                            }
                        }
                    }
                }
    
    
    MPSTensor<Matrix, SymmGroup> ret(out_phys_i, out_left_i, ret_vec.right_basis(), ret_vec, LeftPaired);
    assert( ret.reasonable() );
    return ret;
}

template<class T>
std::vector<T> vectorizer(std::string const & val)
{
    std::string raw = val;
    boost::trim_if(raw, boost::is_any_of("\"'"));
    std::vector<T> ret;
    
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(" ");
    tokenizer tokens(raw, sep);
    BOOST_FOREACH(std::string t, tokens) {
        ret.push_back(boost::lexical_cast<T, std::string>(t));
    }
    return ret;
}

// ******   SIMULATION CLASS   ******
template <class Matrix, class SymmGroup>
class circuit_evolver {
public:
    circuit_evolver(DmrgParameters * parms_, MPS<Matrix, SymmGroup> * mps_,
                    Lattice const& lattice_, Model<Matrix, SymmGroup> const& model_,
                    int init_sweep=0)
    {
        throw std::runtime_error("Circuit evolver only works with U(1) symmetry right now.");
    }
    
    void operator()(unsigned sweep, unsigned nsteps)
    { }
    
    results_collector const& iteration_results() const
    {
        return r;
    }
    
    void prepare_te_terms(unsigned sweep) { }
    
private:
    results_collector r;
};

template<class Matrix>
class circuit_evolver<Matrix, U1> {
    typedef U1 SymmGroup;
    typedef typename operator_selector<Matrix, SymmGroup>::type op_t;
public:
    circuit_evolver(DmrgParameters * parms_, MPS<Matrix, SymmGroup> * mps_,
                    Lattice const& lattice_, Model<Matrix, SymmGroup> const& model_,
                    int init_sweep=0)
    : parms(parms_)
    , mps(mps_)
    , lattice(lattice_) // shallow copy
    , model(model_) // shallow copy
    , L(lattice.size())
    {
        maquis::cout << "Circuit evolution." << std::endl;
        
        std::string input_filename = (*parms)["gatefile"];
        std::ifstream gatefile(input_filename.c_str());
        std::string line;
        
        std::getline(gatefile, line);
        std::vector<int> iconfiguration = vectorizer<int>(line);
        configuration.resize(L);
        for (int i = 0; i < L; ++i) {
            if ((*parms)["MODEL"] == "spinless fermions") {
                if (iconfiguration[i] == 0)
                    configuration[i] = boost::make_tuple(0, 0);
                else if (iconfiguration[i] == 1)
                    configuration[i] = boost::make_tuple(2, 0);
            } else if ((*parms)["MODEL"] == "alternative fermion Hubbard") {
                if (iconfiguration[i] == 0)
                    configuration[i] = boost::make_tuple(0, 0);
                else if (iconfiguration[i] == 1)
                    configuration[i] = boost::make_tuple(0, 1);
            } else
                throw std::runtime_error("Unknown MODEL.");
        }
        
        while (std::getline(gatefile, line))
        {
            std::istringstream iss(line);
            if (static_cast<bool>((*parms)["COMPLEX"])) {
                int pos; double rz, ry;
                if (!(iss >> pos >> rz >> ry)) { break; }
                gates.push_back( std::make_pair(pos, make_cmplx_gate(rz, ry, typename Matrix::value_type())) );
            } else {
                int pos; double phase;
                if (!(iss >> pos >> phase)) { break; }
                gates.push_back( std::make_pair(pos, make_real_gate(phase)) );
            }
        }
        std::reverse(gates.begin(), gates.end());
    }
    
    void operator()(unsigned sweep, unsigned nsteps)
    {
        iteration_results_.clear();
        
        *mps = state_mps<Matrix, SymmGroup>(configuration,
                                            std::vector<Index<SymmGroup> >(L, mps->site_dim(0)),
                                            std::vector<int>(L, 0));
        
        std::size_t Mmax=(*parms)["max_bond_dimension"];
        double cutoff=(*parms)["truncation_final"];
        MPS<Matrix, SymmGroup> const& constmps = *mps;
        
        for (typename std::vector< std::pair<std::size_t, op_t> >::const_iterator it = gates.begin();
             it != gates.end(); ++it)
        {
            std::cout << "Applying gate on sites " << it->first << "," << it->first+1 << std::endl;
            
            mps->canonize(it->first+1);
            constmps[it->first].make_left_paired();
            constmps[it->first+1].make_right_paired();
            
            block_matrix<Matrix, SymmGroup> v0, v1;
            gemm(constmps[it->first].data(), constmps[it->first+1].data(), v0); // outer product of two sites

            v1 = contraction::multiply_with_twosite<Matrix>(v0, it->second,
                                                            constmps[it->first].row_dim(), constmps[it->first+1].col_dim(),
                                                            constmps[it->first].site_dim());
            truncation_results trunc = compression::replace_two_sites_r2l(*mps, Mmax, cutoff, v1, it->first);
            
            // the above leaves the state canonized one to the left, we just need to tell it
            mps->canonization(true);
        }
        
        if (static_cast<bool>((*parms)["project"]))
        {
            std::cout << "Projecting." << std::endl;
            
            MPS<Matrix, SymmGroup> new_mps = *mps;
            
            op_t P;
            Matrix id1(1,1); id1(0,0) = 1;
            P.insert_block(id1, -1, -1);
            P.insert_block(id1, 1, 1);
            
            Index<SymmGroup> new_phys_i = P.right_basis();
            
            for (int p = 0; p < L; ++p)
                new_mps[p] = restrict_with_op((*mps)[p], P, new_phys_i);
            
            std::string projected_path = (*parms)["projfile"];
            save(projected_path, new_mps);
            std::map<std::string, int> status;
            status["sweep"] = 0;
            storage::archive ar(projected_path+"/props.h5", "w");
            ar["/status"] << status;
            ar["/parameters"] << *parms;
        }
    }
    
    results_collector const& iteration_results() const
    {
        return iteration_results_;
    }
    
    void prepare_te_terms(unsigned sweep) { }
    
private:
    op_t make_real_gate(double phase)
    {
        if ((*parms)["MODEL"] == "spinless fermions") {
            op_t gate;
            
            Matrix id1(1,1); id1(0,0) = 1;
            gate.insert_block(id1, 0, 0);
            gate.insert_block(id1, 4, 4);
            
            Matrix rot(2,2);
            rot(0,0) = cos(phase);
            rot(1,0) = -sin(phase);
            rot(0,1) = sin(phase);
            rot(1,1) = cos(phase);
            gate.insert_block(rot, 2, 2);
            
            return gate;
        } else if ((*parms)["MODEL"] == "alternative fermion Hubbard") {
            op_t gate0, gate1, gate;
            
            // 1 4 6 4 1
            Matrix d1(1,1);
            gate0.insert_block(d1, -2, -2);
            gate0.insert_block(d1, 2, 2);
            
            Matrix d4(4,4);
            gate0.insert_block(d4, -1, -1);
            gate0.insert_block(d4, 1, 1);
            
            Matrix d6(6,6);
            gate0.insert_block(d6, 0, 0);
            
            gate1 = gate0;
            gate = gate0;
            
            double theta = phase;
            
            // auto-generated
            gate0( std::make_pair(0,1), std::make_pair(0,1) ) = 1;
            gate0( std::make_pair(-1,0), std::make_pair(-1,0) ) = 1;
            gate0( std::make_pair(1,2), std::make_pair(1,2) ) = cos(theta);
            gate0( std::make_pair(1,2), std::make_pair(1,0) ) = -sin(theta);
            gate0( std::make_pair(0,2), std::make_pair(0,2) ) = cos(theta);
            gate0( std::make_pair(0,2), std::make_pair(0,0) ) = -sin(theta);
            gate0( std::make_pair(-1,2), std::make_pair(-1,2) ) = 1;
            gate0( std::make_pair(-2,0), std::make_pair(-2,0) ) = 1;
            gate0( std::make_pair(0,5), std::make_pair(0,5) ) = cos(theta);
            gate0( std::make_pair(0,5), std::make_pair(0,3) ) = sin(theta);
            gate0( std::make_pair(-1,3), std::make_pair(-1,3) ) = cos(theta);
            gate0( std::make_pair(-1,3), std::make_pair(-1,1) ) = sin(theta);
            gate0( std::make_pair(1,0), std::make_pair(1,2) ) = sin(theta);
            gate0( std::make_pair(1,0), std::make_pair(1,0) ) = cos(theta);
            gate0( std::make_pair(0,0), std::make_pair(0,2) ) = sin(theta);
            gate0( std::make_pair(0,0), std::make_pair(0,0) ) = cos(theta);
            gate0( std::make_pair(2,0), std::make_pair(2,0) ) = 1;
            gate0( std::make_pair(1,1), std::make_pair(1,1) ) = 1;
            gate0( std::make_pair(0,3), std::make_pair(0,5) ) = -sin(theta);
            gate0( std::make_pair(0,3), std::make_pair(0,3) ) = cos(theta);
            gate0( std::make_pair(-1,1), std::make_pair(-1,3) ) = -sin(theta);
            gate0( std::make_pair(-1,1), std::make_pair(-1,1) ) = cos(theta);
            gate0( std::make_pair(1,3), std::make_pair(1,3) ) = 1;
            gate0( std::make_pair(0,4), std::make_pair(0,4) ) = 1;
            
            gate1( std::make_pair(0,1), std::make_pair(0,1) ) = 1;
            gate1( std::make_pair(-1,0), std::make_pair(-1,0) ) = cos(theta);
            gate1( std::make_pair(-1,0), std::make_pair(-1,2) ) = -sin(theta);
            gate1( std::make_pair(1,2), std::make_pair(1,2) ) = 1;
            gate1( std::make_pair(0,2), std::make_pair(0,2) ) = cos(theta);
            gate1( std::make_pair(0,2), std::make_pair(0,5) ) = sin(theta);
            gate1( std::make_pair(-1,2), std::make_pair(-1,0) ) = sin(theta);
            gate1( std::make_pair(-1,2), std::make_pair(-1,2) ) = cos(theta);
            gate1( std::make_pair(-2,0), std::make_pair(-2,0) ) = 1;
            gate1( std::make_pair(0,5), std::make_pair(0,2) ) = -sin(theta);
            gate1( std::make_pair(0,5), std::make_pair(0,5) ) = cos(theta);
            gate1( std::make_pair(-1,3), std::make_pair(-1,3) ) = 1;
            gate1( std::make_pair(1,0), std::make_pair(1,0) ) = 1;
            gate1( std::make_pair(0,0), std::make_pair(0,0) ) = cos(theta);
            gate1( std::make_pair(0,0), std::make_pair(0,3) ) = -sin(theta);
            gate1( std::make_pair(2,0), std::make_pair(2,0) ) = 1;
            gate1( std::make_pair(1,1), std::make_pair(1,1) ) = cos(theta);
            gate1( std::make_pair(1,1), std::make_pair(1,3) ) = sin(theta);
            gate1( std::make_pair(0,3), std::make_pair(0,0) ) = sin(theta);
            gate1( std::make_pair(0,3), std::make_pair(0,3) ) = cos(theta);
            gate1( std::make_pair(-1,1), std::make_pair(-1,1) ) = 1;
            gate1( std::make_pair(1,3), std::make_pair(1,1) ) = -sin(theta);
            gate1( std::make_pair(1,3), std::make_pair(1,3) ) = cos(theta);
            gate1( std::make_pair(0,4), std::make_pair(0,4) ) = 1;

            gemm(gate0, gate1, gate);
            
            return gate;
        } else
            throw std::runtime_error("Unknown MODEL.");
    }

    op_t make_cmplx_gate(double rz, double ry, double) {
        throw std::runtime_error("Complex rotation on real matrices?");
    }
    op_t make_cmplx_gate(double rz, double ry, std::complex<double>)
    {
        if ((*parms)["MODEL"] == "alternative fermion Hubbard") {
            op_t gate0, gate1, gate;
            
            // 1 4 6 4 1
            Matrix d1(1,1);
            gate0.insert_block(d1, -2, -2);
            gate0.insert_block(d1, 2, 2);
            
            Matrix d4(4,4);
            gate0.insert_block(d4, -1, -1);
            gate0.insert_block(d4, 1, 1);
            
            Matrix d6(6,6);
            gate0.insert_block(d6, 0, 0);
            
            gate1 = gate0;
            gate = gate0;
            
            std::complex<double> II(0.0,1.0);
            
            gate0( std::make_pair(0,1), std::make_pair(0,1) ) = 1;
            gate0( std::make_pair(-1,0), std::make_pair(-1,0) ) = 1;
            gate0( std::make_pair(1,2), std::make_pair(1,2) ) = exp(II*rz)*cos(ry);
            gate0( std::make_pair(1,2), std::make_pair(1,0) ) = exp(-II*rz)*sin(ry);
            gate0( std::make_pair(0,2), std::make_pair(0,2) ) = exp(II*rz)*cos(ry);
            gate0( std::make_pair(0,2), std::make_pair(0,0) ) = exp(-II*rz)*sin(ry);
            gate0( std::make_pair(-1,2), std::make_pair(-1,2) ) = 1;
            gate0( std::make_pair(-2,0), std::make_pair(-2,0) ) = 1;
            gate0( std::make_pair(0,5), std::make_pair(0,5) ) = exp(II*rz)*cos(ry);
            gate0( std::make_pair(0,5), std::make_pair(0,3) ) = -exp(-II*rz)*sin(ry);
            gate0( std::make_pair(-1,3), std::make_pair(-1,3) ) = exp(II*rz)*cos(ry);
            gate0( std::make_pair(-1,3), std::make_pair(-1,1) ) = -exp(-II*rz)*sin(ry);
            gate0( std::make_pair(1,0), std::make_pair(1,2) ) = -exp(II*rz)*sin(ry);
            gate0( std::make_pair(1,0), std::make_pair(1,0) ) = exp(-II*rz)*cos(ry);
            gate0( std::make_pair(0,0), std::make_pair(0,2) ) = -exp(II*rz)*sin(ry);
            gate0( std::make_pair(0,0), std::make_pair(0,0) ) = exp(-II*rz)*cos(ry);
            gate0( std::make_pair(2,0), std::make_pair(2,0) ) = 1;
            gate0( std::make_pair(1,1), std::make_pair(1,1) ) = 1;
            gate0( std::make_pair(0,3), std::make_pair(0,5) ) = exp(II*rz)*sin(ry);
            gate0( std::make_pair(0,3), std::make_pair(0,3) ) = exp(-II*rz)*cos(ry);
            gate0( std::make_pair(-1,1), std::make_pair(-1,3) ) = exp(II*rz)*sin(ry);
            gate0( std::make_pair(-1,1), std::make_pair(-1,1) ) = exp(-II*rz)*cos(ry);
            gate0( std::make_pair(1,3), std::make_pair(1,3) ) = 1;
            gate0( std::make_pair(0,4), std::make_pair(0,4) ) = 1;
            
            gate1( std::make_pair(0,1), std::make_pair(0,1) ) = 1;
            gate1( std::make_pair(-1,0), std::make_pair(-1,0) ) = exp(II*rz)*cos(ry);
            gate1( std::make_pair(-1,0), std::make_pair(-1,2) ) = exp(-II*rz)*sin(ry);
            gate1( std::make_pair(1,2), std::make_pair(1,2) ) = 1;
            gate1( std::make_pair(0,2), std::make_pair(0,2) ) = exp(II*rz)*cos(ry);
            gate1( std::make_pair(0,2), std::make_pair(0,5) ) = -exp(-II*rz)*sin(ry);
            gate1( std::make_pair(-1,2), std::make_pair(-1,0) ) = -exp(II*rz)*sin(ry);
            gate1( std::make_pair(-1,2), std::make_pair(-1,2) ) = exp(-II*rz)*cos(ry);
            gate1( std::make_pair(-2,0), std::make_pair(-2,0) ) = 1;
            gate1( std::make_pair(0,5), std::make_pair(0,2) ) = exp(II*rz)*sin(ry);
            gate1( std::make_pair(0,5), std::make_pair(0,5) ) = exp(-II*rz)*cos(ry);
            gate1( std::make_pair(-1,3), std::make_pair(-1,3) ) = 1;
            gate1( std::make_pair(1,0), std::make_pair(1,0) ) = 1;
            gate1( std::make_pair(0,0), std::make_pair(0,0) ) = exp(II*rz)*cos(ry);
            gate1( std::make_pair(0,0), std::make_pair(0,3) ) = exp(-II*rz)*sin(ry);
            gate1( std::make_pair(2,0), std::make_pair(2,0) ) = 1;
            gate1( std::make_pair(1,1), std::make_pair(1,1) ) = exp(II*rz)*cos(ry);
            gate1( std::make_pair(1,1), std::make_pair(1,3) ) = -exp(-II*rz)*sin(ry);
            gate1( std::make_pair(0,3), std::make_pair(0,0) ) = -exp(II*rz)*sin(ry);
            gate1( std::make_pair(0,3), std::make_pair(0,3) ) = exp(-II*rz)*cos(ry);
            gate1( std::make_pair(-1,1), std::make_pair(-1,1) ) = 1;
            gate1( std::make_pair(1,3), std::make_pair(1,1) ) = exp(II*rz)*sin(ry);
            gate1( std::make_pair(1,3), std::make_pair(1,3) ) = exp(-II*rz)*cos(ry);
            gate1( std::make_pair(0,4), std::make_pair(0,4) ) = 1;
            
            gemm(gate0, gate1, gate);
            
            return gate;
        } else
            throw std::runtime_error("Unknown MODEL.");
    }
    
    DmrgParameters * parms;
    MPS<Matrix, SymmGroup> * mps;
    Lattice lattice;
    Model<Matrix, SymmGroup> model;
    size_t L;
    
    results_collector iteration_results_;
    
    std::vector<boost::tuple<typename SymmGroup::charge, size_t> > configuration;
    std::vector< std::pair<std::size_t, op_t> > gates;
};

#endif
