/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Laboratory of Physical Chemistry, ETH Zurich
 *               2014-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2014-2014 by Erik Hedegaard <erik.hedegaard@phys.chem.ethz.ch>
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

#ifdef USE_AMBIENT
#include <mpi.h>
#endif
#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>
#include <string>
#include <sstream>
#include<boost/tokenizer.hpp>
#include<boost/lexical_cast.hpp>
#include <algorithm>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/sim/matrix_types.h"
#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/utils/storage.h"
#include "maquis_dmrg.h"

namespace entanglement_detail {

    std::vector<std::pair<int, int> > get_labels(const std::vector<std::string> & quant_label)
    {
       std::vector<std::pair<int, int> > labels;
       std::pair<int, int> lab;
          for(int j = 0; j < quant_label.size(); j++){
             boost::tokenizer<> tok(quant_label[j]);
             boost::tokenizer<>::iterator beg= tok.begin();
             lab.first = boost::lexical_cast<int>(*beg);
             ++beg;
             lab.second = boost::lexical_cast<int>(*beg);
             labels.push_back(lab);
          }
          return labels;
    }

    std::vector<int> get_labels_vec(const std::vector<std::string> & quant_label)
    {
       std::vector<int> labels;
          for(int j = 0; j < quant_label.size(); j++){
             boost::tokenizer<> tok(quant_label[j]);
             boost::tokenizer<>::iterator beg = tok.begin();
             labels.push_back(boost::lexical_cast<int>(*tok.begin()));
          }
          return labels;
    }

    template <class P>
    bool comp(const P l, const P r)
    {
        return l.first > r.first;
    }

    template <class Matrix>
    Matrix load_vector(storage::archive & ar, std::string path)
    {
        std::vector<std::string> x;
        Matrix y;

        ar[path + "/labels"] >> x;
        std::vector<int> labels = get_labels_vec(x);

        ar[path + "/mean/value"] >> y;
        Matrix ret(num_rows(y), 1);

        for (int i = 0; i < num_rows(y); ++i)
            ret(labels[i], 0) = y(i, 0);

        return ret;
    }

    template <class Matrix>
    Matrix load_vector(const maquis::results_map_type<typename Matrix::value_type>& map, std::string path)
    {
        const maquis::meas_with_results_type<typename Matrix::value_type> & meas = map.at(path);
        const std::vector<std::vector<int> > & labels = meas.first;
        const std::vector<typename Matrix::value_type> & values = meas.second;
        assert(labels.size() > 0);
        assert(labels[0].size() == 1); // Must have only 1 label per entry!

        int L = labels.size();
        Matrix ret(L, 1);

        for (int i = 0; i < L; ++i)
            ret(labels[i][0], 0) = values[i];

        return ret;
    }

    template <class Matrix>
    Matrix assy_hc(const std::vector<std::string> & label_str, Matrix const & values, const int L)
    {
        std::vector<std::pair<int,int> > labels;
        labels = get_labels(label_str);
        Matrix ret(L, L, 0);
        for(int i=0;i<labels.size();i++)
        {
            ret(labels[i].first, labels[i].second) = values(i,0);
            ret(labels[i].second, labels[i].first) = values(i,0);
        }
        return ret;
    }

    template <class Matrix>
    Matrix assy_hc(const maquis::meas_with_results_type<typename Matrix::value_type> & meas, int L)
    {
        const std::vector<std::vector<int> > & labels = meas.first;
        const std::vector<typename Matrix::value_type> & values = meas.second;
        assert(labels.size() > 0);
        assert(labels[0].size() == 2); // Must have only 2 labels per entry!
        Matrix ret(L, L, 0);
        for(int i=0;i<labels.size();i++)
        {
            ret(labels[i][0], labels[i][1]) = values[i];
            ret(labels[i][0], labels[i][1]) = values[i];
        }
        return ret;
    }

    template <class Matrix>
    Matrix load_matrix(storage::archive & ar, std::string path, int L)
    {
        std::vector<std::string> x;
        Matrix y;
        ar[path + "/labels"] >> x;
        ar[path + "/mean/value"] >> y;
        return assy_hc(x, y, L);
    }

    template <class Matrix>
    Matrix load_matrix(const maquis::results_map_type<typename Matrix::value_type> & map, std::string path, int L)
    {
        const maquis::meas_with_results_type<typename Matrix::value_type> & meas = map.at(path);
        return assy_hc<Matrix>(meas, L);
    }

    template <class Matrix>
    Matrix merge_transform(const std::vector<std::string> & lstr1, Matrix const & values1,
                           const std::vector<std::string> & lstr2, Matrix const & values2, const int L)
    {
        assert(lstr1.size() == lstr2.size());

        // output observable = input1 + transpose(input2)
        std::vector<std::pair<int,int> > labels1, labels2;
        labels1 = get_labels(lstr1);
        labels2 = get_labels(lstr2);

        Matrix ret(L, L, 0);

        for(int i = 0; i < labels1.size(); i++)
        {
            ret(labels1[i].first, labels1[i].second) = values1(i,0);
            ret(labels2[i].second, labels2[i].first) = values2(i,0);
        }
        return ret;
    }

    template <class Matrix>
    Matrix load_matrix_pair(storage::archive & ar, std::string path1, std::string path2, int L)
    {
        std::vector<std::string> x1, x2;
        Matrix y1, y2;

        ar[path1 + "/labels"] >> x1;
        ar[path1 + "/mean/value"] >> y1;

        ar[path2 + "/labels"] >> x2;
        ar[path2 + "/mean/value"] >> y2;

        return merge_transform(x1, y1, x2, y2, L);
    }

    template <class Matrix>
    Matrix load_matrix_pair(const maquis::results_map_type<typename Matrix::value_type> & map, std::string path1, std::string path2, int L)
    {
        const std::vector<std::vector<int> > & labels1 = map.at(path1).first;
        const std::vector<typename Matrix::value_type> & values1 = map.at(path1).second;
        const std::vector<std::vector<int> > & labels2 = map.at(path2).first;
        const std::vector<typename Matrix::value_type> & values2 = map.at(path2).second;

        assert(labels1.size() > 0);
        assert(labels1[0].size() == 2); // Must have only 2 labels per entry!
        assert(labels2.size() > 0);
        assert(labels2[0].size() == 2);
        assert(labels1.size() == labels2.size());

        // output observable = input1 + transpose(input2)
        Matrix ret(L, L, 0);

        for(int i = 0; i < labels1.size(); i++)
        {
            ret(labels1[i][0], labels1[i][1]) = values1[i];
            ret(labels2[i][1], labels2[i][0]) = values2[i];
        }

        return ret;
    }

    template <class Matrix>
    struct EntropyData
    {
        Matrix Nup;
        Matrix Ndown;
        Matrix Nupdown;
        Matrix dm_up;
        Matrix dm_down;
        Matrix nupnup;
        Matrix nupndown;
        Matrix ndownnup;
        Matrix ndownndown;
        Matrix doccdocc;
        Matrix transfer_up_while_down;
        Matrix transfer_down_while_up;
        Matrix transfer_up_while_down_at_2;
        Matrix transfer_down_while_up_at_2;
        Matrix transfer_up_while_down_at_1;
        Matrix transfer_down_while_up_at_1;
        Matrix transfer_pair;
        Matrix spinflip;
        Matrix nupdocc;
        Matrix ndowndocc;
        Matrix doccnup;
        Matrix doccndown;
        int L;
    };

#define LOADVEC(path) data.path = load_vector<Matrix>(ar, rpath + #path);
#define LOAD(path) data.path = load_matrix<Matrix>(ar, rpath + #path, L);
#define LOAD_PAIR(path1, path2) data.path1 = load_matrix_pair<Matrix>(ar, rpath + #path1, rpath + #path2, L);

// Loading of matrices is written in macros for both versions
#define LOAD_ALL()  { \
/* load quantities needed for s1 */ \
LOADVEC(Nup) \
LOADVEC(Ndown) \
LOADVEC(Nupdown) \
/* load quantities needed for s2 */ \
LOAD(dm_up) \
LOAD(dm_down) \
LOAD(nupnup) \
LOAD(nupndown) \
LOAD(ndownnup) \
LOAD(ndownndown) \
LOAD(doccdocc) \
LOAD(transfer_up_while_down) \
LOAD(transfer_down_while_up) \
LOAD(transfer_pair) \
LOAD(spinflip) \
LOAD_PAIR(transfer_up_while_down_at_2, transfer_up_while_down_at_1) \
LOAD_PAIR(transfer_up_while_down_at_1, transfer_up_while_down_at_2) \
LOAD_PAIR(transfer_down_while_up_at_2, transfer_down_while_up_at_1) \
LOAD_PAIR(transfer_down_while_up_at_1, transfer_down_while_up_at_2) \
LOAD_PAIR(nupdocc, doccnup) \
LOAD_PAIR(doccnup, nupdocc) \
LOAD_PAIR(ndowndocc, doccndown) \
LOAD_PAIR(doccndown, ndowndocc) \
}

    template<class Matrix>
    EntropyData<Matrix> loadData(storage::archive & ar)
    {
        typedef typename Matrix::value_type value_type;
        DmrgParameters parms;
        ar["/parameters"] >> parms;

        // chain-length
        const std::string rpath = "/spectrum/results/";

        // collect all in EntropyData => get as e.g. data.Nup etc.
        EntropyData<Matrix> data;
        data.L = parms["L"];
        int L = data.L;

        LOAD_ALL()

        return data;
    }

    template<class Matrix>
    EntropyData<Matrix> loadData(const maquis::results_map_type<typename Matrix::value_type> & ar)
    {
        const std::string rpath = "";

        // collect all in EntropyData => get as e.g. data.Nup etc.
        EntropyData<Matrix> data;
        int L = ar.at("Nup").second.size(); // TODO: make sure this exists!
        data.L = L;

        LOAD_ALL()

        return data;
    }
#undef LOADVEC
#undef LOAD
#undef LOAD_PAIR

    template <class Matrix>
    Matrix two_orb_rdm(int p, int q, EntropyData<Matrix> & data)
    {
            Matrix pq_dm_matrix(16,16);
            pq_dm_matrix( 0, 0) = 1 + data.Nupdown(p,0)  + data.Nupdown(q,0) + data.doccdocc(p,q) - data.Ndown(p,0)\
                                  - data.ndowndocc(p,q) - data.Ndown(q,0) - data.doccndown(p,q)\
                                  + data.ndownndown(p,q) - data.Nup(p,0) - data.nupdocc(p,q) + data.nupndown(p,q)\
                                  - data.Nup(q,0) - data.doccnup(p,q) + data.ndownnup(p,q) + data.nupnup(p,q);
            // O(6)/O(1)
            pq_dm_matrix( 1, 1) = - data.Nupdown(p,0) - data.doccdocc(p,q) + data.Ndown(p,0) + data.ndowndocc(p,q) \
                                  + data.doccndown(p,q) - data.ndownndown(p,q) + data.doccnup(p,q) - data.ndownnup(p,q);
            // O(11)O(1)
            pq_dm_matrix( 2, 2) = - data.Nupdown(p,0) - data.doccdocc(p,q) + data.doccndown(p,q) + data.Nup(p,0) \
                                  + data.nupdocc(p,q) - data.nupndown(p,q) + data.doccnup(p,q) - data.nupnup(p,q);
            // O(16)O(1)
            pq_dm_matrix( 3, 3) =  data.Nupdown(p,0)  - data.doccndown(p,q) - data.doccnup(p,q) + data.doccdocc(p,q);
            // O(1)O(6)
            pq_dm_matrix( 4, 4) = -data.Nupdown(q,0)  - data.doccdocc(p,q)  + data.ndowndocc(p,q) + data.Ndown(q,0) \
                                +  data.doccndown(p,q) - data.ndownndown(p,q) + data.nupdocc(p,q) - data.nupndown(p,q);
            // O(6)O(6)
            pq_dm_matrix( 5, 5) = data.ndownndown(p,q) - data.ndowndocc(p,q) - data.doccndown(p,q) + data.doccdocc(p,q);
            // O(11)O(6)
            pq_dm_matrix( 6, 6) = data.nupndown(p,q) - data.doccndown(p,q) - data.nupdocc(p,q) + data.doccdocc(p,q);
            // O(16)O(11)
            pq_dm_matrix( 7, 7) = data.doccndown(p,q) - data.doccdocc(p,q);
            // O(1)O(11)
            pq_dm_matrix( 8, 8) = -data.Nupdown(q,0) - data.doccdocc(p,q) + data.ndowndocc(p,q) + data.nupdocc(p,q)
                                  + data.Nup(q,0) + data.doccnup(p,q) - data.ndownnup(p,q) - data.nupnup(p,q);
            // O(6)O(11)
            pq_dm_matrix( 9, 9) = data.ndownnup(p,q) - data.ndowndocc(p,q) - data.doccnup(p,q) + data.doccdocc(p,q);
            // O(11)O(11)
            pq_dm_matrix(10,10) = data.nupnup(p,q) - data.nupdocc(p,q) - data.doccnup(p,q) + data.doccdocc(p,q);
            // O(16)O(11)
            pq_dm_matrix(11,11) = data.doccnup(p,q) - data.doccdocc(p,q);
            // O(1)O(16)
            pq_dm_matrix(12,12) = data.Nupdown(q,0) - data.nupdocc(p,q) - data.ndowndocc(p,q) + data.doccdocc(p,q);
            // O(1)O(16)
            pq_dm_matrix(13,13) = data.ndowndocc(p,q) - data.doccdocc(p,q);
            // O(11)O(16)
            pq_dm_matrix(14,14) = data.nupdocc(p,q) - data.doccdocc(p,q);
            // O(16)O(16)
            pq_dm_matrix(15,15) = data.doccdocc(p,q);
            // O(9)O(2)
            pq_dm_matrix(1,4)  = data.dm_down(p,q) - data.transfer_down_while_up_at_1(p,q) - data.transfer_down_while_up_at_2(p,q)
                               + data.transfer_down_while_up(p,q);
            pq_dm_matrix(4,1)  = pq_dm_matrix(1,4);
            // O(9)O(2)
            pq_dm_matrix(2,8)  =  data.dm_up(p,q) - data.transfer_up_while_down_at_1(p,q) - data.transfer_up_while_down_at_2(p,q)
                               +  data.transfer_up_while_down(p,q);
            pq_dm_matrix(8,2)  = pq_dm_matrix(2,8);
            // O(15)O(2)
            pq_dm_matrix(3,6) = data.transfer_down_while_up_at_1(p,q) - data.transfer_down_while_up(p,q);
            pq_dm_matrix(6,3) = pq_dm_matrix(3,6);
            // O(14)O(3)
            pq_dm_matrix(3,9) = -data.transfer_up_while_down_at_1(p,q) + data.transfer_up_while_down(p,q);
            pq_dm_matrix(9,3) = pq_dm_matrix(3,9);
            // O(10)O(7)
            pq_dm_matrix(6,9) = data.spinflip(p,q);
            pq_dm_matrix(9,6) = pq_dm_matrix(6,9);
            // O(13)O(4)
            pq_dm_matrix(3,12) = data.transfer_pair(p,q);
            pq_dm_matrix(12,3) = pq_dm_matrix(3,12);
            // O(9)O(8)
            pq_dm_matrix(6,12) = -data.transfer_up_while_down_at_2(p,q) + data.transfer_up_while_down(p,q);
            pq_dm_matrix(12,6) = pq_dm_matrix(6,12);
            // O(8)O(9)
            pq_dm_matrix(9,12) = data.transfer_down_while_up_at_2(p,q) - data.transfer_down_while_up(p,q);
            pq_dm_matrix(12,9) = pq_dm_matrix(9,12);
            // O(14)O(8)
            pq_dm_matrix(7,13) = data.transfer_up_while_down(p,q);
            pq_dm_matrix(13,7) = pq_dm_matrix(7,13);
            // O(15)O(12)
            pq_dm_matrix(11,14) = data.transfer_down_while_up(p,q);
            pq_dm_matrix(14,11) = pq_dm_matrix(11,14);

            return pq_dm_matrix;
    }
} // namespace entanglement_detail

template <class Matrix>
class EntanglementData
{
    typedef typename Matrix::value_type value_type;
    typedef typename maquis::traits::real_type<value_type>::type real_type;

public:
    EntanglementData(std::string rfile)
    {
        // Load results from rfile:
        storage::archive ar(rfile, "r");
        entanglement_detail::EntropyData<Matrix> data = entanglement_detail::loadData<Matrix>(ar);

        calculateData(data);
    }

    EntanglementData(const maquis::results_map_type<value_type> & meas_map)
    {
        entanglement_detail::EntropyData<Matrix> data = entanglement_detail::loadData<Matrix>(meas_map);
        calculateData(data);
    }


    Matrix s1()
    {
        return s1_;
    }

    Matrix s2()
    {
        return s2_;
    }

    Matrix I()
    {
        return I_;
    }

private:
    Matrix s1_, s2_, I_;

    void calculateData(entanglement_detail::EntropyData<Matrix>& data)
    {

        int L = data.L;

        //  calculate s1:
        s1_.resize(1,L);
        Matrix m11(L,1),m22(L,1),m33(L,1),m44(L,1);
        for (int i=0; i<L; ++i)
        {
            m11(i,0) = (data.Nup(i,0) - data.Nupdown(i,0)); // O(11)
            m22(i,0) = (data.Ndown(i,0) - data.Nupdown(i,0)); // O(6)
            m33(i,0) = (1 - data.Nup(i,0) - data.Ndown(i,0) + data.Nupdown(i,0)); //O(1)
            m44(i,0) = data.Nupdown(i,0); // O(16)
            s1_(0,i) = -(m11(i,0)*log( m11(i,0) ) + m22(i,0)*log( m22(i,0) )\
                      + m33(i,0)*log( m33(i,0) ) + m44(i,0)*log( m44(i,0)) );

        }

        // calculate s2
        s2_.resize(L,L);
        for (int p=0; p<L; ++p)
        {
            for (int q=(p+1); q<L; ++q)
            {
                Matrix rdm = entanglement_detail::two_orb_rdm(p,q,data);
                assert (num_rows(rdm) == 16 && num_cols(rdm) == 16);

                Matrix evecs(16,16);
                std::vector<real_type> eval(16);
                alps::numeric::syev(rdm,evecs,eval);
                for (int j=0; j<16; ++j)
                {
                   //avoid negative arguments in log
                   if(eval[j]>0.0){s2_(p,q) = s2_(p,q)  - eval[j]*log(eval[j]);}
                }
                s2_(q,p) = s2_(p,q);
            }
        }

        //  calculate I
        I_.resize(L,L);
        for (int p=0; p<L; ++p)
        {
            for (int q=(p+1); q<L; ++q)
            {
                I_(p,q) = 0.5*(s1_(0,p)+ s1_(0,q) - s2_(p,q));
                I_(q,p) = I_(p,q);
            }
        }

        // get CAS vector and sort HF occ vector
        std::vector<std::pair<real_type, int> > casv_sort(L);
        std::vector<int> casv(L);
        for(int i = 0; i < L; i++){
            casv_sort[i].first = s1_(0,i);
            casv_sort[i].second = i;
        }
        std::sort(casv_sort.begin(), casv_sort.end(), entanglement_detail::comp<std::pair<real_type, int> >);
        for(int i = 0; i<L; i++){
            casv[i] = casv_sort[i].second;
        }
        // print CAS vector and sorted HF vector
        //std::cout << "CAS vector: ";
        //for(int i =0; i<L; i++){std::cout << ","<<casv[i];}
        //std::cout <<std::endl;
    }
};

