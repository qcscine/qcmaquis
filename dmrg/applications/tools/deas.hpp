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

#ifdef USE_AMBIENT
    #include "dmrg/block_matrix/detail/ambient.hpp"
    typedef ambient::numeric::tiles<ambient::numeric::matrix<double> > matrix;
#else
    #include "dmrg/block_matrix/detail/alps.hpp"
    typedef alps::numeric::matrix<double> matrix;
#endif

#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/utils/storage.h"


namespace entanglement_detail {

    //function definition get_labels 
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

    //comparison function
    typedef std::pair<double,int> mpair;
    bool comp(const mpair& l, const mpair& r)
    {
        return l.first > r.first;
    }

    template <class Matrix>
    Matrix assy_hc(const std::vector<std::string> & label_str, Matrix const & values, const int L)
    {
       std::vector<std::pair<int,int> > labels;
       labels = get_labels(label_str);
       Matrix ret(L,L,0);
       //symmetric Matrix build with entries and labels
       for(int i=0;i<labels.size();i++){
          ret(labels[i].first, labels[i].second) = double(values(i,0));
          ret(labels[i].second, labels[i].first) = double(values(i,0));
       }
       return ret;
    }

    template <class Matrix>
    Matrix load_vector(storage::archive & ar, std::string path)
    {
        std::vector<std::string> x;
        Matrix y;

        ar[path + "/mean/labels"] >> x;
        ar[path + "/mean/value"] >> y;

        std::vector<std::pair<int,int> > labels = get_labels(x);

        Matrix ret(num_rows(y));

        return ret;
    }
  
    template <class Matrix>
    struct EntropyData
    {
        Matrix Nup;
        Matrix Ndown; 
        Matrix Nupdown; 
        Matrix dm_up_4x4;
        Matrix dm_down_4x4;
        Matrix nupnup_4x4;
        Matrix nupndown_4x4;
        Matrix ndownnup_4x4;
        Matrix ndownndown_4x4;
        Matrix doccdocc_4x4;
        Matrix transfer_up_while_down_4x4;
        Matrix transfer_down_while_up_4x4;
        Matrix transfer_up_while_down_at_2_4x4;
        Matrix transfer_down_while_up_at_2_4x4;
        Matrix transfer_up_while_down_at_1_4x4;
        Matrix transfer_down_while_up_at_1_4x4;
        Matrix transfer_pair_4x4;
        Matrix spinflip_4x4;
        Matrix nupdocc_4x4;
        Matrix ndowndocc_4x4;
        Matrix doccnup_4x4;
        Matrix doccndown_4x4;
    }; 


    template<class Matrix>
    EntropyData<Matrix> loadData(storage::archive & ar)
    {
        typedef typename Matrix::value_type value_type;
        DmrgParameters parms;
        ar["/parameters"] >> parms;

        // chain-length   
        int L = parms["L"];

        // collect all in EntropyData => get as e.g. data.Nup etc. 
        EntropyData<Matrix> data; 

        // load quantities needed for s1
        ar["/spectrum/results/Nup/mean/value"] >> data.Nup;
        //std::cout << "Nup: " << data.Nup << std::endl;
        ar["/spectrum/results/Ndown/mean/value"] >> data.Ndown;
        //std::cout << "Ndown " << data.Ndown << std::endl;
        ar["/spectrum/results/Nupdown/mean/value"] >> data.Nupdown;
        //std::cout << "docc: " << data.Nupdown << std::endl;

        //Matrix mm = load_vector<Matrix>(ar, "/spectrum/results/Nup");

        // load quantities needed for s2
        Matrix dm_up;
        std::vector<std::string> dm_up_str;
        ar["/spectrum/results/dm_up/labels"] >> dm_up_str;
        ar["/spectrum/results/dm_up/mean/value"] >> dm_up;
        data.dm_up_4x4 = assy_hc(dm_up_str,dm_up,L);
      
        Matrix dm_down;
        std::vector<std::string> dm_down_str;
        ar["/spectrum/results/dm_down/mean/value"] >> dm_down;
        ar["/spectrum/results/dm_down/labels"] >> dm_down_str;
        data.dm_down_4x4 = assy_hc(dm_down_str,dm_down,L);

        Matrix nupnup;
        std::vector<std::string> nupnup_str;
        ar["/spectrum/results/nupnup/mean/value"] >> nupnup;
        ar["/spectrum/results/nupnup/labels"] >> nupnup_str;
        data.nupnup_4x4 = assy_hc(nupnup_str,nupnup,L);

        Matrix nupndown;
        std::vector<std::string> nupndown_str;
        ar["/spectrum/results/nupndown/mean/value"] >> nupndown;
        ar["/spectrum/results/nupndown/labels"] >> nupndown_str;
        data.nupndown_4x4 = assy_hc(nupndown_str, nupndown,L);

        Matrix ndownnup;
        std::vector<std::string> ndownnup_str;
        ar["/spectrum/results/ndownnup/mean/value"] >> ndownnup;
        ar["/spectrum/results/ndownnup/labels"] >> ndownnup_str;
        data.ndownnup_4x4 = assy_hc(ndownnup_str,ndownnup,L);

        Matrix ndownndown;
        std::vector<std::string> ndownndown_str;
        ar["/spectrum/results/ndownndown/mean/value"] >> ndownndown;
        ar["/spectrum/results/ndownndown/labels"] >> ndownndown_str;
        data.ndownndown_4x4 = assy_hc(ndownndown_str,ndownndown,L);

        Matrix doccdocc;
        std::vector<std::string> doccdocc_str;
        ar["/spectrum/results/doccdocc/mean/value"] >> doccdocc;
        ar["/spectrum/results/doccdocc/labels"] >> doccdocc_str;
        data.doccdocc_4x4 = assy_hc(doccdocc_str,doccdocc,L);

        Matrix transfer_up_while_down;
        std::vector<std::string> transfer_up_while_down_str;
        ar["/spectrum/results/transfer_up_while_down/mean/value"] >> transfer_up_while_down;
        ar["/spectrum/results/transfer_up_while_down/labels"] >> transfer_up_while_down_str;
        data.transfer_up_while_down_4x4 = assy_hc(transfer_up_while_down_str,transfer_up_while_down,L);

        Matrix transfer_down_while_up;
        std::vector<std::string> transfer_down_while_up_str;
        ar["/spectrum/results/transfer_down_while_up/mean/value"] >> transfer_down_while_up;
        ar["/spectrum/results/transfer_down_while_up/labels"] >> transfer_down_while_up_str;
        data.transfer_down_while_up_4x4 = assy_hc(transfer_down_while_up_str,transfer_down_while_up,L);

        Matrix transfer_up_while_down_at_2;
        std::vector<std::string> transfer_up_while_down_at_2_str;
        ar["/spectrum/results/transfer_up_while_down_at_2/mean/value"] >> transfer_up_while_down_at_2;
        ar["/spectrum/results/transfer_up_while_down_at_2/labels"] >> transfer_up_while_down_at_2_str;
        data.transfer_up_while_down_at_2_4x4 = assy_hc(transfer_up_while_down_at_2_str,transfer_up_while_down_at_2,L);

        Matrix transfer_up_while_down_at_1;
        std::vector<std::string> transfer_up_while_down_at_1_str;
        ar["/spectrum/results/transfer_up_while_down_at_1/mean/value"] >> transfer_up_while_down_at_1;
        ar["/spectrum/results/transfer_up_while_down_at_1/labels"] >> transfer_up_while_down_at_1_str;
        data.transfer_up_while_down_at_1_4x4 = assy_hc(transfer_up_while_down_at_1_str,transfer_up_while_down_at_1,L);

        Matrix transfer_down_while_up_at_2;
        std::vector<std::string> transfer_down_while_up_at_2_str;
        ar["/spectrum/results/transfer_down_while_up_at_2/mean/value"] >> transfer_down_while_up_at_2;
        ar["/spectrum/results/transfer_down_while_up_at_2/labels"] >> transfer_down_while_up_at_2_str;
        data.transfer_down_while_up_at_2_4x4 = assy_hc(transfer_down_while_up_at_2_str,transfer_down_while_up_at_2,L);

        Matrix transfer_down_while_up_at_1;
        std::vector<std::string> transfer_down_while_up_at_1_str;
        ar["/spectrum/results/transfer_down_while_up_at_1/mean/value"] >> transfer_down_while_up_at_1;
        ar["/spectrum/results/transfer_down_while_up_at_1/labels"] >> transfer_down_while_up_at_1_str;
        data.transfer_down_while_up_at_1_4x4 = assy_hc(transfer_down_while_up_at_1_str,transfer_down_while_up_at_1,L);

        Matrix transfer_pair;
        std::vector<std::string> transfer_pair_str;
        ar["/spectrum/results/transfer_pair/mean/value"] >> transfer_pair;
        ar["/spectrum/results/transfer_pair/labels"] >> transfer_pair_str;
        data.transfer_pair_4x4 = assy_hc(transfer_pair_str,transfer_pair,L);

        Matrix spinflip;
        std::vector<std::string> spinflip_str;
        ar["/spectrum/results/spinflip/mean/value"] >> spinflip;
        ar["/spectrum/results/spinflip/labels"] >> spinflip_str;
        data.spinflip_4x4 = assy_hc(spinflip_str,spinflip,L);

        Matrix nupdocc;
        std::vector<std::string> nupdocc_str;
        ar["/spectrum/results/nupdocc/mean/value"] >> nupdocc;
        ar["/spectrum/results/nupdocc/labels"] >> nupdocc_str;
        data.nupdocc_4x4 = assy_hc(nupdocc_str,nupdocc,L);

        Matrix ndowndocc;
        std::vector<std::string> ndowndocc_str;
        ar["/spectrum/results/ndowndocc/mean/value"] >> ndowndocc;
        ar["/spectrum/results/ndowndocc/labels"] >> ndowndocc_str;
        data.ndowndocc_4x4 = assy_hc(ndowndocc_str,ndowndocc,L);

        Matrix doccnup;
        std::vector<std::string> doccnup_str;
        ar["/spectrum/results/doccnup/mean/value"] >> doccnup;
        ar["/spectrum/results/doccnup/labels"] >> doccnup_str;
        data.doccnup_4x4 = assy_hc(doccnup_str,doccnup,L);

        Matrix doccndown;
        std::vector<std::string> doccndown_str;
        ar["/spectrum/results/doccndown/mean/value"] >> doccndown;
        ar["/spectrum/results/doccndown/labels"] >> doccndown_str;
        data.doccndown_4x4 = assy_hc(doccndown_str,doccndown,L);

        return data;
    }
   

    template <class Matrix>
    Matrix two_orb_rdm(int p, int q, EntropyData<Matrix> & data)
    {
            Matrix pq_dm_matrix(16,16);
            pq_dm_matrix( 0, 0) = 1 + data.Nupdown(p,0)  + data.Nupdown(q,0) + data.doccdocc_4x4(p,q) - data.Ndown(p,0)\
                                  - data.ndowndocc_4x4(p,q) - data.Ndown(q,0) - data.doccndown_4x4(p,q)\
                                  + data.ndownndown_4x4(p,q) - data.Nup(p,0) - data.nupdocc_4x4(p,q) + data.nupndown_4x4(p,q)\
                                  - data.Nup(q,0) - data.doccnup_4x4(p,q) + data.ndownnup_4x4(p,q) + data.nupnup_4x4(p,q);
            // O(6)/O(1)
            pq_dm_matrix( 1, 1) = - data.Nupdown(p,0) - data.doccdocc_4x4(p,q) + data.Ndown(p,0) + data.ndowndocc_4x4(p,q) \
                                  + data.doccndown_4x4(p,q) - data.ndownndown_4x4(p,q) + data.doccnup_4x4(p,q) - data.ndownnup_4x4(p,q);
            // O(11)O(1)
            pq_dm_matrix( 2, 2) = - data.Nupdown(p,0) - data.doccdocc_4x4(p,q) + data.doccndown_4x4(p,q) + data.Nup(p,0) \
                                  + data.nupdocc_4x4(p,q) - data.nupndown_4x4(p,q) + data.doccnup_4x4(p,q) - data.nupnup_4x4(p,q); 
            // O(16)O(1)
            pq_dm_matrix( 3, 3) =  data.Nupdown(p,0)  - data.doccndown_4x4(p,q) - data.doccnup_4x4(p,q) + data.doccdocc_4x4(p,q); 
            // O(1)O(6)
            pq_dm_matrix( 4, 4) = -data.Nupdown(q,0)  - data.doccdocc_4x4(p,q)  + data.ndowndocc_4x4(p,q) + data.Ndown(q,0) \
                                +  data.doccndown_4x4(p,q) - data.ndownndown_4x4(p,q) + data.nupdocc_4x4(p,q) - data.nupndown_4x4(p,q);
            // O(6)O(6) 
            pq_dm_matrix( 5, 5) = data.ndownndown_4x4(p,q) - data.ndowndocc_4x4(p,q) - data.doccndown_4x4(p,q) + data.doccdocc_4x4(p,q);
            // O(11)O(6)
            pq_dm_matrix( 6, 6) = data.nupndown_4x4(p,q) - data.doccndown_4x4(p,q) - data.nupdocc_4x4(p,q) + data.doccdocc_4x4(p,q); 
            // O(16)O(11)
            pq_dm_matrix( 7, 7) = data.doccndown_4x4(p,q) - data.doccdocc_4x4(p,q); 
            // O(1)O(11)
            pq_dm_matrix( 8, 8) = -data.Nupdown(q,0) - data.doccdocc_4x4(p,q) + data.ndowndocc_4x4(p,q) + data.nupdocc_4x4(p,q) 
                                  + data.Nup(q,0) + data.doccnup_4x4(p,q) - data.ndownnup_4x4(p,q) - data.nupnup_4x4(p,q); 
            // O(6)O(11)
            pq_dm_matrix( 9, 9) = data.ndownnup_4x4(p,q) - data.ndowndocc_4x4(p,q) - data.doccnup_4x4(p,q) + data.doccdocc_4x4(p,q);
            // O(11)O(11)
            pq_dm_matrix(10,10) = data.nupnup_4x4(p,q) - data.nupdocc_4x4(p,q) - data.doccnup_4x4(p,q) + data.doccdocc_4x4(p,q); 
            // O(16)O(11)
            pq_dm_matrix(11,11) = data.doccnup_4x4(p,q) - data.doccdocc_4x4(p,q); 
            // O(1)O(16)
            pq_dm_matrix(12,12) = data.Nupdown(q,0) - data.nupdocc_4x4(p,q) - data.ndowndocc_4x4(p,q) + data.doccdocc_4x4(p,q); 
            // O(1)O(16)
            pq_dm_matrix(13,13) = data.ndowndocc_4x4(p,q) - data.doccdocc_4x4(p,q);
            // O(11)O(16)
            pq_dm_matrix(14,14) = data.nupdocc_4x4(p,q) - data.doccdocc_4x4(p,q);
            // O(16)O(16)
            pq_dm_matrix(15,15) = data.doccdocc_4x4(p,q);
            // O(9)O(2)
            pq_dm_matrix(1,4)  = data.dm_down_4x4(p,q) - data.transfer_down_while_up_at_1_4x4(p,q) - data.transfer_down_while_up_at_2_4x4(p,q) 
                               + data.transfer_down_while_up_4x4(p,q);
            pq_dm_matrix(4,1)  = pq_dm_matrix(1,4);
            // O(9)O(2)         
            pq_dm_matrix(2,8)  =  data.dm_up_4x4(p,q) - data.transfer_up_while_down_at_1_4x4(p,q) - data.transfer_up_while_down_at_2_4x4(p,q)
                               +  data.transfer_up_while_down_4x4(p,q);
            pq_dm_matrix(8,2)  = pq_dm_matrix(2,8);                   
            // O(15)O(2)
            pq_dm_matrix(3,6) = data.transfer_down_while_up_at_1_4x4(p,q) - data.transfer_down_while_up_4x4(p,q);
            pq_dm_matrix(6,3) = pq_dm_matrix(3,6);
            // O(14)O(3)
            pq_dm_matrix(3,9) = -data.transfer_up_while_down_at_1_4x4(p,q) + data.transfer_up_while_down_4x4(p,q);
            pq_dm_matrix(9,3) = pq_dm_matrix(3,9);
            // O(10)O(7)
            pq_dm_matrix(6,9) = data.spinflip_4x4(p,q);
            pq_dm_matrix(9,6) = pq_dm_matrix(6,9);
            // O(13)O(4)
            pq_dm_matrix(3,12) = data.transfer_pair_4x4(p,q);
            pq_dm_matrix(12,3) = pq_dm_matrix(3,12);
            // O(9)O(8)
            pq_dm_matrix(6,12) = -data.transfer_up_while_down_at_2_4x4(p,q) + data.transfer_up_while_down_4x4(p,q);
            pq_dm_matrix(12,6) = pq_dm_matrix(6,12);
            // O(8)O(9)
            pq_dm_matrix(9,12) = data.transfer_down_while_up_at_2_4x4(p,q) - data.transfer_down_while_up_4x4(p,q);
            pq_dm_matrix(12,9) = pq_dm_matrix(9,12);
            // O(14)O(8)
            pq_dm_matrix(7,13) = data.transfer_up_while_down_4x4(p,q); 
            pq_dm_matrix(13,7) = pq_dm_matrix(7,13);
            // O(15)O(12)
            pq_dm_matrix(11,14) = data.transfer_down_while_up_4x4(p,q);
            pq_dm_matrix(14,11) = pq_dm_matrix(11,14);

            return pq_dm_matrix;
    }
} // namespace entanglement_detail

template <class Matrix> 
class EntanglementData
{
public:
    EntanglementData(std::string rfile)
    {
        // Load results from rfile:
        storage::archive ar(rfile, "r");
        DmrgParameters parms;
        ar["/parameters"] >> parms;

        // chain-length   
        int L = parms["L"];
        //std::cout << L << std::endl;

        entanglement_detail::EntropyData<Matrix> data = entanglement_detail::loadData<Matrix>(ar);

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
                std::vector<double> eval(16);
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
        std::vector<int> hf_occ = parms["hf_occ"];
        std::vector<int> hf_sort(L);

        // get CAS vector and sort HF occ vector
        std::vector<entanglement_detail::mpair> casv_sort(L);
        std::vector<int> casv(L);
        for(int i = 0; i < L; i++){
            casv_sort[i].first = s1_(0,i);
            casv_sort[i].second = i;
        }
        std::sort(casv_sort.begin(), casv_sort.end(), entanglement_detail::comp); 
        for(int i = 0; i<L; i++){
            casv[i] = casv_sort[i].second;
            hf_sort[i] = hf_occ[casv[i]];
        }
        // print CAS vector and sorted HF vector
        //std::cout << "CAS vector: ";
        //for(int i =0; i<L; i++){std::cout << ","<<casv[i];}
        //std::cout <<std::endl;

        //std::cout << "sorted HF occupation  vector: ";
        //for(int i =0; i<L; i++){std::cout << ","<<hf_sort[i];}
        //std::cout <<std::endl;
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
};

