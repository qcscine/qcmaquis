

/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Yingjin Ma <yingjin.ma@phys.chem.ethz.ch>
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

#ifndef SAMPLING_HPP
#define SAMPLING_HPP

#include <iostream>
#include <iomanip>
#include <boost/random.hpp>
#include <math.h>
#include <string.h>

#ifdef DEBUG
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

#include "ci_encode.hpp"
//#include "ci_deas.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::string;

template<class Matrix, class SymmGroup>
std::vector<typename SymmGroup::charge> 
parse_det(std::vector<typename SymmGroup::charge> const & det, std::vector<Index<SymmGroup> > const & site_dims)
{
    typename SymmGroup::charge A(1), B(0), C(0), D(0);
    B[0]=1; C[1]=1;

    std::vector<typename SymmGroup::charge> tmp;

    for(std::size_t i = 0; i < det.size(); i++) {
        if(det[i][0]==1 && det[i][1]==1)  
          tmp.push_back(site_dims[i][0].first); // doubly occ
        else if(det[i][0]==1 && det[i][1]==0)
          tmp.push_back(site_dims[i][1].first); // up
        else if(det[i][0]==0 && det[i][1]==1)
          tmp.push_back(site_dims[i][2].first); // down 
        else if(det[i][0]==0 && det[i][1]==0)
          tmp.push_back(site_dims[i][3].first); // empty            
        }

    return tmp;
}

// refer from wiki
void quick_sort (double data[], int left,
	int right) {
	int p = (left + right) / 2;
	int pivot = data[p];
	int i = left,j = right;
	for ( ; i < j;) {
		while (! (i>= p || std::abs(pivot) < std::abs(data[i])))
			++i;
		if (i < p) {
			data[p] = data[i];
			p = i;
		}
		while (! (j <= p || std::abs(data[j]) < std::abs(pivot)))
			--j;
		if (j > p) {
			data[p] = data[j];
			p = j;
		}
	}
	data[p] = pivot;
	if (p - left > 1)
		quick_sort(data, left, p - 1);
	if (right - p > 1)
		quick_sort(data, p + 1, right);
}

void quicksort(string dets[], double b[], int left, int right) {

    double pivot = std::abs(b[(left+right)/2]);
    int l = left;
    int r = right;
    while(l <= r) {
     while (std::abs(b[l]) < pivot) l++;
     while (std::abs(b[r]) > pivot) r--;
     if (l <= r) {

        double tmp;
        tmp = b[l];
        b[l]= b[r];
        b[r]= tmp;

        string ctmp;
        ctmp   = dets[l];
        dets[l]= dets[r];
        dets[r]= ctmp;

        l++;
        r--;
     }
    };

    if (left < r)  quicksort(dets, b, left, r);
    if (l < right) quicksort(dets, b, l, right);
}


struct Sampling
{

// Random number generator, work together with boost library
        Sampling(int id)
                :distribution(0.,1.),random_number(generator,distribution)
        {
                unsigned int seed = 123456;
                generator.seed(seed/(id+1)+time(NULL));
        }

        Sampling()
                :distribution(0.,1.),random_number(generator,distribution)
        {
                unsigned int seed = 123456;
                generator.seed(seed+time(NULL));
        }

// template generate_dets 
       template <typename Determinants, typename MPS, typename Hash_value, typename Hash_index, typename SiteSymm, typename Determinant, typename SymmGroup> 
       void generate_dets(Determinants dets,
                          Determinants dets_mclr,
                          MPS mps,
                          Hash_value hash,
                          Hash_index hash_index,
                          SiteSymm site_dims,
                          int norb,
                          int nsample,
                          int nitermax,
                          double CI_threshold,
                          double COM_threshold)
      {
        maquis::cout << "    CI-threshold  : " <<  CI_threshold << std::endl; 
        maquis::cout << "   COM-threshold  : " << COM_threshold << std::endl;

// Set up the initial parameters, fix the symmetry group
        double ci,ci0,ci_ratio,ci_tmp;
        double sum_ci2=0.0,completeness=0.0;	  
        unsigned det_length = dets[0].size();
	unsigned number_of_dets;
        typename Hash_value::iterator iter;
        typename Hash_index::iterator iter_idx;

        typedef typename SymmGroup::charge charge;
        charge target = mps[mps.length()-1].col_dim()[0].first;

// This part will be deleted, since the elements in NU1ChargePG could already be used.
//        std::ifstream dets_file;
//        dets_file.open(file.c_str());      

// This part will be used as the determinants reservoir
        Determinant det;
        for (std::size_t c = 0; c < dets.size(); ++c)
           {
            det=dets[c];
            hash_index[det]=c;
	   }

// determinant spawnning -- preparing part

        Determinant det_queen,det_bee,det_tmp; // determinant "queen bee tmp"
        int iaccept=0,iaccept_queen=0;
 
        for (std::size_t c = 0; c < dets_mclr.size(); ++c)
           {
            det=dets_mclr[c];
            ci0=extract_coefficient(mps, det); 	   
            hash[det]=ci0;
//            maquis::cout << "follow determinant " << c << " with coefficient " << ci0 << std::endl;
//            for(int p = det_length-1; p >= 0; --p)            
//                maquis::cout << det[p][0] << det[p][1] << det[p][2] << std::endl;
	   }

// Get the number of electrons & holes (alpha, beta, total)
        int nele_alpha=0;
        int nele_beta =0;
        int nhole_alpha=0;
        int nhole_beta =0;
        int nele_total=0;
        int igroup_sym=0;  
        for(int p = det_length-1; p >= 0; --p)
           {              
            if(det[p][0]==1)
              nele_alpha++ ;
            else
              nhole_alpha++;
            if(det[p][1]==1) 
              nele_beta++  ;
            else
              nhole_beta++ ;
            }
        nele_total=nele_alpha+nele_beta;

// Initial vectors
        std::vector<int>ele_alpha(nele_alpha); 
        std::vector<int>ele_beta(nele_beta);
        std::vector<int>hole_alpha(nhole_alpha); 
        std::vector<int>hole_beta(nhole_beta); 
        
// determinant spawnning -- doing part

        det_queen=dets[0];
        int nMAX=0;
        do{
          nMAX++;
          // Get the number of excited electrons        
          int nele_excited;
          nele_excited=int(floor(nele_total*random_number())+1);
          maquis::cout << " nele_excited " << nele_excited << " in " << norb << " orbitals" << std::endl;
          for(int isample=0; isample < nsample; isample++ ){
            det_tmp=det_queen;
            int iele_excited=0;
            do{ 
               int i=0,j=0,k=0,l=0;  
               for(int p=0; p<norb; p++ ){
                   if(det_tmp[p][0]==1) 
                     ele_alpha[i++]=p;
                   else
                     hole_alpha[j++]=p;
 
                   if(det_tmp[p][1]==1) 
                     ele_beta[k++]=p;
                   else
                     hole_beta[l++]=p;
                   }
               iele_excited++;

               if(nhole_alpha !=0 && nele_alpha!=0 && random_number()>0.5){
                  int fr =int(floor(random_number()*nele_alpha ));   //fr -> from 
                  int to =int(floor(random_number()*nhole_alpha));    
                  int ifr=ele_alpha[fr] ;
                  int ito=hole_alpha[to];              
//                maquis::cout << " alpha - electron from " << ifr << " goto " << ito << std::endl;              
                  det_tmp[ifr][0]=0;              
                  det_tmp[ito][0]=1;              
                 }
               else if(nhole_beta !=0 && nele_beta!=0){
                  int fr =int(floor(random_number()*nele_beta ));   //fr -> from 
                  int to =int(floor(random_number()*nhole_beta));    
                  int ifr=ele_beta[fr] ;
                  int ito=hole_beta[to];
//                maquis::cout << " beta  - electron from " << ifr << " goto " << ito << std::endl;
                  det_tmp[ifr][1]=0;              
                  det_tmp[ito][1]=1;              
                 }
 
               det_bee=parse_det<matrix, TwoU1PG>(det_tmp,site_dims);
// Check the symmetry group
               charge total = std::accumulate(det_bee.begin(), det_bee.end(), SymmGroup::IdentityCharge);
  
               if(iele_excited >= nele_excited && total == target)           
                 break;
 
              } while(true);

//            maquis::cout << " This det_bee is : " << std::endl;
//            for(int p = det_length-1; p >= 0; --p)            
//                maquis::cout << det_bee[p][0] << det_bee[p][1] << det_bee[p][2]  << std::endl;

//        for(typename Determinants::iterator it = dets.begin();it != dets.end(); ++it)
            ci=extract_coefficient(mps, det_bee); 	    
            sum_ci2=sum_ci2+pow(ci,2.0); 
//            maquis::cout << isample << "-th " << " CI coefficient: " << ci << std::endl;      

// Use Hash(Map) technique, in order to quick access the data
            iter = hash.find(det_bee);
            if(iter == hash.end()){
              if(std::abs(ci) >= CI_threshold){
                hash[det_bee]=ci;
                iaccept++;
                }
              }
            else{
              ci=iter->second;
              } 

            ci_ratio=pow(ci,2.0)/pow(ci0,2);
           
// Whether use this bee-det as the new queen-det 
            if(ci_ratio > random_number()){
              det_queen=det_bee;
              ci0=ci;
              iaccept_queen++;
              }

           }
 
          sum_ci2=0.0;     
          for(iter=hash.begin(); iter!=hash.end(); iter++){
             ci_tmp = iter->second;
             sum_ci2=sum_ci2+pow(ci_tmp,2.0);
//             maquis::cout << "Determinant-" << "  coeff= " << ci_tmp << std::endl;
             }

//        number_of_dets=dets.size();
//        maquis::cout << " Total number of determinants : " << number_of_dets << std::endl;
          maquis::cout << "Determinant-naccept" << iaccept << std::endl;
          maquis::cout << "Determinant-naccept-queen" << iaccept_queen << std::endl;
          completeness=1.0-sum_ci2;
          maquis::cout << " Current completeness (1-\\sum(ci^2)) : " << completeness << std::endl;

          } while(completeness>COM_threshold &&  nMAX< nitermax);

// Back to "2", "u", "d", "0" type electron representation, and print out
          //Determinants  dets_show;
          double    CIs_show[hash.size()]; //value
          uint32_t  CIs_index[hash.size()]; //index in the reservoir
          string   dets_show[hash.size()]; //dets represent
                   
          uint32_t i=0;
          for(iter=hash.begin(); iter!=hash.end(); iter++){
              string ctmp;
              CIs_show[i]=iter->second;
            // idets_show.push_back(iter->first);
              for(int p = 0; p < det_length; p++)  {
                  if(iter->first[p][0]==1 && iter->first[p][1]==1)  
                    ctmp=ctmp+"2";
                  if(iter->first[p][0]==1 && iter->first[p][1]==0)  
                    ctmp=ctmp+"u";
                  if(iter->first[p][0]==0 && iter->first[p][1]==1)  
                    ctmp=ctmp+"d";
                  if(iter->first[p][0]==0 && iter->first[p][1]==0)  
                    ctmp=ctmp+"0";
                 }
         //     maquis::cout << " The ctmp " << ctmp << std::endl; 
              dets_show[i]=ctmp;

              iter_idx = hash_index.find(iter->first);
              if(iter_idx == hash_index.end()){
                CIs_index[i]=0;
                }
              else{
                CIs_index[i]=iter_idx->second+1;
                }
              i++;
              }

//          for(int i=0; i<hash.size();i++) 
//            maquis::cout << " Current dets " << dets_show[i] << std::endl;
          std::ofstream fout;
          fout.open("det_coeff.tmp"); 
          fout << i << std::endl;
          for(uint32_t i=0; i<hash.size() ; i++){
            if(CIs_index[hash.size()-i-1]!=0) 
             {fout << CIs_index[hash.size()-i-1] << "  " << std::fixed << std::setprecision(12) << CIs_show[hash.size()-i-1] << std::endl;}
            }
          
          quicksort(dets_show, CIs_show,0, hash.size()-1); 

          for(int i=0; i<hash.size() ; i++)
            maquis::cout << " The sorted determinants: " << dets_show[hash.size()-i-1] << " " << CIs_show[hash.size()-i-1] << std::endl;            
 
        }    
       
        boost::mt19937 generator;
        boost::uniform_real<> distribution;
        boost::variate_generator<boost::mt19937&, boost::uniform_real<> > random_number;

};



#endif


