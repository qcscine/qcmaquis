/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include <math.h>
#include <algorithm>
#include <alps/numeric/matrix.hpp>

#include "dmrg/models/lattice.h"
#include "../applications/tools/deas.hpp"
#include "../applications/cideas/determinant.hpp"


//function to get orb under consideration
std::vector<int>  el_casv(const int &side, const int &left, std::vector<int> &casv){
   int len = casv.size();
   std::vector<int> casv_el;
   for(int i = 0; i<len; i++){
   //distuinguish left (side=0) and right(=1) part
      if(side==0){
         if(casv[i]>left-1){casv_el.push_back(casv[i]);}
      } else{
         if(casv[i]<left){casv_el.push_back(casv[i]);}
      }
   }
  return casv_el;
}

//function to copy determinants
template <class SymmGroup>
std::vector<Determinant<SymmGroup> > copy_det(std::vector<Determinant<SymmGroup> > dets){
   int end = dets.size();
   for(int i= 0; i<3; i++){    //number of copies
      for(int j =0;j<end;j++){ //number of already existing determinants
         dets.push_back(dets[j]);
      }
   }
   return dets;
}

//function to actually perform deas
template <class SymmGroup>
std::vector<Determinant<SymmGroup> > deas(const int &pos, const int &act_orb, std::vector<Determinant<SymmGroup> > dets)
{
   std::vector<int> orb_list;
   for(int i=0;i<4;i++){orb_list.push_back(i+1);}
   orb_list.erase(orb_list.begin()+dets[0][act_orb]-1);
   for(int i=0; i<3; i++){ 					//loop over possibilities
      for(int j=0; j<pow(4,pos); j++){				//loop over blocks
         dets[pow(4,pos)*(i+1)+j][act_orb] = orb_list[i];
      }
   }
   return dets;
}

//some functions to get info on left part
std::vector<int> get_half_filled(int nelec_left, int left){
   std::vector<int> half_filled;
   if(nelec_left>left){nelec_left=left-(nelec_left-left);}
   if(nelec_left>8){
     if(nelec_left%2==0){nelec_left=8;}
     else{nelec_left=7;}
   }
   for(int i = nelec_left; i>=0; i-=2){
      half_filled.push_back(i);
   }
   return half_filled;
}
//function to reduce the symmetry vector - here actually no pair is required since later only first entry is used
std::vector<std::pair<int,int> > reduce_symvec(const std::vector<int> &symvec_left){
   std::vector<std::pair<int,int> > symvec_red;
   int size  = 0;
   for(int i = 0; i<8; i++){
      size = 0;
      for(int j = 0; j<symvec_left.size(); j++){
         if(symvec_left[j]==i){size++;}
      }
      if(size>0){symvec_red.push_back(std::make_pair(i,size));}
   }
   return symvec_red;
}


//function to determine occupied orbitals in det
std::vector<std::pair<int,int> > get_orb(std::vector<int> hf_occ){
   std::vector<std::pair<int,int> > occ_orb;
   for(int i= 0; i< hf_occ.size(); i++){
      if(hf_occ[i]!=1){occ_orb.push_back(std::make_pair(i, hf_occ[i]));}
   }
   return occ_orb;
}


template <class SymmGroup>
std::vector<Determinant<SymmGroup> > generate_deas(DmrgParameters &parms, EntanglementData<matrix> &em, int &run, std::vector<Determinant<SymmGroup> > &deas_dets)
{
    int L = parms["L"];

    alps::numeric::matrix<int> prd(8,8,0);
    prd(0,1)=prd(2,3)=prd(4,5)=prd(6,7)=1;
    prd(0,2)=prd(1,3)=prd(4,6)=prd(5,7)=2;
    prd(0,3)=prd(1,2)=prd(4,7)=prd(5,6)=3;
    prd(0,4)=prd(1,5)=prd(2,6)=prd(3,7)=4;
    prd(0,5)=prd(1,4)=prd(2,7)=prd(3,6)=5;
    prd(0,6)=prd(1,7)=prd(2,4)=prd(3,5)=6;
    prd(0,7)=prd(1,6)=prd(2,5)=prd(3,4)=7;
    for(int i=0;i<8; i++){
       for(int j=0; j<i;j++){
          prd(i,j)=prd(j,i);
       }
    }

    Determinant<SymmGroup> hf_occ(parms.get<std::vector<int> >("hf_occ"));

    std::vector<mpair> casv_sort(L);
    std::vector<int> casv(L);
    for(int i = 0; i<L; i++){
       casv_sort[i].first = em.s1_(0,i);
       casv_sort[i].second = i;
    }
    std::sort(casv_sort.begin(),casv_sort.end(), comp);
    for(int i = 0; i<L; i++){
       casv[i] = casv_sort[i].second;
     }

    std::cout << "CAS vector: ";
    for(int i =0; i<L; i++){std::cout << casv[i] << " ";}
    std::cout <<std::endl;

//    std::cout << "HF occupation  vector: ";
//    for(int i =0; i<L; i++){std::cout << hf_occ[i];}
//    std::cout <<std::endl;

    Lattice lat(parms);
    std::vector<int> sym_vec(lat.size());
    for (int i = 0; i<lat.size(); i++){
       sym_vec[i] = lat.get_prop<int>("type", i);
    }

    int target_sym, hf_sym;
    target_sym = parms["irrep"];
    hf_sym = hf_occ.sym_check(sym_vec, prd);
    std::cout << "symmetry of HF determinant: " << hf_sym << " vs. target symmetry: " << target_sym << std::endl;

    int nelec = 0;
    nelec = hf_occ.num_el();
    std::cout << "number of electrons: " << nelec <<std::endl;
    int spin = 0;
    spin = hf_occ.spin_check();
    std::cout << "spin of HF determinant: " << spin << std::endl;

    std::vector<int> ci_level(parms.get<std::vector<int> >("ci_level"));
    if(std::find(ci_level.begin(), ci_level.end(), 0) == ci_level.end())
       ci_level.push_back(0);

    std::cout << "excitation levels allowed: ";
    for(int i =0; i<ci_level.size(); i++){std::cout << ci_level[i] << " ";}
    std::cout <<std::endl;

    int act_orb = casv[0];
    if(run == 0){
       deas_dets.push_back(hf_occ);
       int count = 0;
       for(int i = 1; i<5; i++){
          if(hf_occ[act_orb]!=i){
             deas_dets.push_back(hf_occ);
             count++;
             deas_dets[count][act_orb] =i;
          }
       }
    }
    else if (run != 0){

        assert(deas_dets.size() == pow(4,run));
        std::vector<Determinant<SymmGroup> > new_deas;

        deas_dets = copy_det(deas_dets);
        act_orb = casv[run];
        deas_dets = deas(run,act_orb,deas_dets);
    }

//insert ci check here
    return deas_dets;
}
