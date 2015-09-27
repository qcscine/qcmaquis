/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>                 //change that later
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

#include <alps/numeric/matrix.hpp>



template <class SymmGroup>
class Determinant : public std::vector<int>
{
    typedef std::vector<int> base;
public:
   //constructor
    Determinant(base input_det) : base(input_det)  {}
     
    Determinant(int L) : base(L) {}
   //default constructor
    Determinant() {}
   //copy constructor
    Determinant(const Determinant& Copy) {*this = Copy;}

//get number of electrons in a determinant
   int num_el()
   {
      int nelec = 0;
      for(int i=0; i<(*this).size(); i++){
         if((*this)[i] == 4){nelec+=2;}
         else if((*this)[i] == 1){;}
         else{nelec+=1;}
      }
     return nelec;
   }

//check symmetry of a determinant
   int sym_check(const base &sym_vec, const alps::numeric::matrix<int> &prd){
     int sym = 0;
     for(int i=0; i<(*this).size(); i++){
        if((*this)[i]==2||(*this)[i]==3){
           sym = prd(sym,sym_vec[i]);
        }
      }
      return sym;
   }

//check spin of a determinant
   int spin_check(){
      int spin = 0;
      for(int i=0;i<(*this).size();i++){
         if((*this)[i]==2){spin = spin-1;}
         else if((*this)[i]==3){spin = spin+1;}
      }
      return spin;
   }

//function to extract ci determinants of a given level
    bool ci_check(const std::vector<int> &ci_level, const std::vector<std::pair<int,int> > &hf_occ_orb){
      bool wrong_level = false;
      int diff = 0;
    //first check number of changes
      for(int i = 0; i <hf_occ_orb.size(); i++){
         if((*this)[hf_occ_orb[i].first] != hf_occ_orb[i].second){
            if(hf_occ_orb[i].second == 4 && (*this)[hf_occ_orb[i].first] == 3){
               diff++;
            }else if(hf_occ_orb[i].second == 4 && (*this)[hf_occ_orb[i].first] == 2){
               diff++;
            }else if(hf_occ_orb[i].second == 4 && (*this)[hf_occ_orb[i].first] == 1){
               diff = diff + 2;
            }else if(hf_occ_orb[i].second == 3 && (*this)[hf_occ_orb[i].first] == 1){
               diff++;
            }else if(hf_occ_orb[i].second == 2 && (*this)[hf_occ_orb[i].first] == 1){
               diff++;
            }
         }
      }
     //check if number of changes agrees with ci_level 
      for(int i = 0; i<ci_level.size(); i++){
         if(ci_level[i] != diff){
            wrong_level = true;
         }else{
            wrong_level = false;
            break;
         }
      }
      return wrong_level;
   }


};

