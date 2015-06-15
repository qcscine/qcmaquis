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
std::vector<std::vector<int> > copy_det(const int &pos, std::vector<std::vector<int> > dets){
   int end = dets.size();
   for(int i= 0; i<3; i++){    //number of copies
      for(int j =0;j<end;j++){ //number of already existing determinants
         dets.push_back(dets[j]);
      }
   }
   return dets;
}

//function to actually perform deas
std::vector<std::vector<int> > deas(const int &pos, const int &act_orb, std::vector<std::vector<int> > dets){
   std::vector<int> orb_list;
   for(int i=0;i<4;i++){orb_list.push_back(i+1);}
   orb_list.erase(orb_list.begin()+dets[0][act_orb]-1);
   for(int i=0;i<3;i++){ //loop over possibilities
      for(int j=0;j<pow(4,pos);j++){//loop over blocks
         dets[pow(4,pos)*(i+1)+j][act_orb] = orb_list[i];
      }
   }
   return dets;
}

//function to extract ci determinants of a given level
bool ci_check(const std::vector<int> &ci_level, const std::vector<std::pair<int,int> > &hf_occ_orb, const std::vector<int> &det ){
   bool wrong_level = false;
   int diff = 0;
//first check number of changes
   for(int i = 0; i <hf_occ_orb.size(); i++){
      if(det[hf_occ_orb[i].first] != hf_occ_orb[i].second){
         if(hf_occ_orb[i].second == 4 && det[hf_occ_orb[i].first] == 3){
            diff++;
         }else if(hf_occ_orb[i].second == 4 && det[hf_occ_orb[i].first] == 2){
            diff++;
         }else if(hf_occ_orb[i].second == 4 && det[hf_occ_orb[i].first] == 1){
            diff = diff + 2;
         }else if(hf_occ_orb[i].second == 3 && det[hf_occ_orb[i].first] == 1){
            diff++;
         }else if(hf_occ_orb[i].second == 2 && det[hf_occ_orb[i].first] == 1){
            diff++;
         }
      }
   }
//theck if number of changes agrees with ci_level 
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




//get number of electrons in a determinant
int num_el(const std::vector<int> &det){
    int nelec = 0;
    int L = det.size();
    for(int i=0; i<L; i++){
       if(det[i] == 4){nelec+=2;}
       else if(det[i] == 1){;}
       else{nelec+=1;}
    }
   return nelec;
}


//check symmetry of a determinant
int sym_check(const std::vector<int> &det, const std::vector<int> &sym_vec, const alps::numeric::matrix<int> &prd){
  int sym = 0;
  for(int i=0; i<det.size(); i++){
     if(det[i]==2||det[i]==3){
        sym = prd(sym,sym_vec[i]);
     }
  }

   return sym;
}

//check spin of a determinant
int spin_check(const std::vector<int> &det){
   int spin = 0;
   for(int i=0;i<det.size();i++){
      if(det[i]==2){spin = spin-1;}
      else if(det[i]==3){spin = spin+1;}
   }
   return spin;
}


//check if target symmetry can be reached with current determinant

std::vector<std::vector<int> > det_sym_check(const std::vector<std::vector<int> > &dets_left, const std::vector<std::vector<int> > &dets_right, const int &spin, const int &sym, const int &nelec, const std::vector<int> &symvec_left, const std::vector<int> symvec_right, const alps::numeric::matrix<int> &prd){

std::vector<std::vector<int> > deas_dets;
int n_left, n_right, spin_left, spin_right, sym_left, sym_right,tot_sym,tot_spin,tot_nel = 0;
//outer loop over all possible determinants
   for(int i = 0; i<dets_right.size();i++){
      n_right = num_el(dets_right[i]);
      spin_right = spin_check(dets_right[i]);
      sym_right = sym_check(dets_right[i],symvec_right,prd);
   //inner loop over all possible counterparts
      for(int j=0; j<dets_left.size();j++){
         n_left = num_el(dets_left[j]);
         spin_left = spin_check(dets_left[j]);
         sym_left = sym_check(dets_left[j],symvec_left,prd);
         tot_sym = prd(sym_left,sym_right);
         tot_spin = abs(spin_right+spin_left); //up and down spins are indistinguishable -> take absolute value
         tot_nel =n_left + n_right;
         if(tot_spin==abs(spin) && tot_sym == sym && tot_nel==nelec){
            deas_dets.push_back(dets_right[i]);
            break;
         }
         else{;}
      }
   }
   return deas_dets;
}


int main(int argc, char ** argv)
{
    std::cout.precision(12);

    try {
        if (argc != 3) {
            std::cout << "Usage: " << argv[0] << " <result.h5>" << "output.file" << std::endl;
            return 1;
        }

        std::string rfile(argv[1]);
        std::ofstream wfile;
        wfile.open(argv[2], std::ios::out);

        EntanglementData<matrix> em(rfile);
       
        storage::archive ar(rfile, "r");
        DmrgParameters parms;
        ar["/parameters"] >> parms;


        int L = parms["L"];

  //create product table
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

    // get hf_occ vector
       std::vector<int> hf_occ = parms["hf_occ"];

    // get CAS vector
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
    // print CAS vector and HF vector
        std::cout << "CAS vector: ";
        for(int i =0; i<L; i++){std::cout << casv[i] << " ";}
        std::cout <<std::endl;

        std::cout << "HF occupation  vector: ";
        for(int i =0; i<L; i++){std::cout << hf_occ[i];}
        std::cout <<std::endl;

    //get symmetry vector
       Lattice lat(parms);
       std::vector<int> sym_vec(lat.size());
       for (int i = 0; i<lat.size(); i++){
          sym_vec[i] = lat.get_prop<int>("type", i);
       }

     //implement symmetry test for a given part of the HF determinant
       int target_sym,hf_sym;
       target_sym = parms["irrep"];
       hf_sym = sym_check(hf_occ, sym_vec, prd);
       std::cout << "symmetry of HF determinant: " << hf_sym << " vs. target symmetry: " << target_sym << std::endl;      


    // get number of electrons
       int nelec = 0;
       nelec = num_el(hf_occ);
       std::cout << "number of electrons: " << nelec <<std::endl;

    //get spin
       int spin = 0;
       spin = spin_check(hf_occ);
       std::cout << "spin of HF determinant: " << spin << std::endl;

    //determine occupied orbitals in hf_det
       std::vector<std::pair<int,int> > hf_occ_orb;
       int first = 0;
       int second = 0;
       for(int i= 0; i< hf_occ.size(); i++){
          if(hf_occ[i]!=1){
             hf_occ_orb.push_back(std::make_pair(i, hf_occ[i]));
          }
       }


    //initialize stuff
       std::vector<std::vector<int> > dets_left, dets_right, deas_dets, ci_dets;
    //get possible ci_levels (include 0 for HF determinant!!!)
       int ci_level_fill[4] = {0,1,2,4};
       std::vector<int> ci_level;
       for(int i = 0; i<4; i++){ci_level.push_back(ci_level_fill[i]);}

       std::cout << "excitation levels allowed: ";
       for(int i =0; i<ci_level.size(); i++){std::cout << ci_level[i] << " ";}
       std::cout <<std::endl;

    //subtract left part, size is 2 in first microiteration
        int left = 0;//0 not possible here ... would be nice for testing

        if(left!=0){
           int L_env = L-left;
           std::vector<int> hf_left = hf_occ;
           std::vector<int> symvec_left = sym_vec;
           std::vector<int> symvec_right = sym_vec;
           hf_occ.erase(hf_occ.begin(), hf_occ.begin()+left);
           hf_left.erase(hf_left.begin()+left,hf_left.end());
           symvec_right.erase(symvec_right.begin(), symvec_right.begin()+left);
           symvec_left.erase(symvec_left.begin()+left,symvec_left.end());
 
 
           std::cout << "environment only HF occupation  vector: ";
           for(int i =0; i<L_env; i++){std::cout << ","<<hf_occ[i];}
           std::cout <<std::endl;
           std::cout << "system only HF occupation  vector: ";
           for(int i =0; i<left; i++){std::cout << ","<<hf_left[i];}
           std::cout <<std::endl;
 
           std::cout << "environment only symmetry  vector: ";
           for(int i =0; i<L_env; i++){std::cout << ","<<symvec_right[i];}
           std::cout <<std::endl;
           std::cout << "system only symmetry vector: ";
           for(int i =0; i<left; i++){std::cout << ","<<symvec_left[i];}
           std::cout <<std::endl;
 
       //create first set of determinants by DEAS procedure
           std::vector<int> casv_right, casv_left;
           casv_left = el_casv(1,left,casv);
           casv_right = el_casv(0,left,casv);
 
           int act_orb_right = casv_right[0]-left;
           int act_orb_left = casv_left[0];
 
           std::cout << "CAS vector after elimination: ";
           for(int i =0; i<L_env; i++){std::cout << ","<<casv_right[i];}
           std::cout <<std::endl;
 
           //create first four determinants of right part
           dets_right.push_back(hf_occ);
           int count = 0;
           for(int i = 1; i<5; i++){
              if(hf_occ[act_orb_right]!=i){
                 dets_right.push_back(hf_occ);
                 count++;
                 dets_right[count][act_orb_right] =i;
              }
           }
 
           //create first four determinants of left part
           dets_left.push_back(hf_left);
           count = 0;
           for(int i = 1; i<5;i++){
              if(hf_left[act_orb_left]!=i){
                 dets_left.push_back(hf_left);
                 count++;
                 dets_left[count][act_orb_left] =i;
              }
           }
 
 
       //iteratively copy and create new DEAS determinants right
           for(int pos = 0; pos<L_env-1; pos++){
              dets_right = copy_det(pos,dets_right);
              act_orb_right = casv_right[pos+1]-left;
              dets_right = deas(pos+1,act_orb_right,dets_right);
           }
           
       //iteratively copy and create new DEAS determinants right
           for(int pos = 0;pos<left-1;pos++){
              dets_left = copy_det(pos,dets_left);
              act_orb_left = casv_left[pos+1];
              dets_left = deas(pos+1,act_orb_left,dets_left);
           }
 
       //check for spin, symmetry and number of electrons
       deas_dets = det_sym_check(dets_left, dets_right, spin, target_sym, nelec, symvec_left, symvec_right, prd);
     }
     else if(left==0){ //this is only for testing purposes, usually always left>0

        int act_orb = casv[0];
        int deas_pos = 0;
        int ci_pos = 0;
        int length = 0;
        int Mmax = 2000;
        std::vector <int> erase_vec;

           //create first four determinants 
        deas_dets.push_back(hf_occ);
        int count = 0;
        for(int i = 1; i<5;i++){
           if(hf_occ[act_orb]!=i){
              deas_dets.push_back(hf_occ);
              count++;
              deas_dets[count][act_orb] =i;
           }
        }


  //main loop for deas and all checks
  //    first create new set of deas determinantes iteratively
        for(int pos = 0; pos<L-1; pos++){
           deas_dets = copy_det(pos,deas_dets);
           act_orb = casv[pos+1];
           deas_dets = deas(pos+1,act_orb,deas_dets);
        
        //copy deas_dets to ci_dets
           length = deas_dets.size();
           for(int i = deas_pos; i<length; i++){
              ci_dets.push_back(deas_dets[i]);
           }
           deas_pos = length;

        //check ci-dets vector for correct spatial symmetry, spin, number of electrons (and later ci-level)
            int nel_det,sym_det,spin_det = 0;
            erase_vec.clear();
            for(int i = ci_pos; i<ci_dets.size(); i++){
               nel_det = num_el(ci_dets[i]);
               sym_det = sym_check(ci_dets[i], sym_vec, prd); 
               spin_det = spin_check(ci_dets[i]);
        //later also check for suitable ci_level
               if(nel_det!=nelec ||sym_det!=target_sym ||spin_det!=spin ||ci_check(ci_level,hf_occ_orb,ci_dets[i])){
                  erase_vec.push_back(i);
               }
            }
  
         //actual deletion; reverse iteration is a must here
            for(int i = erase_vec.size()-1; i>-1; i--){
               ci_dets.erase(ci_dets.begin()+erase_vec[i]);
            }

            ci_pos = ci_dets.size();

         //import criterium for maximum number of determinants 
           if(ci_pos>Mmax){
              ci_dets.erase(ci_dets.begin()+Mmax,ci_dets.end());
              break;
           } 
        }

     }//end of if case left = 0


           //print determinants
        for(int i = 0; i<ci_dets.size();i++){
           for(int j = 0; j < ci_dets[i].size();j++){
      //        std::cout << ci_dets[i][j] << "  ";
              wfile << ci_dets[i][j] ;
           }
      //     std::cout << std::endl;
           wfile <<std::endl;
        }
        std::cout <<"size of ci_dets: " << ci_dets.size() << std::endl;



    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
