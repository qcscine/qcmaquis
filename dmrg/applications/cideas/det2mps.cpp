/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 * 
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

#include <iostream>
#include <algorithm>
#include <omp.h>
#include <boost/lexical_cast.hpp>

#include "dmrg/sim/matrix_types.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/lattice.h"
#include "alps/numeric/matrix.hpp"
#include "dmrg/models/chem/util.h"
#include "dmrg/utils/DmrgOptions.h"
#include "dmrg/utils/DmrgParameters.h"
//#include "../applications/cideas/determinant.hpp"
#include "../applications/cideas/ci_generator.cpp"
//#include "../applications/tools/deas.hpp"

#ifdef USE_AMBIENT
    #include "dmrg/block_matrix/detail/ambient.hpp"
    typedef ambient::numeric::tiles<ambient::numeric::matrix<double> > matrix;
#else
    #include "dmrg/block_matrix/detail/alps.hpp"
    typedef alps::numeric::matrix<double> matrix;
#endif


#if defined(USE_TWOU1)
typedef TwoU1 grp;
#elif defined(USE_TWOU1PG)
typedef TwoU1PG grp;
#elif defined(USE_SU2U1)
typedef SU2U1 grp;
#elif defined(USE_SU2U1PG)
typedef SU2U1PG grp;
#elif defined(USE_NONE)
typedef TrivialGroup grp;
#elif defined(USE_U1)
typedef U1 grp;
#endif

template<class SymmGroup>
Determinant<SymmGroup> str_from_det(std::vector<typename SymmGroup::charge> const &charge_vec, std::vector<Index<SymmGroup> > const &phys_dims, std::vector<typename SymmGroup::subcharge> const &site_types)
{
   Determinant<SymmGroup> det(charge_vec.size());
   for (int i = 0; i < charge_vec.size(); ++i)
   {
       if (charge_vec[i] == phys_dims[site_types[i]][0].first)
       {   
           det[i] = 4;
       }else if (charge_vec[i] == phys_dims[site_types[i]][1].first)
       {
           det[i] = 3;
       }else if (charge_vec[i] == phys_dims[site_types[i]][2].first) // singly-occ (2)
       {
           det[i] = 2;
       }else if (charge_vec[i] == phys_dims[site_types[i]][3].first)
       {
           det[i] = 1;
       }
   }
   return det; 
}

template<class SymmGroup>
std::vector<std::vector<typename SymmGroup::charge> >  get_charge_determinants(std::vector<Determinant<SymmGroup> > &det_list, std::vector<Determinant<SymmGroup> > &det_list_new,
                                                                                  std::vector<Index<SymmGroup> > const &phys_dims, std::vector<typename SymmGroup::subcharge> site_types, typename SymmGroup::charge right_end)
{
    typedef typename SymmGroup::charge charge;
    std::vector<std::vector<charge> > determinants;
    std::vector< std::vector<std::vector< charge > > > dummy_dets;
     // convert det_list to vec<vec<charge>>
    for (size_t i = 0; i < det_list.size(); ++i)
        dummy_dets.push_back(det_list[i].charge_det(phys_dims, site_types));
    std::cout << "size of dummy_dets: "<< dummy_dets.size() <<std::endl;
    int L = dummy_dets[0].size();
    #pragma omp parallel for
    for (int i = 0; i<dummy_dets.size(); ++i)
    {
        std::vector<std::vector<charge> > single_det(1);
        for (int j = 0; j < det_list[i].size(); ++j)
        {
            if (dummy_dets[i][j].size() != 1)
            {
                int times = single_det.size();
                for (int k = 0; k < times; k++)
                    single_det.push_back(single_det[k]);
                
                for (int k = 0; k < single_det.size(); ++k)
                {
                    if (k < single_det.size()/2) //works only if there are two options, like (1,1) and (1,-1)
                        single_det[k].push_back(dummy_dets[i][j][0]);
                    else
                        single_det[k].push_back(dummy_dets[i][j][1]);
                }        
            }
            else{
                for (int k = 0; k < single_det.size(); ++k)
                    single_det[k].push_back(dummy_dets[i][j][0]); 
            }
        }
        for (int m = 0; m < single_det.size(); ++m)
        {
            bool valid = true;
            charge accumulated_charge = single_det[m][0];
            for (int l = 1; l < L; ++l)
            {
                accumulated_charge = SymmGroup::fuse(accumulated_charge, single_det[m][l]);
                if (!charge_detail::physical<SymmGroup>(accumulated_charge))
                    valid = false;
            }
            charge identity(0);
            if(valid == true  && accumulated_charge == right_end && std::find(determinants.begin(), determinants.end(), single_det[m]) == determinants.end())//letzte Bedingung kann später geloescht werden
            {
		#pragma omp critical
                determinants.push_back(single_det[m]);
                Determinant<SymmGroup> det_str = str_from_det(single_det[m],phys_dims,site_types);
		#pragma omp critical
                det_list_new.push_back(det_str);
            }
        }
        single_det.clear();
    }
    #pragma omp barrier
    return determinants;
}



template<class Matrix, class SymmGroup, class=void>
struct deas_mps_init : public mps_initializer<Matrix,SymmGroup>
{
    deas_mps_init(DmrgParameters parms_,
                EntanglementData<Matrix> em_,
                std::vector<Index<SymmGroup> > const& phys_dims_,
                typename SymmGroup::charge right_end_,
                std::vector<int> const& site_type)
    : parms(parms_)
    , em(em_)
    , phys_dims(phys_dims_)
    , site_types(site_type)
    , di(parms, phys_dims_, right_end_, site_type)
    , right_end(right_end_)
    , det_list()
    , det_list_new()
    , determinants()
   {}

    typedef Lattice::pos_t pos_t;
    typedef std::size_t size_t;
    typedef typename SymmGroup::charge charge;
    typedef std::vector<Index<SymmGroup> > index_vec;
    typedef std::vector<typename SymmGroup::subcharge> site_vec;

    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
        pos_t L = mps.length();

        charge doubly_occ = phys_dims[0].begin()->first, empty = phys_dims[0].rbegin()->first;

        std::vector<int> ci_level(parms.get<std::vector<int> >("ci_level"));
        if(std::find(ci_level.begin(), ci_level.end(), 0) == ci_level.end())
           ci_level.push_back(0);

        Determinant<SymmGroup> hf_occ(parms.get<std::vector<int> >("hf_occ"));
        std::vector<std::pair<int,int> > hf_occ_orb = get_orb(hf_occ);
        int m_value = parms.get<int>("max_bond_dimension");
        if (hf_occ.size() != L)
            throw std::runtime_error("HF occupation vector length != MPS length\n");

        // initialize objects required 
        //idea: include current charges in rows_to_fill, should not affect 2U1
        std::vector<std::vector< int > > rows_to_fill;
        std::vector<std::map<charge, std::map<std::string, int> > > str_to_col_map(L);

   	//another loop has to be inserted here where determinants are generated until m is reached at some site; correct numbers have to be checked here
        std::vector<Determinant<SymmGroup> > deas_dets, new_det_list;
        std::vector<std::vector<charge> >  total_dets;
        std::vector<int> dummy_vec(L);
        std::vector<int> sum_size(L);
        bool keep_running = true;
        int det_nr = 0;
        int num_runs = std::min(L,15);
	//main loop
	for(int run = 0; run < num_runs; ++run){ 
           //generate deas determinants
           deas_dets = generate_deas(parms,em,run,deas_dets);

	   size_t loop_start = pow(4,run)-1;
	   size_t loop_end = pow(4,run+1);
           #pragma omp parallel for
           for(int i = loop_start; i < loop_end; ++i){
              if(!deas_dets[i].ci_check(ci_level,hf_occ_orb))
		 #pragma omp critical
                 new_det_list.push_back(deas_dets[i]);
           }
           #pragma omp barrier
           //convert to charge_vec -> determinants
           std::cout << "size of new_det_list: "<< new_det_list.size() << std::endl;
           determinants = get_charge_determinants(new_det_list, det_list_new, phys_dims, site_types, right_end);
           std::cout << "size of determinants: "<< determinants.size() << std::endl;

           for(int d = 0; d < determinants.size(); ++d)
           {
               rows_to_fill.push_back(dummy_vec);
               charge accumulated_charge = right_end;
               for(int s = L - 1; s > 0; --s)
               {
                   charge site_charge = determinants[d][s];
                   accumulated_charge = SymmGroup::fuse(accumulated_charge, -site_charge);
                   if(charge_detail::physical<SymmGroup>(accumulated_charge))
                   {
         
                       std::string str = det_string(s, det_list_new[det_nr]);
                       std::map<std::string, int> & str_map = str_to_col_map[s-1][accumulated_charge];
         
                       if (str_map[str])
                           rows_to_fill[det_nr][s] = str_map[str] - 1;
         
                       else
                       {
                           //get largest element in map
                           int max_value = str_map.size();
                           str_map[str] = max_value;
                           rows_to_fill[det_nr][s] = max_value - 1;
                           sum_size[s-1] += 1;
                       }
                   }
                   if(sum_size[s-1] >= m_value){
                      keep_running = false;
                      break;
                   }
               }
               det_nr ++;
               total_dets.push_back(determinants[d]);
               if(keep_running == false)
                  break;
           }
           new_det_list.clear();
           if(keep_running == false)
              break;
        }
        //this now calls a function which is part of this structure
        std::cout <<"THE NUMBER OF SUITABLE DETS: " << det_nr <<std::endl;
        std::cout << "sum of sector sizes: " ;
        for (int i = 0; i < L; ++i)
           std::cout << sum_size[i] << " ";
        std::cout << std::endl;

        init_sect(mps, str_to_col_map, true, 0); 

        //this here is absolutely necessary
        for(pos_t i = 0; i < L; ++i)
            mps[i].multiply_by_scalar(0.0);

        //fill loop
        for(int d = 0; d < total_dets.size(); ++d)
        {
            charge accumulated_charge = right_end;
            int prev_row = 0;
            for(int s = L - 1; s > 0; --s)
            {
                charge site_charge = total_dets[d][s];
                charge search_charge = SymmGroup::fuse(accumulated_charge, -site_charge);
                if (charge_detail::physical<SymmGroup>(search_charge) && mps[s].row_dim().has(search_charge))
                {
                   int nrows_fill = mps[s].row_dim().size_of_block(search_charge);
                   //get current matrix
                   size_t max_pos = mps[s].data().left_basis().position(accumulated_charge);
                   Matrix & m_insert = mps[s].data()[max_pos];
                  
                   int nrows = m_insert.num_rows(), off = 0;
                  
                   //get additional offsets for subsectors
                   if(site_charge == phys_dims[site_types[s]][1].first)
                   {
                       if(mps[s].row_dim().has(SymmGroup::fuse(accumulated_charge, -doubly_occ)))
                           off =  mps[s].row_dim().size_of_block(SymmGroup::fuse(accumulated_charge, -doubly_occ));
                   }
                   else if(site_charge == phys_dims[site_types[s]][2].first)
                   {
                       if(mps[s].row_dim().has(accumulated_charge))
                           off = nrows - nrows_fill - mps[s].row_dim().size_of_block(accumulated_charge);
                       else 
                           off = nrows - nrows_fill;
                   }  
                   else if(site_charge == phys_dims[site_types[s]][3].first){
                       off = nrows - nrows_fill;
                   }
                   //actual insertion
                   m_insert(off+rows_to_fill[d][s],prev_row) = 1;
                   prev_row = rows_to_fill[d][s];
                   accumulated_charge = search_charge;
                }
            }
       }
       std::cout << "fill worked" << std::endl;

       //first site needs to be filled as well
       for(int d = 0; d < total_dets.size(); d++){
          charge first_charge = total_dets[d][0];
          if (charge_detail::physical<SymmGroup>(first_charge))
          {
             size_t first_pos = mps[0].data().left_basis().position(first_charge);
             Matrix & m_first = mps[0].data()[first_pos];
             m_first(0,0) = 1;
          }
       }

    }//end of main initialization function 


    //function to get string of left or right part from det
    std::string det_string(int s, Determinant<SymmGroup> det){
       std::string str;
       char c;
       int L = det.size();
       if(s > L/2){
          for(int i = s; i<det.size(); i++){
             c = boost::lexical_cast<char>(det[i]);
             str.push_back(c);
          }
       }else{
          for(int i = 0; i<s; i++){
             c = boost::lexical_cast<char>(det[i]);
             str.push_back(c);
          }
       }
    return str;
    }

    //function to initalize sectors -> copied from mps-initializers

    void init_sect(MPS<Matrix, SymmGroup> & mps,
                   const std::vector<std::map<charge, std::map<std::string, int> > > & str_to_col_map, 
                   bool fillrand = true, 
                   typename Matrix::value_type val = 0)
    {
        parallel::scheduler_balanced scheduler(mps.length());
        std::size_t L = mps.length();

        std::vector<Index<SymmGroup> > allowed = allowed_sectors(site_types, phys_dims, right_end, 5);
        maquis::cout << "allowed sectors created" << std::endl;
        allowed = adapt_allowed(allowed, str_to_col_map);
        maquis::cout << "size of sectors succesfully adapted" << std::endl;

        omp_for(size_t i, parallel::range<size_t>(0,L), {
            parallel::guard proc(scheduler(i));
            mps[i] = MPSTensor<Matrix, SymmGroup>(phys_dims[site_types[i]], allowed[i], allowed[i+1], fillrand, val);
            mps[i].divide_by_scalar(mps[i].scalar_norm());
        });
    }

    std::vector<Index<SymmGroup> > adapt_allowed(std::vector<Index<SymmGroup> > allowed,
                                                     const std::vector<std::map<charge, std::map<std::string, int > > > &str_to_col_map)
    {
        std::vector<Index<SymmGroup> > adapted = allowed;
        pos_t L = str_to_col_map.size();
        for(pos_t i = 1; i < L+1; ++i)
        {
            for(typename Index<SymmGroup>::iterator it = adapted[i].begin();
                     it != adapted[i].end(); ++it)
            {
                if (str_to_col_map[i-1].count(it->first) != 0)
                    it->second = str_to_col_map[i-1].at(it->first).size();
            }
        }
        return adapted;
    }


    DmrgParameters parms;
    EntanglementData<Matrix> em;
    std::vector<Index<SymmGroup> > phys_dims;
    std::vector<typename SymmGroup::subcharge> site_types;
    default_mps_init<Matrix, SymmGroup> di;
    std::vector<Determinant<SymmGroup> > det_list;
    std::vector<Determinant<SymmGroup> > det_list_new;
    std::vector<std::vector<charge> >  determinants;
    charge right_end;
};

//parse determinants
template <class SymmGroup>
std::vector<Determinant<SymmGroup> > dets_from_file(std::string file){
    std::ifstream config_file;
    config_file.open(file.c_str());

    std::vector<Determinant<SymmGroup> > configs;

    for (std::string line; std::getline(config_file, line); ) {
        std::vector<std::string> det_coeff;
        boost::split(det_coeff, line, boost::is_any_of(" "));

        std::string det = det_coeff[0];

        Determinant<SymmGroup> tmp;
        for (std::size_t i = 0; i < det.size(); ++i) {
            int occ = boost::lexical_cast<size_t>(det[i]);
            switch(occ) {
                case 4:
                    tmp.push_back(4); // doubly occ
                    break;
                case 3:
                    tmp.push_back(3); // up
                    break;
                case 2:
                    tmp.push_back(2); // down 
                    break;
                case 1:
                    tmp.push_back(1); // empty
                    break;
            }
        }
        if(std::find(configs.begin(), configs.end(), tmp) == configs.end())
            configs.push_back(tmp);
    }
    return configs;
}




int main(int argc, char ** argv){
    try {
        if (argc != 2) {
            std::cout << "Usage: " << argv[0] << "dmrg-input file" << std::endl;
            return 1;
        }

        DmrgOptions opt(argc, argv);
        DmrgParameters parms = opt.parms;


        std::string rfile(parms.get<std::string>("init_resultfile"));
        EntanglementData<matrix> em(rfile);
       
        Lattice lat(parms);
        Model<matrix,grp> mod(lat,parms);

        /***Create TEST environment not needed in actual implementation***/
        size_t L = parms["L"];

        //create symmetry vector -> site_types
        int max_site_type = 0;
        std::vector<grp::subcharge> site_types(lat.size()), site_types_distinct;
        for (int i = 0; i<lat.size(); i++){
           site_types[i] = lat.get_prop<grp::subcharge>("type", i);
           if(std::find(site_types_distinct.begin(),site_types_distinct.end(),site_types[i]) == site_types_distinct.end())
              site_types_distinct.push_back(site_types[i]);
           max_site_type =std::max(site_types[i],max_site_type);
        }
        std::cout <<"site_types: ";
        for(int i = 0; i<site_types.size(); i++){
           std::cout << site_types[i] <<", "; 
        }
        std::cout <<std::endl;
        maquis::cout << "maximal site type: " << max_site_type << std::endl;

        std::vector<Index<grp> > phys_dims;
        for(int i = 0; i<site_types_distinct.size(); i++)
           phys_dims.push_back(mod.phys_dim(site_types_distinct[i]));

        for(int i = 0; i<max_site_type+1; i++)
           std::cout << "phys_dims["<<i<<"] = " <<phys_dims[i] <<std::endl; 


        //get right end charge
        grp::charge right_end = chem_detail::qn_helper<grp>().total_qn(parms);
        maquis::cout << "Right end: " << right_end <<std::endl;
  
        //create MPS
        MPS<matrix,grp> hf_mps(L);
        deas_mps_init<matrix,grp> hf(parms,em,phys_dims,right_end,site_types);
        hf(hf_mps); 
        maquis::cout << "MPS created" << std::endl;

        std::string chkp = parms["chkpfile"].str();
        save(chkp,hf_mps);
        storage::archive ar2(chkp+"/props.h5", "w");
        ar2["/parameters"] << parms;
//        ar["/version"] << DMRG_VERSION_STRING;
        ar2["/status/sweep"] << -1;
        ar2["/status/site"] << -1;
        
        maquis::cout << "MPS saved" << std::endl; 
  
  
        } catch (std::exception& e) {
            std::cerr << "Error:" << std::endl << e.what() << std::endl;
            return 1;
        }
}
