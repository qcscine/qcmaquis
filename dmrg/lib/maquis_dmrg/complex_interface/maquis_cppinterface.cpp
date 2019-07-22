/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2019         Leon Freitag <lefreita@ethz.ch>
*               2019         Stefan Knecht <stknecht@ethz.ch>
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
#include "maquis_cppinterface.h"
#include <memory>
#include <string>
#include <array>
#include "maquis_dmrg.h"

std::unique_ptr<maquis::DMRGInterface<std::complex<double>> > cpp_interface_ptr;
DmrgParameters cpp_parms;

    typedef std::complex<double> V;
    void qcmaquis_interface_preinit(int nel_, int L_, int irrep_,
                                    int* site_types_,
                                    double conv_thresh_, int m_,
                                    int nsweeps_, int ngrowsweeps_, int nmainsweeps_,
                                    bool meas_2rdm_, bool entropy_, bool magnetism_,
                                    const std::string& init_state_, const std::string& optimization_,
                                    int ietl_jcd_maxiter_,
                                    double ietl_jcd_tol_, double truncation_initial_,
                                    double truncation_final_, double integral_cutoff_,
                                    const std::string& twosite_truncation_, const std::string& orb_order_)

    {
        // relativistic options
        cpp_parms.set("group_id", 8);
        cpp_parms.set("type", 0);
        maquis::prepare_relativistic(cpp_parms);
        // Kramers symmetry?
        if(magnetism_)
            cpp_parms.set("MAGNETIC", 1);

        // optimization conditions
        cpp_parms.set("orbital_order", orb_order_);
        cpp_parms.set("init_state", init_state_);
        cpp_parms.set("optimization", optimization_);
        cpp_parms.set("integral_cutoff", integral_cutoff_);
        cpp_parms.set("truncation_initial", truncation_initial_);
        cpp_parms.set("truncation_final", truncation_final_);
        cpp_parms.set("ietl_jcd_tol", ietl_jcd_tol_);
        cpp_parms.set("ietl_jcd_maxiter", ietl_jcd_maxiter_);
        cpp_parms.set("conv_thresh", conv_thresh_);
        cpp_parms.set("twosite_truncation", twosite_truncation_);
        cpp_parms.set("storagedir", "./tmp/");

        // for now, making things easier
        cpp_parms.set("MEASURE[1rdm]", 1);
        if (meas_2rdm_) cpp_parms.set("MEASURE[2rdm]", 1);

        if(entropy_){
            cpp_parms.set("MEASURE_LOCAL[N]","N");
            cpp_parms.set("MEASURE_HALF_CORRELATIONS[dm]","c_dag:c");
            cpp_parms.set("MEASURE_HALF_CORRELATIONS[doccdocc]","N:N");
        };

        // bond dimension
        cpp_parms.set("max_bond_dimension", m_);

        // sweeps
        cpp_parms.set("nsweeps", nsweeps_);
        cpp_parms.set("ngrowsweeps", ngrowsweeps_);
        cpp_parms.set("nmainsweeps", nmainsweeps_);

        // lattice, #e-, symmetry
        cpp_parms.set("L", L_);
        cpp_parms.set("nelec", nel_);
        cpp_parms.set("irrep", irrep_);

        // site types
        if (site_types_ != nullptr)
        {
            std::string s;
            for (int i = 0; i < L_; i++)
                s += (std::to_string(site_types_[i])+ (i < (L_-1) ? "," : ""));
            cpp_parms.set("site_types", s);
        }
        else
        {
            std::string s;
            for (int i = 0; i < L_; i++)
                s += (std::to_string(1)+ (i < (L_-1) ? "," : ""));
            cpp_parms.set("site_types", s);
        }

        std::cout << " parms are set -> " << std::endl;
        std::cout << cpp_parms << std::endl;
    }

    void qcmaquis_interface_update_integrals(std::vector<std::vector<int>> integral_indices, std::vector<V> integral_values, int integral_size)
    {
        if (cpp_parms.is_set("integral_file")||cpp_parms.is_set("integrals"))
            throw std::runtime_error("updating integrals in the interface not supported yet in the FCIDUMP format");
        // set integrals
//      maquis::integral_map<std::complex<double>> integrals;
        maquis::integral_map<std::complex<double>, maquis::integrals::relativistic> integrals;

        for (int i = 0; i < integral_size; i++)
        {
            std::array<int, 4> idx {integral_indices[i][0], integral_indices[i][1], integral_indices[i][2],
integral_indices[i][3]};
            V value = integral_values.at(i);
            // Print indices in Normal order
            std::cout << "element " ;
            std::copy(idx.begin(), idx.end(), std::ostream_iterator<int>(std::cout,""));
            std::cout << " is " << value << std::endl;
            integrals[idx] = value;
        }

        std::cout << " now from map1 ... " <<  std::endl;
        for (auto&& t: integrals){
            std::cout << "element " << t.first[0] << t.first[1] << t.first[2] << t.first[3] << " is " << t.second << std::endl;
        }

        maquis::integral_map<std::complex<double>, maquis::integrals::relativistic> integrals2{
 { {1  ,1  ,1  ,1},{  0.43604913033777 ,-1.31352163875892e-17   }},
 { {3  ,3  ,1  ,1},{ 0.220879790916985 ,3.56833577884087e-17    }},
 { {1  ,3  ,3  ,1},{ 0.000451512611023166, 6.15109985306372e-20 }},
 { {3  ,1  ,1  ,3},{ 0.000451512611023166, 6.15109985306372e-20 }},
 { {1  ,1  ,3  ,3},{ 0.220879790916985 ,3.56833577884087e-17    }},
 { {3  ,3  ,3  ,3},{ 0.228201307078244 ,-1.04540924734875e-17   }},
 { {1  ,3  ,3  ,2},{ 0.00072840466362894, 0.00153276774662158   }},
 { {1  ,3  ,1  ,4},{-0.000728404663628938, -0.00153276774662158 }},
 { {3  ,2  ,1  ,3},{ 0.00072840466362894, 0.00153276774662158   }},
 { {1  ,4  ,1  ,3},{-0.000728404663628938, -0.00153276774662158 }},
 { {3  ,1  ,4  ,1},{-0.000728404663628938, 0.00153276774662158  }},
 { {3  ,1  ,2  ,3},{ 0.00072840466362894, -0.00153276774662158  }},
 { {1  ,1  ,2  ,2},{ 0.43604913033777, 1.0696431132045e-17      }},
 { {3  ,1  ,4  ,2},{ 0.000451512611023166, -1.35858035557269e-20}},
 { {3  ,3  ,2  ,2},{ 0.220879790916985, 5.25693153440665e-17    }},
 { {1  ,1  ,4  ,4},{ 0.220879790916985, 1.13059686387213e-17    }},
 { {1  ,3  ,2  ,4},{ 0.000451512611023166, 1.18293035589619e-19 }},
 { {3  ,3  ,4  ,4},{ 0.228201307078244 ,-8.34489724655587e-17   }},
 { {3  ,2  ,4  ,1},{-0.00637911362451499 ,-5.80942808457955e-19 }},
 { {1  ,4  ,4  ,1},{ 0.00637911362451499 ,-2.54899693461451e-18 }},
 { {3  ,2  ,2  ,3},{ 0.00637911362451499 ,-1.00834143850827e-18 }},
 { {1  ,4  ,2  ,3},{-0.00637911362451499 ,4.03250794224186e-18  }},
 { {3  ,2  ,4  ,2},{ 0.00072840466362894 ,0.00153276774662158   }},
 { {1  ,4  ,4  ,2},{-0.000728404663628939, -0.00153276774662158 }},
 { {4  ,1  ,3  ,1},{-0.000728404663628938, 0.00153276774662158  }},
 { {2  ,3  ,3  ,1},{ 0.00072840466362894 ,-0.00153276774662158  }},
 { {4  ,1  ,3  ,2},{-0.00637911362451499 ,5.80942808457955e-19  }},
 { {2  ,3  ,3  ,2},{ 0.00637911362451499 ,2.54899693461451e-18  }},
 { {4  ,1  ,1  ,4},{ 0.00637911362451499 ,1.00834143850827e-18  }},
 { {2  ,3  ,1  ,4},{-0.00637911362451499 ,-4.03250794224186e-18 }},
 { {2  ,2  ,1  ,1},{ 0.43604913033777 ,-1.0696431132045e-17     }},
 { {4  ,2  ,3  ,1},{ 0.000451512611023166, 1.35858035557269e-20 }},
 { {4  ,4  ,1  ,1},{ 0.220879790916985, -5.25693153440665e-17   }},
 { {2  ,2  ,3  ,3},{ 0.220879790916985, -1.13059686387213e-17   }},
 { {2  ,4  ,1  ,3},{ 0.000451512611023166, -1.18293035589619e-19}},
 { {4  ,4  ,3  ,3},{ 0.228201307078244 ,8.34489724655587e-17    }},
 { {4  ,2  ,3  ,2},{ 0.00072840466362894 ,0.00153276774662158   }},
 { {4  ,2  ,1  ,4},{-0.000728404663628939, -0.00153276774662158 }},
 { {4  ,1  ,2  ,4},{-0.000728404663628939, 0.00153276774662158  }},
 { {2  ,3  ,2  ,4},{ 0.00072840466362894 ,-0.00153276774662158  }},
 { {2  ,4  ,4  ,1},{-0.000728404663628939, 0.00153276774662158  }},
 { {2  ,4  ,2  ,3},{ 0.00072840466362894, -0.00153276774662158  }},
 { {2  ,2  ,2  ,2},{ 0.43604913033777 ,3.45280786516792e-17     }},
 { {4  ,4  ,2  ,2},{ 0.220879790916985, 2.81919261943354e-17    }},
 { {2  ,4  ,4  ,2},{ 0.000451512611023166, 4.53071983760033e-20 }},
 { {4  ,2  ,2  ,4},{ 0.000451512611023166, 4.53071983760033e-20 }},
 { {2  ,2  ,4  ,4},{ 0.220879790916985, 2.81919261943354e-17    }},
 { {4  ,4  ,4  ,4},{ 0.228201307078244, -1.5644385245702e-16    }},
 { {1  ,1  ,0  ,0},{-1.47144863954766 ,1.13465992495957e-16     }},
 { {3  ,3  ,0  ,0},{-0.594462058458452, -1.27921644083709e-16   }},
 { {2  ,2  ,0  ,0},{-1.47144863954766 ,-2.51035895744511e-16    }},
 { {4  ,4  ,0  ,0},{-0.594462058458451, -3.80166790074138e-17   }},
 { {0  ,0  ,0  ,0},{-97.4693010120879 ,  0.0                    }}
        };

        std::cout << " now from map2 ... " <<  std::endl;
        for (auto&& t: integrals2){
            std::cout << "element " << t.first[0] << t.first[1] << t.first[2] << t.first[3] << " is " << t.second << std::endl;
        }

        if (!cpp_parms.is_set("integrals_binary"))
            cpp_parms.set("integrals_binary", maquis::serialize(integrals));

        // Call an integral update only if the interface has been initialised
        if (cpp_interface_ptr){
            cpp_interface_ptr->update_integrals(integrals);
        }

    }

    void qcmaquis_interface_set_nsweeps(int nsweeps)
    {
        cpp_parms.set("nsweeps", nsweeps);
    }

    // Start a new simulation with stored parameters
    void qcmaquis_interface_reset()
    {
        cpp_interface_ptr.reset(new maquis::DMRGInterface<std::complex<double>>(cpp_parms));
    }

    void qcmaquis_interface_optimize()
    {
        cpp_interface_ptr->optimize();
    }
    std::complex<double> qcmaquis_interface_get_energy()
    {
        return cpp_interface_ptr->energy();
    }

    void qcmaquis_interface_set_state(int state_)
    {
        std::string str;
        for (int i = 0; i < state_; i++)
            str += "chk." + std::to_string(state_-1) + ".h5" + ((i < state_-1) ? ";" : "") ;

        cpp_parms.set("ortho_states", str);
        cpp_parms.set("n_ortho_states", state_);
        cpp_parms.set("chkpfile", "chk." + std::to_string(state_) + ".h5");
        cpp_parms.set("resultfile", "res." + std::to_string(state_) + ".h5");

        qcmaquis_interface_reset();
    }

    void qcmaquis_interface_get_1rdm(std::vector<std::vector<int>> *indices, std::vector<V> *values)
    {
        const typename maquis::DMRGInterface<std::complex<double>>::meas_with_results_type& meas = cpp_interface_ptr->onerdm();
        for (int i = 0; i < meas.first.size(); i++)
        {
            int pos_[2] = {meas.first[i][0], meas.first[i][1]};
            std::vector<int> positions(pos_, pos_ + 2);
            values->push_back(meas.second[i]);
            indices->push_back(positions);
        }
    }

    // hooray for copy-paste
    void qcmaquis_interface_get_2rdm(int* indices, V* values, int size)
    {
        const typename maquis::DMRGInterface<std::complex<double>>::meas_with_results_type& meas = cpp_interface_ptr->twordm();

        assert(size == meas.first.size());
        assert(size == meas.second.size());
        for (int i = 0; i < meas.first.size(); i++)
        {
            values[i] = meas.second[i];
            indices[4*i] = meas.first[i][0];
            indices[4*i+1] = meas.first[i][1];
            indices[4*i+2] = meas.first[i][2];
            indices[4*i+3] = meas.first[i][3];

        }
    }

    void qcmaquis_interface_get_iteration_results(int* nsweeps, std::size_t* m, V* truncated_weight, V* truncated_fraction, V* smallest_ev)
    {
        // Get iteration results from the last sweep
        results_collector& iter = cpp_interface_ptr->get_iteration_results();

        // iter contains results, one element per microiteration
        const std::vector<boost::any>& m_vec = iter["BondDimension"].get();
        const std::vector<boost::any>& tw_vec = iter["TruncatedWeight"].get();
        const std::vector<boost::any>& tf_vec = iter["TruncatedFraction"].get();
        const std::vector<boost::any>& ev_vec = iter["SmallestEV"].get();

        // We return the sum of these values for the last sweep
        // this should be done with transform_reduce

        *m = 0; for (auto&& m_ : m_vec) *m += boost::any_cast<std::size_t>(m_);
        *truncated_weight = 0; for (auto&& tw_ : tw_vec) *truncated_weight += boost::any_cast<V>(tw_);
        *truncated_fraction = 0; for (auto&& tf_ : tf_vec) *truncated_fraction += boost::any_cast<V>(tf_);
        *smallest_ev = 0; for (auto&& ev_ : ev_vec) *smallest_ev += boost::any_cast<V>(ev_);
        *nsweeps = cpp_interface_ptr->get_last_sweep();
    }
