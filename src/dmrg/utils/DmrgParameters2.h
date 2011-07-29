/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef DMRGPARAMETERS2_H
#define DMRGPARAMETERS2_H

#include "BaseParameters.h"

class DmrgParameters : public BaseParameters
{
public:
    DmrgParameters(std::ifstream& param_file)
    : BaseParameters(param_file)
    {
        add_option("truncation_initial", "", "Initial value for the truncation error");
        add_option("truncation_final", "", "Final value for the truncation");
        
        add_option("init_bond_dimension", "5", "");
        add_option("max_bond_dimension", "", "");
        add_option("sweep_bond_dimensions", "", "");
        
        add_option("alpha_initial", "1e-2", "");
        add_option("alpha_main", "1e-6", "");
        add_option("alpha_final", "0", "");
        
        add_option("eigensolver", "IETL_JCD", "");
        add_option("ietl_jcd_tol", "1e-8", "");
        add_option("ietl_jcd_gmres", "5", "");
        add_option("ietl_jcd_maxiter", "100", "");
        
        add_option("nsweeps", "", "");
        add_option("nmainsweeps", "", "");
        add_option("ngrowsweeps", "", "");
        
        add_option("resultfile", "", "", REQUIRED);
        add_option("chkpfile", "", "", REQUIRED);
        add_option("initfile", "", "", EMPTY_VALUE);
        
        add_option("donotsave", "0", "");
        add_option("run_seconds", "0", "");
        add_option("storagedir", "", "", EMPTY_VALUE);
        add_option("use_compressed", "0", "");
        add_option("calc_h2", "0", "");
        add_option("seed", "42", "");
        add_option("always_measure", "", "comma separated list of measurements", EMPTY_VALUE);
        add_option("measure_each", "1", ""); 
        
        add_option("dt", "1e-3", "");
        add_option("nsweeps_img", "0", "");
        
        add_option("ngrainings", "0", "");
        
        add_option("init_state", "default", "");
        
        add_option("model_library", "alps", "");
        
    }
};

class ModelParameters : public BaseParameters
{
public:
    ModelParameters() : BaseParameters() {}
	ModelParameters(std::ifstream& param_file)
    : BaseParameters(param_file)
    {
        
        add_option("model", "", "");
        add_option("lattice", "", "");
        add_option("alps_lattice", "", "");
        
        add_option("L", "", "");
        add_option("W", "", "");
        
        add_option("Jxy", "", "");
        add_option("Jz", "", "");
        add_option("Jxy1", "", "");
        add_option("Jz1", "", "");
        
        add_option("U", "", "");
        add_option("t", "", "");
        add_option("t1", "", "");
        add_option("t2", "", "");
        
        add_option("theta", "", "");
        add_option("h0", "", "");
        add_option("pin", "", "");
        
        add_option("K0", "", "");
        add_option("K1", "", "");
        
        add_option("penalty", "", "");
        add_option("twist", "", "");
        add_option("move", "", "");
        
        add_option("Nmax", "", "");
        add_option("mu", "", "");
        add_option("h", "", "");
        add_option("c", "", "");
        add_option("V0", "", "");
        add_option("k", "", "");
        add_option("a", "", "");
        add_option("Ndiscr", "", "");
        
        add_option("u1_total_charge", "", "");
        add_option("u1_total_charge1", "", "");
        add_option("u1_total_charge2", "", "");
    }
};

#endif
