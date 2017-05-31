/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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

#ifndef DMRGPARAMETERS2_H
#define DMRGPARAMETERS2_H

#include "BaseParameters.h"

class DmrgParameters : public BaseParameters
{
public:
    DmrgParameters() : BaseParameters() { init_options(); }
    DmrgParameters(std::ifstream& param_file)
    : BaseParameters(param_file)
    {
        init_options();
    }
    DmrgParameters(BaseParameters const& p)
    : BaseParameters(p)
    {
        init_options();
    }
    
private:
    
    void init_options()
    {
        using parameters::value;
        
        add_option("COMPLEX", "use complex numbers", value(false));

        add_option("truncation_initial", "Initial value for the truncation error", value(1e-16));
        add_option("truncation_final", "Final value for the truncation", value(1e-16));
        
        add_option("init_bond_dimension", "", value(5));
        add_option("max_bond_dimension", "");
        add_option("sweep_bond_dimensions", "");

        add_option("optimization", "singlesite or twosite", value("twosite"));
        add_option("twosite_truncation", "`svd` on the two-site mps or `heev` on the reduced density matrix (with alpha factor)", value("svd"));
        
        add_option("alpha_initial","", value(1e-2));
        add_option("alpha_main", "", value(1e-4));
        add_option("alpha_final", "", value(1e-8));
        
        add_option("eigensolver", "", value("IETL_JCD"));
        add_option("ietl_jcd_tol", "", value(1e-8));
        add_option("ietl_jcd_gmres", "", value(0));
        add_option("ietl_jcd_maxiter", "", value(8));
        add_option("ietl_moddav_omega", "parameter omega for the modified Davidson algorithm (to compute interior eigenvalues", value(0.));
        
        add_option("nsweeps", "");
        add_option("nmainsweeps", "", 0);
        add_option("ngrowsweeps", "", 0);
        
        add_option("resultfile", "");
        add_option("chkpfile", "");
        add_option("initfile", "", value(""));
        
        add_option("donotsave", "", value(0));
        add_option("run_seconds", "", value(0));
        add_option("storagedir", "", value(""));
        add_option("use_compressed", "", value(0));
        add_option("seed", "", value(42));
        add_option("ALWAYS_MEASURE", "comma separated list of measurements", value(""));
        add_option("measure_each", "", value(1)); 
        add_option("chkp_each", "", value(1)); 
        add_option("update_each", "", value(-1));
        add_option("entanglement_spectra", "", value(0));
        add_option("conv_thresh", "energy convergence threshold to stop the simulation", value(-1));
        
        add_option("expm_method", "algorithm used for exp(-i H dt): heev (default), geev", value("heev"));
        add_option("te_type", "time evolution algorithm: nn (default), mpo", value("nn"));
        add_option("te_optim", "optimized nn time evolution", value(true));
		add_option("te_order", "trotter decomposition: second, fourth (default)", value("fourth"));
        add_option("dt", "time step in time eovlution", value(1e-3));
        add_option("nsweeps_img", "number of imaginary time steps", value(0));
        
        add_option("ngrainings", "", value(0));
        add_option("finegrain_optim", "", value(false));
        
        add_option("init_state", "", value("default"));
        add_option("init_coeff", "coefficients for coherent init", value(""));
        add_option("init_basis_state", "local indexes for basis state init", value(""));
        add_option("ci_level", "number of electrons excited from HF determinant", "1,2,3,4,5,6");

        add_option("symmetry", "none, u1 or 2u1", value("u1"));
        add_option("lattice_library", "", value("alps"));
        add_option("model_library", "", value("alps"));
        add_option("model_file", "path to model parameters", value(""));
        
        add_option("beta_mode", "", value(0));
        
        add_option("NUMBER_EIGENVALUES", "", value(1));
        add_option("n_ortho_states", "", value(0));
        add_option("ortho_states", "comma separated list of filenames", "");
        
        add_option("MEASURE[Energy]", "", value(true));
        add_option("MEASURE[EnergyVariance]", "", value(0));
        add_option("MEASURE[Entropy]", "", value(false));
        add_option("MEASURE[Renyi2]", "", value(false));
    }
    
};

class ModelParameters : public BaseParameters
{
public:
    ModelParameters() : BaseParameters() { init_options(); }
	ModelParameters(std::ifstream& param_file)
    : BaseParameters(param_file)
    {
        init_options();
    }
	ModelParameters(BaseParameters const& p)
    : BaseParameters(p)
    {
        init_options();
    }
    
    
private:
    
    void init_options()
    {
        using parameters::value;
        
        add_option("MODEL", "");
        add_option("LATTICE", "");
        add_option("alps_lattice", "");
        
        add_option("L", "");
        add_option("W", "");
        
        add_option("Jxy", "");
        add_option("Jx", "");
        add_option("Jy", "");
        add_option("Jz", "");
        add_option("Jxy1", "");
        add_option("Jz1", "");
        add_option("J1", "");
        add_option("J2", "");
        
        add_option("U", "");
        add_option("t", "");
        add_option("t1", "");
        add_option("t2", "");
        
        add_option("theta", "");
        add_option("h0", "");
        add_option("pin", "");
        add_option("y", "", value(1));
        add_option("x", "", value(1));
        add_option("z", "", value(1));
        add_option("delta", "");
        
        add_option("K0", "");
        add_option("K1", "");
        
        add_option("penalty", "");
        add_option("twist", "");
        add_option("move", "");
        
        add_option("Nmax", "");
        add_option("mu", "", value(0));
        add_option("mu0", "", value(0));
        add_option("h", "", value(1));
        add_option("c", "", value(0));
        add_option("V0", "", value(0));
        add_option("k", "", value(0));
        add_option("a", "", value(1));
        add_option("Ndiscr", "");
        add_option("omega", "", value(0.));
        add_option("shift", "", value(0.));

        add_option("V"      , "", value(0.));
        add_option("Lambda"  , "", value(0.));
        add_option("Delta"  , "", value(0.));
        add_option("Gamma1a", "", value(0.));
        add_option("Gamma1b", "", value(0.));
        add_option("Gamma2" , "", value(0.));
        add_option("nbar"   , "", value(0.));
        
        add_option("u1_total_charge", "");
        add_option("u1_total_charge1", "");
        add_option("u1_total_charge2", "");

        add_option("orbital_order", "comma separated list of orbital numbers");
        add_option("hf_occ", "comma separated list of orbital occupancies for Hartree Fock initial state");

        add_option("integral_cutoff", "Ignore electron integrals below a certain magnitude", value(1.e-20));
        
        add_option("MEASURE_CONTINUUM[Psi energy]", "", value(false));
        add_option("MEASURE_CONTINUUM[Density]", "", value(true));
        add_option("MEASURE_CONTINUUM[Local density]", "", value(true));
        add_option("MEASURE_CONTINUUM[Onebody density matrix]", "", value(false));

        add_option("MEASURE[Density]", "", value(false));
        add_option("MEASURE[Local density]", "", value(false));
        add_option("MEASURE[Local density^2]", "", value(false));
        add_option("MEASURE[Onebody density matrix]", "", value(false));
        add_option("MEASURE[Density correlation]", "", value(false));

        add_option("RUN_FINITE_T", "", value(false));
   }

};


inline DmrgParameters load_parms_and_model(std::string parms_fname, std::string model_fname="")
{
    /// Load parameters
    std::ifstream param_file(parms_fname.c_str());
    if (!param_file)
        throw std::runtime_error("Could not open parameter file.");
    DmrgParameters parms(param_file);
    
    /// Load model parameters from second input (if needed)
    std::string model_file;
    if (parms.is_set("model_file") && model_fname.empty())
        model_fname = parms["model_file"].str();
    if (!model_fname.empty()) {
        std::ifstream model_ifs(model_fname.c_str());
        if (!model_ifs)
            throw std::runtime_error("Could not open model_parms file.");
        parms << ModelParameters(model_ifs);
    }
    
    return parms;
}


#endif
