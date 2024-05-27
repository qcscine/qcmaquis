/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#ifndef DMRGPARAMETERS2_H
#define DMRGPARAMETERS2_H

#include <fstream>
#include <stdexcept>
#include <string>

#include "BaseParameters.h"

class DmrgParameters : public BaseParameters {
 public:
  DmrgParameters() { init_options(); }
  DmrgParameters(std::ifstream& param_file) : BaseParameters(param_file) {
    init_options();
  }
  DmrgParameters(const BaseParameters& p) : BaseParameters(p) {
    init_options();
  }

 private:
  void init_options() {
    using parameters::value;

    // General settings
    add_option(
        "simulation_type",
        "Simulation type can be optimize, evolve,"
        "ipi,"
        "feast, transcorrelated"
    );
    add_option("verbose", "Verbosity level of the DMRG calculation", value(0));
    add_option(
        "seed",
        "Seed for all random number generators, for instance in the random MPS "
        "initalization, the SRCAS sampling, etc.",
        value(42)
    );
    add_option("COMPLEX", "use complex numbers", value(false));
    add_option("MAGNETIC", "external magnetic field applied", value(false));

    // MPS bond dimension settings
    add_option(
        "init_bond_dimension",
        "Bond dimension m with which the MPS is initialized for init types "
        "containing *const, *default",
        value(5)
    );
    add_option(
        "max_bond_dimension", "Maximum value m of the MPS bond dimension"
    );
    add_option(
        "sweep_bond_dimensions",
        "Comma-seperated list of n bond dimensions to be used in the n first "
        "sweeps"
    );

    // Settings for further truncating than the maximum bond dimension
    add_option(
        "truncation_initial",
        "Initial value for the truncation error, interpolated during "
        "ngrowsweeps",
        value(1e-16)
    );
    add_option(
        "truncation_main",
        "Intermediate value for the truncation error during the nmainsweeps",
        value(1e-16)
    );
    add_option(
        "truncation_final",
        "Final value for the truncation error during the rest of the "
        "nsweeps-nmainsweeps",
        value(1e-16)
    );

    // Settings related to the MPS optimization algorithm
    add_option(
        "optimization", "ALS sweeping modality: singlesite or twosite",
        value("twosite")
    );
    add_option(
        "twosite_truncation",
        "`svd` on the two-site mps or `heev` on the reduced density matrix "
        "(with alpha factor)",
        value("svd")
    );

    // Number of sweeps for different calculation phases
    add_option("nsweeps", "Overall number of sweeps of the calculation", 10);
    add_option(
        "ngrowsweeps",
        "Number of the grow sweeps (used for the initial truncation and noise "
        "parameters)",
        2
    );
    add_option(
        "nmainsweeps",
        "Number of main sweeps (used for the main truncation and noise "
        "parameters)",
        5
    );

    // Setting to terminate DMRG before nsweeps are completed
    add_option(
        "conv_thresh",
        "Energy convergence threshold to stop the simulation (use same units "
        "as integral file is provided in)",
        value(-1)
    );
    add_option(
        "run_seconds",
        "Maximum time in seconds to run the calculation for, activate by "
        "setting a limit > 0",
        value(0)
    );

    // Noise to be added to the MPS during the optimization/evolution
    add_option(
        "alpha_initial",
        "Scaling factor of the noise added to perturb the MPS update during "
        "the optimization/evolution during ngrowsweeps",
        value(1e-2)
    );
    add_option(
        "alpha_main",
        "Scaling factor of the noise added to perturb the MPS update during "
        "the optimization/evolution during nmainweeps",
        value(1e-4)
    );
    add_option(
        "alpha_final",
        "Scaling factor of the noise added to perturb the MPS update during "
        "the optimization/evolution during the remainder of nsweeps",
        value(1e-8)
    );

    // MPS initialization settings
    add_option(
        "init_type",
        "Initialization type of the initial guess MPS. The default is random, "
        "also possible are const, basis_state*/hf, etc.",
        value("default")
    );
    add_option(
        "init_state_type",
        "Type of provided states, can be [csf] for SU2 calcultions, or [det] "
        "by default",
        value("det")
    );
    add_option("init_ckpt", "Checkpoint file name to restore MPS");
    add_option(
        "init_coeffs", "Comma-separated list of coefficients for coherent init",
        value("")
    );
    add_option(
        "init_file", "Filename(s) for coherent initalization", value("")
    );
    add_option(
        "init_basis_state",
        "Local indices (ONV) for basis state init (used if [init_type] is "
        "[basis_state_generic])",
        value("")
    );
    add_option(
        "init_space",
        "Occupation up to which the initial guess MPS should be populated "
        "(used if [init_type] is [basis_state_generic_*])",
        value("")
    );
    add_option(
        "ci_level", "Number of electrons excited from HF determinant",
        "1,2,3,4,5,6"
    );
    add_option(
        "hf_occ",
        "Comma separated list of orbital occupancies for Hartree Fock initial "
        "state"
    );

    // Settings for lattice
    add_option(
        "LATTICE",
        "Definition of the mapping from single-particle basis onto the DMRG "
        "lattice",
        value("orbitals")
    );
    add_option("L", "Length of the DMRG lattice");
    add_option("lattice_library", "", value("coded"));
    add_option("CONSERVED_QUANTUMNUMBERS", "", value("Nup,Ndown"));
    add_option("orbital_order", "Comma separated list of orbital numbers");

    // Settings for model
    add_option("MODEL", "Type of Hamiltonian", value("quantum_chemistry"));
    add_option(
        "symmetry",
        "Total symmetry group of the MPS, e.g. 2u1,2u1pg,su2u1,su2u1pg",
        value("su2u1pg")
    );
    add_option("model_library", "", value("coded"));
    add_option("model_file", "path to model parameters", value(""));

    // Settings for integral read-in
    add_option(
        "integral_file",
        "Path to model parameters, e.g. FCIDUMP-style integral file",
        value("FCIDUMP")
    );
    add_option(
        "integral_cutoff", "Ignore integrals below a certain magnitude",
        value(0)
    );
    add_option("beta_mode", "", value(0));

    // Jacobi-Davidson-related options
    add_option("eigensolver", "Eigensolver to solve site problem", value("IETL_JCD"));
    add_option(
        "ietl_jcd_tol", "Convergence threshold for JCD at each site",
        value(1e-8)
    );
    add_option("ietl_jcd_gmres", "", value(0));
    add_option(
        "ietl_jcd_maxiter", "Maximum number of iterations in JCD at each site",
        value(10)
    );

    // Storing/updating settings
    add_option(
        "resultfile", "Path and name of file in which to store the results"
    );
    add_option("chkpfile", "Path and name of folder in which to store the MPS");
    add_option(
        "chkp_each", "Update the checkpoint every 2*chkp_each sweeps", value(1)
    );
    add_option(
        "storagedir", "Scratch directory for temporary files", value("")
    );
    add_option("donotsave", "", value(0));

    add_option("use_compressed", "", value(0));
    add_option("entanglement_spectra", "", value(0));

    // Parameters related to the Fermi-Hubbard model
    add_option(
        "U_FermiHubbard", "Potential term entering the Fermi-Hubbard model",
        value(0.)
    );
    add_option(
        "t_FermiHubbard", "Hopping term entering the Fermi-Hubbard model",
        value(1.)
    );
    add_option(
        "tx_FermiHubbard",
        "Hopping term for the x dimension entering the two-dimensional "
        "Fermi-Hubbard model",
        value(1.)
    );
    add_option(
        "ty_FermiHubbard",
        "Hopping term for the y dimension entering the two-dimensional "
        "Fermi-Hubbard model",
        value(1.)
    );
    add_option("width_FermiHubbard", "Width of the Fermi-Hubbard lattice");
    add_option("height_FermiHubbard", "Height of the Fermi-Hubbard lattice");

    // Parameters related to the transcorrelation
    add_option(
        "transcorrelated_hamiltonian",
        "If yes, transcorrelates (if possible) the Hamiltonian", value("no")
    );
    add_option(
        "J_Transcorrelated",
        "Scalar factor for the Gutzwiller correlator (used only for "
        "Fermi--Hubbard Hamiltonians)",
        value(0.)
    );
    add_option(
        "transcorrelated_3body",
        "If no, does not add the three-body part of the transcorrelated "
        "Hamiltonian",
        value("yes")
    );
    add_option(
        "transcorrelated_3body_max_coupling",
        "Maximum many-body coupling order for the transcorrelated Hamiltonian",
        value(6)
    );
    add_option(
        "transcorrelated_nsweeps_TI",
        "Number of preliminary TI-DMRG sweeps for a tcDMRG calculation",
        value(5)
    );
    add_option(
        "transcorrelated_nsweeps_TC",
        "Number of iTD-DMRG sweeps for a tcDMRG calculation", value(20)
    );
    add_option(
        "transcorrelated_integral_file",
        "Name of the file storing the transcorrelated integrals"
    );

    // Choice for the format of the Hamiltonian
    add_option(
        "quantum_computing_format",
        "If yes, assumes that the Hamiltonian is in the quantum computing "
        "format",
        value("no")
    );
    add_option(
        "transcorrelated_quantum_computing_format",
        "If yes, assumes that the transcorrelated Hamiltonian is in the "
        "quantum computing format",
        value("no")
    );
    add_option("ngrainings", "", value(0));
    add_option("finegrain_optim", "", value(false));

    // Measurement related settings
    add_option(
        "measure_each",
        "Compute the expectation values every 2*measure_each sweeps", value(1)
    );
    add_option(
        "parallelize_measurements",
        "If true, measurements will be executed in parallel if OPENMP is "
        "enabled",
        value(false)
    );
    add_option("MEASURE[Energy]", "", value(true));
    add_option("MEASURE[EnergyVariance]", "", value(false));
    add_option("MEASURE[Entropy]", "", value(false));
    add_option(
        "MEASURE[ChemEntropy]",
        "Evaluate all expectation valus required for a mututal information "
        "calculation. Only available for 2u1(pg)",
        value(false)
    );
    add_option("MEASURE[Renyi2]", "", value(false));
    add_option("MEASURE[Autocorrelation]", "", value(0));
    add_option(
        "ALWAYS_MEASURE", "comma separated list of measurements", value("")
    );

    // Electronic-structure calculations parameters
    add_option(
        "irrep",
        "Index of the irreducible representation associated with the wave "
        "function",
        value(0)
    );
    add_option(
        "u1_total_charge1", "Number of up electrons in 2u1(pg) calculation"
    );
    add_option(
        "u1_total_charge2", "Number of down electrons in 2u1(pg) calculation"
    );
    add_option(
        "spin",
        "Total spin of target state in su2u1(pg) calculation as int 2*S, so "
        "use 0 for singlet, 1 for doublet, and 2 for triplet"
    );
    add_option("nelec", "Total number of electrons in su2u1(pg) calculation");

    // Canonical vDMRG-related parameters for the Watson Hamiltonian
    add_option(
        "watson_max_coupling",
        "Maximum many-body coupling to be included in the definition of the "
        "PES in canonical quantization",
        6
    );
    add_option(
        "watson_max_coupling_input",
        "Maximum many-body coupling allowed to appear in the input file",
        value(ORDER_NONE)
    );
    add_option(
        "watson_coordinate_type",
        "Type of coordinate used for the Hamiltonian definition",
        value("cartesian")
    );
    add_option(
        "Nmax",
        "Maximum excitation degree for each mode in the canonical "
        "quantization-based vDMRG, either single integer or comma separated "
        "list with the number of basis functions per mode",
        value(6)
    );

    // n-mode vDMRG related parameters
    add_option(
        "nmode_dumpIntegral",
        "If == yes, store the integrals in the result file", value("no")
    );
    add_option(
        "nmode_num_modes",
        "Number of modes of the Hamiltonian expressed in the n-mode "
        "representation"
    );
    add_option(
        "nmode_max_coupling",
        "Maximum many-body coupling order in the potential operator", value(3)
    );
    add_option(
        "nmode_num_basis",
        "Comma separated list with the number of basis functions per mode"
    );

    // Pre-BO
    add_option(
        "PreBO_MaxBondDimVector",
        "Give a maximum bond dimension for each particle type."
    );
    add_option("PreBO_ParticleTypeVector", "Number of particles per type");
    add_option("PreBO_FermionOrBosonVector", "1 if Fermion, 0 if Boson");
    add_option("PreBO_OrbitalVector", "Number of orbitals for each type");
    add_option(
        "PreBO_InitialStateVector",
        "Number of particles in alpha/beta state for each particle."
    );

    // Vibronic settings
    add_option(
        "vibronic_J_coupling",
        "Coulomb coupling defining the excitonic Hamiltonian", value(0.)
    );
    add_option(
        "vibronic_J_excitation",
        "Scaling factor for the single-state component of the Hamiltonian",
        value(0.)
    );
    add_option(
        "vibronic_J_interaction_type",
        "Type of Coulomb coupling. Allowed values: nn (nearest-neighbour) or "
        "all (all molecules are coupled)",
        value("nn")
    );
    add_option(
        "vibronic_num_elestates",
        "Number of the electronic states entering the vibronic Hamiltonian"
    );
    add_option(
        "vibronic_num_vibmodes",
        "Number of vibrational modes per molecule included in the vibronic "
        "Hamiltonian"
    );
    add_option(
        "vibronic_sorting",
        "Method to map the sites onto the DMRG lattice. Can be either equal to "
        "'firstele', or to 'intertwined'",
        "firstele"
    );
    add_option(
        "vibronic_num_molecules",
        "Number of molecules composing the molecular aggregate"
    );
    add_option("vibronic_num_excitons", "Number of excitons", value(1));
    add_option(
        "vibronic_num_connectingmodes",
        "Number of vibrational modes connecting two monomers", value(0)
    );
    add_option(
        "vibronic_max_coupling",
        "Maximum power of operators in the Hamiltonian", value(2)
    );

    add_option(
        "vibronic_max_coupling_nmode",
        "Maximum many body coupling of vibrational modes", value(1)
    );

    // TD-related parameters
    add_option(
        "propagator_accuracy",
        "Accuracy of the iterative approximation of the time-evolution "
        "operator",
        value(1.0E-10)
    );
    add_option(
        "propagator_maxiter",
        "Maximum number of iterations of the iterative approximation of the "
        "propagator",
        value(10)
    );
    add_option("time_step", "Time-step for the TD-DMRG propagation");
    add_option(
        "hamiltonian_units", "Units in which the SQ Hamiltonian is expressed",
        value("Hartree")
    );
    add_option("time_units", "Units in which the time-step is expressed");
    add_option(
        "imaginary_time", "Equal to yes for iTD-DMRG, no for TD-DMRG",
        value("no")
    );
    add_option(
        "TD_backpropagation",
        "Equal to yes if the back-propagation step should be done, no "
        "otherwise",
        value("yes")
    );
    add_option(
        "TD_noise",
        "If set to yes, activates the noise. By default this option is "
        "deactivated.",
        value("no")
    );

    // Tools
    add_option(
        "determinant_file",
        "File where the determinants are stored. Used in the tools."
    );
    add_option(
        "determinant_threshold", "Threshold for the determinant-related tool",
        0.
    );

    // SRCAS settings
    add_option(
        "srcas_veryVerbose",
        "If == yes, SRCAS will print output for each new sample", value("no")
    );
    add_option(
        "srcas_targetCompleteness",
        "Desired completness for SRCAS to terminate sampling", value(0.99)
    );
    add_option(
        "srcas_maxNumIterations",
        "Maximum number of macroiterations until SRCAS sampling is terminated",
        value(10)
    );
    add_option(
        "srcas_numSamples", "Number of samples in each SRCAS macroiteration",
        value(10000)
    );
    add_option(
        "srcas_overlapThreshold",
        "Threshold for overlap coefficient for being added to the determinant "
        "list",
        value(0.001)
    );
    add_option(
        "srcas_samplingSpeed",
        "Controls the speed of sampling across the Hilbert space, as this "
        "parameter determines the number of simultaneously accepted changes.",
        value(0.333)
    );

    // Excited states calculation with ORTHO
    add_option(
        "n_ortho_states",
        "Number of lower-lying MPS to orthogonalize the current calculation to",
        value(0)
    );
    add_option(
        "ortho_states",
        "Comma-separated list of checkpoint names to which to orthogonalize to"
    );
    add_option(
        "completeness_threshold",
        "Threshold for the completeness of the Gram-Schmidt orthogonalization",
        value(1e-10)
    );

    // Solution of linear systems
    add_option(
        "linsystem_precond",
        "If set to [diagonal], applies a diagonal preconditioner to the linear "
        "system solver",
        value("no")
    );
    add_option(
        "linsystem_init",
        "Initial guess for the Krylov basis (either [zero] for a zero MPS or "
        "[last] for the rhs)",
        value("last")
    );
    add_option(
        "linsystem_max_it",
        "Maximum number of times the iterative linear system solver is "
        "repeated (if >1, does basically restarted GMRES",
        value(1)
    );
    add_option(
        "linsystem_tol",
        "Threshold for the error - if the error falls below [linsystem_tol], "
        "the iterative procedure is stopped",
        value(1.0E-5)
    );
    add_option(
        "linsystem_krylov_dim",
        "Maximum dimension of the Krylov subspace for the iterative solution "
        "of the linear system",
        value(50)
    );
    add_option(
        "linsystem_solver",
        "Algorithm to be used to solve the linear system (possible values "
        "[GMRES] and [MINRES])",
        value("GMRES")
    );
    add_option(
        "linsystem_exact_error",
        "If yes, calculates the exact error associated with the solution to "
        "the linear system",
        value("no")
    );
    add_option(
        "linsystem_verbose", "If yes, prints detail of the sweep", value("yes")
    );
    add_option(
        "linsystem_noise",
        "If set to yes, activates the noise. By default this option is "
        "deactivated.",
        value("yes")
    );

    // Parameters related to DMRG[IPI]
    add_option(
        "ipi_sweep_overlap_threshold",
        "If the overlap between the MPSs calculated at two consecutive "
        "iterations is below this threshold, stops",
        value(1.0E-10)
    );
    add_option(
        "ipi_sweep_energy_threshold",
        "Threshold on the energy difference below which the IPI iterations are "
        "defined as converged",
        value(1.0E-10)
    );
    add_option("ipi_shift", "Shift parameter for the DMRG[IPI] algorithm");
    add_option(
        "ipi_iterations",
        "Maximum number of macroiterations for the DMRG[IPI] calculation"
    );

    // Parameters related to the DMRG[FEAST] algorithm
    add_option(
        "feast_num_states", "Number of states to be targeted by DMRG[FEAST]",
        value(1)
    );
    add_option(
        "feast_max_iter", "Maximum number of FEAST iterations", value(1)
    );
    add_option("feast_emin", "Lower bound for the complex contour integration");
    add_option("feast_emax", "Upper bound for the complex contour integration");
    add_option(
        "feast_num_points",
        "Number of quadrature points for approximating the integral", value(8)
    );
    add_option(
        "feast_integral_type",
        "`full' for the complete integration, `half' for the semicircle "
        "integration",
        value("full")
    );
    add_option(
        "feast_truncation_type",
        "`each' for truncating the MPS after each sum, `end' if the truncation "
        "must be done only at the end",
        value("end")
    );
    add_option(
        "feast_overlap_convergence_threshold",
        "Overlap threshold to assess the convergence of the FEAST procedure",
        value(1.0E-5)
    );
    add_option(
        "feast_energy_convergence_threshold",
        "Energy threshold to assess the convergence of the FEAST procedure",
        value(1.0E-5)
    );
    add_option(
        "feast_calculate_standard_deviation",
        "If yes, calculates the standard deviation associated with each FEAST "
        "state.",
        value("no")
    );
    add_option(
        "feast_verbose", "If yes, activate verbose output for FEAST",
        value("no")
    );
    add_option(
        "feast_standard_deviation_threshold",
        "If set, uses this threshold to accept/reject an eigenpair"
    );
    add_option(
        "feast_print_timings",
        "If equal to yes, prints timings spent in each step", value("yes")
    );
  }
};

class ModelParameters : public BaseParameters {
 public:
  ModelParameters() : BaseParameters() { init_options(); }

  explicit ModelParameters(std::ifstream& param_file)
      : BaseParameters(param_file) {
    init_options();
  }

  explicit ModelParameters(const BaseParameters& p) : BaseParameters(p) {
    init_options();
  }

 private:
  void init_options() {
    using parameters::value;

    add_option("MODEL", "quantum_chemistry");
    add_option("LATTICE", "orbitals");
    add_option("CONSERVED_QUANTUMNUMBERS", "Nup,Ndown");
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

    add_option("V", "", value(0.));
    add_option("Lambda", "", value(0.));
    add_option("Delta", "", value(0.));
    add_option("Gamma1a", "", value(0.));
    add_option("Gamma1b", "", value(0.));
    add_option("Gamma2", "", value(0.));
    add_option("nbar", "", value(0.));

    add_option("u1_total_charge1", "");
    add_option("u1_total_charge2", "");

    add_option("orbital_order", "comma separated list of orbital numbers");
    add_option(
        "hf_occ",
        "comma separated list of orbital occupancies for Hartree Fock initial "
        "state"
    );

    add_option("MEASURE[Density]", "", value(false));
    add_option("MEASURE[Local density]", "", value(false));
    add_option("MEASURE[Local density^2]", "", value(false));
    add_option("MEASURE[Onebody density matrix]", "", value(false));
    add_option("MEASURE[Density correlation]", "", value(false));

    add_option("RUN_FINITE_T", "", value(false));
  }
};

inline DmrgParameters load_parms_and_model(
    const std::string& parms_fname, std::string model_fname = ""
) {
  /// Load parameters
  std::ifstream param_file(parms_fname.c_str());
  if (!param_file) {
    throw std::runtime_error("Could not open parameter file " + parms_fname);
  }
  DmrgParameters parms(param_file);

  /// Load model parameters from second input (if needed)
  std::string model_file;
  if (parms.is_set("model_file") && model_fname.empty()) {
    model_fname = parms["model_file"].str();
  }
  if (!model_fname.empty()) {
    std::ifstream model_ifs(model_fname.c_str());
    if (!model_ifs) {
      throw std::runtime_error("Could not open model_file.");
    }
    parms << ModelParameters(model_ifs);
  }

  return parms;
}

#endif