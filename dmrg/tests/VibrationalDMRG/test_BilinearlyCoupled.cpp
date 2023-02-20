/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
  using InterfaceType = maquis::DMRGInterface<double>;
  // Adds the final input parameters
  parametersBilinearly.set("init_type", "basis_state_generic");
  parametersBilinearly.set("init_basis_state", "0,0,0,0,0,0");
  parametersBilinearly.set("optimization", "singlesite");
  parametersBilinearly.set("alpha_initial", 1.0E-8);
  parametersBilinearly.set("alpha_initial", 1.0E-15);
  parametersBilinearly.set("alpha_initial", 0.);
  parametersBilinearly.set("nsweeps", 20);
  parametersBilinearly.set("ngrowsweeps", 2);
  parametersBilinearly.set("nmainsweeps", 2);
  parametersBilinearly.set("max_bond_dimension", 20);
  parametersBilinearly.set("MODEL", "watson");
  // Creates the interface
  InterfaceType interface(parametersBilinearly);
  interface.optimize();
  // The reference value can be calculated based on the Harmonic approximation
  BOOST_CHECK_CLOSE(interface.energy(), 3.8164041, 1.0E-5);
#endif
}

#ifdef HAVE_TrivialGroup

/** @brief As above, but with mode-specific nMax */
BOOST_FIXTURE_TEST_CASE(Test_vDMRG_Calculation_BilinearlyCoupled_nMax, WatsonFixture)
{
  using InterfaceType = maquis::DMRGInterface<double>;
  // Adds the final input parameters
  parametersBilinearly.set("init_type", "basis_state_generic");
  parametersBilinearly.set("init_basis_state", "0,0,0,0,0,0");
  parametersBilinearly.set("optimization", "singlesite");
  parametersBilinearly.set("alpha_initial", 1.0E-8);
  parametersBilinearly.set("alpha_initial", 1.0E-15);
  parametersBilinearly.set("alpha_initial", 0.);
  parametersBilinearly.set("nsweeps", 20);
  parametersBilinearly.set("ngrowsweeps", 2);
  parametersBilinearly.set("nmainsweeps", 2);
  parametersBilinearly.set("max_bond_dimension", 20);
  parametersBilinearly.set("MODEL", "watson");
  parametersBilinearly.set("Nmax", "10,9,8,7,6,7");
  // Creates the interface
  InterfaceType interface(parametersBilinearly);
  interface.optimize();
  // The reference value can be calculated based on the Harmonic approximation
  BOOST_CHECK_CLOSE(interface.energy(), 3.8164041, 1.0E-5);
}

#endif // HAVE_TrivialGroup

#endif // DMRG_VIBRATIONAL
