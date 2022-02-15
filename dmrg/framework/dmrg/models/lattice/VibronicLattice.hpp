/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021- by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef MAQUIS_DMRG_VIBRONIC_LATTICE
#define MAQUIS_DMRG_VIBRONIC_LATTICE

#ifdef DMRG_VIBRONIC

#include "dmrg/models/lattice/lattice.h"

/**
 * @brief Lattice representing a vibronic Hamiltonian.
 * 
 * As for the PreBO case, we have here also a multicomponent lattice.
 * This means that sites are mapped either to electronic or to vibrational
 * degrees of freedom.
 * For the vibrational part, the treatment is the same as for the canonical
 * Watson Hamiltonian. The electronic part is instead treated at a different
 * level of theory.
 */

class VibronicLattice : public lattice_impl
{
public:
    // Types definition
    using pos_t = lattice_impl::pos_t;
    
    /** @brief Class constructor */
    explicit VibronicLattice(BaseParameters & parameters) 
        : L(0), vector_types(0), nElecStates(parameters["vibronic_nstates"].as<int>()), 
          nModes(parameters["vibronic_nmodes"].as<int>()), nParticles(0), eleFirst(false)
    {
        // Checks consistency
        // Determines the number of particles. Note that here by number of particles
        // we mean number of monomer with a manifold of excited states
        if (parameters["MODEL"] == "vibronic") {
            nParticles = 1;
        }
        else if (parameters["MODEL"] == "excitonic") {
            if (nElecStates != 1)
                throw std::runtime_error("Excitonic model currently supports only 1 electronic state");
            nParticles = parameters["n_excitons"].as<int>();
        }
        else {
            throw std::runtime_error("Lattice not coherent with the current MODEL");
        }
        L = parameters["L"];
        // Note that here we assume that, for the excitonic case, all molecules are 
        // described by the same Hamiltonian. This means, in practice, that, for each "particle",
        // we have the same number of modes and electronic states. The number of particles is
        // 1 for the excitonic case.
        if ((nModes+nElecStates)*nParticles != L)
            throw std::runtime_error("Incoherence in lattice size for this vibronic lattice");
        vector_types.resize(L);
        // Sites sorting. If == "firstele", put first all the excited states.
        // Otherwise, intertwine electronic and vibrational DOF (for the excitonic case).
        if (parameters["vibronic_sorting"] == "firstele")
            eleFirst = true;
        // The site type is used to distinguish between electronic and vibrational degrees
        // of freedom. Note that we don't distinguish between different "electronic particles"
        // since it 
        if (eleFirst)
            for (auto idx1 = 0; idx1 < nParticles*nElecStates; idx1++)
                vector_types[idx1] = 1;
        else
            for (auto idx1 = 0; idx1 < nParticles; idx1++)
                for (auto idx2 = 0; idx2 < nElecStates; idx2++)
                    vector_types[(nModes+nElecStates)*idx1+idx2] = 1;
        // Monodimensional chain, the maximum number of vertex is 1.
        maximum_vertex = 1;
    }

    /** @brief Returns the next position in the lattice */
    std::vector<pos_t> forward(pos_t i) const {
        std::vector<pos_t> ret;
        if (i < L-1)
            ret.push_back(i+1);
        return ret;
    }

    /** @brief Returns the neighbors of a given site */    
    std::vector<pos_t> all(pos_t i) const {
        std::vector<pos_t> ret;
        if (i < L-1)
            ret.push_back(i+1);
        if (i > 0)
            ret.push_back(i-1);
        return ret;
    }
    
    /**
     * @brief Property getter.
     * 
     * Note that, compared to other models, we have two additional supported properties:
     * - vibindex: retrieves the index of a given vibrational mode associated with a given 
     *             electronic state
     * - eleindex: retrieves the index of a given electronic state of a given particle.
     * 
     * @param property Property identified
     * @param pos vector with the positions associated with the property to be calculated.
     * @return boost::any requested property.
     */
    boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const
    {
        if (property == "label" && pos.size() == 1)
            return boost::any(site_label(pos[0]));
        else if (property == "label" && pos.size() == 2)
            return boost::any(bond_label(pos[0], pos[1]));
        else if (property == "type" && pos.size() == 1)
            return boost::any(vector_types[pos[0]]);
        else if (property == "ParticleType" && pos.size() == 1)
            return boost::any(vector_types[pos[0]]);
        else if (property == "NumTypes")
            return boost::any(2);
        else if (property == "vibindex" && pos.size() == 2)
            // In this case the first index is the molecule, the second one is the specific
            // mode that molecule.
            return (eleFirst) ? boost::any(nParticles*nElecStates+pos[0]*nModes+pos[1]) :
                                boost::any((nElecStates+nModes)*pos[0]+nElecStates+pos[1]);
        else if (property == "eleindex" && pos.size() == 2)
            // In this case the first index is the molecule, the second one is the specific
            // excited state of that molecule.
            return (eleFirst) ? boost::any(nElecStates*pos[0]+pos[1]) :
                                boost::any((nElecStates+nModes)*pos[0]+pos[1]);
        else {
            std::ostringstream ss;
            ss << "No property '" << property << "' with " << pos.size() << " points implemented.";
            throw std::runtime_error(ss.str());
            return boost::any();
        }
    }

    /** @brief Getter for the lattice size */
    pos_t size() const {return L; }

    /** @brief Getter for the maximum index for the site type */
    int  maximum_vertex_type() const { return maximum_vertex; }

private:
    /** Lattice size */
    pos_t L;
    /** Number of electronic states */
    pos_t nElecStates;
    /** Number of vibrational modes */
    pos_t nModes;
    /** Number of excitons */
    pos_t nParticles;
    /** Maximum number of vertexes */
    int maximum_vertex;
    /** Vector with the type of each site */
    std::vector<int> vector_types;
    /** 
     * If true, put all electronic degrees of freedom at the beginning of the lattice.
     * Otherwise, intertwine electrons and nuclei
     */
    bool eleFirst;

    /** @brief Prints the label of a given site */
    std::string site_label (int i) const
    {
        return "( " + boost::lexical_cast<std::string>(i) + " )";
    }
    
    /** @brief Prints the label of a given bond */
    std::string bond_label (int i, int j) const
    {
        return (  "( " + boost::lexical_cast<std::string>(i) + " )"
       + " -- " + "( " + boost::lexical_cast<std::string>(j) + " )");
    }
} ;

#endif // DMRG_VIBRATIONAL

#endif // MAQUIS_DMRG_VIBRONIC_LATTICE