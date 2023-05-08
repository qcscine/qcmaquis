/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_MPO_TIMES_MPS_HPP
#define MAQUIS_DMRG_MPO_TIMES_MPS_HPP

#include <exception>
#include <map>
#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/generate_mpo/1D_mpo_maker.hpp"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps_sectors.h"

struct MPOTimesMPSException : public std::exception {
    const char* what() const throw() {
        return "MPS times MPO not implemented for this symmetry group"; 
    }
};

template <class Matrix, class SymmGroup, class SymmType = void>
class MPOTimesMPSTraitClass {
public:
    // Types definition 
    using MPSType = MPS<Matrix, SymmGroup>;
    using MPOType = MPO<Matrix, SymmGroup>;
    using ModelType = Model<Matrix, SymmGroup>;
    using ChargeType = typename SymmGroup::charge;

    /**
     * @brief Class constructor
     * @param mps Reference Matrix Product State
     * @param mpo Reference Matrix Product Operator
     * @param model Model class
     * @param lattice DMRG lattice
     */
    MPOTimesMPSTraitClass(const MPSType& mps, const ModelType& model, const Lattice& lattice,
                          ChargeType overallQN, int bondDimension_)
        : mpsRef(mps), modelRef(model), latticeRef(lattice), totalQN(overallQN), bondDimension(bondDimension_)
    {
        // Sets up some data required for mpo_times_mp
        int max_site_type = 0;
        siteTypes.resize(latticeRef.size());
        std::fill(siteTypes.begin(), siteTypes.end(), 0);
        for (int p = 0; p < lattice.size(); ++p) {
          siteTypes[p] = lattice.template get_prop<int>("type", p);
          max_site_type = std::max(siteTypes[p], max_site_type);
        }
        siteBases.resize(max_site_type+1);
        for (int type = 0; type < siteBases.size(); ++type)
          siteBases[type] = modelRef.phys_dim(type);    
    }

    /**
     * @brief Remove an electron from the
     * @param siteToIonize Orbital from where the electron is ionized.
     * @param upOrDown Enum class indicating whether to ionize an alpha or a beta electron.
     * @return MPSType MPS representation of the ionized wave function.
     */
    MPSType ionizeMPS(int siteToIonize, generate_mpo::IonizedOrbital upOrDown) {
        auto ionizedQN = totalQN;
        ionizedQN[0] -= 1;
        auto indexAllowed = allowed_sectors(siteTypes, siteBases, ionizedQN, bondDimension);
        auto destructorOperator = generate_mpo::make_destroy_mpo(latticeRef, modelRef, siteToIonize, upOrDown);
        //
        resetCharges();
        MPS<Matrix, SymmGroup> ionizedMPS(latticeRef.size());
        for (int iMPS = 0; iMPS < ionizedMPS.length(); iMPS++)
          ionizedMPS[iMPS] = mpo_times_mps(destructorOperator, mpsRef, iMPS, charges, indexAllowed);
        return ionizedMPS;
    }

    /**
     * @brief Apply an MPO onto the MPS stored in the class.
     * @param mpo Matrix Product Operator to be applied onto the MPS.
     * @return Result of mpo*mpsRef
     */
    MPSType applyMPO(const MPOType& mpo) {
        auto indexAllowed = allowed_sectors(siteTypes, siteBases, totalQN, bondDimension);
        resetCharges();
        MPS<Matrix, SymmGroup> finalMPS(latticeRef.size());
        for (int iMPS = 0; iMPS < finalMPS.length(); iMPS++)
          finalMPS[iMPS] = mpo_times_mps(mpo, mpsRef, iMPS, charges, indexAllowed);
        return finalMPS;
    }

    /**
     * @brief Resets the charge tracker 
     * To be run before a new MPO is applied.
     */
    void resetCharges() {
        charges = {SymmGroup::IdentityCharge};
        mapTrackingBlocks.clear();
        Index<SymmGroup> tmp;
        tmp.insert(std::make_pair(SymmGroup::IdentityCharge, 1));
        mapTrackingBlocks[0] = tmp;
    }

    /**
     * @brief Multiplication between a MPOTensor and an MPSTensor
     * 
     * This routine performes the following contraction:
     * 
     *   |           |
     * --O--   =   ==O==
     *   |
     * --o--
     * 
     * where the big O represents a mpo tensor and the small o an MPS tensor.
     * The result of the contraction is an MPS with the bond dimension that is larger
     * than the original one by a factor of b_i-1 for the left index and b_i for the
     * right index.
     * Mathematically, the contraction to perform is the following:
     * 
     * +---
     *  \    |\/| (\sigma_i)    \    / (\sigma_i,\sigma_i')  = |\/| (\sigma_i)
     *  /    |  | (a_{i-1},a_i)  \/\/  (b_{i-1},b_i)         = |  | (b_{i-1} a_{i-1} , a_i b_i)
     * +---
     * 
     * where by applying an MPO onto an MPS we increase of its bond dimension.
     * 
     * The code is structured as follows:
     * 
     * 1) all matrices (the input MPSTensor and the output MPSTensor) are first made right paired.
     * 2) each pair of values (b_{i-1},b_i) will contribute to a block of the final MPSTensor.
     *    All pairs must be "stacked" together to give the final, overall tensor.
     *    We must pay attention to the fact that the MPOTensor is stored in a sparse format, therefore only
     *    the pairs (b_{i-1},b_i) giving non-zero elements are actually stored.
     * 3) we first loop over the row (b_{i-1}) and, for each element of the row, we add all blocks
     *    of the MPSTensor obtained for all allowed values of b_i. Blocks are put together with the
     *    [mps_join] method, by exploiting the fact that the row can be recycled. This means that 
     *    the overall structure of the MPS will be:
     * 
     *    +-                                                                       -+
     *    |                       |                     |     |                     |
     *    |   M(a_i, a_{i+1} b_1) | M(a_i, a_{i+1} b_2) | ... | M(a_i, a_{i+1} b_)  |
     *    |                       |                     |     |                     |
     *    +-                                                                       -+
     * 
     *    Clearly, the blocks will then be partitioned based on the symmetry properties.
     * 4) all the resulting blocks are then "stacked" one over the other columnwise.
     *    Again, we can recycle the column index now.
     * 
     * @tparam Matrix Class of the matrix in which the elements of the MPOTensor are stored.
     * @tparam Matrix Class of the matrix in which the elements of the MPSTensor are stored.
     * @tparam SymmGroup Symmetry group of the Hamiltonian.
     * @param mpo MPOTensor object.
     * @param mps MPSTensor object.
     * @param site site for which the MPO/MPS multiplication is applied.
     * @param in_delta Charge difference that is "brought" by the left index of the MPO.
     * @param allowed_sectors Index with the physically allowed symmetry sectors per site
     * (see mps_sectors.h for more details).
     * @return MPSTensor<Matrix, SymmGroup> result of the MPOTensor x MPSTensor operation.
     */
    MPSTensor<Matrix, SymmGroup> mpo_times_mps(MPO<Matrix, SymmGroup> const & mpo, MPS<Matrix, SymmGroup> const & mps,
                                               int site, std::vector< typename SymmGroup::charge> & in_delta,
                                               std::vector<Index<SymmGroup>> const& allowed_sectors)
    {
        // Aliases for derived types
        using MPOTensor_detail::term_descriptor;
        using boost::tuples::get;
        using charge = typename SymmGroup::charge;
        using value_type = typename Matrix::value_type;
        using row_proxy = typename MPOTensor<Matrix, SymmGroup>::row_proxy;
        using col_proxy = typename MPOTensor<Matrix, SymmGroup>::col_proxy;
        // We transform the MPS to the right-paired representation and construct the respective basis.
        mps[site].make_right_paired();
        block_matrix<Matrix, SymmGroup> const & data = mps[site].data();
        Index<SymmGroup> const & right_i = mps[site].col_dim();
        ProductBasis<SymmGroup> right_pb(mps[site].site_dim(), mps[site].col_dim(),
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                        -boost::lambda::_1, boost::lambda::_2));
        assert(in_delta.size() == mpo[site].row_dim());
        std::vector<charge> out_delta(mpo[site].col_dim());
        std::map<int, Index<SymmGroup> > new_right_i_map;
        // == CREATES THE INDEX FOR THE NEW MPS ==
        // We have to separate this step from what happens below because we check that the input indexes
        // are coherent.
        for (int iRow = 0; iRow < mpo[site].row_dim(); iRow++)
        {
            row_proxy row_b2 = mpo[site].row(iRow);
            for (typename row_proxy::const_iterator row_it = row_b2.begin(); row_it != row_b2.end(); ++row_it) {
                // Access to the specific MPO element
                int iCol = row_it.index();
                term_descriptor<Matrix, SymmGroup, true> access = mpo[site].at(iRow, iCol);
                typename operator_selector<Matrix, SymmGroup>::type const & W = access.op();
                // Calculates the difference in symmetry that is "generated" by the operator.
                // Here we must remember that, in an MPS, the rows and the columns will have the same symmetry.
                // This is not true for the MPO, since the local operator can "induce" a change in the symmetry.
                // This also means that, after the application of the MPO, the charges will undergo a "shift".
                charge W_delta = SymmGroup::fuse(W.basis().right_charge(0), -W.basis().left_charge(0));
                out_delta[iCol] = SymmGroup::fuse(in_delta[iRow], W_delta);
                // Loop over the symmetry blocks of the input MPS (that is now right paired).
                // Finds out which block of the final MPS would be populated by the application
                // of that MPO --> note that the final block might be unphysical, in which case
                // we neglect the term.
                for (size_t b = 0; b < data.n_blocks(); ++b)
                {
                    charge lc = data.basis().left_charge(b);
                    charge rc = data.basis().right_charge(b);
                    // Loop over the symmetry block of the MPO (so, over the physical dimensions).
                    for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
                    {
                        // Constructs the final block of the MPS
                        charge phys_in = W.basis().left_charge(w_block);
                        charge phys_out = W.basis().right_charge(w_block);
                        // We extract the right index by difference of the ProductBasis and the input
                        // physical index (times -1.).
                        charge out_l_charge = SymmGroup::fuse(lc, in_delta[iRow]);
                        if (! ChargeDetailClass<SymmGroup>::physical(out_l_charge) || !mapTrackingBlocks[iRow].has(out_l_charge) || !allowed_sectors[site].has(out_l_charge))
                            continue;
                        if (!mps[site].site_dim().has(phys_in))
                            continue;
                        charge in_r_charge = SymmGroup::fuse(rc, phys_in);
                        if (!right_i.has(in_r_charge))
                            continue;
                        charge out_r_charge = SymmGroup::fuse(out_l_charge, phys_out);
                        if (!allowed_sectors[site+1].has(out_r_charge))
                            continue;
                        // Loads the map for the right charges
                        if (new_right_i_map.find(iCol) == new_right_i_map.end())
                            new_right_i_map[iCol] = Index<SymmGroup>();
                        if (!new_right_i_map[iCol].has(out_r_charge)) {
                            new_right_i_map[iCol].insert(std::make_pair(out_r_charge, right_i.size_of_block(in_r_charge)));
                        }
                        else {
                            if (new_right_i_map[iCol].size_of_block(out_r_charge) != right_i.size_of_block(in_r_charge))
                                throw std::runtime_error("Incoherence detected in input data");
                        }
                    }
                }
            }
        }
        // == POPULATES THE FINAL, OVERALL MPS ==
        // Prepares the data structure
        Index<SymmGroup> finalLeft = mapTrackingBlocks[0], 
                         finalRight = new_right_i_map[0],
                         finalPhys = mps[site].site_dim();
        for (int iRow = 1; iRow < mpo[site].row_dim(); iRow++)
            for (int iCharge = 0; iCharge < mapTrackingBlocks[iRow].size(); iCharge++)
                if (finalLeft.has(mapTrackingBlocks[iRow][iCharge].first))
                    finalLeft[finalLeft.position(mapTrackingBlocks[iRow][iCharge].first)].second += mapTrackingBlocks[iRow][iCharge].second;
                else
                    finalLeft.insert(std::make_pair(mapTrackingBlocks[iRow][iCharge].first, mapTrackingBlocks[iRow][iCharge].second));

        for (int iCol = 1; iCol < mpo[site].col_dim(); iCol++)
            for (int iCharge = 0; iCharge < new_right_i_map[iCol].size(); iCharge++)
                if (finalRight.has(new_right_i_map[iCol][iCharge].first))
                    finalRight[finalRight.position(new_right_i_map[iCol][iCharge].first)].second += new_right_i_map[iCol][iCharge].second;
                else
                    finalRight.insert(std::make_pair(new_right_i_map[iCol][iCharge].first, new_right_i_map[iCol][iCharge].second));

        auto finalMPS = MPSTensor<Matrix, SymmGroup>(finalPhys, finalLeft, finalRight, false, 0.);
        finalMPS.make_right_paired();
        Index<SymmGroup> finalPhysAndRight = adjoin(finalPhys)*finalRight;
        common_subset(finalLeft, finalPhysAndRight);

        // Calculates the internal thresholds for the MPS
        std::map< charge, std::vector<int> > thresholdLeft;
        std::map< std::pair<charge, charge>, std::vector<int> > thresholdRight;

        for (int iLeftCharge = 0; iLeftCharge < finalLeft.size(); iLeftCharge++)
            thresholdLeft[finalLeft[iLeftCharge].first] = std::vector<int>(mpo[site].row_dim(), 0);
        for (int iRightCharge = 0; iRightCharge < finalRight.size(); iRightCharge++)
            for (int iPhysCharge = 0; iPhysCharge < finalPhys.size(); iPhysCharge++)
                thresholdRight[std::make_pair(finalPhys[iPhysCharge].first, finalRight[iRightCharge].first)] = std::vector<int>(mpo[site].col_dim(), 0);

        for (int iRow = 0; iRow < mpo[site].row_dim(); iRow++) {
            for (auto& iCharge: thresholdLeft) {
                iCharge.second[iRow] = (iRow == 0) ? 0 : iCharge.second[iRow-1];
                if (iRow > 0)
                    if (mapTrackingBlocks[iRow-1].has(iCharge.first))
                        iCharge.second[iRow] += mapTrackingBlocks[iRow-1].size_of_block(iCharge.first);
            }
        }

        for (int iCol = 0; iCol < mpo[site].col_dim(); iCol++) {
            for (auto& iCharge: thresholdRight) {
                iCharge.second[iCol] = (iCol == 0) ? 0 : iCharge.second[iCol-1];
                if (iCol > 0 && finalPhys.has(iCharge.first.first) && new_right_i_map[iCol-1].has(iCharge.first.second))
                    iCharge.second[iCol] += /*finalPhys.size_of_block(iCharge.first.first) *  */new_right_i_map[iCol-1].size_of_block(iCharge.first.second);
            }
        }

        // Load the data inside the finalMPS MPSTensor
        ProductBasis<SymmGroup> out_right_pb(finalPhys, finalRight,
                                             boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                            -boost::lambda::_1, boost::lambda::_2));
        // Loop over the columns of the MPO
        for (int iCol = 0; iCol < mpo[site].col_dim(); iCol++)
        {
            // This will store the MPS obtained by blocking all values of b_i.
            col_proxy col_b2 = mpo[site].column(iCol);
            // Load the available (b_{i-1}, b_i) pairs.
            std::vector<int> availableRows;
            for (typename col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it)
                availableRows.push_back(col_it.index());
            // Loop over the rows of the MPO
            for (int iRow = 0; iRow < mpo[site].row_dim(); iRow++) 
            {
                block_matrix<Matrix, SymmGroup>& prod = finalMPS.data();
                // Check if the pair is available
                if (std::find(availableRows.begin(), availableRows.end(), iRow) != availableRows.end()) {
                    term_descriptor<Matrix, SymmGroup, true> access = mpo[site].at(iRow, iCol);
                    for (size_t oi = 0; oi < access.size(); ++oi)
                    {
                        typename operator_selector<Matrix, SymmGroup>::type const & W = access.op(oi);
                        // Loop over the MPS blocks (remember that the MPS is right paired)
                        for (size_t b = 0; b < data.n_blocks(); ++b)
                        {
                            auto lc = data.basis().left_charge(b);
                            auto rc = data.basis().right_charge(b);
                            auto out_l_charge = SymmGroup::fuse(lc, in_delta[iRow]);
                            for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
                            {
                                auto phys_in = W.basis().left_charge(w_block);
                                auto phys_out = W.basis().right_charge(w_block);
                                if (!ChargeDetailClass<SymmGroup>::physical(out_l_charge) || 
                                    !mapTrackingBlocks[iRow].has(out_l_charge) ||
                                    !allowed_sectors[site].has(out_l_charge))
                                    continue;
                                if (!mps[site].site_dim().has(phys_in))
                                    continue;
                                auto in_r_charge = SymmGroup::fuse(rc, phys_in);
                                if (!right_i.has(in_r_charge))
                                    continue;
                                auto out_r_charge = SymmGroup::fuse(out_l_charge, phys_out);
                                if (!allowed_sectors[site+1].has(out_r_charge))
                                    continue;
                                size_t in_right_offset  = right_pb(phys_in, in_r_charge);
                                size_t out_right_offset = out_right_pb(phys_out, out_r_charge);
                                size_t l_size = data.basis().left_size(b);
                                size_t r_size = right_i.size_of_block(in_r_charge);
                                Matrix const & iblock = data[b];
                                size_t o = prod.find_block(out_l_charge, out_l_charge);
                                if (o == prod.n_blocks())
                                    throw std::runtime_error("Block not found in the MPS");
                                Matrix & oblock = prod[o];
                                /*
                                for(size_t rr = 0; rr < r_size; ++rr) {
                                    maquis::dmrg::detail::iterator_axpy(&iblock(0, in_right_offset + rr),
                                                                        &iblock(0, in_right_offset + rr) + l_size,
                                                                        &oblock(thresholdLeft[out_l_charge][iRow], thresholdRight[std::make_pair(phys_out, out_r_charge)][iCol] + out_right_offset + rr),
                                                                        alfa);
                                }
                                */
                                for (int iRowPhys = 0; iRowPhys < W.basis().left_size(w_block); iRowPhys++) {
                                    for (int iColPhys = 0; iColPhys < W.basis().right_size(w_block); iColPhys++) {
                                        value_type alfa = access.scale(oi) * W[w_block](iRowPhys, iColPhys);
                                        auto thresholdRightElement = thresholdRight[std::make_pair(phys_out, out_r_charge)][iCol];
                                        for(int rr = 0; rr < r_size; ++rr) {
                                            maquis::dmrg::detail::iterator_axpy(&iblock(0, in_right_offset + iRowPhys*r_size + rr),
                                                                                &iblock(0, in_right_offset + iRowPhys*r_size + rr) + l_size,
                                                                                &oblock(thresholdLeft[out_l_charge][iRow],
                                                                                        out_right_offset                                  // Offset given by the product basis
                                                                                      + iColPhys*finalRight.size_of_block(out_r_charge)   // (sigma*m) values for all preceding sigma values
                                                                                      + thresholdRightElement                             // Threshold induced by b
                                                                                      + rr),                                              // Different m values for the columns of the tensor
                                                                                alfa);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        // Modifies the input delta/map for the indexes.
        std::swap(in_delta, out_delta);
        std::swap(mapTrackingBlocks, new_right_i_map);
        return finalMPS;
    }

    /**
     * @brief Old version of [mpo_times_mps], kept for back-compatibility
     * 
     * This version of [mpo_times_mps] works only for elementary operators - in other words,
     * for MPOs with bond dimension b=1.
     *
     * @param mpo Input Matrix Product operator
     * @param mps Input Matrix Product State
     * @param in_delta Charge difference "accumulated" by the MPO
     * @return MPSTensor<Matrix, SymmGroup> 
     */
    static MPSTensor<Matrix, SymmGroup> mpo_times_mps_singleop(MPOTensor<Matrix, SymmGroup> const & mpo,
                                                                  MPSTensor<Matrix, SymmGroup> const & mps,
                                                                  typename SymmGroup::charge & in_delta)
    {
        using MPOTensor_detail::term_descriptor;
        using boost::tuples::get;
        using charge = typename SymmGroup::charge;
        using value_type = typename Matrix::value_type;
    
        mps.make_right_paired();
        block_matrix<Matrix, SymmGroup> const & data = mps.data();
    
        Index<SymmGroup> const & right_i = mps.col_dim();
        //maquis::cout << "      mps.site_dim: " << mps.site_dim() << std::endl;
        //maquis::cout << "      mps.col_dim : " << mps.col_dim() << std::endl;
        ProductBasis<SymmGroup> right_pb(mps.site_dim(), mps.col_dim(),
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                            -boost::lambda::_1, boost::lambda::_2));
    
        term_descriptor<Matrix, SymmGroup, true> access = mpo.at(0,0);
        typename operator_selector<Matrix, SymmGroup>::type const & W = access.op();
    
        charge W_delta = SymmGroup::fuse(W.basis().right_charge(0), -W.basis().left_charge(0));
        charge out_delta = SymmGroup::fuse(in_delta, W_delta);
    
        Index<SymmGroup> new_left_i, new_right_i, new_phys_i;
        for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
        {
            charge phys_in = W.basis().left_charge(w_block);
            if (! mps.site_dim().has(phys_in) ) continue;
            charge phys_out = W.basis().right_charge(w_block);
            new_phys_i.insert(std::make_pair(phys_out, W.basis().right_size(w_block)));
        }
    
        for (size_t b = 0; b < data.n_blocks(); ++b)
        {
            charge lc = data.basis().left_charge(b);
            charge rc = data.basis().right_charge(b); 
    
            for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
            {
                charge phys_in = W.basis().left_charge(w_block);
                charge phys_out = W.basis().right_charge(w_block);
    
                // add the operator deltas from previous sites to the left charge
                charge out_l_charge = SymmGroup::fuse(lc, in_delta); // unpaired
                if (! ChargeDetailClass<SymmGroup>::physical(out_l_charge)) continue;
    
                charge in_r_charge = SymmGroup::fuse(rc, phys_in); // unpaired
                if (!right_i.has(in_r_charge)) continue; // do we have phys_in in block b?
    
                charge out_r_charge = SymmGroup::fuse(out_l_charge, phys_out); // unpaired
    
                if (!new_left_i.has(out_l_charge)) new_left_i.insert(std::make_pair(out_l_charge, data.basis().left_size(b)));
                if (!new_right_i.has(out_r_charge)) new_right_i.insert(std::make_pair(out_r_charge, right_i.size_of_block(in_r_charge)));
            }
        }
    
        //maquis::cout << "      new_left_i: " << new_left_i << std::endl;
        //maquis::cout << "      new_right_i: " << new_right_i << std::endl;
        //maquis::cout << "      new_phys_i: " << new_phys_i << std::endl;
    
        ProductBasis<SymmGroup> out_right_pb(new_phys_i, new_right_i,
                                             boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                                     -boost::lambda::_1, boost::lambda::_2));
        block_matrix<Matrix, SymmGroup> prod;
    
        for (size_t b = 0; b < data.n_blocks(); ++b)
        {
            charge lc = data.basis().left_charge(b);
            charge rc = data.basis().right_charge(b); 
    
            for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
            {
                charge phys_in = W.basis().left_charge(w_block);
                charge phys_out = W.basis().right_charge(w_block);
    
                charge out_l_charge = SymmGroup::fuse(lc, in_delta); // unpaired
                if (! ChargeDetailClass<SymmGroup>::physical(out_l_charge)) continue;
    
                charge in_r_charge = SymmGroup::fuse(rc, phys_in); // unpaired
                if (!mps.site_dim().has(phys_in)) continue; // do we have phys_in in block b ?
                if (!right_i.has(in_r_charge)) continue; // do we have phys_in in block b?
    
                charge out_r_charge = SymmGroup::fuse(out_l_charge, phys_out); // unpaired
    
                // source     -> data[b](·, in_right_offset + 1:rsize)  
                // destination -> prod[o](·, out_right_offset + 1:rsize)
                size_t in_right_offset  = right_pb(phys_in,  in_r_charge); 
                size_t out_right_offset = out_right_pb(phys_out, out_r_charge); 
                size_t l_size = data.basis().left_size(b);
                size_t r_size = right_i.size_of_block(in_r_charge);
    
                Matrix const & iblock = data[b];
    
                size_t o = prod.find_block(out_l_charge, out_l_charge); // out_l_charge = out_r_charge right-paired
                if (o == prod.n_blocks())
                    o = prod.insert_block(Matrix(l_size, out_right_pb.size(-phys_out, out_r_charge)), out_l_charge, out_l_charge);
    
                Matrix & oblock = prod[o];
    
                value_type alfa = access.scale() * W[w_block](0,0);
                //maquis::cout << " access.scale()  ... " << alfa << std::endl;
                //maquis::cout << " W[w_block](0,0) ... " << W[w_block](0,0) << "for block " << w_block << std::endl;
                //maquis::cout << " alfa            ... " << access.scale() << std::endl;
                for(size_t rr = 0; rr < r_size; ++rr)
                    maquis::dmrg::detail::iterator_axpy(&iblock(0, in_right_offset + rr),
                                                        &iblock(0, in_right_offset + rr) + l_size,
                                                        &oblock(0, out_right_offset + rr),
                                                        alfa);
            }
        } 
        std::swap(in_delta, out_delta);
    
        MPSTensor<Matrix, SymmGroup> ret;
        ret.make_right_paired();
        ret.left_i = new_left_i;
        ret.right_i = new_right_i;
        ret.phys_i = new_phys_i;
        swap(ret.data(), prod);
    
        return ret;
    }

private:
    // Class members
    const MPSType& mpsRef;
    const ModelType& modelRef;
    const Lattice& latticeRef;
    std::vector<Index<SymmGroup> > siteBases;           // Physical basis per site
    std::vector<int> siteTypes;                         // Type of each site
    std::vector<ChargeType> charges;                    // Charges of the MPS that is currently constructed
    std::map<int, Index<SymmGroup> > mapTrackingBlocks; // Symmetry tracker
    ChargeType totalQN;
    int bondDimension;
};

/** @brief Overload for the SU2U1 class, yet to be implemented */
template <class Matrix, class SymmGroup>
class MPOTimesMPSTraitClass<Matrix, SymmGroup, symm_traits::enable_if_su2_t<SymmGroup>> {
public:

    // Types definition 
    using MPSType = MPS<Matrix, SymmGroup>;
    using MPOType = MPO<Matrix, SymmGroup>;
    using ModelType = Model<Matrix, SymmGroup>;
    using ChargeType = typename SymmGroup::charge;

    /** @brief Class constructor */
    MPOTimesMPSTraitClass(const MPSType& mps, const ModelType& model, const Lattice& lattice,
                          ChargeType overallQN, int bondDimension_) 
    {
        throw MPOTimesMPSException();
        //throw std::runtime_error("MPOTimesMPSTraitClass not available for spin-adapted Hamiltonians");
    }

    /** @brief General overload */
    MPSTensor<Matrix, SymmGroup> mpo_times_mps(MPO<Matrix, SymmGroup> const & mpo, MPS<Matrix, SymmGroup> const & mps,
                                               int site, std::vector< typename SymmGroup::charge> & in_delta,
                                               std::vector<Index<SymmGroup>> const& allowed_sectors)
    {
        throw MPOTimesMPSException();
        //throw std::runtime_error("[mpo_times_mps] not yet implemented for SU2U1 symmetry");
    }

    /** @brief Single-operator specialization */
    static MPSTensor<Matrix, SymmGroup> mpo_times_mps_singleop(MPOTensor<Matrix, SymmGroup> const & mpo,
                                                           MPSTensor<Matrix, SymmGroup> const & mps,
                                                           typename SymmGroup::charge & in_delta)
    {
        throw MPOTimesMPSException();
        //throw std::runtime_error("[mpo_times_mps_singleop] not yet implemented for SU2U1 symmetry");
    }

    /** @brief AppplyOp method */
    MPSType applyMPO(const MPOType& mpo) 
    {
        throw MPOTimesMPSException();
    }
};

#endif
