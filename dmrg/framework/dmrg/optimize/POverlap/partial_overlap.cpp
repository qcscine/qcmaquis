/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2017-2017 by Alberto Baiardi <alberto.baiardi@sns.it>
 *
 * This software is part of the ALPS libraries, published under the ALPS
 * Library License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Library License along with
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
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

// +------------+
//  CONSTRUCTORS
// +------------+

// -- Empty constroctor --

template<class Matrix, class SymmGroup>
partial_overlap<Matrix,SymmGroup>::partial_overlap() :
        lattice_L_(0),
        data_left_(),
        data_right_(),
        MPS_reference_(),
        ortho_MPS_()
{ };

// -- Constructor that takes the MPS and a basis vector --

template<class Matrix, class SymmGroup>
partial_overlap<Matrix,SymmGroup>::partial_overlap(const MPSWave& MPS,
                                                   const MPSWave& MPS_reference) :
    lattice_L_(0),
    data_left_(),
    data_right_(),
    MPS_reference_(MPS_reference),
    ortho_MPS_()
{
    // Builds the MPS associated to the ONV provided in input
    assert (MPS.size() == MPS_reference.size()) ;
    lattice_L_ = MPS.size() ;
    // First site
    data_left_.resize(lattice_L_+1) ;
    data_right_.resize(lattice_L_+1) ;
    data_left_[0] = MPS.left_boundary()[0] ;
    data_right_[lattice_L_] = MPS.right_boundary()[0] ;
    // Prepares left boundary data
    for (int i = lattice_L_-1; i >= 0; --i)
        data_right_[i] = contr::overlap_right_step(MPS[i], MPS_reference[i], data_right_[i+1]);
    for (int i = 1; i < lattice_L_; ++i)
        data_left_[i+1] = contr::overlap_left_step(MPS[i], MPS_reference[i], data_left_[i]);
};

// +-------+
//  METHODS
// +-------+

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::update(const MPSWave& MPS,
                                                const dim_type& l,
                                                const int& direction)
{
    // Set the sites where the
    if (direction == 1)
        data_left_[l+1] = contr::overlap_left_step(MPS[l], MPS_reference_[l], data_left_[l]);
    else if (direction == -1)
        data_right_[l] = contr::overlap_right_step(MPS[l], MPS_reference_[l], data_right_[l+1]);
    else
        std::cout << " Direction parameter not recognized" << std::endl ;
};


template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::update(const MPSWave& MPS,
                                                const dim_type& l1,
                                                const dim_type& l2,
                                                const int& direction)
{

    // Set the sites where the
    if (direction == 1) {
        data_left_[l1+1] = contr::overlap_left_step(MPS[l1], MPS_reference_[l1], data_left_[l1]);
        data_left_[l2+1] = contr::overlap_left_step(MPS[l2], MPS_reference_[l2], data_left_[l2]);
    } else if (direction == -1) {
        data_right_[l2] = contr::overlap_right_step(MPS[l2], MPS_reference_[l2], data_right_[l2+1]);
        data_right_[l1] = contr::overlap_right_step(MPS[l1], MPS_reference_[l1], data_right_[l1+1]);
    } else {
        std::cout << " Direction parameter not recognized" << std::endl ;
    }
};

// +-------------------------------+
//  METHODS TO COMPUTE THE OVERLAPS
// +-------------------------------+

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::prepare(const MPSTensor& MPS,
                                                 const dim_type& idx)
{
    ortho_MPS_ = contraction::site_ortho_boundaries(MPS, MPS_reference_[idx], data_left_[idx], data_right_[idx+1]) ;
};

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::prepare(const MPSTensor& TS,
                                                 const dim_type& idx1,
                                                 const dim_type& idx2)
{
    TSTensor ts_ortho(MPS_reference_[idx1], MPS_reference_[idx2]);
    ortho_MPS_ = contraction::site_ortho_boundaries(TS, ts_ortho.make_mps(), data_left_[idx1], data_right_[idx2+1]);
}

template<class Matrix, class SymmGroup>
typename partial_overlap<Matrix, SymmGroup>::value_type partial_overlap<Matrix, SymmGroup>::overlap(const dim_type &i)
{
    // Check data consistency and declaration
    assert (i >= 0 && i < lattice_L_) ;
    throw std::runtime_error("NYI") ;
};

template<class Matrix, class SymmGroup>
typename partial_overlap<Matrix, SymmGroup>::value_type partial_overlap<Matrix, SymmGroup>::overlap(const MPSTensor& MPSTns)
{
    return ietl::dot(MPSTns, ortho_MPS_) ;
};


