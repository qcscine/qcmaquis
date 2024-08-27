/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2021 Laboratory for Physical Chemistry, ETH Zurich
 *               2021-2022 by Alberto Baiardi <abaiardi@ethz.ch>
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

/**
 * @brief Class representing a Two-Dimensional Real-Space Fermi-Hubbard model.
 *
 * The Hamiltonian is written as:
 *
 * H = - t * \sum_{<i,j>} \sum_{\sigma} a_{i,\sigma}^\dagger a_{j,\sigma}
 *     + U \sum_i \sum_{\sigma} n_{i,\sigma} n_{i,\bar{\sigma}}
 *
 * where <i,j> denotes a sum over all nearest-neighbour pairs of a square
 * spin lattice.
 */

#ifndef REALSPACE_FERMIHUBBARD_MODEL
#define REALSPACE_FERMIHUBBARD_MODEL

template <class Matrix>
class FermiHubbardRealTwoU1 : public model_impl<Matrix, TwoU1> {
 public:
  // Types definition
  using base = model_impl<Matrix, TwoU1>;
  using tag_type = typename base::tag_type;
  using op_t = typename base::op_t;
  using measurements_type = typename base::measurements_type;
  using value_type = typename Matrix::value_type;
  using table_ptr = typename base::table_ptr;
  using pos_t = typename Lattice::pos_t;

  /** @brief Class constructor */
  FermiHubbardRealTwoU1(
      const Lattice& lat_, BaseParameters& parms_,
      bool isTranscorrelated = false
  )
      : lat(lat_),
        parms(parms_),
        tag_handler(new TagHandler<Matrix, TwoU1>()),
        isTranscorrelated_(isTranscorrelated) {
    // Definition of the charges (i.e., the relevant QN)
    TwoU1::charge A(0), B(0), C(0), D(1);
    B[0] = 1;
    C[1] = 1;
    phys.insert(std::make_pair(A, 1));
    phys.insert(std::make_pair(B, 1));
    phys.insert(std::make_pair(C, 1));
    phys.insert(std::make_pair(D, 1));
    // Elementary operators
    op_t create_up_op, create_down_op, destroy_up_op, destroy_down_op, fill_op,
        ident_op;
    ident_op.insert_block(Matrix(1, 1, 1), A, A);
    ident_op.insert_block(Matrix(1, 1, 1), B, B);
    ident_op.insert_block(Matrix(1, 1, 1), C, C);
    ident_op.insert_block(Matrix(1, 1, 1), D, D);
    create_up_op.insert_block(Matrix(1, 1, 1), A, B);
    create_up_op.insert_block(Matrix(1, 1, 1), C, D);
    create_down_op.insert_block(Matrix(1, 1, 1), A, C);
    create_down_op.insert_block(Matrix(1, 1, 1), B, D);
    destroy_up_op.insert_block(Matrix(1, 1, 1), B, A);
    destroy_up_op.insert_block(Matrix(1, 1, 1), D, C);
    destroy_down_op.insert_block(Matrix(1, 1, 1), C, A);
    destroy_down_op.insert_block(Matrix(1, 1, 1), D, B);
    fill_op.insert_block(Matrix(1, 1, 1), A, A);
    fill_op.insert_block(Matrix(1, 1, -1), B, B);
    fill_op.insert_block(Matrix(1, 1, -1), C, C);
    fill_op.insert_block(Matrix(1, 1, 1), D, D);
    /* Create operator tag table */
    ident = tag_handler->register_op(ident_op, tag_detail::bosonic);
    fill = tag_handler->register_op(fill_op, tag_detail::bosonic);
    create_up = tag_handler->register_op(create_up_op, tag_detail::fermionic);
    create_down =
        tag_handler->register_op(create_down_op, tag_detail::fermionic);
    destroy_up = tag_handler->register_op(destroy_up_op, tag_detail::fermionic);
    destroy_down =
        tag_handler->register_op(destroy_down_op, tag_detail::fermionic);
  }

  void create_terms() override {
    // Definition of the parameters of the Fermi-Hubbard lattice
    value_type U = parms["U_FermiHubbard"];
    value_type tx, ty;
    if (!parms.is_set("tx_FermiHubbard") && !parms.is_set("tx_FermiHubbard")) {
      tx = parms["t_FermiHubbard"].as<value_type>();
      ty = parms["t_FermiHubbard"].as<value_type>();
    } else if (parms.is_set("tx_FermiHubbard") && parms.is_set("tx_FermiHubbard")) {
      tx = parms["tx_FermiHubbard"].as<value_type>();
      ty = parms["ty_FermiHubbard"].as<value_type>();
    } else {
      throw std::runtime_error("Please set *both* tx and ty");
    }
    int width = parms["width_FermiHubbard"];
    int height = parms["height_FermiHubbard"];
    if (isTranscorrelated_ && !parms.is_set("J_Transcorrelated"))
      throw std::runtime_error("Please set the transcorrelation parameter");
    value_type J =
        parms.is_set("J_Transcorrelated") ? parms["J_Transcorrelated"] : 0.;
    // == HAMILTONIAN CREATION ==
    double hamiltonianNorm = 0.;
    double nonHermitianNorm = 0.;
    auto jw = JordanWignerHandler<Matrix, TwoU1>(
        lat, fill, create_up, create_down, destroy_up, destroy_down
    );
    // Potential energy term
    for (int iSite = 0; iSite < height * width; iSite++) {
      std::vector<OperatorType> opVector = {
          OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
          OperatorType::CreateBeta, OperatorType::DestroyBeta};
      std::vector<pos_t> positions = {iSite, iSite, iSite, iSite};
      this->terms_.push_back(
          jw.getTerm(positions, opVector, tag_handler, true, U)
      );
      hamiltonianNorm += std::norm(U);
    }
    // Hopping term
    for (int iSite = 0; iSite < width * height; iSite++) {
      // Note that this is specific of the Two-dimensional Fermi-Hubbard model.
      bool hasLower = ((iSite + 1) % height != 0);
      bool hasRight = (iSite + height) < height * width;
      std::vector<pos_t> posVectorDown, posVectorDownHerm;
      std::vector<pos_t> posVectorLeft, posVectorLeftHerm;
      if (hasLower) {
        posVectorDown = {iSite, iSite + 1};
        posVectorDownHerm = {iSite + 1, iSite};
      } else {
        posVectorDown = {iSite, iSite - height + 1};
        posVectorDownHerm = {iSite - height + 1, iSite};
      }

      if (hasRight) {
        posVectorLeft = {iSite, iSite + height};
        posVectorLeftHerm = {iSite + height, iSite};
      } else {
        posVectorLeft = {iSite, iSite % height};
        posVectorLeftHerm = {iSite % height, iSite};
      }
      // Adds the hopping terms to the term vector
      value_type coeffx = -tx;
      value_type coeffy = -ty;
      std::vector<OperatorType> opVector = {
          OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
      this->terms_.push_back(
          jw.getTerm(posVectorLeft, opVector, tag_handler, true, coeffx)
      );
      hamiltonianNorm += std::norm(coeffx);
      this->terms_.push_back(
          jw.getTerm(posVectorLeftHerm, opVector, tag_handler, true, coeffx)
      );
      hamiltonianNorm += std::norm(coeffx);
      this->terms_.push_back(
          jw.getTerm(posVectorDown, opVector, tag_handler, true, coeffy)
      );
      hamiltonianNorm += std::norm(coeffy);
      this->terms_.push_back(
          jw.getTerm(posVectorDownHerm, opVector, tag_handler, true, coeffy)
      );
      hamiltonianNorm += std::norm(coeffy);
      opVector = {OperatorType::CreateBeta, OperatorType::DestroyBeta};
      this->terms_.push_back(
          jw.getTerm(posVectorLeft, opVector, tag_handler, true, coeffx)
      );
      hamiltonianNorm += std::norm(coeffx);
      this->terms_.push_back(
          jw.getTerm(posVectorLeftHerm, opVector, tag_handler, true, coeffx)
      );
      hamiltonianNorm += std::norm(coeffx);
      this->terms_.push_back(
          jw.getTerm(posVectorDown, opVector, tag_handler, true, coeffy)
      );
      hamiltonianNorm += std::norm(coeffy);
      this->terms_.push_back(
          jw.getTerm(posVectorDownHerm, opVector, tag_handler, true, coeffy)
      );
      hamiltonianNorm += std::norm(coeffy);
      // Transcorrelated-specific contributions
      if (isTranscorrelated_) {
        std::vector<pos_t> posVectorDownTC_1, posVectorDownHermTC_1,
            posVectorDownTC_2, posVectorDownHermTC_2, posVectorDownTC_3,
            posVectorDownHermTC_3;
        std::vector<pos_t> posVectorLeftTC_1, posVectorLeftHermTC_1,
            posVectorLeftTC_2, posVectorLeftHermTC_2, posVectorLeftTC_3,
            posVectorLeftHermTC_3;
        if (hasLower) {
          posVectorDownTC_1 = {iSite, iSite + 1, iSite + 1, iSite + 1};
          posVectorDownHermTC_1 = {iSite + 1, iSite, iSite, iSite};
          posVectorDownTC_2 = {iSite, iSite + 1, iSite, iSite};
          posVectorDownHermTC_2 = {iSite + 1, iSite, iSite + 1, iSite + 1};
          posVectorDownTC_3 = {iSite, iSite + 1, iSite,
                               iSite, iSite + 1, iSite + 1};
          posVectorDownHermTC_3 = {iSite + 1, iSite, iSite + 1,
                                   iSite + 1, iSite, iSite};
        } else {
          posVectorDownTC_1 = {
              iSite, iSite - height + 1, iSite - height + 1,
              iSite - height + 1};
          posVectorDownHermTC_1 = {iSite - height + 1, iSite, iSite, iSite};
          posVectorDownTC_2 = {iSite, iSite - height + 1, iSite, iSite};
          posVectorDownHermTC_2 = {
              iSite - height + 1, iSite, iSite - height + 1,
              iSite - height + 1};
          posVectorDownTC_3 = {iSite, iSite - height + 1, iSite,
                               iSite, iSite - height + 1, iSite - height + 1};
          posVectorDownHermTC_3 = {
              iSite - height + 1, iSite, iSite - height + 1,
              iSite - height + 1, iSite, iSite};
        }
        if (hasRight) {
          posVectorLeftTC_1 = {
              iSite, iSite + height, iSite + height, iSite + height};
          posVectorLeftHermTC_1 = {iSite + height, iSite, iSite, iSite};
          posVectorLeftTC_2 = {iSite, iSite + height, iSite, iSite};
          posVectorLeftHermTC_2 = {
              iSite + height, iSite, iSite + height, iSite + height};
          posVectorLeftTC_3 = {iSite, iSite + height, iSite,
                               iSite, iSite + height, iSite + height};
          posVectorLeftHermTC_3 = {iSite + height, iSite, iSite + height,
                                   iSite + height, iSite, iSite};
        } else {
          posVectorLeftTC_1 = {
              iSite, iSite % height, iSite % height, iSite % height};
          posVectorLeftHermTC_1 = {iSite % height, iSite, iSite, iSite};
          posVectorLeftTC_2 = {iSite, iSite % height, iSite, iSite};
          posVectorLeftHermTC_2 = {
              iSite % height, iSite, iSite % height, iSite % height};
          posVectorLeftTC_3 = {iSite, iSite % height, iSite,
                               iSite, iSite % height, iSite % height};
          posVectorLeftHermTC_3 = {iSite % height, iSite, iSite % height,
                                   iSite % height, iSite, iSite};
        }
        // TC operator 1
        coeffx = -tx * (std::exp(J) - 1.);
        coeffy = -ty * (std::exp(J) - 1.);
        if (std::abs(coeffx) > 1.0E-10) {
          opVector = {
              OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
              OperatorType::CreateBeta, OperatorType::DestroyBeta};
          this->terms_.push_back(
              jw.getTerm(posVectorLeftTC_1, opVector, tag_handler, true, coeffx)
          );
          hamiltonianNorm += std::norm(coeffx);
          nonHermitianNorm += std::norm(coeffx);
          this->terms_.push_back(jw.getTerm(
              posVectorLeftHermTC_1, opVector, tag_handler, true, coeffx
          ));
          hamiltonianNorm += std::norm(coeffx);
          nonHermitianNorm += std::norm(coeffx);
          opVector = {
              OperatorType::CreateBeta, OperatorType::DestroyBeta,
              OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
          this->terms_.push_back(
              jw.getTerm(posVectorLeftTC_1, opVector, tag_handler, true, coeffx)
          );
          hamiltonianNorm += std::norm(coeffx);
          nonHermitianNorm += std::norm(coeffx);
          this->terms_.push_back(jw.getTerm(
              posVectorLeftHermTC_1, opVector, tag_handler, true, coeffx
          ));
          hamiltonianNorm += std::norm(coeffx);
          nonHermitianNorm += std::norm(coeffx);
        }
        if (std::abs(coeffy) > 1.0E-10) {
          opVector = {
              OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
              OperatorType::CreateBeta, OperatorType::DestroyBeta};
          this->terms_.push_back(
              jw.getTerm(posVectorDownTC_1, opVector, tag_handler, true, coeffy)
          );
          hamiltonianNorm += std::norm(coeffy);
          nonHermitianNorm += std::norm(coeffy);
          this->terms_.push_back(jw.getTerm(
              posVectorDownHermTC_1, opVector, tag_handler, true, coeffy
          ));
          hamiltonianNorm += std::norm(coeffy);
          nonHermitianNorm += std::norm(coeffy);
          opVector = {
              OperatorType::CreateBeta, OperatorType::DestroyBeta,
              OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
          this->terms_.push_back(
              jw.getTerm(posVectorDownTC_1, opVector, tag_handler, true, coeffy)
          );
          hamiltonianNorm += std::norm(coeffy);
          nonHermitianNorm += std::norm(coeffy);
          this->terms_.push_back(jw.getTerm(
              posVectorDownHermTC_1, opVector, tag_handler, true, coeffy
          ));
          hamiltonianNorm += std::norm(coeffy);
          nonHermitianNorm += std::norm(coeffy);
        }
        // TC operator 2
        coeffx = -tx * (std::exp(-J) - 1.);
        coeffy = -ty * (std::exp(-J) - 1.);
        if (std::abs(coeffx) > 1.0E-10) {
          opVector = {
              OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
              OperatorType::CreateBeta, OperatorType::DestroyBeta};
          this->terms_.push_back(
              jw.getTerm(posVectorLeftTC_2, opVector, tag_handler, true, coeffx)
          );
          hamiltonianNorm += std::norm(coeffx);
          nonHermitianNorm += std::norm(coeffx);
          this->terms_.push_back(jw.getTerm(
              posVectorLeftHermTC_2, opVector, tag_handler, true, coeffx
          ));
          hamiltonianNorm += std::norm(coeffx);
          nonHermitianNorm += std::norm(coeffx);
          opVector = {
              OperatorType::CreateBeta, OperatorType::DestroyBeta,
              OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
          this->terms_.push_back(
              jw.getTerm(posVectorLeftTC_2, opVector, tag_handler, true, coeffx)
          );
          hamiltonianNorm += std::norm(coeffx);
          nonHermitianNorm += std::norm(coeffx);
          this->terms_.push_back(jw.getTerm(
              posVectorLeftHermTC_2, opVector, tag_handler, true, coeffx
          ));
          hamiltonianNorm += std::norm(coeffx);
          nonHermitianNorm += std::norm(coeffx);
        }
        if (std::abs(coeffy) > 1.0E-10) {
          opVector = {
              OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
              OperatorType::CreateBeta, OperatorType::DestroyBeta};
          this->terms_.push_back(
              jw.getTerm(posVectorDownTC_2, opVector, tag_handler, true, coeffy)
          );
          hamiltonianNorm += std::norm(coeffy);
          nonHermitianNorm += std::norm(coeffy);
          this->terms_.push_back(jw.getTerm(
              posVectorDownHermTC_2, opVector, tag_handler, true, coeffy
          ));
          hamiltonianNorm += std::norm(coeffy);
          nonHermitianNorm += std::norm(coeffy);
          opVector = {
              OperatorType::CreateBeta, OperatorType::DestroyBeta,
              OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
          this->terms_.push_back(
              jw.getTerm(posVectorDownTC_2, opVector, tag_handler, true, coeffy)
          );
          hamiltonianNorm += std::norm(coeffy);
          nonHermitianNorm += std::norm(coeffy);
          this->terms_.push_back(jw.getTerm(
              posVectorDownHermTC_2, opVector, tag_handler, true, coeffy
          ));
          hamiltonianNorm += std::norm(coeffy);
          nonHermitianNorm += std::norm(coeffy);
        }
        // TC operator 3
        if (parms["transcorrelated_3body"] == true) {
          coeffx = 2. * tx * (std::cosh(J) - 1.);
          coeffy = 2. * ty * (std::cosh(J) - 1.);
          if (std::abs(coeffx) > 1.0E-10) {
            opVector = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
                        OperatorType::CreateBeta,  OperatorType::DestroyBeta,
                        OperatorType::CreateBeta,  OperatorType::DestroyBeta};
            this->terms_.push_back(jw.getTerm(
                posVectorLeftTC_3, opVector, tag_handler, true, coeffx
            ));
            hamiltonianNorm += std::norm(coeffx);
            this->terms_.push_back(jw.getTerm(
                posVectorLeftHermTC_3, opVector, tag_handler, true, coeffx
            ));
            hamiltonianNorm += std::norm(coeffx);
            opVector = {OperatorType::CreateBeta,  OperatorType::DestroyBeta,
                        OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
                        OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
            this->terms_.push_back(jw.getTerm(
                posVectorLeftTC_3, opVector, tag_handler, true, coeffx
            ));
            hamiltonianNorm += std::norm(coeffx);
            this->terms_.push_back(jw.getTerm(
                posVectorLeftHermTC_3, opVector, tag_handler, true, coeffx
            ));
            hamiltonianNorm += std::norm(coeffx);
          }
          if (std::abs(coeffy) > 1.0E-10) {
            opVector = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
                        OperatorType::CreateBeta,  OperatorType::DestroyBeta,
                        OperatorType::CreateBeta,  OperatorType::DestroyBeta};
            this->terms_.push_back(jw.getTerm(
                posVectorDownTC_3, opVector, tag_handler, true, coeffy
            ));
            hamiltonianNorm += std::norm(coeffy);
            this->terms_.push_back(jw.getTerm(
                posVectorDownHermTC_3, opVector, tag_handler, true, coeffy
            ));
            hamiltonianNorm += std::norm(coeffy);
            opVector = {OperatorType::CreateBeta,  OperatorType::DestroyBeta,
                        OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
                        OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
            this->terms_.push_back(jw.getTerm(
                posVectorDownTC_3, opVector, tag_handler, true, coeffy
            ));
            hamiltonianNorm += std::norm(coeffy);
            this->terms_.push_back(jw.getTerm(
                posVectorDownHermTC_3, opVector, tag_handler, true, coeffy
            ));
            hamiltonianNorm += std::norm(coeffy);
          }
        }
      }
    }
    maquis::cout << " Hamiltonian Norm = " << std::sqrt(hamiltonianNorm)
                 << std::endl;
    maquis::cout << " Non-Hermtian Hamiltonian Norm = "
                 << std::sqrt(nonHermitianNorm) << std::endl;
    // ALB TODO USE ADD TERM FROM model.hpp (to be put in the utils)
    for (auto it1 = this->terms_.begin(); it1 != this->terms_.end(); it1++) {
      for (auto it2 = it1 + 1; it2 != this->terms_.end();) {
        auto term1 = *it1;
        auto term2 = *it2;
        bool isEqual = true;
        if (term1.size() == term2.size()) {
          for (int idx = 0; idx < term1.size(); idx++)
            if (term1[idx] != term2[idx]) isEqual = false;
        } else {
          isEqual = false;
        }
        if (isEqual) {
          it1->coeff += it2->coeff;
          it2 = this->terms_.erase(it2);
        } else {
          it2++;
        }
      }
    }
  }

  /** @brief Updates the parameters underlying the model */
  void update(BaseParameters const& p) override {
    throw std::runtime_error("update() not yet implemented for this model.");
    return;
  }

  /** @brief Getter for the physical dimensions */
  Index<TwoU1> const& phys_dim(size_t type) const override { return phys; }

  /**
   * @brief Getter for the measurements
   * Note that, so far, no measurements are implemented.
   */
  measurements_type measurements() const override {
    measurements_type meas;
    return meas;
  }

  /** @brief Getter for the identity operator */
  tag_type identity_matrix_tag(size_t type) const override { return ident; }

  /** @brief Getter for the filling operator */
  tag_type filling_matrix_tag(size_t type) const override { return fill; }

  /** @brief Getter for the total quantum number */
  typename TwoU1::charge total_quantum_numbers(BaseParameters& parms
  ) const override {
    typename TwoU1::charge ret(0);
    ret[0] = static_cast<int>(parms["u1_total_charge1"]);
    ret[1] = static_cast<int>(parms["u1_total_charge2"]);
    return ret;
  }

  /** @brief Getter of operators from a string */
  tag_type get_operator_tag(const std::string& name, std::size_t type)
      const override {
    if (name == "create_up")
      return create_up;
    else if (name == "create_down")
      return create_down;
    else if (name == "destroy_up")
      return destroy_up;
    else if (name == "destroy_down")
      return destroy_down;
    else
      throw std::runtime_error("Operator not valid for this model.");
    return 0;
  }

  /** @brief Returns the tag handler */
  table_ptr operators_table() const override { return tag_handler; }

 private:
  // Class members
  Index<TwoU1> phys;
  Lattice const& lat;
  BaseParameters& parms;
  std::shared_ptr<TagHandler<Matrix, TwoU1> > tag_handler;
  tag_type create_up, create_down, destroy_up, destroy_down, ident, fill;
  bool isTranscorrelated_;
};

#endif  // REALSPACE_FERMIHUBBARD_MODEL
