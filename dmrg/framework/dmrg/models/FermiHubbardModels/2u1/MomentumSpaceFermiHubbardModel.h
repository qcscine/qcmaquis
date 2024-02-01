/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2021 Laboratory for Physical Chemistry, ETH Zurich
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

/**
 * @brief Class representing a Two-Dimensional Real-Space Fermi-Hubbard model.
 * The Hamiltonian is obtained as the Fourier transformation of the real-space
 * Hamiltonian, and has the same structure as the electronic-structure one.
 */

#ifndef MOMENTUMSPACE_FERMIHUBBARD_MODEL
#define MOMENTUMSPACE_FERMIHUBBARD_MODEL

template <class Matrix>
class FermiHubbardMomentumTwoU1 : public model_impl<Matrix, TwoU1> {
 public:
  // Types definition
  using base = model_impl<Matrix, TwoU1>;
  using table_type = typename base::table_type;
  using table_ptr = typename base::table_ptr;
  using tag_type = typename base::tag_type;
  using term_descriptor = typename base::term_descriptor;
  using terms_type = typename base::terms_type;
  using op_t = typename base::op_t;
  using measurements_type = typename base::measurements_type;
  using value_type = typename Matrix::value_type;
  using pos_t = typename Lattice::pos_t;
  using MapOfOperatorsType = std::unordered_map<
      std::vector<std::pair<int, unsigned int> >, value_type,
      boost::hash<std::vector<std::pair<int, unsigned int> > > >;

  /**
   * @brief Class constructor
   */
  FermiHubbardMomentumTwoU1(
      const Lattice& lat_, BaseParameters& parms_,
      bool isTranscorrelated = false
  )
      : lat(lat_),
        parms(parms_),
        tag_handler(new TagHandler<Matrix, TwoU1>()),
        order(lat_.size()),
        reverseOrder(lat_.size()) {
    // Definition of the charges (i.e., the relevant QN)
    TwoU1::charge A(0), B(0), C(0), D(1);
    B[0] = 1;
    C[1] = 1;
    phys.insert(std::make_pair(A, 1));
    phys.insert(std::make_pair(B, 1));
    phys.insert(std::make_pair(C, 1));
    phys.insert(std::make_pair(D, 1));
    // Elementary operators
    // (note that we define also these that are needed for the entropy
    // measurements)
    op_t create_up_op, create_down_op, destroy_up_op, destroy_down_op, fill_op,
        ident_op;
    op_t count_up_op, count_down_op, count_up_down_op, docc_op, d2e_op, e2d_op;
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
    count_up_op.insert_block(Matrix(1, 1, 1), B, B);
    count_up_op.insert_block(Matrix(1, 1, 1), D, D);
    count_down_op.insert_block(Matrix(1, 1, 1), C, C);
    count_down_op.insert_block(Matrix(1, 1, 1), D, D);
    count_up_down_op.insert_block(Matrix(1, 1, 1), B, B);
    count_up_down_op.insert_block(Matrix(1, 1, 1), C, C);
    count_up_down_op.insert_block(Matrix(1, 1, 2), D, D);
    docc_op.insert_block(Matrix(1, 1, 1), D, D);
    e2d_op.insert_block(Matrix(1, 1, 1), A, D);
    d2e_op.insert_block(Matrix(1, 1, 1), D, A);
    /* Create operator tag table */
#define REGISTER(op, kind) op = tag_handler->register_op(op##_op, kind);
    REGISTER(ident, tag_detail::bosonic)
    REGISTER(fill, tag_detail::bosonic)
    REGISTER(create_up, tag_detail::fermionic)
    REGISTER(create_down, tag_detail::fermionic)
    REGISTER(destroy_up, tag_detail::fermionic)
    REGISTER(destroy_down, tag_detail::fermionic)
    REGISTER(count_up, tag_detail::bosonic)
    REGISTER(count_down, tag_detail::bosonic)
    REGISTER(count_up_down, tag_detail::bosonic)
    REGISTER(docc, tag_detail::bosonic)
    REGISTER(e2d, tag_detail::bosonic)
    REGISTER(d2e, tag_detail::bosonic)
#undef REGISTER
    std::pair<tag_type, value_type> cutf =
        tag_handler->get_product_tag(create_up, fill);
    std::pair<tag_type, value_type> cdtf =
        tag_handler->get_product_tag(create_down, fill);
    std::pair<tag_type, value_type> ftdu =
        tag_handler->get_product_tag(fill, destroy_up);
    std::pair<tag_type, value_type> ftdd =
        tag_handler->get_product_tag(fill, destroy_down);
    std::pair<tag_type, value_type> cund =
        tag_handler->get_product_tag(create_up, count_down);
    std::pair<tag_type, value_type> dund =
        tag_handler->get_product_tag(destroy_up, count_down);
    std::pair<tag_type, value_type> cdnu =
        tag_handler->get_product_tag(create_down, count_up);
    std::pair<tag_type, value_type> ddnu =
        tag_handler->get_product_tag(destroy_down, count_up);
    std::pair<tag_type, value_type> cundtf =
        tag_handler->get_product_tag(cund.first, fill);
    std::pair<tag_type, value_type> ftdund =
        tag_handler->get_product_tag(fill, dund.first);
    std::pair<tag_type, value_type> cdnutf =
        tag_handler->get_product_tag(cdnu.first, fill);
    std::pair<tag_type, value_type> ftddnu =
        tag_handler->get_product_tag(fill, ddnu.first);
    std::pair<tag_type, value_type> ddcu =
        tag_handler->get_product_tag(destroy_down, create_up);
    std::pair<tag_type, value_type> ducd =
        tag_handler->get_product_tag(destroy_up, create_down);
    //#define HERMITIAN(op1, op2) tag_handler->hermitian_pair(op1, op2);
    //        HERMITIAN(create_up, destroy_up)
    //        HERMITIAN(create_down, destroy_down)
    //        HERMITIAN(e2d, d2e);
    //        HERMITIAN(cutf.first, ftdu.first)
    //        HERMITIAN(cdtf.first, ftdd.first)
    //        HERMITIAN(cund.first, dund.first)
    //        HERMITIAN(cdnu.first, ddnu.first)
    //        HERMITIAN(cundtf.first, ftdund.first)
    //        HERMITIAN(cdnutf.first, ftddnu.first)
    //        HERMITIAN(ddcu.first, ducd.first)
    //#undef HERMITIAN
    // General variables
    value_type U = parms["U_FermiHubbard"];
    value_type t = parms["t_FermiHubbard"];
    int width = parms["width_FermiHubbard"];
    int height = parms["height_FermiHubbard"];
    int kxMin = (height % 2 == 0) ? -height / 2 + 1 : -height / 2;
    int kxMax = (height % 2 == 0) ? height / 2 : height / 2;
    int kyMin = (width % 2 == 0) ? -width / 2 + 1 : -width / 2;
    int kyMax = (width % 2 == 0) ? width / 2 : width / 2;
    if (isTranscorrelated && !parms.is_set("J_Transcorrelated"))
      throw std::runtime_error("Please set the transcorrelation parameter");
    value_type J =
        parms.is_set("J_Transcorrelated") ? parms["J_Transcorrelated"] : 0.;
    // Orbital sorting
    if (!parms.is_set("orbital_order")) {
      for (pos_t p = 0; p < width * height; ++p) {
        order[p] = p;
        reverseOrder[p] = p;
      }
    } else {
      order = parms["orbital_order"].as<std::vector<int> >();
      for (pos_t p = 0; p < width * height; ++p) {
        auto it = std::find(order.begin(), order.end(), p);
        if (it != order.end())
          reverseOrder[p] = std::distance(order.begin(), it);
        else
          throw std::runtime_error("Incorrect sorting");
      }
    }
    // == CONSTRUCTION OF THE OPERATOR ==
    auto jw = JordanWignerHandler<Matrix, TwoU1>(
        lat, fill, create_up, create_down, destroy_up, destroy_down
    );
    MapOfOperatorsType mapOfOperators;
    for (int iSite = 0; iSite < width * height; iSite++) {
      auto kVector = getK(iSite, width, height);
      std::vector<pos_t> posVector = {reverseOrder[iSite], reverseOrder[iSite]};
      std::vector<OperatorType> opVector = {
          OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
      value_type coeff = -t * getEnergy(kVector, width, height, t);
      term_descriptor term =
          jw.getTerm(posVector, opVector, tag_handler, true, coeff);
      // this->terms_.push_back(term);
      addTerm(mapOfOperators, term);
      opVector = {OperatorType::CreateBeta, OperatorType::DestroyBeta};
      term = jw.getTerm(posVector, opVector, tag_handler, true, coeff);
      addTerm(mapOfOperators, term);
    }
    // Two-body operator
    for (int iSite = 0; iSite < height * width; iSite++) {
      for (int jSite = 0; jSite < height * width; jSite++) {
        for (int deltaSite = 0; deltaSite < height * width; deltaSite++) {
          auto iPair = getK(iSite, width, height);
          auto jPair = getK(jSite, width, height);
          auto deltaPair = getK(deltaSite, width, height);
          int iRow1 = iPair.first - deltaPair.first,
              iRow2 = jPair.first + deltaPair.first;
          int iCol1 = iPair.second - deltaPair.second,
              iCol2 = jPair.second + deltaPair.second;
          if (iRow1 < kxMin) iRow1 += height;
          if (iRow1 > kxMax) iRow1 -= height;
          if (iRow2 < kxMin) iRow2 += height;
          if (iRow2 > kxMax) iRow2 -= height;
          if (iCol1 < kyMin) iCol1 += width;
          if (iCol1 > kyMax) iCol1 -= width;
          if (iCol2 < kyMin) iCol2 += width;
          if (iCol2 > kyMax) iCol2 -= width;
          iRow1 -= kxMin;
          iRow2 -= kxMin;
          iCol1 -= kyMin;
          iCol2 -= kyMin;
          int index1 = iRow1 + iCol1 * height;
          int index2 = iRow2 + iCol2 * height;
          auto pair1 = getK(index1, width, height);
          auto pair2 = getK(index2, width, height);
          std::vector<OperatorType> opVector1 = {
              OperatorType::CreateAlpha, OperatorType::CreateBeta,
              OperatorType::DestroyBeta, OperatorType::DestroyAlpha};
          std::vector<OperatorType> opVector2 = {
              OperatorType::CreateBeta, OperatorType::CreateAlpha,
              OperatorType::DestroyAlpha, OperatorType::DestroyBeta};
          std::vector<pos_t> positions = {
              reverseOrder[index1], reverseOrder[index2], reverseOrder[jSite],
              reverseOrder[iSite]};
          value_type coeff =
              (!isTranscorrelated)
                  ? U / static_cast<value_type>(2 * height * width)
                  : (U / 2. - t * ((std::exp(J) - 1.) *
                                       getEnergy(pair1, width, height, t) +
                                   (std::exp(-J) - 1.) *
                                       getEnergy(iPair, width, height, t))) /
                        static_cast<value_type>(height * width);
          auto term1 =
              jw.getTerm(positions, opVector1, tag_handler, true, coeff);
          // this->terms_.push_back(term1);
          addTerm(mapOfOperators, term1);
          auto term2 =
              jw.getTerm(positions, opVector2, tag_handler, true, coeff);
          // this->terms_.push_back(term2);
          addTerm(mapOfOperators, term2);
        }
      }
    }
    // Three-body operator
    if (isTranscorrelated && parms["transcorrelated_3body"] == "yes") {
      for (int iSite = 0; iSite < height * width; iSite++) {
        for (int jSite = 0; jSite < height * width; jSite++) {
          for (int kSite = 0; kSite < height * width; kSite++) {
            for (int deltaSite1 = 0; deltaSite1 < height * width;
                 deltaSite1++) {
              for (int deltaSite2 = 0; deltaSite2 < height * width;
                   deltaSite2++) {
                auto iPair = getK(iSite, width, height);
                auto jPair = getK(jSite, width, height);
                auto kPair = getK(kSite, width, height);
                auto deltaPair1 = getK(deltaSite1, width, height);
                auto deltaPair2 = getK(deltaSite2, width, height);
                int iRow1 = iPair.first - deltaPair1.first,
                    iRow2 = jPair.first + deltaPair2.first,
                    iRow3 = kPair.first + deltaPair1.first - deltaPair2.first,
                    iRow4 = iPair.first - deltaPair1.first + deltaPair2.first;
                int iCol1 = iPair.second - deltaPair1.second,
                    iCol2 = jPair.second + deltaPair2.second,
                    iCol3 =
                        kPair.second + deltaPair1.second - deltaPair2.second,
                    iCol4 =
                        iPair.second - deltaPair1.second + deltaPair2.second;
                while (iRow1 < kxMin) iRow1 += height;
                while (iRow1 > kxMax) iRow1 -= height;
                while (iRow2 < kxMin) iRow2 += height;
                while (iRow2 > kxMax) iRow2 -= height;
                while (iRow3 < kxMin) iRow3 += height;
                while (iRow3 > kyMax) iRow3 -= height;
                while (iRow4 < kxMin) iRow4 += height;
                while (iRow4 > kxMax) iRow4 -= height;
                while (iCol1 < kyMin) iCol1 += width;
                while (iCol1 > kyMax) iCol1 -= width;
                while (iCol2 < kyMin) iCol2 += width;
                while (iCol2 > kyMax) iCol2 -= width;
                while (iCol3 < kyMin) iCol3 += width;
                while (iCol3 > kyMax) iCol3 -= width;
                while (iCol4 < kyMin) iCol4 += width;
                while (iCol4 > kyMax) iCol4 -= width;

                iRow1 -= kxMin;
                iRow2 -= kxMin;
                iRow3 -= kxMin;
                iRow4 -= kxMin;

                iCol1 -= kyMin;
                iCol2 -= kyMin;
                iCol3 -= kyMin;
                iCol4 -= kyMin;

                int index1 = iRow1 + iCol1 * height;
                int index2 = iRow2 + iCol2 * height;
                int index3 = iRow3 + iCol3 * height;
                int index4 = iRow4 + iCol4 * height;

                auto pair4 = getK(index4, width, height);

                std::vector<pos_t> positions = {
                    reverseOrder[index1], reverseOrder[index2],
                    reverseOrder[index3], reverseOrder[kSite],
                    reverseOrder[jSite],  reverseOrder[iSite]};
                std::vector<OperatorType> opVector1 = {
                    OperatorType::CreateAlpha, OperatorType::CreateBeta,
                    OperatorType::CreateBeta,  OperatorType::DestroyBeta,
                    OperatorType::DestroyBeta, OperatorType::DestroyAlpha};
                std::vector<OperatorType> opVector2 = {
                    OperatorType::CreateBeta,   OperatorType::CreateAlpha,
                    OperatorType::CreateAlpha,  OperatorType::DestroyAlpha,
                    OperatorType::DestroyAlpha, OperatorType::DestroyBeta};

                value_type coeff =
                    2. * t * (std::cosh(J) - 1.) *
                    getEnergy(pair4, width, height, t) /
                    static_cast<value_type>(height * height * width * width);
                if (std::abs(coeff) > 1.0E-10 && !(index2 == index3) &&
                    !(kSite == jSite)) {
                  auto term1 = jw.getTerm(
                      positions, opVector1, tag_handler, true, coeff
                  );
                  // this->terms_.push_back(term1);
                  addTerm(mapOfOperators, term1);
                  auto term2 = jw.getTerm(
                      positions, opVector2, tag_handler, true, coeff
                  );
                  // this->terms_.push_back(term2);
                  addTerm(mapOfOperators, term2);
                }
              }
            }
          }
        }
      }
    }
    this->terms_.reserve(mapOfOperators.size());
    for (const auto& idx : mapOfOperators)
      this->terms_.push_back(term_descriptor(idx.first, idx.second, true));
  }

  void addTerm(MapOfOperatorsType& mapOfOperators, const term_descriptor& term)
      const {
    if (mapOfOperators.find(term.getBase()) == mapOfOperators.end())
      mapOfOperators.insert({term.getBase(), term.coeff});
    else
      mapOfOperators[term.getBase()] += term.coeff;
  }

  /** @brief Gets the k vector of a given point in the real-space lattice */
  std::pair<int, int> getK(int iSite, int width, int height) const {
    int iRow = (height % 2 == 0) ? (iSite % height - height / 2 + 1)
                                 : (iSite % height - height / 2);
    int iCol = (height % 2 == 0) ? (iSite / height - width / 2 + 1)
                                 : (iSite / height) % width - width / 2;
    return std::make_pair(iRow, iCol);
  }

  /** @brief Gets the energy of a given k vector value */
  value_type getEnergy(
      std::pair<int, int> kVector, int width, int height, value_type t
  ) const {
    value_type kx = 2. * M_PI * (static_cast<value_type>(kVector.first)) /
                    static_cast<value_type>(height);
    value_type ky = 2. * M_PI * (static_cast<value_type>(kVector.second)) /
                    static_cast<value_type>(width);
    return 2. * cos(kx) + 2. * cos(ky);
  }

  /** @brief Update method */
  void update(BaseParameters const& p) {
    // TODO: update this->terms_ with the new parameters
    throw std::runtime_error("update() not yet implemented for this model.");
    return;
  }

  /** @brief Getter for the physical dimensions */
  Index<TwoU1> const& phys_dim(size_t type) const { return phys; }

  /**
   * @brief Getter for the measurements
   * This implements the entanglement measurements
   */
  measurements_type measurements() const {
    measurements_type meas;
    /*
    // TODO COMMENTED OUT ONLY FOR THE MOMENT
    auto jw = JordanWignerHandler<Matrix, TwoU1>(lat, fill, create_up,
    create_down, destroy_up, destroy_down);
    //
    std::vector<tag_type> ident_ops = std::vector<tag_type>(1, ident);
    std::vector<tag_type> fill_ops = std::vector<tag_type>(1, fill);
    //
    if (parms.is_set("MEASURE[ChemEntropy]")) {
    for (int p1 = 0; p1 < lat.size(); p1++) {
        // == SINGLE-ORBITAL ENTROPY ==
        // Nup
        std::vector< pos_t > posVector = { reverseOrder[p1], reverseOrder[p1] };
        std::vector< OperatorType > opVector = {OperatorType::CreateAlpha,
    OperatorType::DestroyAlpha}; terms_type termsVector = {jw.getTerm(posVector,
    opVector, tag_handler, true, 1.)}; generate_mpo::TaggedMPOMaker<Matrix,
    TwoU1> mpoLocal1(lat, ident_ops, ident_ops, fill_ops, tag_handler,
    termsVector); MPO<Matrix, TwoU1> mpo1 = mpoLocal1.create_mpo();
        meas.push_back( new measurements::expvalMeas<Matrix,
    TwoU1>("Nup_"+std::to_string(p1), mpo1) );
        // Ndown
        opVector = {OperatorType::CreateBeta, OperatorType::DestroyBeta};
        termsVector = {jw.getTerm(posVector, opVector, tag_handler, true, 1.)};
        generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2(lat, ident_ops,
    ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2.create_mpo(); meas.push_back( new measurements::expvalMeas<Matrix,
    TwoU1>("Ndown_"+std::to_string(p1), mpo1) );
        // NupNdown
        opVector = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
    OperatorType::CreateBeta, OperatorType::DestroyBeta}; posVector = {
    reverseOrder[p1], reverseOrder[p1], reverseOrder[p1], reverseOrder[p1] };
        termsVector = {jw.getTerm(posVector, opVector, tag_handler, true, 1.)};
        generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal3(lat, ident_ops,
    ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal3.create_mpo(); meas.push_back( new measurements::expvalMeas<Matrix,
    TwoU1>("Nupdown_"+std::to_string(p1), mpo1) ); for (int p2 = 0; p2 <
    lat.size(); p2++) {
            // == TWO-ORBITAL ENTROPY ==
            // dm_up
            opVector = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
            posVector = { reverseOrder[p1], reverseOrder[p2] };
            termsVector = {jw.getTerm(posVector, opVector, tag_handler,
    true, 1.)}; generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_0(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_0.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("dm_up_"+std::to_string(p1)+"_"+std::to_string(p2), mpo1) );
            // dm_down
            opVector = {OperatorType::CreateBeta, OperatorType::DestroyBeta};
            posVector = { reverseOrder[p1], reverseOrder[p2] };
            termsVector = {jw.getTerm(posVector, opVector, tag_handler,
    true, 1.)}; generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_00(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_00.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("dm_down_"+std::to_string(p1)+"_"+std::to_string(p2), mpo1) );
            // NupNup
            opVector = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
    OperatorType::CreateAlpha, OperatorType::DestroyAlpha}; posVector = {
    reverseOrder[p1], reverseOrder[p1], reverseOrder[p2], reverseOrder[p2] };
            termsVector = {jw.getTerm(posVector, opVector, tag_handler,
    true, 1.)}; generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_1(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_1.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("NupNup_"+std::to_string(p1)+"_"+std::to_string(p2), mpo1) );
            // NupNdown
            opVector = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
    OperatorType::CreateBeta, OperatorType::DestroyBeta}; termsVector =
    {jw.getTerm(posVector, opVector, tag_handler, true, 1.)};
            generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_2(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_2.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("NupNdown_"+std::to_string(p1)+"_"+std::to_string(p2), mpo1) );
            // NdownNup
            opVector = {OperatorType::CreateBeta, OperatorType::DestroyBeta,
    OperatorType::CreateAlpha, OperatorType::DestroyAlpha}; termsVector =
    {jw.getTerm(posVector, opVector, tag_handler, true, 1.)};
            generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_3(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_3.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("NdownNup_"+std::to_string(p1)+"_"+std::to_string(p2), mpo1) );
            // NdownNdown
            opVector = {OperatorType::CreateBeta, OperatorType::DestroyBeta,
    OperatorType::CreateBeta, OperatorType::DestroyBeta}; termsVector =
    {jw.getTerm(posVector, opVector, tag_handler, true, 1.)};
            generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_4(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_4.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("NdownNdown_"+std::to_string(p1)+"_"+std::to_string(p2), mpo1) );
            // NupdownNupdown
            opVector = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
    OperatorType::CreateBeta, OperatorType::DestroyBeta,
                        OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
    OperatorType::CreateBeta, OperatorType::DestroyBeta}; posVector = {
    reverseOrder[p1], reverseOrder[p1], reverseOrder[p1], reverseOrder[p1],
    reverseOrder[p2], reverseOrder[p2], reverseOrder[p2], reverseOrder[p2] };
            termsVector = {jw.getTerm(posVector, opVector, tag_handler,
    true, 1.)}; generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_5(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_5.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("NupdownNupdown_"+std::to_string(p1)+"_"+std::to_string(p2), mpo1) );
            // cdag_up*Ndown:c_up*Ndown
            opVector = {OperatorType::CreateAlpha, OperatorType::CreateBeta,
    OperatorType::DestroyBeta, OperatorType::DestroyAlpha,
    OperatorType::CreateBeta, OperatorType::DestroyBeta}; posVector = {
    reverseOrder[p1], reverseOrder[p1], reverseOrder[p1], reverseOrder[p2],
    reverseOrder[p2], reverseOrder[p2] }; termsVector = {jw.getTerm(posVector,
    opVector, tag_handler, true, 1.)}; generate_mpo::TaggedMPOMaker<Matrix,
    TwoU1> mpoLocal2_6(lat, ident_ops, ident_ops, fill_ops, tag_handler,
    termsVector); mpo1 = mpoLocal2_6.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("transfer_up_while_down_"+std::to_string(p1)+"_"+std::to_string(p2),
    mpo1) );
            // cdag_down*Nup:c_down*Nup
            opVector = {OperatorType::CreateBeta, OperatorType::CreateAlpha,
    OperatorType::DestroyAlpha, OperatorType::DestroyBeta,
    OperatorType::CreateAlpha, OperatorType::DestroyAlpha}; posVector = {
    reverseOrder[p1], reverseOrder[p1], reverseOrder[p1], reverseOrder[p2],
    reverseOrder[p2], reverseOrder[p2] }; termsVector = {jw.getTerm(posVector,
    opVector, tag_handler, true, 1.)}; generate_mpo::TaggedMPOMaker<Matrix,
    TwoU1> mpoLocal2_7(lat, ident_ops, ident_ops, fill_ops, tag_handler,
    termsVector); mpo1 = mpoLocal2_7.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("transfer_down_while_up_"+std::to_string(p1)+"_"+std::to_string(p2),
    mpo1) );
            // cdag_up:c_up*Ndown
            opVector = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
    OperatorType::CreateBeta, OperatorType::DestroyBeta}; posVector = {
    reverseOrder[p1], reverseOrder[p2], reverseOrder[p2], reverseOrder[p2] };
            termsVector = {jw.getTerm(posVector, opVector, tag_handler,
    true, 1.)}; generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_8(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_8.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("transfer_up_while_down_at_2_"+std::to_string(p1)+"_"+std::to_string(p2),
    mpo1) );
            // cdag_up*Ndown:c_up
            opVector = {OperatorType::CreateAlpha, OperatorType::CreateBeta,
    OperatorType::DestroyBeta, OperatorType::DestroyAlpha}; posVector = {
    reverseOrder[p1], reverseOrder[p1], reverseOrder[p1], reverseOrder[p2] };
            termsVector = {jw.getTerm(posVector, opVector, tag_handler,
    true, 1.)}; generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_9(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_9.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("transfer_up_while_down_at_1_"+std::to_string(p1)+"_"+std::to_string(p2),
    mpo1) );
            // cdag_down:c_down*Nup
            opVector = {OperatorType::CreateBeta, OperatorType::DestroyBeta,
    OperatorType::CreateAlpha, OperatorType::DestroyAlpha}; posVector = {
    reverseOrder[p1], reverseOrder[p2], reverseOrder[p2], reverseOrder[p2] };
            termsVector = {jw.getTerm(posVector, opVector, tag_handler,
    true, 1.)}; generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_10(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_10.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("transfer_down_while_up_at_2_"+std::to_string(p1)+"_"+std::to_string(p2),
    mpo1) );
            // cdag_down*Nup:c_down
            opVector = {OperatorType::CreateBeta, OperatorType::CreateAlpha,
    OperatorType::DestroyAlpha, OperatorType::DestroyBeta}; posVector = {
    reverseOrder[p1], reverseOrder[p1], reverseOrder[p1], reverseOrder[p2] };
            termsVector = {jw.getTerm(posVector, opVector, tag_handler,
    true, 1.)}; generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_11(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_11.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("transfer_down_while_up_at_1_"+std::to_string(p1)+"_"+std::to_string(p2),
    mpo1) );
            // cdag_up*cdag_down:c_up*c_down
            opVector = {OperatorType::CreateAlpha, OperatorType::CreateBeta,
    OperatorType::DestroyAlpha, OperatorType::DestroyBeta}; posVector = {
    reverseOrder[p1], reverseOrder[p1], reverseOrder[p2], reverseOrder[p2] };
            termsVector = {jw.getTerm(posVector, opVector, tag_handler,
    true, 1.)}; generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_12(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_12.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("transfer_pair_"+std::to_string(p1)+"_"+std::to_string(p2), mpo1) );
            // cdag_up*c_down:cdag_down*c_up
            opVector = {OperatorType::CreateAlpha, OperatorType::DestroyBeta,
    OperatorType::CreateBeta, OperatorType::DestroyAlpha}; posVector = {
    reverseOrder[p1], reverseOrder[p1], reverseOrder[p2], reverseOrder[p2] };
            termsVector = {jw.getTerm(posVector, opVector, tag_handler,
    true, 1.)}; generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_13(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_13.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("spinflip_"+std::to_string(p1)+"_"+std::to_string(p2), mpo1) );
            // Nup:Nup*Ndown
            opVector = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
    OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
                        OperatorType::CreateBeta, OperatorType::DestroyBeta};
            posVector = { reverseOrder[p1], reverseOrder[p1], reverseOrder[p2],
    reverseOrder[p2], reverseOrder[p2], reverseOrder[p2] }; termsVector =
    {jw.getTerm(posVector, opVector, tag_handler, true, 1.)};
            generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_14(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_14.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("nupdocc_"+std::to_string(p1)+"_"+std::to_string(p2), mpo1) );
            // Ndown:Nup*Ndown
            opVector = {OperatorType::CreateBeta, OperatorType::DestroyBeta,
    OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
                        OperatorType::CreateBeta, OperatorType::DestroyBeta};
            posVector = { reverseOrder[p1], reverseOrder[p1], reverseOrder[p2],
    reverseOrder[p2], reverseOrder[p2], reverseOrder[p2] }; termsVector =
    {jw.getTerm(posVector, opVector, tag_handler, true, 1.)};
            generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_15(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_15.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("ndowndocc_"+std::to_string(p1)+"_"+std::to_string(p2), mpo1) );
            // Nup*Ndown:Nup
            opVector = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
    OperatorType::CreateBeta, OperatorType::DestroyBeta,
                        OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
            posVector = { reverseOrder[p1], reverseOrder[p1], reverseOrder[p1],
    reverseOrder[p1], reverseOrder[p2], reverseOrder[p2] }; termsVector =
    {jw.getTerm(posVector, opVector, tag_handler, true, 1.)};
            generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_16(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_16.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("doccnup_"+std::to_string(p1)+"_"+std::to_string(p2), mpo1) );
            // Nup*Ndown:Ndown
            opVector = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
    OperatorType::CreateBeta, OperatorType::DestroyBeta,
                        OperatorType::CreateBeta, OperatorType::DestroyBeta};
            posVector = { reverseOrder[p1], reverseOrder[p1], reverseOrder[p1],
    reverseOrder[p1], reverseOrder[p2], reverseOrder[p2] }; termsVector =
    {jw.getTerm(posVector, opVector, tag_handler, true, 1.)};
            generate_mpo::TaggedMPOMaker<Matrix, TwoU1> mpoLocal2_17(lat,
    ident_ops, ident_ops, fill_ops, tag_handler, termsVector); mpo1 =
    mpoLocal2_17.create_mpo(); meas.push_back( new
    measurements::expvalMeas<Matrix,
    TwoU1>("doccndown_"+std::to_string(p1)+"_"+std::to_string(p2), mpo1) );
        }
    }
    }
    */
    return meas;
  }

  /** @brief Getter for the identity matrix tag */
  tag_type identity_matrix_tag(size_t type) const { return ident; }

  /** @brief Getter for the filling matrix tag */
  tag_type filling_matrix_tag(size_t type) const { return fill; }

  /** @brief Getter for the charge that is conserved by the MPS */
  typename TwoU1::charge total_quantum_numbers(BaseParameters& parms) const {
    typename TwoU1::charge ret(0);
    ret[0] = static_cast<int>(parms["u1_total_charge1"]);
    ret[1] = static_cast<int>(parms["u1_total_charge2"]);
    return ret;
  }

  /** @brief Getter for the operator given its string representation */
  tag_type get_operator_tag(std::string const& name, size_t type) const {
    if (name == "create_up")
      return create_up;
    else if (name == "create_down")
      return create_down;
    else if (name == "destroy_up")
      return destroy_up;
    else if (name == "destroy_down")
      return destroy_down;
    else if (name == "count_up")
      return count_up;
    else if (name == "count_down")
      return count_down;
    else
      throw std::runtime_error("Operator not valid for this model.");
    return 0;
  }

  /** @brief Getter for the operator map */
  table_ptr operators_table() const { return tag_handler; }

 private:
  Index<TwoU1> phys;
  Lattice const& lat;
  BaseParameters& parms;
  std::shared_ptr<TagHandler<Matrix, TwoU1> > tag_handler;
  tag_type create_up, create_down, destroy_up, destroy_down, ident, fill,
      count_up, count_down, count_up_down, docc, d2e, e2d;
  std::vector<int> order, reverseOrder;
};

#endif
