/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2020 Laboratory for Physical Chemistry, ETH Zurich
 *               2020- by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef JORDAN_WIGNER_MANAGER_H
#define JORDAN_WIGNER_MANAGER_H

enum class OperatorType {
  Normal,
  CreateAlpha,
  CreateBeta,
  DestroyAlpha,
  DestroyBeta,
  Filling
};

namespace OperatorTypeUtils {

inline bool isAlpha(OperatorType tmp) {
  return tmp == OperatorType::CreateAlpha || tmp == OperatorType::DestroyAlpha;
}
inline bool isBeta(OperatorType tmp) {
  return tmp == OperatorType::CreateBeta || tmp == OperatorType::DestroyBeta;
}
inline bool isCreate(OperatorType tmp) {
  return tmp == OperatorType::CreateAlpha || tmp == OperatorType::CreateBeta;
}
inline bool isDestroy(OperatorType tmp) {
  return tmp == OperatorType::DestroyAlpha || tmp == OperatorType::DestroyBeta;
}

}  // namespace OperatorTypeUtils

template <class TagType, class PosType>
class OperatorAndPosition {
 public:
  /* Constructor */
  OperatorAndPosition(OperatorType opType, TagType tagType, PosType posType)
      : opType_(opType), tagType_(tagType), posType_(posType) {}
  /* Members */
  OperatorType opType_;
  TagType tagType_;
  PosType posType_;
};

/**
 * @brief Class containing the methods to implement the JW transformation.
 */
template <class Matrix, class SymmGroup>
class JordanWignerHandler {
 public:
  // Types declaration
  using value_type = typename Matrix::value_type;
  using term_descriptor = ::term_descriptor<value_type>;
  using tag_type = typename TagHandler<Matrix, SymmGroup>::tag_type;
  using pos_t = Lattice::pos_t;
  using sc_t = typename SymmGroup::subcharge;
  using PairType = std::pair<pos_t, OperatorType>;

  /**
   * @brief Base constructor
   */
  JordanWignerHandler(
      const Lattice& lat, const std::vector<tag_type>& fill_op,
      const std::vector<tag_type>& create_up_op,
      const std::vector<tag_type>& create_down_op,
      const std::vector<tag_type>& destroy_up_op,
      const std::vector<tag_type>& destroy_down_op
  )
      : lat_(lat),
        fillOp_(fill_op),
        createUpOp_(create_up_op),
        createDownOp_(create_down_op),
        destroyUpOp_(destroy_up_op),
        destroyDownOp_(destroy_down_op) {}

  /**
   * @brief Construct from a single set of tags
   */
  JordanWignerHandler(
      const Lattice& lat, tag_type fill_op, tag_type create_up_op,
      tag_type create_down_op, tag_type destroy_up_op, tag_type destroy_down_op
  )
      : lat_(lat) {
    fillOp_ = {fill_op};
    createUpOp_ = {create_up_op};
    createDownOp_ = {create_down_op};
    destroyUpOp_ = {destroy_up_op};
    destroyDownOp_ = {destroy_down_op};
  }

  /** @brief Method to generate the term_descriptor of a given SQ operator
   * string. */
  term_descriptor getTerm(
      std::vector<pos_t> positions, std::vector<OperatorType> tags,
      std::shared_ptr<TagHandler<Matrix, SymmGroup> > op_table, bool sign,
      value_type scale
  ) const {
    // Variables declaration
    assert(positions.size() == tags.size());
    int size = positions.size();
    std::vector<OperatorAndPosition<pos_t, tag_type> > localOperatorBuffer;
    // Sorts and calculate the sign
    std::vector<PairType> tmpVector;
    tmpVector.reserve(size);
    for (int idx = 0; idx < size; idx++)
      tmpVector.push_back(std::make_pair(positions[idx], tags[idx]));
    double scalarSign = this->calculatePermutationSign(tmpVector);
    auto compareLambda = [](const PairType& a, const PairType& b) {
      return a.first < b.first;
    };
    std::stable_sort(tmpVector.begin(), tmpVector.end(), compareLambda);
    // Generation of the overall operator
    for (int idx = 0; idx < size; idx++) {
      auto pos = tmpVector[idx].first;
      auto opTag = tmpVector[idx].second;
      // -- Creation operators --
      if (OperatorTypeUtils::isCreate(opTag)) {
        int initCycle = (OperatorTypeUtils::isAlpha(opTag)) ? pos - 1 : pos;
        auto tag = (OperatorTypeUtils::isAlpha(opTag))
                       ? createUpOp_[lat_.get_prop<sc_t>("type", pos)]
                       : createDownOp_[lat_.get_prop<sc_t>("type", pos)];
        localOperatorBuffer.push_back(
            OperatorAndPosition<pos_t, tag_type>(OperatorType::Normal, tag, pos)
        );
        for (int iFill = initCycle; iFill >= 0; iFill--)
          if (std::find(positions.begin(), positions.end(), iFill) !=
              positions.end())
            localOperatorBuffer.push_back(OperatorAndPosition<pos_t, tag_type>(
                OperatorType::Filling,
                fillOp_[lat_.get_prop<sc_t>("type", iFill)], iFill
            ));
      }
      // -- Annihilation operators --
      else if (OperatorTypeUtils::isDestroy(opTag)) {
        int endCycle = (OperatorTypeUtils::isAlpha(opTag)) ? pos - 1 : pos;
        auto tag = (OperatorTypeUtils::isAlpha(opTag))
                       ? destroyUpOp_[lat_.get_prop<sc_t>("type", pos)]
                       : destroyDownOp_[lat_.get_prop<sc_t>("type", pos)];
        for (int iFill = 0; iFill <= endCycle; iFill++)
          if (std::find(positions.begin(), positions.end(), iFill) !=
              positions.end())
            localOperatorBuffer.push_back(OperatorAndPosition<pos_t, tag_type>(
                OperatorType::Filling,
                fillOp_[lat_.get_prop<sc_t>("type", iFill)], iFill
            ));
        localOperatorBuffer.push_back(
            OperatorAndPosition<pos_t, tag_type>(OperatorType::Normal, tag, pos)
        );
      }
    }

    std::stable_sort(
        localOperatorBuffer.begin(), localOperatorBuffer.end(),
        [](const OperatorAndPosition<pos_t, tag_type>& a,
           const OperatorAndPosition<pos_t, tag_type>& b) {
          return a.posType_ < b.posType_;
        }
    );

    // Cleaning of the JW string.
    for (auto i = localOperatorBuffer.begin();
         i != localOperatorBuffer.end();) {
      auto n = std::next(i);
      if (n == localOperatorBuffer.end()) break;
      if (i->posType_ == n->posType_ && i->opType_ == n->opType_ &&
          i->opType_ == OperatorType::Filling) {
        i = localOperatorBuffer.erase(i);
        i = localOperatorBuffer.erase(i);
      } else {
        i++;
      }
    }

    // Generation of the term
    term_descriptor term;

    /*
    std::vector<pos_t> positions_final;
    std::vector<tag_type> operators_final;
    for (int idx = 0; idx < localOperatorBuffer.size(); idx++) {
      positions_final.push_back(localOperatorBuffer[idx].posType_);
      tag_type localTag = (localOperatorBuffer[idx].opType_ ==
    OperatorType::Filling) ? fillOp_[lat_.get_prop<sc_t>("type",
    localOperatorBuffer[idx].posType_)] : localOperatorBuffer[idx].tagType_;
      operators_final.push_back(localTag);
    }

    term = arrange_operators(positions_final, operators_final, op_table);
    */

    term.is_fermionic = sign;
    term.coeff = scale * scalarSign;
    bool isNull = false;

    for (int idx = 0; idx < localOperatorBuffer.size();) {
      tag_type localTag =
          (localOperatorBuffer[idx].opType_ == OperatorType::Filling)
              ? fillOp_[lat_.get_prop<sc_t>(
                    "type", localOperatorBuffer[idx].posType_
                )]
              : localOperatorBuffer[idx].tagType_;
      std::pair<tag_type, value_type> ptag = std::make_pair(localTag, 1.);
      int refPos = localOperatorBuffer[idx].posType_;
      if (idx != localOperatorBuffer.size() - 1) {
        bool isNextSame = (localOperatorBuffer[idx + 1].posType_ == refPos);
        while (isNextSame) {
          if (localOperatorBuffer[idx + 1].opType_ == OperatorType::Filling) {
            if (op_table->product_is_null(
                    fillOp_[lat_.template get_prop<sc_t>(
                        "type", localOperatorBuffer[idx + 1].posType_
                    )],
                    localTag
                )) {
              isNull = true;
              break;
            }
            ptag = op_table->get_product_tag(
                fillOp_[lat_.template get_prop<sc_t>(
                    "type", localOperatorBuffer[idx + 1].posType_
                )],
                localTag
            );
          } else {
            if (op_table->product_is_null(
                    localOperatorBuffer[idx + 1].tagType_, localTag
                )) {
              isNull = true;
              break;
            }
            ptag = op_table->get_product_tag(
                localOperatorBuffer[idx + 1].tagType_, localTag
            );
          }
          localTag = ptag.first;
          term.coeff *= ptag.second;
          idx++;
          if (idx == localOperatorBuffer.size() - 1)
            isNextSame = false;
          else
            isNextSame = (localOperatorBuffer[idx + 1].posType_ == refPos);
        }
      }
      if (!isNull)
        term.push_back(std::make_pair(refPos, localTag));
      else
        break;
      idx++;
    }
    // If a null operator has been found, just returns a zero-sized term
    if (isNull) term = term_descriptor();
    return term;
  }

 private:
  /**
   * @brief Simple function to check the parity of a permutation.
   */
  double calculatePermutationSign(const std::vector<PairType>& positionVector
  ) const {
    if (positionVector.size() == 1) return 1.;
    int transitionsCount = 0;
    for (int idx1 = 0; idx1 < positionVector.size() - 1; idx1++) {
      for (int idx2 = idx1 + 1; idx2 < positionVector.size(); idx2++) {
        if (positionVector[idx1].first > positionVector[idx2].first)
          transitionsCount += 1;
      }
    }
    return (transitionsCount % 2 == 0) ? 1. : -1.;
  }

  /* Class members */
  const Lattice& lat_;
  std::vector<tag_type> fillOp_, createUpOp_, createDownOp_, destroyUpOp_,
      destroyDownOp_;
};

#endif
