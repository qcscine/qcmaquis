#ifndef MAQUIS_SRCAS_DETERMINANT_H
#define MAQUIS_SRCAS_DETERMINANT_H

#include "onv.h"
#include <set>

namespace maquis {
namespace srcas {

/** @brief enum to specify the spin of an orbital */
enum class Spin { alpha, beta };

/**
 * @class Determinant
 * @brief Handle excitations from a queen.
 */
class Determinant : public virtual ONV {
 public:
  /** @brief default constructor */
  Determinant() = default;
  /**
   * @brief Construct a determiant from a vector filled with 4, 3, 2, 1
   * for doubly, alpha, beta, unoccupied orbitals, respectively.
   *
   * @param onvVec vector<int> filled with 4, 3, 2, 1
   */
  Determinant(const std::vector<int>& onvVec);
  /**
   * @brief return the vector representation of a determinant
   */
  std::vector<int> vector() const override;
  /**
   * @brief Make two determiants comparable
   *
   * @param rhs other determinants
   * @return true if both determiants have same alpha and beta occupation
   */
  bool operator==(const Determinant& rhs) const; //  override;
  /**
   * @brief excite electron from one orbital to another
   *
   * @param from occupied orbital index
   * @param to virtual orbital index
   * @param spin alpha or beta
   */
  void excite_electron(int from, int to, Spin spin);
  /**
   * @brief get the occupied orbital indeces for corresponding spin
   *
   * @param spin alpha or beta
   */
  const std::vector<int>& occupied(Spin spin) const;
  /**
   * @brief get the unoccupied orbital indeces for corresponding spin
   *
   * @param spin alpha or beta
   */
  const std::vector<int>& unoccupied(Spin spin) const;

 private:
  /** @brief store occupied alpha orbital indices */
  std::vector<int> alpha_occupied_;
  /** @brief store occupied beta orbital indices */
  std::vector<int> beta_occupied_;
  /** @brief store unoccupied alpha orbital indices */
  std::vector<int> alpha_unoccupied_;
  /** @brief store unoccupied alpha orbital indices */
  std::vector<int> beta_unoccupied_;
};

} // namespace srcas
} // namespace maquis

#endif
