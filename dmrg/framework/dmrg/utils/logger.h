/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#ifndef LOGGER_H
#define LOGGER_H

template <class Archive>
class Logger {
 public:
  template <class T>
  void operator<<(const T& p) {}

  template <class OtherArchive>
  void save(OtherArchive& ar) const {}
};

#endif
