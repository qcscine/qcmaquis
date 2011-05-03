
#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <cstddef>

struct cmp_with_prefactor {
	static double prefactor;
	bool operator() (std::size_t i, std::size_t j) {
		bool ret = (i < j);
		if (ret) prefactor *= -1.;
		return ret;
	}
};


#endif /* UTILS_HPP_ */
