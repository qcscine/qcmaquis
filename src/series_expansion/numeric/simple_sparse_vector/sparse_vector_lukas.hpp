#ifndef SPARSE_VECTOR_HPP
#define SPARSE_VECTOR_HPP

#include <vector>
#include <iostream>

#include <boost/cstdint.hpp>
#include <boost/shared_array.hpp>

namespace blas
{

template <typename T, std::size_t L>
class sparse_vector {
	static const boost::uint64_t one = 0x01;
	
	class proxy {
		public:
	
			proxy(sparse_vector<T, L> & o, std::size_t i)
				: obj(o)
				, index(i)
			{}
	
            inline T& get_ref() const
            {
                return obj[index];
            }

			operator T() const {
				return get_ref();
			}
	
			T operator=(T const & v) {
				obj.set(v, index);
				return v;
			}

            T operator += (T const& v) {
                obj.set(get_ref() += v, index);
                return get_ref();
            }

            T operator -= (T const& v) {
                obj.set(get_ref() -= v, index);
                return get_ref();
            }
            
            T operator *= (T const& v) {
                obj.set(get_ref() *= v, index);
                return get_ref();
            }

            T operator /= (T const& v) {
                obj.set(get_ref() /= v, index);
                return get_ref();
            } 


		private:
			sparse_vector<T, L> & obj;
			std::size_t index;
	};
	
	public:

		typedef T           value_type;
        typedef proxy       reference_type;
        typedef std::size_t size_type;

		static std::size_t size() {
			return L;
		}

		sparse_vector()
			: last(0)
			, data()
			, index(new boost::uint64_t[(L + 63) / 64])
		{
			std::memset(index.get(), 0, static_cast<std::size_t>((L + 7) / 8));
		}

		sparse_vector(size_type t)
			: last(0)
			, data()
			, index(new boost::uint64_t[(L + 63) / 64])
		{
            assert( t <=L );
			std::memset(index.get(), 0, static_cast<std::size_t>((L + 7) / 8));
		}
	
		sparse_vector(sparse_vector<T, L> const & rhs)
			: last(rhs.last)
			, data(rhs.data)
			, index(new boost::uint64_t[(L + 63) / 64])
		{
			std::memcpy(index.get(), rhs.index.get(), static_cast<std::size_t>((L + 7) / 8));
		}
	
		T operator[] (std::size_t i) const {
			if (exists(i))
				return data[count_to(i)];
			else
				return T();
		}

		proxy operator[] (std::size_t i) {
			return proxy(*this, i);
		}

		sparse_vector<T, L> operator+(sparse_vector<T, L> const & v) const {
			sparse_vector<T, L> r;
			std::size_t p = 0, q = 0;
			for (std::size_t i = 0; i < (L + 63) / 64; ++i)
				if ((r.index[i] = index[i] | v.index[i]))
					for (std::size_t j = 0; j < (i == (L + 63) / 64 ? L & 0x3F : 64); ++j)
						if (r.index[i] & one << j) {
							if (index[i] & one << j) {
								if (v.index[i] & one << j)
									r.data.push_back(data[p++] + v.data[q++]);
								else
									r.data.push_back(data[p++]);
							} else
								r.data.push_back(v.data[q++]);
							r.last = (i << 6) + j;
						}
			return r;
		}

		T operator*(sparse_vector<T, L> const & v) const {
			T r = T();
			std::size_t p = 0;
			for (std::size_t i = 0; i < (L + 63) / 64; ++i)
				if (index[i] & v.index[i])
					for (std::size_t j = 0; j < (i == (L + 63) / 64 ? L & 0x3F : 64); ++j)
						if (index[i] & v.index[i] & one << j) {
							r += data[p] * v.data[p];
							++p;
						}
			return r;
		}

		void output(std::ostream & os) const {
			std::size_t p = 0;
			for (std::size_t i = 0; i < (L + 63) / 64; ++i)
				if (index[i])
					for (std::size_t j = 0; j < (i == (L + 63) / 64 ? L & 0x3F : 64); ++j)
						if (index[i] & (one << j)) {
							std::cout << (p ? ", " : "") << ((i << 6) + j) << ":" << data[p];
							++p;
						}
		}

        size_type index_of_next_nonzro(size_type i) const {

            if( i > last) return L+1;

            for(size_type j = i>>6; j > 0; --j)
                bit_count(j-1);
        }
	private:

		void set(T const & v, std::size_t i) {
			if (exists(i))
				data[count_to(i)] = v;
			else {
				index[i >> 6] |= (one << (i & 0x3F));
				if (last < i) {
					data.push_back(v);
					last = i;
				} else {
					std::size_t j = count_to(i);
					data.resize(data.size() + 1);
					std::copy(data.begin() + j, data.end() - 1, data.begin() + j + 1);
					data[j] = v;
				}
			}
		}

		bool exists(std::size_t i) {
			return i < last && index[i >> 6] & (one << (i & 0x3F));
		}

		std::size_t count_to(std::size_t i) const {
			std::size_t c = 0;
			for (std::size_t j = 0; j < i >> 6; ++j)
				c += bit_count(j);
			for (std::size_t j = 0; j < i & 0x3F; ++j)
				c += index[i >> 6] & (one << (i & 0x3F));
			return c;
		}

		std::size_t bit_count(std::size_t i) const {
			// see: http://graphics.stanford.edu/~seander/bithacks.html
			boost::uint64_t c = 0;
			c = index[i] - ((index[i] >> 1) & 0x5555555555555555);
			c = ((c >> 2) & 0x3333333333333333) + (c & 0x3333333333333333);
			c = ((c >> 4) + c) & 0x0F0F0F0F0F0F0F0F;
			c = ((c >> 8) + c) & 0x00FF00FF00FF00FF;
			c = ((c >> 16) + c) & 0x0000FFFF0000FFFF;
			c = ((c >> 32) + c) & 0x00000000FFFFFFFF;
			return c;
		}

		std::size_t last;
		std::vector<T> data;
		boost::shared_array<boost::uint64_t> index;
};

template <typename T, std::size_t L>
std::ostream & operator<<(std::ostream & os, sparse_vector<T, L> const & v) {
	v.output(os);
	return os;
}


}

#endif //SPARSE_VECTOR_HPP
