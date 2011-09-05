#ifndef HP2C__SIMPLE_SPARSE_VECTOR_HPP
#define HP2C__SIMPLE_SPARSE_VECTOR_HPP

#include "simple_sparse_vector_iterator.hpp"
#include <map>
#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <iostream>

namespace series_expansion
{

template <typename T>
class simple_sparse_vector
{
    public:
        typedef T value_type;
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
    private:
        typedef typename std::map<size_type,T>::iterator         map_iterator;
        typedef typename std::map<size_type,T>::const_iterator   const_map_iterator;
    public:
        typedef simple_sparse_vector_iterator<T,map_iterator,size_type> iterator;
        typedef simple_sparse_vector_iterator<const T,const_map_iterator,size_type> const_iterator;
	
    class proxy {
		public:
	
			proxy(simple_sparse_vector<T> & o, size_type i)
				: obj(o)
				, index(i)
			{}
	
            inline T& get_ref()
            {
                map_iterator it = obj.elements.find(index);
                if( it != obj.elements.end())
                    return it->second;
                else
                    return obj.elements.insert(std::make_pair(index,obj.init_value)).first->second;
                /*
                    // Astonisingly doing everything in a single STL call isn't really faster:
                    return obj.elements.insert(std::pair<size_type,T&>(index,obj.init_value)).first->second;
                */
            }

			operator T() const {
                map_iterator it = obj.elements.find(index);
                if( it != obj.elements.end())
                    return it->second;
                else
				    return obj.init_value;
			}
	
			T operator=(T const & v) {
				return get_ref() = v;
			}

            T operator += (T const& v) {
                return get_ref() += v;
            }

            T operator -= (T const& v) {
                return get_ref() -= v;
            }
            
            T operator *= (T const& v) {
                return get_ref() *= v;
            }

            T operator /= (T const& v) {
                return get_ref() /= v;
            } 
		private:
			simple_sparse_vector<T> & obj;
			size_type index;
	};
    public:
        simple_sparse_vector(size_type size = 0, T const& t = T())
            : init_value(t)
        {
        }

        iterator begin()
        {
            return iterator(elements.begin());
        }
        
        iterator end()
        {
            return iterator(elements.end());
        }
        
        const_iterator begin() const
        {
            return const_iterator(elements.begin());
        }
        
        const_iterator end() const
        {
            return const_iterator(elements.end());
        }

        const T operator()(size_type i) const
        {
            map_iterator it =elements.find(i);
            if(it != elements.end())
                return it->second;
            else
                return init_value;
        }

        proxy operator()(size_type i)
        {
            return proxy(*this,i);
        }

        size_type size() const
        {
//            std::cerr<<"WARNING Size returns number of non-zero elements, not total elements."<<std::endl;
            return elements.size();
        }
        
        inline const size_type get_index_by_iterator(iterator it)
        { return it->first; }

        void print(std::ostream& o)
        {
            for(map_iterator it(elements.begin()); it != elements.end(); ++it)
                o<<it->first<<": "<<it->second<<std::endl;
        }

        simple_sparse_vector& operator += (simple_sparse_vector const& v)
        {
            for(const_map_iterator it(v.elements.begin()); it != v.elements.end(); ++it)
            {
                map_iterator lookup = elements.find(it->first);
                if( lookup != elements.end() )
                {
                    // A similar entry has been found in our simple_sparse_vector
                    // => add elements
                    lookup->second += it->second;
                    if( lookup->second == value_type(0))
                        elements.erase(lookup);
                }
                else
                {
                    // There's no similar entry in our simple_sparse_vector
                    // => create and add one
                    elements.insert(lookup,*it);
                }
            }
            return *this;
        }

        simple_sparse_vector& operator -= (simple_sparse_vector const& v)
        {
	        *this += T(-1) * v;
            return *this;
        }

        simple_sparse_vector& operator *= (T const& t)
        {
            if(t == T(0))
                elements.clear();
            else
                std::for_each(begin(), end(), boost::lambda::_1 *= t);
            return *this;
        }
        
        friend inline T operator * (simple_sparse_vector<T> const& v1, simple_sparse_vector<T> const& v2)
        {
            // The first state should be smaller than the second one, for best performance
            if(v1.elements.size() > v2.elements.size())
            {
                return v2*v1;
            }

            T result(0);
            // Go through all bras in the smaller state and look if we find a corresponding ket in the larger state
            for(typename simple_sparse_vector<T>::const_map_iterator it(v1.elements.begin()); it != v1.elements.end(); ++it)
            {
                typename simple_sparse_vector<T>::const_map_iterator found = v2.elements.find(it->first);
                if(found != v2.elements.end())
                {	
                    // It really found something
                    // so we have a ket and the corresponding bra
                    // => multiply amplitudes and add it to our result
                    result += found->second * it->second;
                }
            }
            return result;
        }

    private:
        T init_value;
        std::map<size_type,T>   elements;
        friend class simple_sparse_vector_iterator<T,size_type,map_iterator>;
        friend class simple_sparse_vector_iterator<const T,size_type,const_map_iterator>;
};

template <typename T>
inline const simple_sparse_vector<T> operator + (simple_sparse_vector<T> v1, simple_sparse_vector<T> const& v2)
{
    v1 += v2;
    return v1;
}

template <typename T>
inline const simple_sparse_vector<T> operator - (simple_sparse_vector<T> v1, simple_sparse_vector<T> const& v2)
{
    v1 -= v2;
    return v1;
}

template <typename T>
inline const simple_sparse_vector<T> operator * (simple_sparse_vector<T> v, T const& t)
{
    v *= t;
    return v;
}

template <typename T>
inline const simple_sparse_vector<T> operator * (T const& t, simple_sparse_vector<T> v)
{
    v *= t;
    return v;
}

} //namespace series_expansion

#endif // HP2C__SIMPLE_SPARSE_VECTOR_HPP
