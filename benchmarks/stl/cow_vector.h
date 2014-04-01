#ifndef COW_VECTOR_H
#define COW_VECTOR_H

#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>

template<class T, class Allocator = std::allocator<T> >
class copy_on_write_vector
{
protected:
    typedef std::vector<T, Allocator> data_t;
    boost::shared_ptr<data_t> data_;
    
    void make_unique()
    {
        // std::cout << "Making unique" << std::endl;
        if (!data_.unique())
            data_.reset(new data_t(*data_.get()));
    }
    
public:
    // let's forward some interface
    // basically copy&pasted this from the C++ std
    
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef typename data_t::iterator iterator;
    typedef typename data_t::const_iterator const_iterator;
    typedef typename data_t::size_type size_type;
    typedef typename data_t::difference_type difference_type;
    typedef T value_type;
    typedef Allocator allocator_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    
    explicit copy_on_write_vector(const Allocator& a = Allocator())
    : data_(new data_t(a)) { }
    
    explicit copy_on_write_vector(size_type n, const T& value = T(), const Allocator& a = Allocator())
    : data_(new data_t(n, value, a)) { }
    
    template <class InputIterator>
    copy_on_write_vector(InputIterator first, InputIterator last, const Allocator& a = Allocator())
    : data_(new data_t(first, last, a)) { }
    
    copy_on_write_vector(const data_t& x)
    : data_(new data_t(x)) { }
    
    copy_on_write_vector(const copy_on_write_vector& x)
    : data_(x.data_) { }
    
    copy_on_write_vector<T,Allocator>& operator=(copy_on_write_vector<T,Allocator> x)
    {
        swap(x, *this);
    }
    
    copy_on_write_vector<T,Allocator>& operator=(const data_t &x)
    {
        data_.reset(new data_t(x));
    }
    
    // template <class InputIterator> void assign(InputIterator first, InputIterator last);
    // void assign(size_type n, const T& u);
    // allocator_type get_allocator() const;
    
    iterator begin()
    {
        make_unique();
        return data_->begin();
    }
    iterator end()
    {
        make_unique();
        return data_->end();
    }
    
    const_iterator begin() const
    {
        return data_->begin();
    }
    const_iterator end() const
    {
        return data_->end();
    }
    
    // reverse_iterator rbegin();
    // const_reverse_iterator rbegin() const;
    // reverse_iterator rend();
    // const_reverse_iterator rend() const;
    
    size_type size() const
    {
        return data_->size();
    }
    // size_type max_size() const;
    void resize(size_type sz, T c = T())
    {
        make_unique();
        data_->resize(sz, c);
    }
    // size_type capacity() const;
    // bool empty() const;
    // void reserve(size_type n);
    
    // element access:
    reference operator[](size_type n)
    {
        make_unique();
        return (*data_)[n];
    }
    const_reference operator[](size_type n) const
    {
        return (*data_)[n];
    }
    
    // yet to be done...
    // const_reference at(size_type n) const;
    // reference at(size_type n);
    // reference front();
    // const_reference front() const;
    // reference back();
    // const_reference back() const;
    
    // 23.2.4.3 modifiers:
    // void push_back(const T& x);
    // void pop_back();
    // iterator insert(iterator position, const T& x);
    // void insert(iterator position, size_type n, const T& x);
    // template <class InputIterator>
    // void insert(iterator position,
    // InputIterator first, InputIterator last);
    // iterator erase(iterator position);
    // iterator erase(iterator first, iterator last);
    // void swap(vector<T,Allocator>&);
    // void clear();
    
    // specialized algorithms:
    friend void swap(copy_on_write_vector<T,Allocator>& x, copy_on_write_vector<T,Allocator>& y)
    {
        swap(x.data_, y.data_);
    }
};

// template <class T, class Allocator>
// bool operator==(const vector<T,Allocator>& x,
// const vector<T,Allocator>& y);
// template <class T, class Allocator>
// bool operator< (const vector<T,Allocator>& x,
// const vector<T,Allocator>& y);
// template <class T, class Allocator>
// bool operator!=(const vector<T,Allocator>& x,
// const vector<T,Allocator>& y);
// template <class T, class Allocator>
// bool operator> (const vector<T,Allocator>& x,
// const vector<T,Allocator>& y);
// template <class T, class Allocator>
// bool operator>=(const vector<T,Allocator>& x,
// const vector<T,Allocator>& y);
// template <class T, class Allocator>
// bool operator<=(const vector<T,Allocator>& x,
// const vector<T,Allocator>& y);

#endif
