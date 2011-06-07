// nested inside ambient.hpp in ambient namespace
using namespace blas;

#include "ambient/core/operation/operation.pp.hpp"
#define MAX_NUM_CHAR_LEN 10
#define scope_select(...) scope_select(std::string(std::string() + __VA_ARGS__).c_str());

void copy_l(p_dense_matrix<double>& ac, pinned const p_dense_matrix<double>& a);
void copy_c(p_dense_matrix<double>& ac, pinned const p_dense_matrix<double>& a);

class void_pt: public p_profile 
{ 
public: 
    template <typename T> void_pt(const T* ptr) : p_profile()
    { breakdown_model(this, ptr); }
    ~void_pt(){ };  // use with caution
private: 
    friend class T; // the container can delete its profile
    template<class T> friend inline void boost::checked_delete(T * x); // so as boost >_<
};

template<typename T>
void null_deleter(T* object){ }

template<typename T>
void poke_deleter(T* object){
//    printf("Poking... %d - %d ---- proxy: %d\n", object->use_count, ambient::rank(), object->breakdown()->id);
    if(--object->use_count) return;
    delete object;
}

template<class T, policy P>
class livelong
{
public:
   ~livelong(){
        if(this->is_loose_copied()){                  // I have a copy
            if(this->is_loose_copy()){                // I am a copy
                self->duplicant->self->original = self->original;
                self->original->self->duplicant = self->duplicant;
            }else{ 
                if(!this->is_loose()){
                    //this->bind(); // works
                    resize_bind_model((T*)this, self->duplicant->num_rows(), 
                                                self->duplicant->num_cols());
                    self->meta.loose_copied = self->duplicant->self->meta.loose_copied;
                    std::swap(this->handle, self->duplicant->handle);
                    self->duplicant->self = this->self;
                }else{ 
                    self->duplicant->self->meta.loose_copy = false;
                    self->duplicant->self->meta.loose      = true;
                }
            }              // I am original
        }else if(this->is_loose_copy()){              // I am a copy
            self->original->self->meta.loose_copied = false;
        }

        //if(this->p != ANY){
        //    this->breakdown()->deallocate(); // deallocate profile
        //}
    }
    
    void bind() const {
        assert(self != NULL);
        self->modifier = NULL;
        if(this->is_loose_copy()){
            self->original->bind();
        }else if(this->is_loose_copied()){
            self->meta.loose_copied = false;
            self->duplicant->self->meta.loose_copy = false; // replica
            bool loose_copied = self->duplicant->self->meta.loose_copied;
            self->duplicant->self->meta.loose_copied = false; // preventing false stacking
            resize_bind_model(self->duplicant);
            ambient::push(ambient::copy_l, ambient::copy_c, 
                         *self->duplicant, *(const T*)this);
            self->duplicant->self->meta.loose_copied = loose_copied;
        }else if(this->is_loose()){
            self->meta.loose = false;
            bind_model((T*)this);
        }
    }

    T* snapshot(T* holder) const {
        if(this->is_loose_copied()){ // somebody loose copied me
            return self->duplicant->snapshot(holder);
        }
        self->duplicant = holder; 
        self->meta.loose_copied = true;
        T* snapshot = (T*)(new typename T::replica(*self));
        snapshot->original = (T*)this;
        return snapshot;
    }

    livelong(const livelong& o){ // copy constructor
        this->init();
        if(P != REPLICA){
            if(P == ANY) handle.reset( self = o.snapshot((T*)this) );
            else if(P == MANUAL) handle.reset( self = o.snapshot((T*)this), null_deleter<T> );
            else if(P == WEAK) self = o.snapshot((T*)this);
            self->meta.loose = false;         // avoiding binding
            copy_bind_model((T*)this);
            self->meta.loose_copy = true;     // avoiding binding
        }else if(P == REPLICA){
            self = (T*)this;
        }
    }

    void init(){
        this->p                 = P;
        this->use_count         = 0;
        this->meta.loose        = true;
        this->meta.loose_copy   = false;
        this->meta.loose_copied = false;
        this->self              = NULL;
        this->modifier          = NULL;
    }

    livelong(){
        this->init();
        if(P == ANY) handle.reset( self = (T*)(new typename T::replica()) );
        else if(P == MANUAL) handle.reset( self = (T*)this, null_deleter<T> );
        else if(P == REPLICA || P == WEAK) self = (T*)this;
        thyself = self;
    }
    livelong(void_pt* p):profile(p){
        this->init();
        if(P == ANY) handle.reset( self = (T*)(new typename T::replica(p)) );
        else if(P == MANUAL) handle.reset( self = (T*)this, null_deleter<T> );
        else if(P == REPLICA || P == WEAK) self = (T*)this;
        thyself = self;
    }
    /*template<typename A1>
    livelong(A1 a1){
        printf("Args const of livelong!\n");
        this->init();
        if(P == ANY) handle.reset( self = (T*)(new typename T::replica(a1)) );
        else if(P == MANUAL) handle.reset( self = (T*)this, null_deleter<T> );
        else if(P == REPLICA || P == WEAK) self = (T*)this;
        thyself = self;
    }*/
    template<typename A1, typename A2>
    livelong(A1 a1, A2 a2){
        this->init();
        if(P == ANY) handle.reset( self = (T*)(new typename T::replica(a1, a2)) );
        else if(P == MANUAL) handle.reset( self = (T*)this, null_deleter<T> );
        else if(P == REPLICA || P == WEAK) self = (T*)this;
        thyself = self;
    }
    template<typename A1, typename A2, typename A3>
    livelong(A1 a1, A2 a2, A3 a3){
        this->init();
        if(P == ANY) handle.reset( self = (T*)(new typename T::replica(a1, a2, a3)) );
        else if(P == MANUAL) handle.reset( self = (T*)this, null_deleter<T> );
        else if(P == REPLICA || P == WEAK) self = (T*)this;
        thyself = self;
    }
    boost::shared_ptr<T> get_handle() const {
        this->bind();
        if(this->p == ANY) return this->handle;
        else if(this->p == MANUAL) return this->handle;
        else if(this->p == WEAK){ ((livelong*)this)->use_count++; return boost::shared_ptr<T>(this->self, &poke_deleter<T>); } // for one-touch only
        else if(this->p == REPLICA) return boost::shared_ptr<T>(this->self, null_deleter<T>); // should never occur except copyed objects (not deleting in init :)) <- can lead to errors
        return boost::shared_ptr<T>();
    }
    void_pt*& breakdown() const {
        return self->profile;
    }
    void set_breakdown() {
        if(this->p == ANY) return;
        self->profile = new void_pt(self);
    }

    void(*init_fp)(T&);
    void(*get_init_fp() const)(T&){
        return (void(*)(T&))self->init_fp;
    }

    template<typename O>
    void set_init(void(*fp)(O&)){ // incomplete typename is not allowed
        this->latch_meta();
        self->init_fp = (void(*)(T&))fp;
        self->profile->set_init(new core::operation(fp, (T*)this)); // T casting in order to avoid copy-construction!
        this->unlatch_meta();
    }
    bool is_loose()        const { return  self->meta.loose;            } 
    bool is_loose_copied() const { return  self->meta.loose_copied;     }
    bool is_loose_copy()   const { return  self->meta.loose_copy;       }
    bool is_abstract()     const { return (this->is_loose()            || 
                                           this->is_loose_copy());      }
    void latch_meta(){ 
        self->meta_latch        = self->meta;
        self->meta.loose        = false;
        self->meta.loose_copy   = false;
        self->meta.loose_copied = false;
    }
    void unlatch_meta(){
        self->meta = self->meta_latch;
    }
    void decouple_copy(){
        if(this->is_loose_copy()){
            self->meta.loose_copy = false;
            self->original->self->meta.loose_copied = false;
            self->meta.loose = true;
        }
    }
    struct loose_state{
        bool loose;
        bool loose_copy;
        bool loose_copied;
    } meta, meta_latch;

public:
    mutable std::queue< std::vector< std::pair<std::pair<size_t,size_t>,void*> > > modifiers;
    mutable std::vector< std::pair<std::pair<size_t,size_t>,void*> >* modifier; // currently used modifier

    T* self;
    T* thyself;
    size_t use_count;
    void_pt* profile;
    policy p; // the same as P (avoiding casting collisions)
    boost::shared_ptr<T> handle;
    mutable T* duplicant;
private:
    T* original;
};

template<typename T>
boost::shared_ptr<T> get_handle(T& a){
    return boost::shared_ptr<T>(new T(a));
}

template<typename T>
boost::shared_ptr<T> get_handle(const T& a){
    return boost::shared_ptr<T>(new T(a));
}

template<typename T, policy P>
boost::shared_ptr< p_dense_matrix<T> > get_handle(p_dense_matrix<T,P>& a){
    return boost::static_pointer_cast< p_dense_matrix<T> >(boost::static_pointer_cast<void>(a.get_handle()));
}

template<typename T, policy P>
boost::shared_ptr< p_dense_matrix<T> > get_handle(const p_dense_matrix<T,P>& a){
    return boost::static_pointer_cast< p_dense_matrix<T> >(boost::static_pointer_cast<void>(a.get_handle()));
}

// user types //
template <typename T> 
void_pt& breakdown(p_dense_matrix<T>& obj){
    void_pt** profile_ptr = const_cast<void_pt**>(&obj.breakdown());
    *profile_ptr = (void_pt*)(*profile_ptr)->dereference();
    (*profile_ptr)->inconstant();
    return **profile_ptr;
}

template <typename T> 
void_pt& breakdown(const p_dense_matrix<T>& obj){
    void_pt** profile_ptr = const_cast<void_pt**>(&obj.breakdown());
    *profile_ptr = (void_pt*)(*profile_ptr)->dereference();
    (*profile_ptr)->constant();
    return **profile_ptr;
}
// user types //

// aliases for computational kernels
template <typename T>
void_pt& current(T& obj){
    return breakdown(obj);
}

template <typename T>
void plus_reduce(memblock* grp, void* update){
    assert(false); // only partially specialized
}

template <char R, typename T>
void_pt& reduced(T& obj){
    if(breakdown(obj).associated_proxy == NULL){
        if(R == '+') breakdown(obj).associate_proxy(new void_pt((T*)NULL), plus_reduce<T>);
        breakdown_proxy_model((void_pt*)breakdown(obj).associated_proxy, &(breakdown(obj)), &obj);
    }
    return *((void_pt*)breakdown(obj).associated_proxy);
}

std::string & operator+(std::string & lhs, double rhs){
    char* rhs_str = (char*)malloc(sizeof(char)*MAX_NUM_CHAR_LEN);
    if(rhs - (int)rhs) sprintf(rhs_str,"%.1f",rhs);
    else sprintf(rhs_str,"%d",(int)rhs);
    lhs += rhs_str;
    free(rhs_str);
    return lhs;
}
std::string & operator+(std::string & lhs, int rhs){
    char* rhs_str = (char*)malloc(sizeof(char)*MAX_NUM_CHAR_LEN);
    sprintf(rhs_str,"%d", rhs);
    lhs += rhs_str;
    free(rhs_str);
    return lhs;
}
std::string & operator+(const std::string & lhs, double rhs){
    return const_cast<std::string&>(lhs)+rhs;
}
std::string & operator+(const std::string & lhs, int rhs){
    return const_cast<std::string&>(lhs)+rhs;
}
std::string & operator+(std::string & lhs, const char* rhs){
    return lhs += rhs;
}
std::string & operator+(const std::string & lhs, const char* rhs){
    return const_cast<std::string&>(lhs) += rhs;
}
std::string & operator+(std::string & lhs, std::pair<unsigned int*,size_t> rhs){
    lhs + (int)(*rhs.first);
    lhs += "-";
    lhs + (int)(rhs.second);
    return lhs;
}
std::string & operator+(const std::string & lhs, std::pair<unsigned int*,size_t> rhs){
    const_cast<std::string&>(lhs) + (int)(*rhs.first);
    const_cast<std::string&>(lhs) += "-";
    return const_cast<std::string&>(lhs) + (int)rhs.second;
}
