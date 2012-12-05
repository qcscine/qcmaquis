#ifndef AMBIENT_CONTROLLERS_VELVET_ITERATABLE
#define AMBIENT_CONTROLLERS_VELVET_ITERATABLE

namespace ambient { namespace controllers { namespace velvet {

    using ambient::models::velvet::revision;

    class c_revision { public: template<typename T> operator T* (); }; // check
    class p_revision { public: template<typename T> operator T* (); }; // purge
    class w_revision { public: template<typename T> operator T* (); }; // weak
    class s_revision { public: template<typename T> operator T* (); }; // refresh

} } }

#endif
