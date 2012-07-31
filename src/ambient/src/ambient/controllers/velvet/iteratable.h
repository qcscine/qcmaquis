#ifndef AMBIENT_CONTROLLERS_VELVET_ITERATABLE
#define AMBIENT_CONTROLLERS_VELVET_ITERATABLE

namespace ambient { namespace controllers { namespace velvet {

    using ambient::models::velvet::revision;

    struct revision_sub {
        inline revision* get_parent();
    };

    class c_revision : public revision_sub { public: template<typename T> inline operator T* (); }; // check
    class w_revision : public revision_sub { public: template<typename T> inline operator T* (); }; // weak
    class p_revision : public revision_sub { public: template<typename T> inline operator T* (); }; // purge
    class r_revision : public revision_sub { public: template<typename T> inline operator T* (); }; // reuse

} } }

#endif
