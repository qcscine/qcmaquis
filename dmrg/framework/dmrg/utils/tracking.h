#ifndef TRACKING_H
#define TRACKING_H

#ifdef AMBIENT_TRACKING
    #define ambient_track(variable) ambient_track_as(variable, #variable)
    #define ambient_track_array(variable, i) ambient_track_as( variable[i], (std::string(#variable) + "[" + std::to_string(i) + "]") )

    template<class Matrix, class SymmGroup> class Boundary;
    template<class Matrix, class SymmGroup> class MPSTensor;
    template<class Matrix, class SymmGroup> class block_matrix;

    template<class Matrix, class SymmGroup> 
    void track(block_matrix<Matrix, SymmGroup>& t, const std::string& label){
        for(int i = 0; i < t.n_blocks(); ++i) 
        ambient::track(t[i], label);
        t.label = label;
    }

    template<class Matrix, class SymmGroup> 
    void track(MPSTensor<Matrix, SymmGroup>& t, const std::string& label){
        track(t.data(), label);
    }

    template<typename T>
    void ambient_track_as(T& object, const std::string& label){
        track(object, label);
    }
#else
    #define ambient_track(variable)
    #define ambient_track_array(variable, i)
    #define ambient_track_as(variable, label)
#endif

#endif
