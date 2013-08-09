#ifndef TRACKING_H
#define TRACKING_H

#ifdef AMBIENT_TRACKING
    #define __ambient_track(variable) __ambient_track_as(variable, #variable)
    #define __ambient_track_array(variable, i) __ambient_track_as( variable[i], (std::string(#variable) + "[" + std::to_string(i) + "]") )

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
    void __ambient_track_as(T& object, const std::string& label){
        track(object, label);
    }
#else
    #define __ambient_track(variable)
    #define __ambient_track_array(variable, i)
    #define __ambient_track_as(variable, label)
#endif

#endif
