#ifndef TENSOR_SYMMETRY_H
#define TENSOR_SYMMETRY_H

namespace tensor {
    // Class representing symmetries
    // This implements the charge calculus of arXiv:1010.3595
    class ExampleSymmetry
    {
    public:
        // Charge
        typedef /* implementation defined */ charge;
        // this could also be implementation-defined
        typedef std::vector<charge> charge_v;
        // mapping from charge to an integer index
        // why in the symmetry class? this is often trivial for symmetries like Z_q
        typedef /* implementation defined */ charge_map;
        // the only method required in charge_map is:
        // std::size_t operator[](charge) const;
        // more or less a convenience function, could be replaced by
        // a constructor in the charge_map
        static charge_map get_map(charge_v);
    
        // this is called \mathbb{I} in the paper
        static const charge SingletCharge;
    
        static charge fuse(charge, charge);
        static template<int L> charge fuse(boost::array<charge, L>);
    };
    ExampleSymmetry::charge operator-(ExampleSymmetry::charge);
} // namespace tensor

#endif
