#ifndef TENSOR_SYMMETRY_H
#define TENSOR_SYMMETRY_H

namespace tensor {
    // Class representing symmetries
    // This implements the charge calculus
    class ExampleSymmetry
    {
    public:
        // Charge
        typedef /* implementation defined */ charge;
        // this could also be implementation-defined
        typedef std::vector<charge> charge_v;
        // mapping from integer indices 0, ..., n to charge
        // why in the symmetry class? this is often trivial for symmetries like Z_q
        typedef /* implementation defined */ charge_map;
        static charge_map get_map(charge_v);
    
        static const charge SingletCharge;
    
        static charge fuse(charge, charge);
        static template<int L> charge fuse(boost::array<charge, L>);
    };
    ExampleSymmetry::charge operator-(ExampleSymmetry::charge);
} // namespace tensor

#endif
