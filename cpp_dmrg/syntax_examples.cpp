using namespace tensor;

// there will most likely be a global typedef to decide
// whether we work with real or complex matrices
typedef double T;

template<typename T, class SymmGroup>
tensor<T, 2, SymmGroup> overlap_left_step(
    tensor<T, 3, SymmGroup> & bra_tensor,
    tensor<T, 3, SymmGroup> & ket_tensor,
    tensor<T, 2, SymmGroup> & c,
    tensor<T, 2, SymmGroup> * local_operator = NULL)
{
    // error handling skipped
    
    if (local_operator == NULL) {
        tensor<T, 3, SymmGroup> temp3;
        tensor<T, 2, SymmGroup> temp2;
        // contract c with bra_tensor, where the alpha index of bra_tensor is
        // contracted with the up index of c
        mult(c(down, up),
            bra_tensor(alpha, (beta, sigma)),
            temp3);
        temp3.rename(beta, up);
        // make beta (from bra_tensor) the new up
        mult(temp3(beta, (sigma, low)),
            conj(ket_tensor((sigma, alpha), beta)),
            temp2);
        temp2.rename(beta, low);
        return temp2;
    } else { /* ... */ }
}

enum DM { qr, svd };

template<typename T, class SymmGroup>
tensor<T, 2, SymmGroup> MPSTensor::normalize_left(DM decomposition_method, double truncation = 0, DIndex<SymmGroup> bond_dim)
{
    // error handling skipped
    
    tensor<T, 3, SymmGroup> q;
    tensor<T, 2, SymmGroup> r;
    
    if (decomposition_method == qr)
        qr(me((alpha, sigma), beta), q, r);
        // or, see below
        // qr(left_join(me), q, r);
    else
        svd(me((alpha, sigma), beta), q, r, NULL, j, truncation, bond_dim, SV_TO_RIGHT);
    me = q;
    return r;
}

// we could define convenience functions, like:
template<typename T, class SymmGroup>
tensor_matrix_wrapper<T, 2, 1, SymmGroup, NoOp> left_join(tensor<T, 3, SymmGroup> &t)
{
    return t((alpha, sigma), beta);
}
