using namespace tensor;

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
        mult(c(_(down), _(up)),
            bra_tensor(_(alpha), (beta, sigma)),
            temp3);
        temp3.rename(beta, up);
        mult(temp3(_(beta), (sigma, low)),
            conj(ket_tensor((sigma, alpha), _(beta))),
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
    
    if (decomposition_method == qr) {
        tensor<T, 3, SymmGroup> q;
        tensor<T, 2, SymmGroup> r;
        qr(me((alpha, sigma), _(beta)), q, r);
        me = q;
        return r;
    } else {
        tensor<T, 3, SymmGroup> newme;
        tensor<T, 2, SymmGroup> rubbish;
        svd(me((alpha, sigma), _(beta)), newme, rubbish, NULL, j, truncation, bond_dim, SV_TO_RIGHT);
        me = newme;
        return rubbish;
    }
}
