/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2015-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef QC_OPERATOR_BUNDLE_SU2_H
#define QC_OPERATOR_BUNDLE_SU2_H

    template <class Matrix, class SymmGroup>
    class OperatorBundle_
    {
        // This class manages the different fill / spin_input / spin_output variants of creators and destructors
        typedef typename TagHandler<Matrix, SymmGroup>::tag_type tag_type;
        typedef typename TagHandler<Matrix, SymmGroup>::op_t op_t;
        typedef typename SpinDescriptor<symm_traits::SU2Tag>:: spin_t spin_t;
        
    public:
        OperatorBundle_(op_t bop_, op_t bfop_, boost::shared_ptr<TagHandler<Matrix, SymmGroup> > th)
            : base_op(bop_)
            , base_fill_op(bfop_)
            , tag_handler(th)
            , spin(base_op.spin.get())
        {
                 variants[std::make_pair(base_op.spin.input(), base_op.spin.output())] = base_op;
            fill_variants[std::make_pair(base_fill_op.spin.input(), base_fill_op.spin.output())] = base_fill_op;
        }

        tag_type operator()(bool fill, spin_t spin_in, spin_t spin_out)
        {
            try {
#if defined(__xlC__) || defined(__FCC_VERSION)
                if (fill) {
                    if (fill_variants.count(std::make_pair(spin_in, spin_out)) == 0)
                        throw std::out_of_range("");

                    return fill_variants[std::make_pair(spin_in, spin_out)];
                }
                else  {
                    if (variants.count(std::make_pair(spin_in, spin_out)) == 0)
                        throw std::out_of_range("");

                    return variants[std::make_pair(spin_in, spin_out)];
                }
#else
                if (fill)
                    return fill_variants.at(std::make_pair(spin_in, spin_out));
                else 
                    return variants.at(std::make_pair(spin_in, spin_out));
#endif
            }
            catch(const std::out_of_range& e) {

                // triangle condition
                assert (spin_out >= std::abs(spin_in - spin) && spin_out <= std::abs(spin_in + spin));

                SpinDescriptor<symm_traits::SU2Tag> new_descriptor(spin, spin_in, spin_out);

                op_t new_variant = (fill) ? base_fill_op : base_op;
                new_variant.spin = new_descriptor;

                tag_type new_tag = tag_handler->register_op(new_variant);

                if (fill)
                    fill_variants[std::make_pair(spin_in, spin_out)] = new_tag;                          
                else
                    variants[std::make_pair(spin_in, spin_out)] = new_tag;                          
            }
        }

    private:
        boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
        std::map<std::pair<spin_t, spin_t>, tag_type> variants;
        std::map<std::pair<spin_t, spin_t>, tag_type> fill_variants;

        spin_t spin;
        op_t base_op;
        op_t base_fill_op;
    };

#endif
