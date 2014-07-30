/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include "libpscan/scheduler.hpp"

#include <alps/utility/copyright.hpp>
#include <iostream>

#include "dmrg/version.h"

int main(int argc, char ** argv)
{
    try {
        std::cout << "ALPS/MPS version " DMRG_VERSION_STRING " (2013-2014)\n"
                  << "  Density Matrix Renormalization Group algorithm\n"
                  << "  available from http://alps.comp-phys.org/\n"
                  << "  copyright (c) 2013 Institute for Theoretical Physics, ETH Zurich\n"
                  << "  copyright (c) 2010-2011 by Bela Bauer\n"
                  << "  copyright (c) 2011-2013 by Michele Dolfi\n"
                  << "  for details see the publication: \n"
                  << "  M. Dolfi et al, in preparation\n"
                  << std::endl;
        alps::print_copyright(std::cout);
        
        Options opt(argc,argv);
        if (opt.valid) {
            Scheduler pscan(opt);
            pscan.run();
        }
    } catch (std::exception & e) {
        std::cerr << "Exception thrown:" << std::endl;
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}

