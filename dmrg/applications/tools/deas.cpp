/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Laboratory of Physical Chemistry, ETH Zurich
 *               2014-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2014-2014 by Erik Hedegaard <erik.hedegaard@phys.chem.ethz.ch>
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
#include "dmrg/models/measurements/chementropy.h"


int main(int argc, char ** argv)
{
    std::cout.precision(12);

    try {
        if (argc != 2) {
            std::cout << "Usage: " << argv[0] << " <result.h5>" << std::endl;
            return 1;
        }

        std::string rfile(argv[1]);

        EntanglementData<matrix> em(rfile);
        std::cout << "single orbital entropy:" << std::endl << em.s1() << std::endl;
        std::cout << "two orbital entropy:" << std::endl << em.s2() << std::endl;
        std::cout << "mutual information:" << std::endl << em.I() << std::endl;

    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
