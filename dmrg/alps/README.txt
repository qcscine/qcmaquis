The ALPS distribution in this subdirectory has been modified to provide a minimal ALPS library to work flawlessly with the QCMaquis distribution.
Only utility/ expression/ parameter/ parser/ and hdf5/ subdirectories are used, the remaining source code has been removed. In addition, the src/boost subdirectory
contains the boost::numeric bindings that compile with ALPS 2.3.0.

The full version of this library is available from
http://alps.comp-phys.org/static/software/releases/alps-2.3.0-src.tar.gz
or
http://alps.comp-phys.org/static/software/releases/alps-2.3.0-src-with-boost.tar.gz

Below you can find the original ALPS readme file.
-----------------------------------------------------------------------------------

The ALPS project (Algorithms and Libraries for Physics Simulations) aims at providing generic parallel algorithms for classical and quantum lattice models and provides utility classes and algorithm for many other problems. It strives to increase software reuse in the physics community.

The ALPS Libraries are published under the ALPS Library License; you can use, redistribute it and/or modify it under the terms of the license, either version 1 or (at your option) any later version.

You should have received a copy of the ALPS Library License along with the ALPS Libraries; see the file LICENSE.txt. If not, the license is also available from http://alps.comp-phys.org/.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,  FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT  SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE  FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER  DEALINGS IN THE SOFTWARE.

Any publication for which one of the following libraries are used has to
acknowledge the use of the ALPS libraries, and the papers listed below:

When alps/model.h or any header in alps/model was used:
reference the web page http://alps.comp-phys.org/ and cite the publication:
A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007)
B. Bauer et al., J. Stat. Mech. (2011) P05001

When alps/lattice.h or any header in alps/lattice was used:
reference the web page http://alps.comp-phys.org/ and cite the publication:
A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007)
B. Bauer et al., J. Stat. Mech. (2011) P05001

When alps/alea.h or any header in alps/alea was used:
reference the web page http://alps.comp-phys.org/ and cite the publications:
A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007)
B. Bauer et al., J. Stat. Mech. (2011) P05001

When alps/scheduler.h or any header in alps/scheduler was used:
reference the web page http://alps.comp-phys.org/ and cite the publications:
A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007)
B. Bauer et al., J. Stat. Mech. (2011) P05001
M. Troyer et al., Lecture Notes in Computer Science, Vol. 1505, p. 191 (1998).

The use of any other library, in particular those with headers in the
subdirectories alps/parser, alps/osiris, alps/random do not carry any
citation requirement but acknowledgment of the ALPS project is encouraged.

* When the SSE quantum Monte Carlo program sse or sse_mpi was used:
  - reference the ALPS web page http://alps.comp-phys.org/
  - and cite the publications:
    A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).
    B. Bauer et al., J. Stat. Mech. (2011) P05001

* When the loop quantum Monte Carlo program loop or loop_mpi was used:
  - reference the ALPS and ALPS/looper web pages
      http://alps.comp-phys.org/
      http://wistaria.comp-phys.org/alps-looper/
  - and cite the publications:
      S. Todo and K. Kato, Phys. Rev. Lett. 87 047203 (2001).
      A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).
      B. Bauer et al., J. Stat. Mech. (2011) P05001

* When the classical Monte Carlo program spinmc or spinmc_mpi was
  used:
  - reference the ALPS web page http://alps.comp-phys.org/
  - and cite the publications:
      A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).
      B. Bauer et al., J. Stat. Mech. (2011) P05001

* When the diagonalization programs fulldiag, sparsediag or fulldiag_mpi was
  used:
  - reference the ALPS web page http://alps.comp-phys.org/
  - and cite the publications:
    A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).
    B. Bauer et al., J. Stat. Mech. (2011) P05001

* When the worm quantum Monte Carlo program is used:
  - reference the ALPS web page http://alps.comp-phys.org/
  - and cite the publications:
    A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).
    B. Bauer et al., J. Stat. Mech. (2011) P05001

* When the quantum Wang-Landau program is used:
  - reference the ALPS web page http://alps.comp-phys.org/
  - and cite the publications:
    M. Troyer, S. Wessel, and F. Alet, Phys. Rev. Lett. 90, 120201 (2003).
    A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).
    B. Bauer et al., J. Stat. Mech. (2011) P05001

* When the Continuous-Time quantum Monte Carlo impurity solver programs or the DMFT framework are used:
 - cite the publication:
   A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).
   B. Bauer et al., J. Stat. Mech. (2011) P05001
 - cite the ALPS DMFT publication:
   E. Gull, P. Werner, S. Fuchs, B. Surer, T. Pruschke, and M. Troyer,
   Computer Physics Communications 182, 1078 (2011).

* When the Matrix Product States applications mps_optim, mps_evolve, mps_meas or mps_overlap are used:
 -cite the publcations:
  B. Bauer et al., J. Stat. Mech. (2011) P05001.
  M. Dolfi, B. Bauer, S. Keller, et al., Computer Physics Communications 185, 3430 (2014).

  Copyright ALPS collaboration 2002 - 2010
  Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE_1_0.txt or copy at
          http://www.boost.org/LICENSE_1_0.txt)

