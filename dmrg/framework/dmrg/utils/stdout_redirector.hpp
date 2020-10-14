/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2020         Leon Freitag <lefreita@ethz.ch>
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
#ifndef STDOUT_REDIRECTOR_HPP
#define STDOUT_REDIRECTOR_HPP

#include <iostream>
#include <fstream>
#include <streambuf>

namespace maquis
{
    // Class to help redirect std::cout to a file
    class StdoutRedirector
    {
        public:
            StdoutRedirector() : coutbuf_(std::cout.rdbuf()) {}

            StdoutRedirector(const std::string & filename) : coutbuf_(std::cout.rdbuf()), out_(filename)
            {
                if (!out_)
                    throw std::runtime_error("stdout redirection to " + filename + " failed.");
                // redirect output to filename
                std::cout.rdbuf(out_.rdbuf());
            }

            ~StdoutRedirector()
            {
                // restore the standard cout stream on destruction
                std::cout.rdbuf(coutbuf_);
            }

            // restore the standard cout stream without destroying the object
            void restore()
            {
                std::cout.rdbuf(coutbuf_);
                out_.close();
            }

            // change output filename (close the old file and open a new one)
            void set_filename(const std::string & filename)
            {
                out_.close();
                out_.clear();
                out_.open(filename);
                std::cout.rdbuf(out_.rdbuf());
            }
        private:
            std::streambuf* coutbuf_; // cout buffer to restore cout
            std::ofstream out_;
    };
}
#endif