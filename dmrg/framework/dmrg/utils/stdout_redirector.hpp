/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
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