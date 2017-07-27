/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2017 by Alberto Baiardi <alberto.baiardi@sns.it>
 *
 * This software is part of the ALPS libraries, published under the ALPS
 * Library License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Library License along with
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

#ifndef PRINTER_H
#define PRINTER_H

struct printer {
    // Forward declaration of all the functions
    void print_endline_simple(void) ;
    void print_header_table_simple(void) ;
    void print_newline_table_simple(const int& iter, const int& dim, const float& error, const float& en) ;
    void print_newline_table_simple_onlyenergy(const float& error, const float& en) ;
};

void printer::print_endline_simple(void) {
    std::cout << "-----------+-----------+-------------+-------------" << std::endl ;
} ;

void printer::print_header_table_simple(void) {
    print_endline_simple() ;
    std::cout << " Iteration | Sub. Dim. |    Error    |    Energy    " << std::endl ;
    print_endline_simple() ;
} ;

void printer::print_newline_table_simple(const int& iter, const int& dim, const float& error, const float& en){
    char buf[100];
    int n = sprintf(buf, "%5d      | %7d   | %1.4E  | %6.5F ", iter , dim, error, en);
    std::cout << buf << std::endl ;
}

void printer::print_newline_table_simple_onlyenergy(const float& error, const float& en){
    char buf[100];
    int n = sprintf(buf, "           |           | %1.4E  | %6.5F ", error, en);
    std::cout << buf << std::endl;
}

#endif
