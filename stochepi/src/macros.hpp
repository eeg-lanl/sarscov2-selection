/* mlev-hiv-model - A program for running multi-level HIV-1 simulations
 * Copyright (C) 2017 Christiaan H. van Dorp <chvandorp@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MACROS_HPP_
#define MACROS_HPP_

#include <string>
#include <iostream>

#define RIGHT_HERE (std::string(" (in function ") + __PRETTY_FUNCTION__ + ") ")

#define WARN_UNTESTED_FUN \
    static bool fun_called(false); \
    if ( !fun_called ) { \
        std::cerr << "WARNING: called untested function " << __PRETTY_FUNCTION__ << std::endl; \
        fun_called = true; \
    }
// end define WARN_UNTESTED_FUN

#define WARN_NOT_IMPLEMENTED \
    static bool fun_called(false); \
    if ( !fun_called ) { \
        std::cerr << "WARNING: not implemented " << __PRETTY_FUNCTION__ << std::endl; \
        fun_called = true; \
    }
// end define WARN_NOT_IMPLEMENTED

#define WARN_DEPRECATED_FUN \
    static bool fun_called(false); \
    if ( !fun_called ) { \
        std::cerr << "WARNING: called deprecated function " << __PRETTY_FUNCTION__ << std::endl; \
        fun_called = true; \
    }
// end define WARN_DEPRECATED_FUN

#endif
