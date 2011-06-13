/*******************************************************************************
* nemesis. an experimental finite element code.                                *
* Copyright (C) 2004-2011 F.E.Karaoulanis [http://www.nemesis-project.org]     *
*                                                                              *
* This program is free software; you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License version 3, as           *
* published by the Free Software Foundation.                                   *
*                                                                              *
* This program is distributed in the hope that it will be useful,              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see < http://www.gnu.org/licenses/>.       *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#ifndef SRC_NUMERIC_NUMERIC_H_
#define SRC_NUMERIC_NUMERIC_H_

// C++ system files
#include <cmath>
#include <cstdio>
#include <limits>
#include "exception/sexception.h"

namespace num {
  const double pi = 4.*std::atan(1.0);
  const double sq2 = std::sqrt(2.);
  const double sq3 = std::sqrt(3.);
  const double sq6 = std::sqrt(6.);
  const double d13 = 1./3.;
  const double d23 = 2./3.;
  const double d43 = 4./3.;
  const double d16 = 1./6.;
  const double eps = std::numeric_limits <double>::epsilon();
  inline bool equal(double d1, double d2) {
    return std::fabs(d2-d1) < eps ? true : false;
  }
  inline bool smaller(double d1, double e) {
    return std::fabs(d1)    < e   ? true : false;
  }
  inline bool smaller(double d1, double d2, double e) {
    return std::fabs(d1-d2) < e   ? true : false;
  }
  inline bool tiny(double d) {
    return std::fabs(d)     < eps ? true : false;
  }
  inline double deg2rad(double d) {
    return d*pi/180.;
  }
  inline double rad2deg(double d) {
    return d*180./pi;
  }

  inline void print_d(double d, int total, int decimal) {
    // (total)-(decimal)-(sign)-(period)
    int integral = total-decimal-1-1;
    // find number of integral digits on d
    int digits = 0;
    int temp = d>0 ? static_cast<int>(d) : -static_cast<int>(d);
    while (temp > 0) {
      digits++;
      temp /= 10;
    }
    // now we are ready to print
    if (integral < 0) {
      for (int i = 0; i < total; i++) printf("*");
    } else if (digits > integral) {
      printf(" ");
      for (int i = 0; i < integral; i++) printf("*");
      printf(".");
      for (int i = 0; i < decimal; i++) printf("*");
    } else {
      printf("%*.*f", total, decimal, d);
    }
  }

  inline void print_i(int n, int total) {
    // (total)-(decimal)-(sign)-(period)
    int integral = total-1;
    // find number of integral digits on d
    int digits = 0;
    int temp = n>0 ? n : -n;
    while (temp > 0) {
      digits++;
      temp /= 10;
    }
    // now we are ready to print
    if (integral < 0) {
      for (int i = 0; i < total; i++) printf("*");
    } else if (digits > integral) {
      printf(" ");
      for (int i = 0; i < integral; i++) printf("*");
    } else {
      printf("%*d", total, n);
    }
  }


  inline double sign(double d) {
    return d < 0  ? -1 : 1;
  }
  
  inline bool is_zero(int i) {
    return i == 0 ? true : false;  
  }

  inline bool is_negative(int i) {
    return i < 0 ? true : false;  
  }
  
  namespace check {
    inline void array_range(int i, int n) {
      if ((i < 0)||(i > n))
        throw SException("[nemesis:%d] %s", 1002, "Out of array bounds.");
    }
    inline void array_range(int i, int j, int n, int m) {
      if ((i < 0)||(i > n)||(j < 0)||(j > m))
        throw SException("[nemesis:%d] %s", 1003, "Out of array bounds.");
    }
    inline void array_size(int n, int size) {
      if (n != size)
        throw SException("[nemesis:%d] %s", 1004, "Sizes do not match.");
    }
    inline void array_size(int n, int m, int size1, int size2) {
      if (n != size1||m != size2)
        throw SException("[nemesis:%d] %s", 1005, "Sizes do not match.");
    }
  }
}
#endif  // SRC_NUMERIC_NUMERIC_H_
