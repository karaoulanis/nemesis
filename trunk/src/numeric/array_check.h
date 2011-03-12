/*******************************************************************************
* nemesis. an experimental finite element code.                                *
* Copyright (C) 2004-2010 F.E.Karaoulanis [http://www.nemesis-project.org]     *
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

#ifndef SRC_NUMERIC_ARRAY_CHECK_H_
#define SRC_NUMERIC_ARRAY_CHECK_H_

#include "exception/sexception.h"

inline void array_range_check(int i, int n) {
  if ((i < 0)||(i > n))
    throw SException("[nemesis:%d] %s", 1002, "Out of array bounds.");
}
inline void array_range_check(int i, int j, int n, int m) {
  if ((i < 0)||(i > n)||(j < 0)||(j > m))
    throw SException("[nemesis:%d] %s", 1003, "Out of array bounds.");
}
inline void array_size_check(int n, int size) {
  if (n != size)
    throw SException("[nemesis:%d] %s", 1004, "Sizes do not match.");
}
inline void array_size_check(int n, int m, int size1, int size2) {
  if (n != size1||m != size2)
    throw SException("[nemesis:%d] %s", 1005, "Sizes do not match.");
}
#endif  // SRC_NUMERIC_ARRAY_CHECK_H_
