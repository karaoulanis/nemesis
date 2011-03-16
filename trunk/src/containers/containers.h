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

#ifndef SRC_CONTAINERS_CONTAINERS_H_
#define SRC_CONTAINERS_CONTAINERS_H_

#include <algorithm>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
typedef std::vector<int> IDContainer;

namespace Containers {
  template<class TC> void map_delete(TC& c) {
    typename TC::const_iterator p;
    for (p = c.begin(); p != c.end(); p++) delete p->second;
    c.clear();
  }
  template<class TC> void map_print(const TC& c) {
    typename TC::const_iterator p;
    cout << "Container : ";
    for (p = c.begin();p != c.end();p++) cout << p->first << ' ';
    cout << endl;
  }
  template<class TC> void vector_delete(TC& c) {
    int n = c.size();
    for (int i = 0;i < n;i++) delete c[i];
    c.clear();
  }
  template<class TC> void vector_print(const TC& c) {
    int n = c.size();
    for (int i = 0;i < n;i++) cout << c[i] <<' ';
    cout << endl;
  }
  template<class TE> inline int index_find(const std::vector<TE>& c, TE e) {
    typename std::vector<TE>::const_iterator i = find(c.begin(), c.end(), e);
    if (i == c.end())  return -1;
    else      return static_cast < int>(std::distance(c.begin(), i));
  }
  template<class TE> inline bool all_positive(const std::vector<TE>& c) {
    int n = c.size();
    for (int i = 0;i < n;i++) if (c[i] < 0) return false;
    return true;
  }
}
#endif  // SRC_CONTAINERS_CONTAINERS_H_
