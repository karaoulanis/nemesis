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

#include "reorderer/forward_cuthill_mckee.h"
#include < boost/graph/cuthill_mckee_ordering.hpp>

ForwardCuthillMckee::ForwardCuthillMckee() {
  myTag = TAG_REORDERER_FORWARD_CUTHILL_MCKEE;
}
ForwardCuthillMckee::~ForwardCuthillMckee() {
}
int ForwardCuthillMckee::get_perm(std::vector < int>& perm) {
  // Create the Graph and additional vectors
  UndirectedGraph G(pA->get_model()->get_num_eqns());
  pA->get_model()->get_undirected_graph(G);
  property_map<UndirectedGraph, vertex_index_t>::type
    index_map = get(vertex_index, G);
  std::vector<int> inv_perm(num_vertices(G));
  perm.resize(num_vertices(G));

  // Find old bandwidth and profile
  int oldBandwidth = bandwidth(G);
  int oldProfile = profile(G);

  // Call Cuthill McKee algorithm
  cuthill_mckee_ordering(G, inv_perm.begin(), get(vertex_color, G),
                         make_degree_map(G));
  for (unsigned k = 0; k != inv_perm.size(); k++)
    perm[index_map[inv_perm[k]]]=k;

  // Find new bandwidth and profile
  int newBandwidth = bandwidth(G,
                      make_iterator_property_map(&perm[0], index_map, perm[0]));
  int newProfile = profile(G, make_iterator_property_map(&perm[0], index_map));

  // Print optimized sizes
  cout  << "reo: Optimized (original) bandwidth : "
        << newBandwidth  << " (" << oldBandwidth << ")" << endl;
  cout  << "reo: Optimized (original) profile   : "
        <<newProfile    << " (" << oldProfile  << ")"  << endl;

  // Check if optimization is needed
  if (newBandwidth > oldBandwidth) return -1;
  return 1;
}
