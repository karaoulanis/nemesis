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

#include "reorderer/reverse_sloan.h"
#include < boost/graph/sloan_ordering.hpp>

ReverseSloan::ReverseSloan(double w1, double w2) {
  myTag = TAG_REORDERER_REVERSE_SLOAN;
  weight1 = w1;
  weight2 = w2;
}
ReverseSloan::~ReverseSloan() {
}
int ReverseSloan::get_perm(std::vector<int>& perm) {
  // Create the Graph and additional vectors
  UndirectedGraph G(pA->get_model()->get_num_eqns());
  pA->get_model()->get_undirected_graph(G);
  property_map<UndirectedGraph, vertex_index_t>::type
    index_map = get(vertex_index, G);
  std::vector<int> inv_perm(num_vertices(G));
  perm.resize(num_vertices(G));

  // Find old bandwidth, profile and wavefronts
  int oldBandwidth = bandwidth(G);
  int oldProfile = profile(G);
  /// @todo check for warnings in wavefront()
  // double oldMaxWavefront = max_wavefront(G);
  // double oldAverWavefront = aver_wavefront(G);
  // double oldRmsWavefront = rms_wavefront(G);

  // Call Sloan algorithm
  sloan_ordering(G, inv_perm.rbegin(), get(vertex_color, G),
    make_degree_map(G), get(vertex_priority, G), weight1, weight2);
  for (unsigned k = 0; k != inv_perm.size(); k++)
    perm[index_map[inv_perm[k]]]=k;

  // Find new bandwidth and profile
  int newBandwidth = bandwidth(G, make_iterator_property_map(&perm[0],
                               index_map, perm[0]));
  int newProfile = profile(G, make_iterator_property_map(&perm[0], index_map));
  // double newMaxWavefront = max_wavefront(G,
  //                          make_iterator_property_map(&perm[0], index_map));
  // double newAverWavefront = aver_wavefront(G,
  //                          make_iterator_property_map(&perm[0], index_map));
  // double newRmsWavefront = rms_wavefront(G,
  //                          make_iterator_property_map(&perm[0], index_map));

  // Print optimized sizes
  cout  << "reo: Optimized (original) bandwidth : "
        << newBandwidth  << " (" << oldBandwidth << ")" << endl;
  cout  << "reo: Optimized (original) profile   : "
        <<newProfile    << " (" << oldProfile  << ")"  << endl;
  // cout   << "rSloan: Optimized (original) maximum wavefront   : "
  //        <<newMaxWavefront   << " (" << oldMaxWavefront   << ")" << endl;
  // cout   << "rSloan: Optimized (original) average wavefront   : "
  //        <<newAverWavefront  << " (" << oldAverWavefront  << ")" << endl;
  // cout   << "rSloan: Optimized (original) rms wavefront       : "
  //        << newRmsWavefront  << " (" << oldRmsWavefront   << ")" << endl;

  // Check if optimization is needed
  if (newBandwidth > oldBandwidth) return -1;
  return 1;
}
