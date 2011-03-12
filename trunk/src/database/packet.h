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

#ifndef SRC_DATABASE_PACKET_H_
#define SRC_DATABASE_PACKET_H_

#include <cstring>

struct Packet {
  const static int dblSize = 64;
  const static int intSize = 32;
  const static int chrSize = 32;
  const static int size = 1+1+dblSize+chrSize+1;
  const static double dblDefault;
  const static int    intDefault;

  int   tag;
  int   id;
  double  dblArray[dblSize];
  int   intArray[intSize];
  char  chrArray[chrSize];

  void zero();
  void print();
};
#endif  // SRC_DATABASE_PACKET_H_
