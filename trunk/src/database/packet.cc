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
* along with this program.  If not, see < http://www.gnu.org/licenses/>.        *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

// Corresponding header file
#include "database/packet.h"

// C++ system files
#include < iostream>

// Namespaces
using std::cout;
using std::endl;

// Constant initializations
const double Packet::dblDefault=-9999999999.99;
const int    Packet::intDefault=-999999;

void Packet::zero() {
  tag = 0;
  id = 0;
  int i;
  for (i = 0;i < dblSize;i++) dblArray[i]=dblDefault;
  for (i = 0;i < intSize;i++) intArray[i]=intDefault;
  strcpy(chrArray, "empty");
}
void Packet::print() {
  cout << "Tag      : " <<tag   <<endl;
  cout << "Id       : " <<id    <<endl;
  int i;
  cout << "Doubles  : ";
  for (i = 0;i < dblSize; i++) cout << dblArray[i]  <<'\t'; cout << endl;
  cout << "Integers : ";
  for (i = 0; i < intSize; i++) cout << intArray[i]  <<'\t'; cout << endl;
  cout << "Chars    : " <<chrArray  <<endl;
}
