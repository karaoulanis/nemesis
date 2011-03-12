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

#ifndef SRC_DOMAIN_DOMAIN_OBJECT_H_
#define SRC_DOMAIN_DOMAIN_OBJECT_H_

#include <iostream>
#include "feobject/fe_object.h"

// Forward declarations
class Domain;

/**
 * The DomainObject Class.
 */
class DomainObject: public FEObject {
 protected:
  int myID;
  static Packet thePacket;
  static Domain* pD;
  public:
  // Constructors
  DomainObject();
  DomainObject(int ID);
  virtual ~DomainObject();

  virtual int get_id();

  virtual const Packet& get_packet();
  virtual void set_packet(const Packet& p);
  virtual void save(std::ostream& /*s*/)  {}
  virtual void load(std::istream& /*s*/)  {}

  void set_domain(Domain* pDomain);
};
#endif  // SRC_DOMAIN_DOMAIN_OBJECT_H_
