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

#ifndef SRC_IMPOSER_IMPOSER_H_
#define SRC_IMPOSER_IMPOSER_H_

#include <map>
#include "containers/containers.h"

// Forward declarations
class Element;
class Model;
class Node;
class Constraint;

class Imposer {
 public:
  Imposer();
  Imposer(const std::map<int, Node*>& nodes,
          const std::map<int, Element*>& elements,
          const std::map<int, Constraint*>& constraints);
  virtual ~Imposer();
  virtual int impose(Model* model)=0;

 protected:
  const std::map<int, Node*>*       nodes_;
  const std::map<int, Element*>*    elements_;
  const std::map<int, Constraint*>* constraints_;

 private:
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  Imposer(const Imposer&);
  void operator=(const Imposer&);
};
#endif  // SRC_IMPOSER_IMPOSER_H_
