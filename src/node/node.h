/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2010 F.E.Karaoulanis [http://www.nemesis-project.org]  *
*                                                                             *
*   This program is free software; you can redistribute it and/or modify      *
*   it under the terms of the GNU General Public License version 3, as        *
*   published by the Free Software Foundation.                                *
*                                                                             *
*   This program is distributed in the hope that it will be useful,           *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License         *
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
******************************************************************************/

//*****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
//*****************************************************************************

#ifndef _NODE_H
#define _NODE_H

#include <sstream>

#include "containers/containers.h"
#include "domain/domain.h"
#include "domain/domain_object.h"
#include "elements/element.h"
#include "numeric/matrix.h"
#include "numeric/vector.h"
#include "tracker/tracker.h"

class Domain;
class Element;
const int MAX_NUMBER_OF_DOFS=16;

/**
 * The Node class.                                                
 */
class Node: public DomainObject
{
private:
	double x1,x2,x3;					// The nodal coordinates			
	IDContainer myConnectedElements;	// Array of ID's of connected elements
	IDContainer myActivatedDofs;		// Array of activated dofs
	IDContainer myConstrainedDofs;		// Array of constrained dofs
	int nActivatedDofs;

	Vector P;
	Vector dispTrial;
	Vector dispConvg;
	Vector velcTrial;
	Vector velcConvg;
	Vector acclTrial;
	Vector acclConvg;
	Matrix dispSensi;
	Matrix eigenVecs;

	Tracker* myTracker;

	int avgStress; ///@todo: remove this when full nodal recovering is implemented.
	Vector stress;
	Vector strain;
	bool isLoadApplied;

public:
	// Constructors and destructor
	Node();
	Node(int ID,double xc1,double xc2=0,double xc3=0);
	~Node();

	// Get nodal coordinates
	double getx1();
	double getx2();
	double getx3();

	// Connected elements
	int addEleToNode(Element *pElement);
	const IDContainer& getConnectedElements() const;
	bool isActive();

	// Activated Dofs
	int addDofToNode(int dof);
	const IDContainer& getActivatedDofs() const;
	int getActivatedDof(int localDof) const;
	int getnActivatedDofs() const;

	// Loads
	void addLoad(int dof,double value,double factor=1.0);
	void addInitialDisp(int dof,double disp);
	void addInitialVelc(int dof,double velc);
	bool existsLoad();
	void zeroLoad();
	const Vector& getR();


	// Functions that handle the trial and converged states
	void incTrialDisp(const Vector& du);
	void addTrialVelc(const Vector& dv);
	void addTrialAccl(const Vector& da);
	void setTrialDisp(const Vector& u);
	void setTrialVelc(const Vector& v);
	void setTrialAccl(const Vector& a);
	void commit();
	void rollback();

	const Vector& getDispTrial();
	const Vector& getVelcTrial();
	const Vector& getAcclTrial();
	const Vector& getDispConvg();
	const Vector& getVelcConvg();
	const Vector& getAcclConvg();
	double getDispTrialAtDof(int dof);
	double getVelcTrialAtDof(int dof);
	double getAcclTrialAtDof(int dof);
	double getDispConvgAtDof(int dof);
	double getVelcConvgAtDof(int dof);
	double getAcclConvgAtDof(int dof);

	void zeroStress();
	void addStress(const Vector& s);
	void multDisp(double facD);
	
	const Packet& getPacket();
	void setPacket(const Packet& p);
	void save(std::ostream& s);

	// Tracker member functions
	void addTracker();
	Tracker* getTracker();
	void track();

	// Sensitivity functions
	void initSensitivityMatrix(int nGrads);
	void commitSens(const Vector& v,int param);

	// Enrichment functions
	void evalLevelSets();

};

#endif
