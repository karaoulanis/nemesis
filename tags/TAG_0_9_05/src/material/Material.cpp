/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2007 F.E.Karaoulanis [http://www.nemesis-project.org]  *
*                                                                             *
*   This program is free software; you can redistribute it and/or modify      *
*   it under the terms of the GNU General Public License version 2, as        *
*   published by the Free Software Foundation.                                *
*                                                                             *
*   This program is distributed in the hope that it will be useful,           *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License along   *
*   with this program; if not, write to the Free Software Foundation, Inc.,   *
*   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.               *
******************************************************************************/

//*****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
//*****************************************************************************

// Included files
#include <Material.h>

// Variable initialization
int Material::counter=0;

/**
 * Default constructor.
 */
Material::Material()
{
	myTracker=0;
}
/**
 * Constructor.
 * Creates a Material and increases Material's index by one.
 * @param ID Material id.
 * @param rho Density.
 * @param aT Thermal coefficient.
 */
Material::Material(int ID,double rho,double aT)
:DomainObject(ID)
{
	index=counter++;
	MatParams.resize(32);
	MatParams[30]=rho;
	MatParams[31]=aT;
	myTracker=0;
}
/**
 * Destructor.
 * Deletes only \a myTracker if exists. No other objects created here.
 */
Material::~Material()
{
	if(myTracker!=0) delete myTracker;
}
/**
 * Add a Tracker to a Material.
 * \a myTracker pointer should be \a null up to this point. If not this means
 * that a Tracker is already added and nothing is changed.
 * The Material deconstructor should take the responsibility to delete the Tracker.
 * Then track is called, in order to initialize the tracker.
 * @todo Check tracker initialization.
 */
void Material::addTracker()
{
	if(myTracker!=0) return;
    myTracker=new Tracker();
	this->track();
}
/**
 * Get the Tracker.
 * Return \a myTracker pointer either it is null or not. 
 * @todo Change this to a constant pointer.
 * @return  pointer to \a myTracker.
 */
Tracker* Material::getTracker()
{
	return myTracker;
}
/**
 * Add a record to the tracker.
 * This function should be overwritten in derived classes.
 */
void Material::track()
{
}
