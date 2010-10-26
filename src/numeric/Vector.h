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
#ifndef __VECTOR_H
#define __VECTOR_H

#include <iostream>
#include <cmath>
#include "numeric/ArrayCheck.h"
#include "numeric/Numeric.h"
#include "numeric/Lapack.h"

class Vector
{
protected: 
	int size_;
	double* data_;
public: 
	/**
	 * Default constructor.
	 * Initializes everything to zero.
	 */
	Vector()
		:size_(0),data_(0)
	{
	}
	/**
	 * Constructor.
	 * Allocates a size n double. 
	 * Exception handling for bad allocation is provided.
	 * @param n The size of the vector.
	 */
	Vector(int n)
		:size_(n)
	{
		try
		{
			data_=new double[size_];
		}
		catch(std::bad_alloc)
		{
			throw SException("[nemesis:%d] %s",1001,"Run out of memory.\n");
		}
	}
	/**
	 * Constructor.
	 * Allocates a size n double and initializes entries to c. 
	 * Exception handling for bad allocation is provided.
	 * @param n The size of the vector.
	 * @param c Initial value for all entries.
	 */
	Vector(int n,double c)
		:size_(n)
	{
		try
		{
			data_=new double[size_];
		}
		catch(std::bad_alloc)
		{
			throw SException("[nemesis:%d] %s",1001,"Run out of memory.\n");
		}
		for(int i=0;i<size_;i++) data_[i]=c;
	}
	/**
	 * Copy constructor.
	 * Allocates a v.size double and copies entries from Vector v2.
	 * Exception handling for bad allocation is provided.
	 * @param v The Vector that is copied.
	 */
	Vector(const Vector& v)
		:size_(v.size_)
	{
		if(size_!=0)
		{
			try
			{
				data_=new double[size_];
			}
			catch(std::bad_alloc)
			{
			throw SException("[nemesis:%d] %s",1001,"Run out of memory.\n");
			}
			for(int i=0;i<size_;i++) data_[i]=v.data_[i];
		}
		else data_=0;
	}	
	/**
	 * Destructor.
	 */
	~Vector()
	{
		if(data_!=0) delete[] data_;
	}
	/**
	 * Copy assignment.
	 * Size checking is provided for the debug version.
	 */
	inline Vector& operator=(const Vector& v)
	{
		#ifdef _DEBUG
		array_size_check(v.size_,size_);
		#endif
		if (this!=&v) 
			for (int i=0;i<size_;i++) data_[i]=v.data_[i];
		return *this;
	}
	/**
	 * Implements [i] operator: v[i]
	 * Range checking is provided for the debug version.
	 */
	inline double& operator[](int i)
	{
		#ifdef _DEBUG
		array_range_check(i,size_);
		#endif
		return data_[i];
	}	
	/**
	 * Implements [i] operator: v[i]
	 * Range checking is provided for the debug version.
	 */
	inline double operator[](int i) const
	{
		#ifdef _DEBUG
		array_range_check(i,size_);
		#endif
		return data_[i];
	}	
	/**
	 * Returns the size of the Vector.
	 */
	inline const int size() const
	{
		return size_;
	}
	/**
	 * Returns a pointer to Vector data.
	 */
	inline double* data()
	{
		return data_;
	}
	/**
	 * Returns a const pointer to Vector data.
	 */
	inline double* data() const
	{
		return data_;
	}
	/**
	 * Resizes the Vector.
	 * Does not preserves data and does not initialize entries.
	 * If new Vector of the same size no allocation takes place.
	 * Exception handling for bad allocation is provided.
	 * @param n The size of the vector.
	 */
	inline void resize(int n)
	{
		if(size_!=n)
		{
			size_=n;
			if(data_!=0) delete[] data_;
			try
			{
				data_=new double[size_];
			}
			catch(std::bad_alloc)
			{
				throw SException("[nemesis:%d] %s",1001,"Run out of memory.\n");
			}
		}
	}
	/**
	 * Resizes the Vector.
	 * Does not preserves data but initializes all entries to c.
	 * If new Vector of the same size no allocation takes place.
	 * Exception handling for bad allocation is provided.
	 * @param n The size of the vector.
	 * @param c Initial value for all entries.
	 */
	inline void resize(int n,double c)
	{
		if(size_!=n)
		{
			size_=n;
			if(data_!=0) delete[] data_;
			try
			{
				data_=new double[size_];
			}
			catch(std::bad_alloc)
			{
				throw SException("[nemesis:%d] %s",1001,"Run out of memory.\n");
			}
		}
		for(int i=0;i<size_;i++) data_[i]=c;
	}
	/**
	 * Clears the contents of a Vector.
	 */
	inline void clear()
	{
		for(int i=0;i<size_;i++) data_[i]=0.;
	}
	/**
	 * Implements += operator: this+=c.
	 */
	inline Vector& operator+=(const double c)
	{
		for(int i=0;i<size_;i++) data_[i]+=c;
		return *this;
	}
	/**
	 * Implements -= operator: v1-=c.
	 */
	inline Vector& operator-=(const double c)
	{
		for(int i=0;i<size_;i++) data_[i]-=c;
		return *this;
	}
	/**
	 * Implements *= operator: v1*=c.
	 */
	inline Vector& operator*=(const double c)
	{
		for(int i=0;i<size_;i++) data_[i]*=c;
		return *this;
	}
	/**
	 * Implements /= operator: v1/=c.
	 */
	inline Vector& operator/=(const double c)
	{
		for(int i=0;i<size_;i++) data_[i]/=c;
		return *this;
	}
	/**
	 * Implements * operator: v=this*c.
	 * A temporary object is created.
	 */
	inline Vector operator*(const double c) const
	{
		Vector res(*this);
		res*=c;
		return res;
	}
 	/**
	 * Implements * operator: v=c*this, using member * operator.
	 * A temporary object is created.
	 */
    inline friend Vector operator*(const double c,const Vector& v)
	{
		return v*c;
	}	
	/**
	 * Implements / operator: v=this/c.
	 * A temporary object is created.
	 */
	inline Vector operator/(const double c) const
	{
		Vector res(*this);
		res/=c;
		return res;	
	}
	/**
	 * Implements += operator: this+=v.
	 * Size checking is provided for the debug version.
	 */
	inline Vector& operator+=(const Vector& v)
	{
		#ifdef _DEBUG
		array_size_check(v.size_,size_);
		#endif
		for(int i=0;i<size_;i++) data_[i]+=v.data_[i];
		return *this;
	}
	/**
	 * Implements -= operator: this-=v.
	 * Size checking is provided for the debug version.
	 */
	inline Vector& operator-=(const Vector& v)
	{
		#ifdef _DEBUG
		array_size_check(v.size_,size_);
		#endif
		for(int i=0;i<size_;i++) data_[i]-=v.data_[i];
		return *this;
	}
	/**
	 * Implements c0*this+=c*v.
	 * Size checking is provided for the debug version.
	 */
	inline void add_cV(double c,const Vector& v,double c0=1.0)
	{
		#ifdef _DEBUG
		array_size_check(v.size_,size_);
		#endif
		if(c==0.0&&c0!=0)	return;
		else if(c0==0.0)	for(int i=0;i<size_;i++) data_[i]=0;
		else if(c0!=1.0)	for(int i=0;i<size_;i++) data_[i]*=c0;
		for(int i=0;i<size_;i++) data_[i]+=c*v.data_[i];
	}
	/**
	 * Implements + operator: this+v.
	 * A temporary object is created.
	 */
	inline Vector operator+(const Vector& v)
	{
		Vector res(*this);
		res.add_cV(1.0,v,1.0);
		return res;
	}
	/**
	 * Implements - operator: this-v.
	 * A temporary object is created.
	 * Size checking is provided for the debug version.
	 */
	inline Vector operator-(const Vector& v)
	{
		Vector res(*this);
		res.add_cV(-1.0,v,1.0);
		return res;
	}
	/**
	 * Implements + operator: +this.
	 */
	inline Vector& operator+()
	{
		return *this;
	}
	/**
	 * Implements - operator: -this.
	 * A temporary object is created.
	 */
	inline Vector operator-()
	{
		Vector res(size_);
		for(int i=0;i<size_;i++) res.data_[i]=-data_[i];
		return res;
	}
	/**
	 * Implements dot product.
	 * Size checking is provided for the debug version.
	 */
	inline friend double operator*(const Vector& v1,const Vector& v2)
	{
		#ifdef _DEBUG
		array_size_check(v1.size_,v2.size_);
		#endif
		if(v1.size_==0) return 0;
		double res=v1.data_[0]*v2.data_[0];
		for(int i=1;i<v1.size_;i++) res+=v1.data_[i]*v2.data_[i];
		return res;
	}
	/**
	 * Implements cross product.
	 * Size checking is provided for the debug version.
	 */
	inline friend Vector cross(const Vector& v1,const Vector& v2)
	{
		#ifdef _DEBUG
		array_size_check(v1.size_,3);
		array_size_check(v2.size_,3);
		#endif
		Vector res(3);
		res.data_[0]=v1.data_[1]*v2.data_[2]-v1.data_[2]*v2.data_[1];
		res.data_[1]=v1.data_[2]*v2.data_[0]-v1.data_[0]*v2.data_[2];
		res.data_[2]=v1.data_[0]*v2.data_[1]-v1.data_[1]*v2.data_[0];
		return res;
	}
	/**
	 * Implements Euclidean norm.
	 * The implementation provides overflow protection.
	 */
	inline double twonorm() const
	{
		if(size_==0) return 0;
		double norm=fabs(data_[0]);
		for (int i=1;i<size_;i++) 
		{
			double vi=fabs(data_[i]);
			if (norm<100) 
			{
				norm=sqrt(norm*norm+vi*vi);
			}
			else 
			{  
				double temp=vi/norm;
				norm*=sqrt(1.0+temp*temp);
			}
		} 
		return norm;
	}
	/**
	 * Implements maximum norm.
	 */
	inline double maxnorm() const 
	{
		if(size_==0) return 0;
		double norm=fabs(data_[0]);
		for (int i=1;i<size_;i++) if(norm<fabs(data_[i])) norm=fabs(data_[i]);
		return norm;
	}
	/**
	 * Returns the normal of a Vector.
	 * Provides check for zero division.
	 */
	inline Vector& normalize() 
	{
		double norm=twonorm();
		if(num::tiny(norm))
			throw SException("[nemesis:%d] %s",9999,"Zero vector length encountered.");
		for(int i=0;i<size_;i++) data_[i]/=norm;
		return *this;
	}
	/**
	 * Appends the entries of a Vector v.
	 * Range checking is provided for the debug version.
	 * @param v The Vector to be copied.
	 * @param row The row from which copy starts.
	 * @param c A factor to be multiplied with the appended entries.
	 * @param c0 A factor to be multiplied with the existing entries.
	 */
	inline Vector& append(const Vector& v,int row,double c=1.0,double c0=0.)
	{
		#ifdef _DEBUG
		array_range_check(row+v.size_,size_);
		#endif
		if(c0==0.0)			for (int i=0;i<v.size_;i++) data_[row+i]=0.;
		else if(c0!=1.0)	for (int i=0;i<v.size_;i++) data_[row+i]*=c0;
		for(int i=0;i<v.size_;i++) data_[row+i]+=c*v.data_[i];
		return *this;
	}
	inline friend double angle(const Vector& v1,const Vector& v2)
	{
		double denom=v1.twonorm()*v2.twonorm();
		if(denom==0.) return -360.;
		double beta=(v1*v2)/denom;
		if(beta>= 1.)	beta= 0.999999999999999999;
		if(beta<=-1.)	beta=-0.999999999999999999;
		return acos(beta)*180/num::pi;
	}
	/**
	 * Return mean stress.
	 * sb=1/3*(sxx+syy+szz)
	 */
	inline double sb() const
	{
		#ifdef _DEBUG
		array_size_check(size_,6);
		#endif
		return num::d13*(data_[0]+data_[1]+data_[2]);
	}
	/**
	 * Return invariant I1.
	 * I1=sxx+syy+szz
	 */
	inline double I1() const
	{
		#ifdef _DEBUG
		array_size_check(size_,6);
		#endif
		return data_[0]+data_[1]+data_[2];
	}
	/**
	 * Return invariant I2.
	 * @todo this
	 */
	inline double I2() const
	{
		#ifdef _DEBUG
		array_size_check(size_,6);
		#endif
		return 0.;
	}
	/**
	 * Return invariant I3.
	 * I3=sxx*syy*szz+2.*txy*tzx*tyz-sxx*tyz*tyz-syy*tzx*tzx-szz*txy*txy
	 */
	inline double I3() const
	{
		#ifdef _DEBUG
		array_size_check(size_,6);
		#endif
		return   data_[0]*data_[1]*data_[2]
			+2.0*data_[3]*data_[4]*data_[5]
				-data_[0]*data_[4]*data_[4]
				-data_[1]*data_[5]*data_[5]
				-data_[2]*data_[3]*data_[3];
	}
	/**
	 * Return invariant J1.
	 * J1=0
	 */
	inline double J1() const
	{
		#ifdef _DEBUG
		array_size_check(size_,6);
		#endif
		return 0;
	}
	/**
	 * Return invariant J2.
	 * J2=0.5*(sx*sx+sy*sy+sz*sz)+txy*txy+tyz*tyz+tzx*tzx
	 */
	inline double J2() const
	{
		#ifdef _DEBUG
		array_size_check(size_,6);
		#endif
		double sm=(data_[0]+data_[1]+data_[2])/3.;
		double sx=data_[0]-sm;
		double sy=data_[1]-sm;
		double sz=data_[2]-sm;
		return 0.5*(sx*sx+sy*sy+sz*sz)
				+data_[3]*data_[3]+data_[4]*data_[4]+data_[5]*data_[5];
	}
	/**
	 * Return invariant J3.
	 * J3=sx*sy*sz+2.*txy*tyz*tzx-sx*tyz*tyz-sy*tzx*tzx-sz*txy*txy
	 */
	inline double J3() const
	{
		#ifdef _DEBUG
		array_size_check(size_,6);
		#endif
		double sm=(data_[0]+data_[1]+data_[2])/3.;
		double sx=data_[0]-sm;
		double sy=data_[1]-sm;
		double sz=data_[2]-sm;
		return sx*sy*sz+2.*data_[3]*data_[4]*data_[5]
				-sx*data_[4]*data_[4]
				-sy*data_[5]*data_[5]
				-sz*data_[3]*data_[3];
	}
	/**
	 * Return theta.
	 * theta=-1.5*sqrt(3.)*J3/(sqJ2*sqJ2*sqJ2)/3.0
	 * @todo check this
	 */
	inline double theta() const
	{
		#ifdef _DEBUG
		array_size_check(size_,6);
		#endif
		double facJ3=this->J3();
		double facSqJ2=sqrt(this->J2());
		if(abs(facJ3)>1e10)
		{
			double fac=.001;
			double sm=(data_[0]+data_[1]+data_[2])/3.0;
			double sx=fac*(data_[0]-sm);
			double sy=fac*(data_[1]-sm);
			double sz=fac*(data_[2]-sm);
			double txy=fac*data_[3];
			double tyz=fac*data_[4];
			double tzx=fac*data_[5];
			facSqJ2=sqrt(0.5*(sx*sx+sy*sy+sz*sz)+txy*txy+tyz*tyz+tzx*tzx);
			facJ3=sx*sy*sz+2.*txy*tyz*tzx-sx*tyz*tyz-sy*tzx*tzx-sz*txy*txy;
		}
		double arg;
		if(facSqJ2<1e-12)	arg= 0.0;
		else				arg= (-1.5*sqrt(3.)*facJ3/(facSqJ2*facSqJ2*facSqJ2));
		if(arg>1.0)			arg= 1.0;
		else if(arg<-1.)	arg=-1.0;
		return asin(arg)/3.0;
	}
	/**
	 * Return p (hydrostatic).
	 * p=-I1/3.0 (sign convention).
	 */
	inline double p() const
	{
		#ifdef _DEBUG
		array_size_check(size_,6);
		#endif
		return -I1()/3.0;
	}
	/**
	 * Return q (deviatoric).
	 * p=sqrt(3*J2)
	 */
	inline double q() const
	{
		#ifdef _DEBUG
		array_size_check(size_,6);
		#endif
		return sqrt(3.*J2());
	}
	/**
	 * Return dp/ds.
	 * @todo remove?
	 */
	inline Vector dpds() const
	{
		#ifdef _DEBUG
		array_size_check(size_,6);
		#endif
		Vector ret(6,0.);
		ret[0]=-num::d13;
		ret[1]=-num::d13;
		ret[2]=-num::d13;
		return ret;
	}
	/**
	 * Return dq/ds. 
	 * @todo remove?
	 */
	inline Vector dqds() const
	{
		#ifdef _DEBUG
		array_size_check(size_,6);
		#endif
		double sm =sb();
		double sx =data_[0]-sm;
		double sy =data_[1]-sm;
		double sz =data_[2]-sm;
		double txy=data_[3];
		double tyz=data_[4];
		double tzx=data_[5];

		Vector ret(6,0.);
		ret[0]=sx;
		ret[1]=sy;
		ret[2]=sz;
		ret[0]=2*txy;
		ret[1]=2*tyz;
		ret[2]=2*tzx;
		ret*=num::sq3/(2*J2());
		return ret;
	}
	/**
	 * Return eigenvalues.
	 * Eigenvalues are returned in descending order.
	 */
	const Vector& eigenvalues()
	{
		#ifdef _DEBUG
		array_size_check(size_,6);
		#endif
		int N=3;
		char JOBZ='N';
		char UPLO='L';
		int LDA=3;
		int LWORK=102;
		int INFO;
		static Vector res(3);
		static Vector A(9);
		static Vector WORK(LWORK);
		res.clear();

		A[0]=data_[0];	A[1]=data_[3];	A[2]=data_[5];
		A[3]=0.;		A[4]=data_[1];	A[5]=data_[4];
		A[6]=0.;		A[7]=0.;		A[8]=data_[2];

		dsyev(&JOBZ,&UPLO,&N,A.data(),&LDA,res.data(),WORK.data(),&LWORK,&INFO,1,1);
	//	dsyev(&JOBZ,&UPLO,&N,A.data(),&LDA,res.data(),WORK.data(),&LWORK,&INFO);
	//	cout<<"Optimal LWORK : "<<WORK[0]<<endl;
		double d=res[0];
		res[0]=res[2];
		res[2]=d;
		return res;
	}
	inline friend std::ostream& operator<<(std::ostream& s, const Vector& v)
	{
		s<<1100<<' '<<v.size_<<' ';
		for(int i=0;i<v.size_;i++) s<<v.data_[i]<<' ';
		return s;
	}

};
#endif
