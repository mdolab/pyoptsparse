/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search - version 3.7.1        */
/*                                                                                     */
/*  Copyright (C) 2001-2015  Mark Abramson        - the Boeing Company, Seattle        */
/*                           Charles Audet        - Ecole Polytechnique, Montreal      */
/*                           Gilles Couture       - Ecole Polytechnique, Montreal      */
/*                           John Dennis          - Rice University, Houston           */
/*                           Sebastien Le Digabel - Ecole Polytechnique, Montreal      */
/*                           Christophe Tribes    - Ecole Polytechnique, Montreal      */
/*                                                                                     */
/*  funded in part by AFOSR and Exxon Mobil                                            */
/*                                                                                     */
/*  Author: Sebastien Le Digabel                                                       */
/*                                                                                     */
/*  Contact information:                                                               */
/*    Ecole Polytechnique de Montreal - GERAD                                          */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada                  */
/*    e-mail: nomad@gerad.ca                                                           */
/*    phone : 1-514-340-6053 #6928                                                     */
/*    fax   : 1-514-340-5665                                                           */
/*                                                                                     */
/*  This program is free software: you can redistribute it and/or modify it under the  */
/*  terms of the GNU Lesser General Public License as published by the Free Software   */
/*  Foundation, either version 3 of the License, or (at your option) any later         */
/*  version.                                                                           */
/*                                                                                     */
/*  This program is distributed in the hope that it will be useful, but WITHOUT ANY    */
/*  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    */
/*  PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.   */
/*                                                                                     */
/*  You should have received a copy of the GNU Lesser General Public License along     */
/*  with this program. If not, see <http://www.gnu.org/licenses/>.                     */
/*                                                                                     */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad               */
/*-------------------------------------------------------------------------------------*/
/**
 \file   Point.cpp
 \brief  Custom class for points (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-05
 \see    Point.hpp
 */
#include "Point.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
#ifdef MEMORY_DEBUG
int NOMAD::Point::_cardinality     = 0;
int NOMAD::Point::_max_cardinality = 0;
#endif

int NOMAD::Point::_display_limit   = NOMAD::DEFAULT_POINT_DISPLAY_LIMIT;

/*-----------------------------------------------------------*/
/*                         constructor                       */
/*-----------------------------------------------------------*/
NOMAD::Point::Point ( int n , const NOMAD::Double & d ) : _n (n) , _coords (NULL)
{
#ifdef MEMORY_DEBUG
	++NOMAD::Point::_cardinality;
	if ( NOMAD::Point::_cardinality > NOMAD::Point::_max_cardinality )
		++NOMAD::Point::_max_cardinality;
#endif
	if (_n > 0)
    {
		_coords = new NOMAD::Double [_n];  
		if ( d.is_defined() )
			std::fill ( _coords , _coords+_n , d );
	}
	else
		_n = 0;
}

/*-----------------------------------------------------------*/
/*                        copy constructor                   */
/*-----------------------------------------------------------*/
NOMAD::Point::Point ( const NOMAD::Point & p ) : _n (p._n) , _coords (NULL)
{
#ifdef MEMORY_DEBUG
	++NOMAD::Point::_cardinality;
	if ( NOMAD::Point::_cardinality >= NOMAD::Point::_max_cardinality )
		++NOMAD::Point::_max_cardinality;
#endif
	if ( _n > 0 )
    {
		NOMAD::Double       * p1 =   _coords = new NOMAD::Double [_n];
		const NOMAD::Double * p2 = p._coords;
		for ( int k = 0 ; k < _n ; ++k , ++p1 , ++p2 )
			*p1 = *p2;
	}
}

/*-----------------------------------------------*/
/*                    destructor                 */
/*-----------------------------------------------*/
NOMAD::Point::~Point ( void )
{
#ifdef MEMORY_DEBUG
	--NOMAD::Point::_cardinality;
#endif
	delete [] _coords;
}

/*-----------------------------------------------*/
/*   This method changes the point's dimension   */
/*   and sets all values to d                    */
/*-----------------------------------------------*/
void NOMAD::Point::reset ( int n , const NOMAD::Double & d )
{
	if ( n <= 0 )
    {
		_n = 0; 
		delete [] _coords;
		_coords = NULL;
	}
	else
    {
		if ( _n != n )
        {
			delete [] _coords;
			_n      = n;
			_coords = new NOMAD::Double [_n];
		}
		if ( d.is_defined() )
			std::fill ( _coords , _coords+_n , d );
	}
}

/*----------------------------------------------------------------*/
/*  This method changes the array's dimension (and keeps values)  */
/*----------------------------------------------------------------*/
void NOMAD::Point::resize ( int n )
{
	if ( n == _n )
		return;
    
	if ( n <= 0 )
    {
		_n = 0; 
		delete [] _coords;
		_coords = NULL;
		return;
	}
	NOMAD::Double * new_coords = new NOMAD::Double [n];
	if ( _coords )
    {
		int min = ( n < _n ) ? n : _n;
		
		NOMAD::Double       * p1 = new_coords;
		const NOMAD::Double * p2 = _coords;
		
		for ( int i = 0 ; i < min ; ++i , ++p1 , ++p2 )
			*p1 = *p2;
		
		delete [] _coords;
	}
	_coords = new_coords;
	_n      = n;
}

/*-----------------------------------------------------------*/
/*                       '[]' operators                      */
/*-----------------------------------------------------------*/

// const version:
const NOMAD::Double & NOMAD::Point::operator [] ( int i ) const
{
	if ( !_coords )
		throw NOMAD::Point::Not_Defined ( "Point.cpp" , __LINE__ ,
										 "operator x[i] (const): 'x' not defined" );
	if ( i < 0 || i >= _n )
		throw NOMAD::Point::Bad_Access ( "Point.cpp" , __LINE__ ,
										"operator x[i] (const): 'i' outside the array's bounds." );
	return _coords[i];
}

// non-const version:
NOMAD::Double & NOMAD::Point::operator [] ( int i )
{
	if ( !_coords )
		throw NOMAD::Point::Not_Defined ( "Point.cpp" , __LINE__ ,
										 "operator x[i]: 'x' not defined" );
	if ( i < 0 || i >= _n )
		throw NOMAD::Point::Bad_Access ( "Point.cpp" , __LINE__ ,
										"operator x[i] (const): 'i' outside the array's bounds." );
	return _coords[i];
}

/*-----------------------------------------------------------*/
/*                     affectation operator                  */
/*-----------------------------------------------------------*/
const NOMAD::Point & NOMAD::Point::operator = ( const NOMAD::Point & p )
{
	if ( this == &p )
		return *this;
	
	if ( _n != p._n )
    {
		delete [] _coords;
		_n = p._n;
		if (_n > 0)
			_coords = new NOMAD::Double [_n];
		else
			_coords = NULL;
	}
	
	NOMAD::Double       * p1 =   _coords;
	const NOMAD::Double * p2 = p._coords;
	for ( int k = 0 ; k < _n ; ++k , ++p1 , ++p2 )
		*p1 = *p2;
	
	return *this;
}

/*------------------------------------*/
/*  projection to mesh of size delta  */
/*      (*this = ref + k * delta)     */
/*------------------------------------*/
void NOMAD::Point::project_to_mesh ( const NOMAD::Point & ref   ,
									const NOMAD::Point & delta ,
									const NOMAD::Point & lb    ,
									const NOMAD::Point & ub      )
{
	if ( delta._n != _n               ||
		ref._n   != _n               ||
		( lb._n > 0 && lb._n != _n ) ||
		( ub._n > 0 && ub._n != _n )    )
		throw NOMAD::Point::Bad_Operation ( "Point.cpp" , __LINE__ ,
										   "Point::project_to_mesh(): invalid Point sizes" );
	
	NOMAD::Double * p  = _coords       ,
	* pr = ref._coords   ,
	* pd = delta._coords ,
	* pl = lb._coords    ,
	* pu = ub._coords;
	int k;
	
	if ( lb._n == 0 && ub._n == 0 )
		for ( k = 0 ; k < _n ; ++k , ++pr , ++p , ++pd )
			p->project_to_mesh ( *pr , *pd );
	else if ( lb._n == 0 )
		for ( k = 0 ; k < _n ; ++k , ++pr , ++p , ++pd , ++pu )
			p->project_to_mesh ( *pr , *pd , NOMAD::Double() , *pu );
	else if ( ub._n == 0 )
		for ( k = 0 ; k < _n ; ++k , ++pr , ++p , ++pd , ++pl )
			p->project_to_mesh ( *pr , *pd , *pl );
	else
		for ( k = 0 ; k < _n ; ++k , ++pr , ++p , ++pd , ++pl , ++pu )
			p->project_to_mesh ( *pr , *pd , *pl , *pu );
}


/*-----------------------------------------------------------*/
/*                             display                       */
/*-----------------------------------------------------------*/
void NOMAD::Point::display ( const NOMAD::Display & out ,
							const std::string    & sep ,
							int                    w   ,
							int                    lim   ) const
{
	int nm1 = _n-1;
	
	// for a limited display of maximum lim elements:
	if ( lim > 0 && lim < _n )
	{
		
		int l1 = (lim + 1) / 2 , l2 = lim / 2 , i;
		
		// first coordinates:
		for ( i = 0 ; i < l1 ; ++i )
			out << std::setw ( w ) << _coords[i] << sep;
		
		// separator:
		out << "..." << sep;
		
		// last coordinates:
		for ( i = _n - l2 ; i < nm1 ; ++i )
			out << std::setw ( w ) << _coords[i] << sep;
	}
	
	// normal display (lim <= 0 or lim >= _n):
	else 
	{
		const NOMAD::Double * p = _coords;
		for ( int i = 0 ; i < nm1 ; ++i , ++p )
			out << std::setw ( w ) << *p << sep;
	}
	
	// last coordinate (different because there is no separator after that):
	if ( _n > 0 )
		out << std::setw ( w ) << _coords[nm1];
}

/*-----------------------------------------------------------*/
/*              input (can read undefined coordinates)       */
/*-----------------------------------------------------------*/
std::istream & NOMAD::operator >> ( std::istream & in , NOMAD::Point & p )
{
	int n = p.size();
	for ( int k = 0 ; k < n ; ++k )
		in >> p[k];
	if ( in.fail() )
		throw NOMAD::Point::Bad_Input ( "Point.cpp" , __LINE__ , "in >> x: bad input" );
	return in;
}

/*-----------------------------------------------------------*/
/*  set the point's coordinate with the array 'a' of size n  */
/*  the Point's dimension is changed to n                    */
/*-----------------------------------------------------------*/
void NOMAD::Point::set ( int n , const NOMAD::Double * a )
{
	if ( n <= 0 || !a )
		return;
	
	if ( _n != n )
    {
		delete [] _coords;
		_n      = n;
		_coords = new NOMAD::Double [_n];
	}
	
	NOMAD::Double * p = _coords;
	for ( int k = 0 ; k < _n ; ++k , ++p , ++a )
		*p = *a;
}

/*-----------------------------------------------------------*/
/*        check if all the point values are defined          */
/*-----------------------------------------------------------*/
bool NOMAD::Point::is_complete ( void ) const
{
	if ( _n <= 0 )
		return false;
	const NOMAD::Double * p = _coords;
	for ( int i = 0 ; i < _n ; ++i , ++p )
		if ( !p->is_defined() )
			return false;
	return true;
}

/*---------------------------------------------------------------*/
/*  check if at least one value is defined in the _coords array  */
/*---------------------------------------------------------------*/
bool NOMAD::Point::is_defined ( void ) const
{
	if ( _n <= 0 )
		return false;
	const NOMAD::Double * p = _coords;
	for ( int i = 0 ; i < _n ; ++i , ++p )
		if ( p->is_defined() )
			return true;
	return false;
}

/*---------------------------------------------------------------*/
/*          count the number of values that are defined          */
/*---------------------------------------------------------------*/
int NOMAD::Point::nb_defined ( void ) const
{
	const NOMAD::Double * p = _coords;
	int                   k = 0;
	for ( int i = 0 ; i < _n ; ++i , ++p )
		if ( p->is_defined() )
			++k;
	return k;
} 

/*-----------------------------------------------------------*/
/*                           negation                        */
/*-----------------------------------------------------------*/
const NOMAD::Point NOMAD::Point::operator - ( void ) const
{
	NOMAD::Point          tmp (_n);
	NOMAD::Double       * p1 = tmp._coords;
	const NOMAD::Double *  p2 = _coords;
	for ( int k = 0 ; k < _n ; ++k , ++p1 , ++p2 )
		*p1 = - *p2;
	return tmp;
}

/*----------------------------------------------------------*/
/*                    scalar multiplication                 */
/*----------------------------------------------------------*/
const NOMAD::Point & NOMAD::Point::operator *= ( const NOMAD::Double & d )
{
	NOMAD::Double * p = _coords;
	for ( int k = 0 ; k < _n ; ++k , ++p )
		*p *= d;
	return *this;
}

/*----------------------------------------------------------*/
/*               multiplication of two points               */
/*----------------------------------------------------------*/
const NOMAD::Point NOMAD::Point::operator * ( const NOMAD::Point & p ) const
{
	if ( p._n != _n )
		throw NOMAD::Point::Bad_Operation ( "Point.cpp" , __LINE__ ,
										   "x * y: x.size != y.size" );
	NOMAD::Point          tmp ( _n );
	NOMAD::Double       * p1 = tmp._coords;
	const NOMAD::Double * p2 =     _coords;
	const NOMAD::Double * p3 =   p._coords;
    
	for ( int k = 0 ; k < _n ; ++k , ++p1 , ++p2 , ++p3 )
		*p1 = *p2 * *p3;
	
	return tmp;
}

/*----------------------------------------------------------*/
/*                     division of two points               */
/*----------------------------------------------------------*/
const NOMAD::Point NOMAD::Point::operator / ( const NOMAD::Point & p ) const
{
	if ( p._n != _n )
		throw NOMAD::Point::Bad_Operation ( "Point.cpp" , __LINE__ ,
										   "x / y: x.size != y.size" );
	NOMAD::Point          tmp ( _n );
	NOMAD::Double       * p1 = tmp._coords;
	const NOMAD::Double * p2 =     _coords;
	const NOMAD::Double * p3 =   p._coords;
    
	for ( int k = 0 ; k < _n ; ++k , ++p1 , ++p2 , ++p3 )
		*p1 = *p2 / *p3;
	
	return tmp;
}

/*----------------------------------------------------------*/
/*                   addition of two points                 */
/*----------------------------------------------------------*/
const NOMAD::Point NOMAD::Point::operator + ( const NOMAD::Point & p ) const
{
	if ( p._n != _n )
		throw NOMAD::Point::Bad_Operation ( "Point.cpp" , __LINE__ ,
										   "x + y: x.size != y.size" );
	NOMAD::Point          tmp ( _n );
	NOMAD::Double       * p1 = tmp._coords;
	const NOMAD::Double * p2 =     _coords;
	const NOMAD::Double * p3 =   p._coords;
    
	for ( int k = 0 ; k < _n ; ++k , ++p1 , ++p2 , ++p3 )
		*p1 = *p2 + *p3;
	
	return tmp;
}

/*----------------------------------------------------------*/
/*                  subtraction of two points               */
/*----------------------------------------------------------*/
const NOMAD::Point NOMAD::Point::operator - ( const NOMAD::Point & p ) const
{
	if ( p._n != _n )
		throw NOMAD::Point::Bad_Operation ( "Point.cpp" , __LINE__ ,
										   "x - y: x.size != y.size" );
	NOMAD::Point          tmp ( _n );
	NOMAD::Double       * p1 = tmp._coords;
	const NOMAD::Double * p2 =     _coords;
	const NOMAD::Double * p3 =   p._coords;
    
	for ( int k = 0 ; k < _n ; ++k , ++p1 , ++p2 , ++p3 )
		*p1 = *p2 - *p3;
	
	return tmp;
}

/*--------------------------------------------------------------------------*/
/* comparison operator '<': it is used to find and store the points         */
/*                           in the cache                                   */
/*--------------------------------------------------------------------------*/
bool NOMAD::Point::operator < ( const NOMAD::Point & p ) const
{
	if ( this == &p )
		return false;
	
	if ( _n < p._n )
		return true;
	if ( _n > p._n )
		return false;
	
	const NOMAD::Double * p1 =   _coords;
	const NOMAD::Double * p2 = p._coords;
	
	for ( int k = 0 ; k < _n ; ++k , ++p1 , ++p2 )
    {
		
		if ( *p1 < *p2 )
			return true;
		
		if ( *p1 > *p2 )
			return false;
	}
	
	return false;
}

/*---------------------------------------------------------------------*/
/*  the same as operator < but with consideration of undefined values  */
/*---------------------------------------------------------------------*/
bool NOMAD::Point::comp_with_undef ( const NOMAD::Point & p ) const
{
	if ( this == &p )
		return false;
	
	if ( _n < p._n )
		return true;
	if ( _n > p._n )
		return false;
	
	const NOMAD::Double * p1 =   _coords;
	const NOMAD::Double * p2 = p._coords;
	
	bool p1d , p2d;
	
	for ( int k = 0 ; k < _n ; ++k , ++p1 , ++p2 )
    {
		p1d = p1->is_defined();
		p2d = p2->is_defined();
		
		if ( !p1d && !p2d )
			continue;
		
		if ( !p1d )
			return true;
		
		if ( !p2d )
			return false;
		
		if ( *p1 < *p2 )
			return true;
		
		if ( *p1 > *p2 )
			return false;
	}
	return false;
}

/*-----------------------------------------------------------*/
/*                           operator ==                     */
/*-----------------------------------------------------------*/
bool NOMAD::Point::operator == ( const NOMAD::Point & p ) const
{
	if ( this == &p )
		return true;
	if ( p._n != _n )
		return false;
	
	const NOMAD::Double * p1 =   _coords;
	const NOMAD::Double * p2 = p._coords;
	for ( int k = 0 ; k < _n ; ++k , ++p1 , ++p2 )
		if ( *p1 != *p2 )
			return false;
	
	return true;
}

/*-----------------------------------------------------------*/
/*        computation of the angle with another point        */
/*-----------------------------------------------------------*/
const NOMAD::Double NOMAD::Point::get_angle ( const NOMAD::Point & x ) const
{
	if ( _n != x._n )
		return NOMAD::Double();
	
	NOMAD::Double inner_product = 0.0 , norm_1 = 0.0 , norm_2 = 0.0;
	
	const NOMAD::Double * p1 =   _coords;
	const NOMAD::Double * p2 = x._coords;
	
	for ( int i = 0 ; i < _n ; ++i , ++p1 , ++p2 )
    {
		norm_1        += *p1 * *p1;
		norm_2        += *p2 * *p2;
		inner_product += *p1 * *p2;
	}
	
	if ( norm_1 == 0.0 || norm_2 == 0.0 )
		return NOMAD::Double();
	
	return acos ( ( inner_product / ( norm_1.sqrt() * norm_2.sqrt() ) ).value() );
}
