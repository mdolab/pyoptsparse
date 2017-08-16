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
  \file   Cache_File_Point.cpp
  \brief  Class for points in binary files (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-06
  \see    Cache_File_Point.hpp
*/
#include "Cache_File_Point.hpp"
#include "Eval_Point.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
#ifdef MEMORY_DEBUG
int NOMAD::Cache_File_Point::_cardinality     = 0;
int NOMAD::Cache_File_Point::_max_cardinality = 0;
#endif

/*---------------------------------------------------------------------*/
/*                             constructor 1                           */
/*---------------------------------------------------------------------*/
NOMAD::Cache_File_Point::Cache_File_Point ( void )
  : _n           ( 0    ) ,
    _m           ( 0    ) ,
    _m_def       ( 0    ) ,
    _eval_status ( 3    ) ,
    _coords      ( NULL ) ,
    _bbo_def     ( NULL ) ,
    _bbo_index   ( NULL )
{
#ifdef MEMORY_DEBUG
  ++NOMAD::Cache_File_Point::_cardinality;
  if ( NOMAD::Cache_File_Point::_cardinality >
       NOMAD::Cache_File_Point::_max_cardinality )
    ++NOMAD::Cache_File_Point::_max_cardinality;
#endif
}

/*---------------------------------------------------------------------*/
/*                             constructor 2                           */
/*---------------------------------------------------------------------*/
NOMAD::Cache_File_Point::Cache_File_Point ( const NOMAD::Eval_Point & x )
  : _n         ( x.size() ) ,
    _m         ( 0        ) ,
    _m_def     ( 0        ) ,
    _coords    ( NULL     ) ,
    _bbo_def   ( NULL     ) ,
    _bbo_index ( NULL     )
{
	int i;
	
	// eval_status:
	switch ( x.get_eval_status() ) {
		case NOMAD::EVAL_FAIL:
			_eval_status = 0;
			break;
		case NOMAD::EVAL_OK:
			_eval_status = 1;
			break;
		case NOMAD::EVAL_IN_PROGRESS:
			_eval_status = 2;
			break;
		case NOMAD::UNDEFINED_STATUS:
			_eval_status = 3;
			break;
		case NOMAD::EVAL_USER_REJECT:
			_eval_status = 3;
			break;
	}
	
	// inputs:
	if ( _n > 0 ) {
		_coords = new double [_n];
		for ( i = 0 ; i < _n ; ++i )
			_coords[i] = x[i].value();
	}
	else
		_n = 0;
	
	// outputs:
	const NOMAD::Point & bbo = x.get_bb_outputs();
	_m = bbo.size();
	if ( _m > 0 ) {
		
		std::vector<double> vd;
		std::vector<int>    vi;
		
		for ( i = 0 ; i < _m ; ++i )
			if ( bbo[i].is_defined() ) {
				vd.push_back ( bbo[i].value() );
				vi.push_back ( i              );
			}
		
		_m_def = static_cast<int> ( vd.size() );
		if ( _m_def > 0 ) {
			_bbo_def   = new double [_m_def];
			_bbo_index = new int    [_m_def];
			for ( i = 0 ; i < _m_def ; ++i ) {
				_bbo_def  [i] = vd[i];
				_bbo_index[i] = vi[i];
			}
		}
	}
	else
		_m = 0;
	
#ifdef MEMORY_DEBUG
	++NOMAD::Cache_File_Point::_cardinality;
	if ( NOMAD::Cache_File_Point::_cardinality >
		NOMAD::Cache_File_Point::_max_cardinality )
		++NOMAD::Cache_File_Point::_max_cardinality;
#endif
}

/*---------------------------------------------------------------------*/
/*                              destructor                             */
/*---------------------------------------------------------------------*/
NOMAD::Cache_File_Point::~Cache_File_Point ( void ) 
{
  delete [] _coords;
  delete [] _bbo_def;
  delete [] _bbo_index;
#ifdef MEMORY_DEBUG
  --NOMAD::Cache_File_Point::_cardinality;
#endif
}

/*---------------------------------------------------------------------*/
/*                             get_coord(i)                            */
/*---------------------------------------------------------------------*/
double NOMAD::Cache_File_Point::get_coord ( int i ) const
{
  if ( !_coords || i < 0 || i >= _n )
    throw NOMAD::Exception ( "Cache_File_Point.cpp" , __LINE__ ,
			     "bad access in Cache_File_Point::get_coord()" );
  return _coords[i];
}

/*---------------------------------------------------------------------*/
/*                               get_bb_outputs                        */
/*---------------------------------------------------------------------*/
const NOMAD::Point NOMAD::Cache_File_Point::get_bb_outputs ( void  ) const
{
  NOMAD::Point bbo ( _m );
  for ( int i = 0 ; i < _m_def ; ++i )
    bbo [ _bbo_index[i] ] = _bbo_def[i];
  return bbo;
}

/*---------------------------------------------------------*/
/*                       reset (private)                   */
/*---------------------------------------------------------*/
void NOMAD::Cache_File_Point::reset ( void )
{
  _n = _m = _m_def = 0;
  _eval_status = 3;

  delete [] _coords;
  delete [] _bbo_def;
  delete [] _bbo_index;

  _coords    = NULL;
  _bbo_def   = NULL;
  _bbo_index = NULL;
}

/*---------------------------------------------------------*/
/*                   write in binary file                  */
/*---------------------------------------------------------*/
bool NOMAD::Cache_File_Point::write ( std::ofstream & fout ) const
{
  // do nothing if no point:
  if ( _n <= 0 )
    return true;

  // 1. _eval_status:
  fout.write ( (char *) &_eval_status , sizeof(_eval_status)  );

  // 2. _n:
  fout.write ( (char *) &_n           , sizeof(_n)            );

  // 3. _m:
  fout.write ( (char *) &_m           , sizeof(_m)            );

  // 4. _m_def:
  fout.write ( (char *) &_m_def       , sizeof(_m_def)        );

  // 5. _coords:
  fout.write ( (char *)  _coords      , _n*sizeof(double)     );

  if ( _m_def > 0 ) {

    // 6. _bbo_def:
    fout.write ( (char *)  _bbo_def   , _m_def*sizeof(double) );

    // 7. _bbo_index:
    fout.write ( (char *)  _bbo_index , _m_def*sizeof(int)    );
  }

  return !fout.fail();
}

/*---------------------------------------------------------*/
/*                    read in binary file                  */
/*---------------------------------------------------------*/
bool NOMAD::Cache_File_Point::read ( std::ifstream & fin )
{
  reset();

  // 1. _eval_status:
  fin.read ( (char *) &_eval_status , sizeof(_eval_status) );
  if ( fin.fail() || _eval_status > 3 )
    return false;

  // 2. _n:
  fin.read ( (char *) &_n , sizeof(_n) );
  if ( fin.fail() || _n <= 0 ) {
    _n = 0;
    return false;
  }

  // 3. _m:
  fin.read ( (char *) &_m , sizeof(_m) );
  if ( fin.fail() || _m < 0 ) {
    _n = _m = 0;
    return false;
  }

  // 4. _m_def:
  fin.read ( (char *) &_m_def , sizeof(_m_def) );
  if ( fin.fail() || _m_def < 0 ) {
    _m_def = _n = _m = 0;
    return false;
  }

  // 5. _coords:
  _coords = new double [_n];
  fin.read ( (char *) _coords , _n*sizeof(double) );
  if ( fin.fail() ) {
    reset();
    return false;
  }

  if ( _m_def > 0 ) {

    // 6. _bb_def:
    _bbo_def = new double [_m_def];
    fin.read ( (char *) _bbo_def , _m_def*sizeof(double) );
    if ( fin.fail() ) {
      reset();
      return false;
    }
  
    // 7. _bbo_index:
    _bbo_index = new int [_m_def];
    fin.read ( (char *) _bbo_index , _m_def*sizeof(int) );
    if ( fin.fail() ) {
      reset();
      return false;
    }
  }

  return true;
}

/*---------------------------------------------------------*/
/*                         display                         */
/*---------------------------------------------------------*/
void NOMAD::Cache_File_Point::display ( const NOMAD::Display & out ) const
{
  out << "n      : " << _n     << std::endl
      << "m      : " << _m     << std::endl
      << "m_def  : " << _m_def << std::endl;

  int i;
  if ( _n > 0 ) {
    out << "coords    : ( ";
    for ( i = 0 ; i < _n ; ++i )
      out << _coords[i] << " ";
    out << ")" << std::endl;
  }
  if ( _m_def > 0 ) {
    out << "bbo_def   : [ ";
    for ( i = 0 ; i < _m_def ; ++i )
      out << _bbo_def[i] << " ";
    out << "]" << std::endl
	<< "bbo_index : [ ";
    for ( i = 0 ; i < _m_def ; ++i )
      out << _bbo_index[i] << " ";
    out << "]" << std::endl;
  }
}
