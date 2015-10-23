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
  \file   Direction.cpp
  \brief  Polling direction (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-05
  \see    Direction.hpp
*/
#include "Direction.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
#ifdef MEMORY_DEBUG
int NOMAD::Direction::_cardinality     = 0;
int NOMAD::Direction::_max_cardinality = 0;
#endif

/*---------------------------------------------------------*/
/*                       constructor 1                     */
/*---------------------------------------------------------*/
NOMAD::Direction::Direction ( void )
  : NOMAD::Point  (                            ) ,
    _type         ( NOMAD::UNDEFINED_DIRECTION ) ,
_index        ( -1                         ),
_dir_group_index (-1)
{
#ifdef MEMORY_DEBUG
  ++NOMAD::Direction::_cardinality;
  if ( NOMAD::Direction::_cardinality > NOMAD::Direction::_max_cardinality )
    ++NOMAD::Direction::_max_cardinality;
#endif
}

/*---------------------------------------------------------*/
/*                       constructor 2                     */
/*---------------------------------------------------------*/
NOMAD::Direction::Direction ( int                     n    ,
			      const NOMAD::Double   & v    ,
			      NOMAD::direction_type   type ,		 
					int dir_group_index		 )
  : NOMAD::Point  ( n , v ) ,
    _type         ( type  ) ,
    _index        ( -1   ) ,
    _dir_group_index (dir_group_index)
{
#ifdef MEMORY_DEBUG
  ++NOMAD::Direction::_cardinality;
  if ( NOMAD::Direction::_cardinality > NOMAD::Direction::_max_cardinality )
    ++NOMAD::Direction::_max_cardinality;
#endif
}

/*---------------------------------------------------------*/
/*                       constructor 2b                    */
/*---------------------------------------------------------*/
NOMAD::Direction::Direction ( int                     n    ,
							 const NOMAD::Double   & v    ,
							 NOMAD::direction_type   type   )
: NOMAD::Point  ( n , v ) ,
_type         ( type  ) ,
_index        ( -1    ),
_dir_group_index (-1)
{
#ifdef MEMORY_DEBUG
  ++NOMAD::Direction::_cardinality;
  if ( NOMAD::Direction::_cardinality > NOMAD::Direction::_max_cardinality )
    ++NOMAD::Direction::_max_cardinality;
#endif
}

/*---------------------------------------------------------*/
/*                       constructor 3                     */
/*---------------------------------------------------------*/
NOMAD::Direction::Direction ( const NOMAD::Point    & x    , 
			      NOMAD::direction_type   type   )
  : NOMAD::Point  ( x    ) ,
    _type         ( type ) ,
_index        ( -1   ),
_dir_group_index (-1)
{
#ifdef MEMORY_DEBUG
  ++NOMAD::Direction::_cardinality;
  if ( NOMAD::Direction::_cardinality > NOMAD::Direction::_max_cardinality )
    ++NOMAD::Direction::_max_cardinality;
#endif
}

/*---------------------------------------------------------*/
/*                      copy constructor                   */
/*---------------------------------------------------------*/
NOMAD::Direction::Direction ( const Direction & d )
  : NOMAD::Point  ( d        ) ,
    _type         ( d._type  ) ,
    _index        ( d._index ),
_dir_group_index (d._dir_group_index)
{
#ifdef MEMORY_DEBUG
  ++NOMAD::Direction::_cardinality;
  if ( NOMAD::Direction::_cardinality > NOMAD::Direction::_max_cardinality )
    ++NOMAD::Direction::_max_cardinality;
#endif
}

/*---------------------------------------------------------*/
/*                         destructor                      */
/*---------------------------------------------------------*/
NOMAD::Direction::~Direction ( void )
{
#ifdef MEMORY_DEBUG
  --NOMAD::Direction::_cardinality;
#endif
}

/*---------------------------------------------------------*/
/*                  affectation operator                   */
/*---------------------------------------------------------*/
NOMAD::Direction & NOMAD::Direction::operator = ( const NOMAD::Direction & d )
{
  if ( this == &d )
    return *this;

  NOMAD::Point::operator = ( d );
  
  _type  = d._type;
  _index = d._index;

  return *this;
}

/*---------------------------------------------------------*/
/*                            clear                        */
/*---------------------------------------------------------*/
void NOMAD::Direction::clear ( void )
{
  NOMAD::Point::clear();
  _type  = NOMAD::UNDEFINED_DIRECTION;
  _index = -1;
}

/*-----------------------------------------------------------*/
/*                           negation                        */
/*-----------------------------------------------------------*/
const NOMAD::Direction NOMAD::Direction::operator - ( void ) const
{
  return NOMAD::Direction ( this->NOMAD::Point::operator-() , _type );
}

/*---------------------------------------------------------*/
/*                          display                        */
/*---------------------------------------------------------*/
void NOMAD::Direction::display ( const NOMAD::Display & out ,
				 const std::string    & sep ,
				 int                    w   ,
				 int                    lim   ) const
{
  if ( is_defined() ) {
    out << "( ";
    NOMAD::Point::display ( out , sep , w , lim );
    out << " ) " << _type;
  }
  else
    out << "undefined";
}
