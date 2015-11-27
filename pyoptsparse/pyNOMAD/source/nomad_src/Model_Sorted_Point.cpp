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
  \file   Model_Sorted_Point.cpp
  \brief  Interpolation point with distance to model center (implementation)
  \author Sebastien Le Digabel
  \date   2010-11-15
  \see    Model_Sorted_Point.hpp
*/
#include "Model_Sorted_Point.hpp"

/*---------------------------------------------------------*/
/*                        constructor                      */
/*---------------------------------------------------------*/
NOMAD::Model_Sorted_Point::Model_Sorted_Point
( NOMAD::Point * x , const NOMAD::Point & center ) : _x(x)
{
  int i , n = center.size();
  if ( x && x->size() == n ) {
    _dist = 0.0;
    for ( i = 0 ; i < n ; ++i )
      if ( (*x)[i].is_defined() && center[i].is_defined() ) {
	_dist += ( (*x)[i] - center[i] ).pow2();
      }
      else {
	_dist.clear();
	break;
      }
  }
}

/*---------------------------------------------------------*/
/*                   affectation operator                  */
/*---------------------------------------------------------*/
NOMAD::Model_Sorted_Point & NOMAD::Model_Sorted_Point::operator =
( const NOMAD::Model_Sorted_Point & x )
{
  _x    = x._x;
  _dist = x._dist;
  return *this;
}

/*---------------------------------------------------------*/
/*                    comparison operator                  */
/*---------------------------------------------------------*/
bool NOMAD::Model_Sorted_Point::operator <
( const Model_Sorted_Point & x ) const
{
  if ( _dist.is_defined() ) {
    if ( !x._dist.is_defined() )
      return true;
    return _dist < x._dist;
  }
  return false;
}
