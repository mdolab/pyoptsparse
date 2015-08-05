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
  \file   L_Curve.cpp
  \brief  L_CURVE_TARGET stopping criterion (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-09
  \see    L_Curve.hpp
*/
#include "L_Curve.hpp"

/*-----------------------------------------------*/
/*          insertion of a pair bbe/f            */
/*-----------------------------------------------*/
void NOMAD::L_Curve::insert ( int bbe , const NOMAD::Double & f )
{
  if ( _f.empty() ) {   
    _f.push_back   ( f );
    _bbe.push_back (bbe);
  }
  else {
    size_t nm1 = _bbe.size()-1;
    if ( _bbe[nm1] == bbe )
      _f[nm1] = f;
    else {
      _f.push_back   ( f );
      _bbe.push_back (bbe);
    }
  }
}

/*---------------------------------------------------------------*/
/*         check the L_CURVE_TARGET stopping criterion           */
/*  returns true if it detects that the target won't be reached  */
/*  after bbe evaluations)                                       */
/*---------------------------------------------------------------*/
bool NOMAD::L_Curve::check_stop ( int bbe ) const
{
  // we check the p last successes and approximate the L-curve
  // with a line joining the extremities:
  const size_t p = 7;

  if ( _f.size() >= p ) {
	
    size_t n = _f.size();
    
    NOMAD::Double f2 = _f[n-1];
    if ( f2 <= _target )
      return false;
    
    size_t       nmp = n-p;
    int         bbe1 = _bbe [ nmp ];
    NOMAD::Double f1 = _f   [ nmp ];
    NOMAD::Double  a = ( f2 - f1 ) / ( bbe - bbe1 );
    NOMAD::Double  b = f1 - a * bbe1;
    int   bbe_target = static_cast<int> ( ceil ( ( ( _target - b ) / a ).value() ) );
    
    // test: if ( bbe_target > bbe+(bbe-bbe1) )
    return ( bbe_target > 2*bbe - bbe1 );
  }
  return false;
}
