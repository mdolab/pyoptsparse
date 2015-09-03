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
  \file   Pareto_Point.cpp
  \brief  Pareto point (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-09
  \see    Pareto_Point.hpp
*/
#include "Pareto_Point.hpp"

/*--------------------------------------------------------*/
/* comparison operator:                                   */
/*  . supposes that argument fp is a Pareto_Point         */
/*  . used for the insertion in a set (the Pareto front)  */
/*  . we compare f1(x) and f1(y)                          */
/*--------------------------------------------------------*/
bool NOMAD::Pareto_Point::operator <
( const NOMAD::Set_Element<NOMAD::Eval_Point> & fp ) const
{
  if ( this == &fp || get_element() == fp.get_element() )
    return false;
  
  int i1 = NOMAD::Multi_Obj_Evaluator::get_i1();

  return get_element()->get_bb_outputs()[i1].value() <
         fp.get_element()->get_bb_outputs()[i1].value();
}

/*---------------------------------------------------------------*/
/* dominance notion:                                             */
/*  . used for the comparison (dominance) of two Pareto points,  */
/*    before they are inserted into the Pareto front             */
/*---------------------------------------------------------------*/
bool NOMAD::Pareto_Point::dominates ( const NOMAD::Pareto_Point & pp ) const
{
  if ( this == &pp || get_element() == pp.get_element() )
    return false;

  int i1 = NOMAD::Multi_Obj_Evaluator::get_i1();
  int i2 = NOMAD::Multi_Obj_Evaluator::get_i2();

  // we compare F(x)=[f1(x),f2(x)] and F(y)=[f1(y),f2(y)]:
  double f1x  = get_element()->get_bb_outputs   ()[i1].value();
  double f2x  = get_element()->get_bb_outputs   ()[i2].value();
  double f1y  = pp.get_element()->get_bb_outputs()[i1].value();
  double f2y  = pp.get_element()->get_bb_outputs()[i2].value();
  
  if ( f1x < f1y )
    return f2x <= f2y;

  if ( f1x == f1y )
    return ( f2x < f2y );
  
  return false;
}

/*---------------------------------------------------------------*/
/*                           display                             */
/*---------------------------------------------------------------*/
void NOMAD::Pareto_Point::display ( const NOMAD::Display & out ) const
{
  const NOMAD::Point & bbo = get_element()->get_bb_outputs();
  int                  w   = 13;

  out << "x=( ";
  get_element()->NOMAD::Point::display ( out , " " , w , -1 );
  out << " ) F(x)=[ ";
  bbo.Point::display ( out , " " , w , -1 );
  out << " ] [ f1(x) f2(x) ]=[ "
      << std::setw(w) << bbo[NOMAD::Multi_Obj_Evaluator::get_i1()] << " "
      << std::setw(w) << bbo[NOMAD::Multi_Obj_Evaluator::get_i2()] << " ]";
}
