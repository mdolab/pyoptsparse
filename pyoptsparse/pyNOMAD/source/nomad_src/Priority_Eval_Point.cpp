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
  \file   Priority_Eval_Point.cpp
  \brief  Evaluation point with a priority (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-22
  \see    Priority_Eval_Point.hpp
*/
#include "Priority_Eval_Point.hpp"

bool NOMAD::Priority_Eval_Point::_lexicographic_order=false;

/*------------------------------------------------*/
/*                comparison operator             */
/*------------------------------------------------*/
/*  . x1.dominates(x2) returns true if x1 should  */
/*    be evaluated before x2                      */
/*  . x is a Priority_Eval_Point                  */
/*------------------------------------------------*/
bool NOMAD::Priority_Eval_Point::dominates
( const NOMAD::Set_Element<NOMAD::Eval_Point> & x ) const
{
  if ( this == &x )
    return false;
     const NOMAD::Eval_Point * x1 = get_element();
  const NOMAD::Eval_Point * x2 = x.get_element();
  
  // criterion 0: lexicographic order
	if (_lexicographic_order)
		return NOMAD::Point(*x1) < NOMAD::Point(*x2);

  
  // criterion 1: user criterion:
  // ------------
  const NOMAD::Double uep1 = x1->get_user_eval_priority();
  if ( uep1.is_defined() ) 
  {
    const NOMAD::Double uep2 = x2->get_user_eval_priority();
    if ( uep2.is_defined() )
	{
      if ( uep1 > uep2 )
	return true;
      if ( uep2 > uep1 )
	return false;
    }
  }

  // specific Priority_Eval_Point elements of comparison:
  NOMAD::Double x_f_sgte;
  NOMAD::Double x_h_sgte;
  NOMAD::Double x_f_model;
  NOMAD::Double x_h_model;
  NOMAD::Double x_angle_success_dir;
  NOMAD::Double x_angle_simplex_grad;

  x.get_priority_criteria ( x_f_sgte             ,
			    x_h_sgte             ,
			    x_f_model            ,
			    x_h_model            ,
			    x_angle_success_dir  ,
			    x_angle_simplex_grad   );

  // criterion 2: give priority to already evaluated cache points:
  // ------------
  if ( x1->is_in_cache() && !x2->is_in_cache() )
    return true;
  if ( x2->is_in_cache() && !x1->is_in_cache() )
    return false;

  // criterion 3: give priority to already evaluated points
  // ------------ that are eval_ok:
  if ( x1->is_eval_ok() && !x2->is_eval_ok() )
    return true;
  if ( x2->is_eval_ok() && !x1->is_eval_ok() )
    return false;

  // criterion 4: true f and h values:
  // -----------
  int flag = compare_hf_values ( x1->get_h() ,
				 x1->get_f() ,
				 x2->get_h() ,
				 x2->get_f()   );
  if ( flag )
    return ( flag > 0 );

  // criterion 5: surrogate f and h values:
  // ------------
  flag = compare_hf_values ( _h_sgte , _f_sgte , x_h_sgte , x_f_sgte );
  if ( flag )
      return ( flag > 0 );  

  // criterion 6: model f and h values:
  // ------------
  flag = compare_hf_values ( _h_model , _f_model , x_h_model , x_f_model );
  if ( flag )
    return ( flag > 0 );
	
	

  // criterion 7: check the angle with the last successful direction:
  // ------------
  if ( _angle_success_dir.is_defined() && x_angle_success_dir.is_defined() ) 
  {
    if ( _angle_success_dir < x_angle_success_dir )
      return true;
    if ( x_angle_success_dir < _angle_success_dir )
      return false;
  }

	
  // criterion 8: take the point with the best h value:
  // ------------
  flag = compare_h_values ( x1->get_h() , x2->get_h() );
  if ( flag )
    return ( flag > 0 );

  flag = compare_h_values ( _h_sgte , x_h_sgte );
  if ( flag )
    return ( flag > 0 );

  flag = compare_h_values ( _h_model , x_h_model );
  if ( flag )
    return ( flag > 0 );

  // criterion 9: random criterion for randomly generated directions:
  // -------------
  const NOMAD::Double rep1 = x1->get_rand_eval_priority();
  if ( rep1.is_defined() ) 
  {
    const NOMAD::Double rep2 = x2->get_rand_eval_priority();
    if ( rep2.is_defined() )
	{
      if ( rep1 < rep2 )
	return true;
      if ( rep2 < rep1 )
	return false;
    }
  }

	// criterion 10: compare the tags:
  // -------------
  return x1->get_tag() < x2->get_tag();
	
}

/*-----------------------------------------------*/
/*       compare the h values of two points      */
/*-----------------------------------------------*/
/*  . return h(x1) < h(x2) with the format:      */
/*                      1 : x1 better than x2    */
/*                     -1 : x2 better than x1    */
/*                      0 : undetermined         */
/*  . private method                             */
/*-----------------------------------------------*/
int NOMAD::Priority_Eval_Point::compare_h_values ( const NOMAD::Double & hx1 ,
						   const NOMAD::Double & hx2   ) const
{
  if ( hx1.is_defined() && hx2.is_defined() )
  {
    if ( hx1 < hx2 )
      return 1;
    if ( hx2 < hx1 )
      return -1;
  }
  return 0;
}

/*-----------------------------------------------*/
/*    compare the h and f values of two points   */
/*-----------------------------------------------*/
/*  . return ( h(x1),f(x1) ) < ( h(x2),f(x2) )   */
/*    with the following format:                 */
/*                      1 : x1 better than x2    */
/*                     -1 : x2 better than x1    */
/*                      0 : undetermined         */
/*  . private method                             */
/*-----------------------------------------------*/
int NOMAD::Priority_Eval_Point::compare_hf_values ( const NOMAD::Double & hx1 ,
						    const NOMAD::Double & fx1 ,
						    const NOMAD::Double & hx2 ,
						    const NOMAD::Double & fx2   ) const
{
	if ( fx1.is_defined() && fx2.is_defined() ) 
	{
		
		if ( hx1.is_defined() && hx2.is_defined() )
		{
			// x1 is feasible:
			if ( hx1 <= _h_min ) 
			{
				// both points are feasible:
				if ( hx2 <= _h_min  ) 
				{
					if ( fx1 < fx2 )
						return 1;
					if ( fx2 < fx1 )
						return -1;
				}
				
				// x1 feasible and x2 infeasible:
				else
					return 1;
			}
			
			// x1 is infeasible:
			else 
			{
				// x2 is feasible:
				if ( hx2 <= _h_min  )
					return -1;
				
				// both points are infeasible:
				if ( ( hx1  < hx2 && fx1  < fx2 ) ||
					( hx1 == hx2 && fx1  < fx2 ) ||
					( hx1  < hx2 && fx1 == fx2 )    )
					return 1;
				
				if ( ( hx2  < hx1 && fx2  < fx1 ) ||
					( hx2 == hx1 && fx2  < fx1 ) ||
					( hx2  < hx1 && fx2 == fx1 )    )
					return -1; 
			}
		}
		
		// we only have f values:
		else 
		{
			if ( fx1 < fx2 )
				return 1;
			if ( fx2 < fx1 )
				return -1;
		}
	}
	return 0;
}
