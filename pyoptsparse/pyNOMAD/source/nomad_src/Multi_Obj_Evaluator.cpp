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
  \file   Multi_Obj_Evaluator.cpp
  \brief  NOMAD::Evaluator subclass for multiobjective optimization (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-09
  \see    Multi_Obj_Evaluator.hpp
*/
#include "Multi_Obj_Evaluator.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
int NOMAD::Multi_Obj_Evaluator::_i1 = -1;
int NOMAD::Multi_Obj_Evaluator::_i2 = -1;

/*---------------------------------------*/
/*  initialization of objective indexes  */
/*                (static)               */
/*---------------------------------------*/
void NOMAD::Multi_Obj_Evaluator::set_obj_indexes ( const std::list<int> & index_obj )
{
  if ( index_obj.size() != 2 )
    throw NOMAD::Exception ( "Multi_Obj_Evaluator.cpp" , __LINE__ ,
	  "Multi_Obj_Evaluator defined with a number of indexes different than two" );

  std::list<int>::const_iterator it = index_obj.begin();

  NOMAD::Multi_Obj_Evaluator::_i1 = *it;
  it++;
  NOMAD::Multi_Obj_Evaluator::_i2 = *it;
}

/*-------------------------------------------------------------*/
/*  computation of f, taking into account the two objectives,  */
/*  with weights _w1 and _w2                                   */
/*-------------------------------------------------------------*/
void NOMAD::Multi_Obj_Evaluator::compute_f ( NOMAD::Eval_Point & x ) const
{
  if ( _i1 < 0 || _i2 < 0 )
    throw NOMAD::Exception ( "Multi_Obj_Evaluator.cpp" , __LINE__ ,
	  "Multi_Obj_Evaluator::compute_f(): no objective indexes defined" );

  int obj_index [2];
  obj_index[0] = _i1;
  obj_index[1] = _i2;

  const NOMAD::Point & bbo = x.get_bb_outputs();

  // a reference is available:
  if ( _ref ) 
  {
      
	  NOMAD::multi_formulation_type mft = _p.get_multi_formulation();
	  
	  if ( mft == NOMAD::UNDEFINED_FORMULATION )
		  throw NOMAD::Exception ( "Multi_Obj_Evaluator.cpp" , __LINE__ ,
								  "Multi_Obj_Evaluator::compute_f(): no formulation type is defined" );
	  
	  // normalized formulation:
	  if ( mft == NOMAD::NORMALIZED || mft == NOMAD::DIST_LINF ) {
		  
		  // f1 - r1:
		  NOMAD::Double d = bbo[obj_index[0]] - (*_ref)[0];
		  
		  // f2 - r2:
		  NOMAD::Double f2mr2 = bbo[obj_index[1]] - (*_ref)[1];
		  
		  // we take the max:
		  if ( f2mr2 > d )
			  d = f2mr2;
		  
		  x.set_f ( d );
	  }
	  
	  // product formulation:
	  else if ( mft == NOMAD::PRODUCT ) {
		  
		  NOMAD::Double prod = 1.0 , ri , fi;
		  
		  for ( int i = 0 ; i < 2 ; ++i ) {
			  
			  ri = (*_ref)[i];
			  fi = bbo[obj_index[i]];
			  
			  if ( fi > ri ) {
				  prod = 0.0;
				  break;
			  }
			  prod = prod * (ri-fi).pow2();
		  }
		  
		  x.set_f ( -prod );
	  }
	  
	  // distance formulation:
	  else {
		  
		  NOMAD::Double d;  
		  NOMAD::Double r1mf1 = (*_ref)[0] - bbo[obj_index[0]];
		  NOMAD::Double r2mf2 = (*_ref)[1] - bbo[obj_index[1]];
		  
		  if ( r1mf1 >= 0.0 && r2mf2 >= 0.0 ) {
			  d = r1mf1.pow2();
			  NOMAD::Double tmp = r2mf2.pow2();
			  if ( tmp < d )
				  d = tmp;
			  d = -d;
		  }
		  else if ( r1mf1 <= 0.0 && r2mf2 <= 0.0 ) {
			  
			  // with L2 norm:
			  if ( mft == NOMAD::DIST_L2 )
				  d = r1mf1.pow2() + r2mf2.pow2();
			  
			  // with L1 norm:
			  else
				  d = (r1mf1.abs() + r2mf2.abs()).pow2();
			  
			  // Linf norm: treated as NORMALIZED
		  }
		  else if ( r1mf1 > 0.0 )
			  d = r2mf2.pow2();
		  else
			  d = r1mf1.pow2();
		  
		  x.set_f ( d );
	  }
  }
    
  // no reference is available (use weights):
  else
    x.set_f ( _w1 * bbo[obj_index[0]] + _w2 * bbo[obj_index[1]] );
}
