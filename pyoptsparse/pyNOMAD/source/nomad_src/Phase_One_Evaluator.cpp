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
 \file   Phase_One_Evaluator.cpp
 \brief  NOMAD::Evaluator subclass for the phase one (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-09
 \see    Phase_One_Evaluator.hpp
 */
#include "Phase_One_Evaluator.hpp"

/*------------------------------------------------------------------*/
/*         compute f(x) from the blackbox outputs of a point        */
/*              (special objective for MADS phase 1)                */
/*------------------------------------------------------------------*/
void NOMAD::Phase_One_Evaluator::compute_f ( NOMAD::Eval_Point & x ) const
{
	if ( x.get_bb_outputs().size() != _p.get_bb_nb_outputs() )
	{
		std::ostringstream err;
		err << "Phase_One_Evaluator::compute_f(x): "
		<< "x has a wrong number of blackbox outputs ("
		<< x.get_bb_outputs().size() <<  " != " << _p.get_bb_nb_outputs() << ")";
		throw NOMAD::Exception ( "Phase_One_Evaluator.cpp" , __LINE__ , err.str() );
	}
	
	// objective value for MADS phase 1: the squared sum of all EB constraint violations
	// (each EB constraint has been previously transformed into OBJ values):
	const std::list<int>               & index_obj = _p.get_index_obj();
	const std::list<int>::const_iterator end       = index_obj.end();
	const NOMAD::Point                 & bbo       = x.get_bb_outputs();
	NOMAD::Double                        h_min     = _p.get_h_min();
	NOMAD::Double                        sum       = 0.0;
	NOMAD::Double                        v;
	
	for ( std::list<int>::const_iterator it = index_obj.begin() ; it != end ; ++it )
	{
		v = bbo[*it];
		if ( v > h_min )
			sum += v.pow2();
	}
	
	x.set_f ( sum );
}
