/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search - version (3.7.0.beta)*/
/*                                                                                     */
/*  Copyright (C) 2001-2015 Mark Abramson        - the Boeing Company, Seattle        */
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
  \file   Quad_Model_Evaluator.cpp
  \brief  NOMAD::Evaluator subclass for quadratic model optimization (implementation)
  \author Sebastien Le Digabel
  \date   2010-08-31
  \see    Quad_Mopel_Evaluator.hpp
*/
#include "Quad_Model_Evaluator.hpp"

/*-----------------------------*/
/*         constructor         */
/*-----------------------------*/
NOMAD::Quad_Model_Evaluator::Quad_Model_Evaluator
( const NOMAD::Parameters & p     ,
  const NOMAD::Quad_Model & model   )
  :  
	_n					( model.get_n()         ) ,
    _nm1				( _n-1                  ) , 
    _m					( p.get_bb_nb_outputs() ) ,
    _x					( NULL                  ) ,
    _alpha				( NULL                  ) ,
    _model_ready		( model.check()         )
{
	if ( _model_ready )
	{
		
		int i , j , k , k2 , nalpha = (_n+1)*(_n+2)/2 , nfree = model.get_nfree();
		NOMAD::Point ** model_alpha = model.get_alpha();
		
		_x     = new double   [_n];
		_alpha = new double * [_m];
		
		for ( int io = 0 ; io < _m ; ++io )
		{
			_alpha[io] = NULL;
			if ( model_alpha[io] )
			{				
				_alpha[io] = new double[nalpha];
				_alpha[io][0] = (*model_alpha[io])[0].value();
				
				for ( i = 1 ; i < nalpha ; ++i ) 
					_alpha[io][i] = 0.0;
				
				k = 0;
				
				for ( i = 0 ; i < _n ; ++i )
				{
					if ( !model.variable_is_fixed(i) ) 
					{
						++k;
						_alpha[io][i+1   ] = (*model_alpha[io])[k      ].value();
						_alpha[io][i+1+_n] = (*model_alpha[io])[k+nfree].value();
					}
				}
				
				k += nfree;
				k2 = 2*_n;
				
				for ( i = 0 ; i < _nm1 ; ++i )
				{
					if ( !model.variable_is_fixed(i) ) 
					{
						for ( j = i+1 ; j < _n ; ++j )
						{
							++k2;
							if ( !model.variable_is_fixed(j) )
								_alpha[io][k2] = (*model_alpha[io])[++k].value();
						}
					}
					else
						for ( j = i+1 ; j < _n ; ++j )
							++k2;
				}
			}
		}
	}
}

/*-----------------------------*/
/*          destructor         */
/*-----------------------------*/
NOMAD::Quad_Model_Evaluator::~Quad_Model_Evaluator ( void )
{
  if ( _model_ready ) 
  {
    for ( int i = 0 ; i < _m ; ++i )
      if ( _alpha[i] )
	delete [] _alpha[i];
    delete [] _alpha;
    delete [] _x;
  }
}

/*------------------------------------------------------------------------*/
/*      evaluate the blackboxes quad model at a given trial point         */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* x is in [-1000;1000] and must be scaled to [-1;1] for the evalutation  */
/*                                                                        */
/*------------------------------------------------------------------------*/
bool NOMAD::Quad_Model_Evaluator::eval_x
( NOMAD::Eval_Point   & x          ,
  const NOMAD::Double & h_max      ,
  bool                & count_eval   ) const
{
	
	count_eval=false;
	if ( !_model_ready )
		return false;
	
	int    i , j , k;
	double z , * alpha , * p;
	
	for ( i = 0 ; i < _n ; ++i )
		_x[i] = x[i].value() / 1000.0;
	
	for ( int oi = 0 ; oi < _m ; ++oi )
	{
		
		alpha = _alpha[oi];
		
		if ( alpha )
		{
			
			z = alpha[0];
			p = _x;
			
			for ( k = 1 ; k <= _n ; ++k , ++p )
				z += *p * ( alpha[k] + 0.5 * alpha[k+_n] * *p );
			
			k += _n-1;
			
			for ( i = 0 ; i < _nm1 ; ++i )
				for ( j = i+1 ; j < _n ; ++j )
					z += alpha[++k] * _x[i] * _x[j];
			
			x.set_bb_output ( oi , z );
		}
		
		else
			x.set_bb_output ( oi , 0.0 );
	}
	
	count_eval = true;
	return true;
}


/*------------------------------------------------------------------------*/
/* evaluate the gradient of a blackbox quad model at a given trial point  */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* x is in [-1000;1000] and must be scaled to [-1;1] for the evalutation  */
/*                                                                        */
/*------------------------------------------------------------------------*/
bool NOMAD::Quad_Model_Evaluator::evalGrad_x
(const NOMAD::Point   & x   ,
 NOMAD::Point   & g         ,
 const int & output_index   ,
 bool                & count_eval   ) const
{
	
	if ( !_model_ready )
		return false;
	
	double * alpha ;
	
	for (int  i = 0 ; i < _n ; ++i )
		_x[i] = x[i].value() / 1000.0;
	
	alpha = _alpha[output_index];
	
	if ( !alpha )
		return false;
		
	int i , j , k ;
		
	for ( k = 1 ; k <= _n ; ++k )
		g[k-1] = alpha[k] + alpha[k+_n] * _x[k-1];
	
	
	k += _n-1;
	
	for ( i = 0 ; i < _nm1 ; ++i )
		for ( j = i+1 ; j < _n ; ++j )
			g[i] += alpha[++k] * _x[j];
	
	count_eval = true;
	return true;
}
