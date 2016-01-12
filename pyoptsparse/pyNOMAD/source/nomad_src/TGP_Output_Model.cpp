/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search - version 3.7.0.beta        */
/*                                                                                     */
/*  Copyright (C) 2001-2014  Mark Abramson        - the Boeing Company, Seattle        */
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
  \file   TGP_Output_Model.cpp
  \brief  TGP (Bayesian treed Gaussian process) model for one output (implementation)
  \author Sebastien Le Digabel
  \date   2011-02-07
  \see    TGP_Output_Model.hpp
*/

#ifndef USE_TGP

int TGP_OUTPUT_MODEL_DUMMY; // avoids that TGP_Output_Model.o has no symbols with ranlib

#else

#include "TGP_Output_Model.hpp"

/*---------------------------------------------------------*/
/*  NOMAD-TGP callback function (called regularly by TGP)  */
/*---------------------------------------------------------*/
// SLD -- 2014-09-04
// void NOMAD::TGP_callback ( bool & TGP_interrupt )
// {
//   if ( NOMAD::TGP_Output_Model::get_force_quit() )
//     TGP_interrupt = true;
// }

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
double NOMAD::TGP_Output_Model::_ditemps[7] =
  { 1.0 , 0.0 , 0.0 , 1.0 , 1.0 , 0.0 , 1.0 };

bool NOMAD::TGP_Output_Model::_force_quit = false;

/*-----------------------------------*/
/*            constructor            */
/*-----------------------------------*/
NOMAD::TGP_Output_Model::TGP_Output_Model
( const std::list<const NOMAD::Eval_Point *> & X_pts     ,
  int                                          bbo_index ,
  int                                          seed      ,
  const NOMAD::Display                       & out         )
  : _out         ( out            ) ,
    _p           ( X_pts.size()   ) ,
    _Z           ( new double[_p] ) ,
    _Z_is_scaled ( false          ) ,
    _is_binary   ( true           ) ,
    _bin_values  ( 2              ) ,
    _is_fixed    ( false          ) ,
    _tgp_state   ( NULL           ) ,
    _tgp_model   ( NULL           ) ,
    _tgp_its     ( NULL           )
{
  NOMAD::TGP_Output_Model::_force_quit = false;

  _Z_scaling[0] = _Z_scaling[1] = 0.0;

  std::list<const NOMAD::Eval_Point *>::const_iterator it , end = X_pts.end();

  NOMAD::Double tmp , Zmin , Zmax , Zsum = 0.0;
  int           j = 0;

  for ( it = X_pts.begin() ; it != end ; ++it ) {

    // the output value:
    tmp = (*it)->get_bb_outputs()[bbo_index];
    _Z[j++] = tmp.value();
    
    // Z scaling parameters (1/2):
    Zsum += tmp;
    if ( !Zmin.is_defined() || tmp < Zmin )
      Zmin = tmp;
    if ( !Zmax.is_defined() || tmp > Zmax )
      Zmax = tmp;

    // check if the output is binary:
    if ( _is_binary ) {
      if ( !_bin_values[0].is_defined() )
	_bin_values[0] = tmp;
      else if ( !_bin_values[1].is_defined() && tmp != _bin_values[0] )
	_bin_values[1] = tmp;
      else {
	if ( tmp != _bin_values[0] && tmp != _bin_values[1] )
	  _is_binary = false;
      }
    }
  }

  // Z scaling parameters (1/2):
  _Z_scaling[0] = (Zmax - Zmin).value();

  // the output is fixed:
  if ( _Z_scaling[0] == 0.0 )
    _is_fixed = true;

  else {

    _Z_scaling[1] = (Zsum/_p).value() / _Z_scaling[0];

    if ( !_is_binary )
      _bin_values = NOMAD::Point(2);

    // RNG (random number generator):
    int state[] = { 896 , 265 , 371 };

    // if seed==0, the model will be the same as the one constructed in R,
    // with values taken from the R tgp functions for:
    //   set.seed(0)
    //   state <- sample(seq(0, 999), 3)

    // otherwise use rand() to get three different integers in [0;999]:
    if ( seed != 0 ) {
      state[0] = rand()%1000;
      while ( state[1] == state[0] )
	state[1] = rand()%1000;
      while ( state[2] == state[0] || state[2] == state[1] )
	state[2] = rand()%1000;
    }
    _tgp_state = newRNGstate ( three2lstate ( state ) );

    // importance tempering:
    _tgp_its = new Temper ( NOMAD::TGP_Output_Model::_ditemps );
  }
}

/*--------------------------------------------*/
/*                  destructor                */
/*--------------------------------------------*/
NOMAD::TGP_Output_Model::~TGP_Output_Model ( void )
{
  if ( _Z )
    delete [] _Z;
  if ( _tgp_model )
    delete _tgp_model;
  if ( _tgp_its )
    delete _tgp_its;
  if ( _tgp_state )
    deleteRNGstate ( _tgp_state );
}

/*--------------------------------------------*/
/*              compute the model             */
/*--------------------------------------------*/
void NOMAD::TGP_Output_Model::compute ( double ** X              ,
					double ** XX             ,
					double ** Xsplit         ,
					int       n              ,
					int       n_XX           ,
					int       nsplit         ,
					Params  * tgp_params     ,
					double ** tgp_rect       ,
					int     * tgp_BTE        ,
					bool      tgp_linburn    ,
					bool      tgp_verb       ,
					double  * ZZ             ,   // OUT
					double  * Ds2x           ,   // OUT, may be NULL
					int     * improv           ) // OUT, may be NULL
{
  bool compute_Ds2x   = ( Ds2x   != NULL );
  bool compute_improv = ( improv != NULL );

  // the output is fixed:
  if ( _is_fixed ) {
    for ( int j = 0 ; j < n_XX ; ++j ) {
      ZZ[j] = _bin_values[0].value();
      if ( compute_Ds2x )
	Ds2x[j] = 0.0;
    }
    return;
  }

  // scale Z:
  scale_Z();

  // construct the TGP model:
  _tgp_model = new Model ( tgp_params ,
			   n          ,
			   tgp_rect   ,
			   0          , // Id=0
			   false      , // trace=false
			   _tgp_state   );

  _tgp_model->Init ( X        ,
		     _p       ,
		     n        ,
		     _Z       ,
		     _tgp_its ,
		     NULL     ,    // dtree=NULL
		     0        ,    // ncol=0
		     NULL       ); // dhier=NULL

  // set the NOMAD-TGP callback function:
  // _tgp_model->set_callback ( NOMAD::TGP_callback ); // SLD -- 2014-09-04

  // TGP verbim (display):
#ifdef TGP_DEBUG
  _tgp_model->Outfile ( stdout , 1 ); // set 10 for full display
#else
  if ( tgp_verb )
    _tgp_model->Outfile ( stdout , 1 );
  else
    _tgp_model->Outfile ( NULL , 0 );
#endif

  // set the splitting locations (Xsplit):
  _tgp_model->set_Xsplit ( Xsplit , nsplit , n );

  // linear model initialization rounds -B thru 1:
  if ( tgp_linburn )
    _tgp_model->Linburn ( tgp_BTE[0] , _tgp_state );

  // do model rounds 1 thru B (burn in):
  _tgp_model->Burnin ( tgp_BTE[0] , _tgp_state );

  // do the MCMC rounds B,...,T:
  Preds * tgp_preds = new_preds ( XX                      ,
 				  n_XX                    ,
 				  0                       ,
 				  n                       ,
 				  tgp_rect                ,
 				  tgp_BTE[1]-tgp_BTE[0]   ,
 				  false                   , // pred_n
 				  true                    , // krige
 				  _tgp_its->IT_ST_or_IS() ,
				  compute_Ds2x            ,
 				  compute_improv          ,
 				  false                   , // sens
 				  tgp_BTE[2]                );

  _tgp_model->Sample ( tgp_preds , tgp_BTE[1]-tgp_BTE[0] , _tgp_state );

  // if importance tempering, then update the pseudo-prior
  // based on the observation counts:
  if ( _tgp_its->Numit() > 1 ) 
    _tgp_its->UpdatePrior ( _tgp_model->update_tprobs() , _tgp_its->Numit() );

  // copy back the itemps:
  _tgp_model->DupItemps ( _tgp_its );

  // kriging mean:
  wmean_of_columns ( ZZ             ,
		     tgp_preds->ZZm ,
		     tgp_preds->R   ,
		     n_XX           ,
		     (_tgp_its->IT_ST_or_IS()) ? tgp_preds->w : NULL );

  int i;

  // expected reduction in predictive variance (Ds2x):
  if ( compute_Ds2x ) {
    for ( i = 0 ; i < n_XX ; ++i )
      Ds2x[i] = 0.0;
    if ( tgp_preds->Ds2x )
      wmean_of_columns ( Ds2x            ,
			 tgp_preds->Ds2x ,
			 tgp_preds->R    ,
			 tgp_preds->nn   ,
			 (_tgp_its->IT_ST_or_IS()) ? tgp_preds->w : NULL );
  }

  // expected improvement of objective (improv):
  if ( compute_improv ) {

    // double * improvec = new double [n_XX];
    // wmean_of_columns ( improvec ,
    // 		       tgp_preds->improv,
    // 		       tgp_preds->R,
    // 		       tgp_preds->nn,
    // 		       (_tgp_its->IT_ST_or_IS()) ? tgp_preds->w : NULL );
    // for ( i = 0 ; i < n_XX ; ++i )
    //   _out << "IMPROVEC[" << i<< "] = " << improvec[i] << std::endl;   
    // delete [] improvec;

    int *ir = (int*) GetImprovRank ( tgp_preds->R      ,
 				     tgp_preds->nn     ,
 				     tgp_preds->improv ,
 				     true              , // improv=true
 				     tgp_preds->nn     , // numirank = n_XX
 				     (_tgp_its->IT_ST_or_IS()) ? tgp_preds->w : NULL );
    
    for ( i = 0 ; i < n_XX ; ++i ) {
      improv[i] = ir[i];
      if ( improv[i] == 0 )
	improv[i] = n_XX;
    }

    free ( ir );

    // for ( i = 0 ; i < n_XX ; ++i )
    //  _out << "RANK[" << i<< "] = " << improv[i] << std::endl;
  }

  delete_preds ( tgp_preds );

#ifdef TGP_DEBUG
  _tgp_model->Outfile ( NULL , 0 );
#endif

  // unscale Z and ZZ:
  unscale_Z ( ZZ , n_XX );
  unscale_Z();

  // treat binary output:
  if ( _is_binary )
    treat_binary_output ( ZZ , n_XX );

  // disable TGP display:
  _tgp_model->Outfile ( NULL , 0 );
}

/*--------------------------------------------*/
/*           prediction at one point          */
/*--------------------------------------------*/
bool NOMAD::TGP_Output_Model::predict ( double  * XX       ,
					int       n        ,
					double  & ZZ       ,
					double ** tgp_rect   ) const
{
  if ( _is_fixed ) {
    ZZ = _bin_values[0].value();
    return true;
  }

  // do the MCMC rounds B,...,T:
  Preds * tgp_preds = new_preds ( &XX                     ,
				  1                       ,
				  0                       ,
				  n                       ,
				  tgp_rect                ,
				  1                       ,    // instead of T-B
				  false                   ,    // pred_n
				  true                    ,    // krige
				  _tgp_its->IT_ST_or_IS() ,
				  false                   ,    // delta_s2 (flag for Ds2x)
				  false                   ,    // improv
				  false                   ,    // sens
				  1                         ); // instead of E
 
  // new TGP function for the one point prediction:
  _tgp_model->MAPreplace();

  // prediction:
  _tgp_model->Predict ( tgp_preds  ,
			1          , // instead of T-B
			_tgp_state   );

  // kriging mean:
  ZZ = *tgp_preds->ZZm[0];
  // no need to do:
  //   wmean_of_columns ( &ZZ            ,
  // 		          tgp_preds->ZZm ,
  // 		          tgp_preds->R   ,
  // 		          1              ,
  // 		          (_tgp_its->IT_ST_or_IS()) ? tgp_preds->w : NULL );

  delete_preds ( tgp_preds );

  // unscale:
  unscale_Z ( &ZZ , 1 );

  // treat binary output:
  if ( _is_binary )
    treat_binary_output ( &ZZ , 1 );

  return true;
}

/*--------------------------------------------*/
/*            scale Z or ZZ (private)         */
/*--------------------------------------------*/
void NOMAD::TGP_Output_Model::scale_Z ( void )
{
  if ( _Z_is_scaled || _is_fixed )
    return;
  scale_Z ( _Z , _p );
  _Z_is_scaled = true;
}

void NOMAD::TGP_Output_Model::scale_Z ( double * Z , int n ) const
{
  if ( _is_fixed )
    return;
  for ( int i = 0 ; i < n ; ++i )
    Z[i] = Z[i] / _Z_scaling[0] - _Z_scaling[1];
}

/*--------------------------------------------*/
/*          unscale Z or ZZ (private)         */
/*--------------------------------------------*/
void NOMAD::TGP_Output_Model::unscale_Z ( void )
{
  if ( !_Z_is_scaled || _is_fixed )
    return;
  unscale_Z ( _Z , _p );
  _Z_is_scaled = false;
}

void NOMAD::TGP_Output_Model::unscale_Z ( double * Z , int n ) const
{
  if ( _is_fixed )
    return;
  for ( int i = 0 ; i < n ; ++i )
    Z[i] = ( Z[i] + _Z_scaling[1] ) * _Z_scaling[0];
}

/*--------------------------------------------*/
/*         treat binary output (private)      */
/*--------------------------------------------*/
void NOMAD::TGP_Output_Model::treat_binary_output ( double * ZZ , int nout ) const
{
  NOMAD::Double d0 , d1;
  for ( int j = 0 ; j < nout ; ++j ) {     
    d0 = (_bin_values[0] - ZZ[j]).abs();
    d1 = (_bin_values[1] - ZZ[j]).abs();
    if ( d0 < d1 )
      ZZ[j] = _bin_values[0].value();
    else
      ZZ[j] = _bin_values[1].value();
  }
}

/*--------------------------------------------*/
/*                    display                 */
/*--------------------------------------------*/
void NOMAD::TGP_Output_Model::display ( const NOMAD::Display & out ) const
{
  out << "Z (";
  if ( !_Z_is_scaled )
    out << "un";
  out << "scaled)" << " = [ ";
  for ( int i = 0 ; i < _p ; ++i ) {
    out << std::setw(15) << _Z[i];
    if ( (i+1)%5 == 0 )
      out << " ;" << std::endl << "      ";
    else
      out << " ";
  }
  out << "]" << std::endl
      << "size(Z)=" << _p << std::endl;
  if ( _is_fixed )
    out << "fixed output";
  else if ( _is_binary )
    out << "binary output: values in {"
	<< _bin_values[0] << "," << _bin_values[1]
	<< "}";
  out << std::endl
      << "Z_scaling=[" << _Z_scaling[0] << "," << _Z_scaling[1]
      << "]" << std::endl;
}

#endif
