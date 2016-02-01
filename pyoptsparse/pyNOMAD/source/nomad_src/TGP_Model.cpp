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
  \file   TGP_Model.cpp
  \brief  TGP (Bayesian treed Gaussian process) model for all outputs (implementation)
  \author Sebastien Le Digabel
  \date   2011-02-07
  \see    TGP_Model.hpp
*/

#ifndef USE_TGP

int TGP_DUMMY; // avoids that TGP_Model.o has no symbols with ranlib

#else

#include "TGP_Model.hpp"

/*-----------------------------------*/
/*          constructor 1/2          */
/*-----------------------------------*/
NOMAD::TGP_Model::TGP_Model ( int                                        n0   ,
			      const std::vector<NOMAD::bb_output_type> & bbot ,
			      const NOMAD::Display                     & out  ,
			      NOMAD::TGP_mode_type                       mode   )
  : _out            ( out     ) ,
    _bbot           ( bbot    ) ,
    _p              ( 0       ) ,
    _n0             ( n0      ) ,
    _n              ( 0       ) ,
    _n_XX           ( 0       ) ,
    _av_index       ( NULL    ) ,
    _fv_index       ( NULL    ) ,
    _X              ( NULL    ) ,
    _XX             ( NULL    ) ,
    _Ds2x           ( NULL    ) ,
    _improv         ( NULL    ) ,
    _model_computed ( false   ) ,
    _nep_flag       ( false   ) ,
    _tgp_mode       ( mode    ) ,
    _tgp_models     ( NULL    ) ,
    _tgp_rect       ( NULL    ) ,
    _tgp_linburn    ( true    ) ,
    _usr_pmax       ( -1      ) ,
    _tgp_verb       ( false   )
{
  _usr_BTE[0] = _usr_BTE[1] = _usr_BTE[2] = -1;

  // user mode: this is not the good constructor to call:
  if ( mode == NOMAD::TGP_USER )
    throw NOMAD::Exception ( "TGP_Model.cpp" , __LINE__ ,
	  "this constructor only accepts fast or precise TGP modes" );

  _error_str = "NOMAD::TGP_Model::set_X() has not been called";
}

/*--------------------------------------------------------*/
/*  constructor 2/2 (with user BTE and pmax and no bbot)  */
/*--------------------------------------------------------*/
NOMAD::TGP_Model::TGP_Model ( int                    n0      ,
			      const NOMAD::Display & out     ,
			      int                    BTE[3]  ,
			      bool                   linburn ,
			      int                    pmax    ,
			      bool                   verb      )
  : _out            ( out             ) ,
    _p              ( 0               ) ,
    _n0             ( n0              ) ,
    _n              ( 0               ) ,
    _n_XX           ( 0               ) ,
    _av_index       ( NULL            ) ,
    _fv_index       ( NULL            ) ,
    _X              ( NULL            ) ,
    _XX             ( NULL            ) ,
    _Ds2x           ( NULL            ) ,
    _improv         ( NULL            ) ,
    _model_computed ( false           ) ,
    _nep_flag       ( false           ) ,
    _tgp_mode       ( NOMAD::TGP_USER ) ,
    _tgp_models     ( NULL            ) ,
    _tgp_rect       ( NULL            ) ,
    _tgp_linburn    ( linburn         ) ,
    _usr_pmax       ( pmax            ) ,
    _tgp_verb       ( verb            )
{

  // user BTE parameters:
  _usr_BTE[0] = BTE[0];
  _usr_BTE[1] = BTE[1];
  _usr_BTE[2] = BTE[2];

  // default bbot: just one objective output:
  _bbot.push_back ( NOMAD::OBJ );

  _error_str  = "NOMAD::TGP_Model::set_X() has not been called";
}

/*--------------------------------------------*/
/*                 destructor                 */
/*--------------------------------------------*/
NOMAD::TGP_Model::~TGP_Model ( void )
{
  clear();

  if ( _Ds2x ) {
    int i , m = _bbot.size();
    for ( i = 0 ; i < m ; ++i )
      if ( _Ds2x[i] )
	delete [] _Ds2x[i];
    delete [] _Ds2x;
  }

  if ( _improv )
    delete [] _improv;
}

/*--------------------------------------------*/
/*           clear memory (private)           */
/*--------------------------------------------*/
/*  members that depend only on constructor   */
/*  arguments are not reseted here            */
/*--------------------------------------------*/
void NOMAD::TGP_Model::clear ( void )
{
  int i;

  if ( _av_index ) {
    delete [] _av_index;
    _av_index = NULL;
  }

  if ( _fv_index ) {
    delete [] _fv_index;
    _fv_index = NULL;
  }

  if ( _X ) {
    for ( i = 0 ; i < _p ; ++i )
      delete [] _X[i];
    delete [] _X;
    _X = NULL;
  }

  if ( _XX ) {
    for ( i = 0 ; i < _n_XX ; ++i )
      delete [] _XX[i];
    delete [] _XX;
    _XX = NULL;
  }

  if ( _tgp_models ) {
    int m = _bbot.size();
    for ( i = 0 ; i < m ; ++i )
      if ( _tgp_models[i] ) {
	delete _tgp_models[i];
      }
    delete [] _tgp_models;
    _tgp_models = NULL;
  }

  if ( _tgp_rect ) {
    delete_matrix ( _tgp_rect );
    _tgp_rect = NULL;
  }

  _error_str.clear();
  _lb.clear();
  _ub.clear();
  _fv.clear();
  
  _p = _n = _n_XX = 0;

  _model_computed = false;
}

/*--------------------------------------------*/
/*      set the interpolation set X (1/2)     */
/*--------------------------------------------*/
bool NOMAD::TGP_Model::set_X ( const NOMAD::Cache & cache     , // only truth evals
			       const NOMAD::Point * center    , // the incumbent
			       int                  seed      ,
			       bool                 remove_fv   )
{
  std::vector<const NOMAD::Eval_Point *> X;
  
  const NOMAD::Eval_Point * cur = cache.begin();
  while ( cur ) {
    X.push_back ( cur );
    cur = cache.next();
  }

  return set_X ( X , center , seed , remove_fv );
}

/*--------------------------------------------*/
/*      set the interpolation set X (2/2)     */
/*--------------------------------------------*/
bool NOMAD::TGP_Model::set_X
( const std::vector<const NOMAD::Eval_Point *> & X         ,
  const NOMAD::Point                           * center    ,
  int                                            seed      ,
  bool                                           remove_fv   )
{
  clear();
  _nep_flag = false;

  if ( X.size() <= 1 ) {
    _error_str = ( X.empty() ) ?
      "interpolation set X is empty" : "only one point in interpolation set X";
    _nep_flag = true;
    return false;
  }

  // filter and sort X:
  std::list<const NOMAD::Eval_Point *> filtered_X;
  _p = filter_and_sort_X ( X , center , filtered_X );

  if ( _p == 0 ) {
    _error_str = "no valid interpolation point";
    _nep_flag  = true;
    return false;
  }

  int  pmin , pmax;
  bool valid = false;

  while ( !valid ) {

    _n = check_fixed_variables ( filtered_X , remove_fv );

    if ( _n <= 1 ) {
      _error_str = ( _n == 0 ) ?
	"only fixed variables" : "interpolation set of dimension 1";
      return false;
    }

    // get the limits pmin <= p <= pmax:
    if ( !get_p_limits ( _n , pmin , pmax ) )
      return false;

    if ( _p < pmin ) {
      std::ostringstream oss;
      oss << "not enough interpolation points ("
	  << _p << "<=" << pmin << ")";
      _error_str = oss.str();
      _nep_flag  = true;
      return false;
    }

    // reduce the number of interpolation points:
    if ( _p > pmax ) {
      filtered_X.resize ( pmax );
      _p    = pmax;
      valid = false;
    }

    // the interpolation set has a valid size:
    else
      valid = true;
  }

  // display limits on p:
#ifdef TGP_DEBUG
  _out << std::endl
       << "number of interpolation points: "
       << pmin << " <= p=" << _p << " <= " << pmax
       << std::endl;
#endif

  // create interpolation matrix X:
  int i = 0 , j;
  _X = new double * [_p];

  std::list<const NOMAD::Eval_Point *>::const_iterator
    it , end = filtered_X.end();
  for ( it = filtered_X.begin() ; it != end ; ++it , ++i ) {
    _X[i] = new double[_n];
    for ( j = 0 ; j < _n ; ++j )
      _X[i][j] = (**it)[_fv_index[j]].value();
  }

  // check that X is of full rank:
  if ( !NOMAD::TGP_Model::check_full_rank ( _X , _p , _n ) ) {
    _error_str = "X is not of full rank";
    return false;
  }
   
  // construct individual TGP models (one for each output):
  int m = _bbot.size();
  _tgp_models = new NOMAD::TGP_Output_Model * [m];
  for ( i = 0 ; i < m ; ++i ) {
    _tgp_models[i] = NULL;
    if ( _bbot[i] == NOMAD::OBJ || NOMAD::bbot_is_constraint(_bbot[i]) )
      _tgp_models[i] = new NOMAD::TGP_Output_Model ( filtered_X ,
						     i          ,
						     seed       ,
						     _out         );
  }

  return _error_str.empty();
}

/*----------------------------------------------------*/
/*     get the limits pmin <= p <= pmax (private)     */
/*----------------------------------------------------*/
bool NOMAD::TGP_Model::get_p_limits ( int n , int & pmin , int & pmax )
{
  if ( n <= 0 ) {
    _error_str = "bad value of n for the computation of pmin and pmax";
    return false;
  }

  // TGP fast mode:
  // --------------
  //
  // example of some values for pmin and pmax:
  //
  //   n pmin pmax
  //
  //   2    7   30
  //   3    8   30
  //   4    9   30
  //   5   10   30
  //   6   11   32
  //   7   12   35
  //   8   13   37
  //   9   14   40
  //  10   15   42
  //  20   25   60
  //  50   55   94
  //  55   60   99
  //  56   61  100
  //  57   62  100
  //  95  100  100
  //  96  101  101
  // 100  105  105

  if ( _tgp_mode == NOMAD::TGP_FAST ) {
    pmin = n+5;
    
    pmax = static_cast<int> ( floor ( sqrt(180.0*n) ) );

    if ( pmax > 100 )
      pmax = 100;
    
    if ( pmax < 30 )
      pmax = 30;

    if ( pmax < pmin )
      pmax = pmin;
  
    return true;
  }

  // TGP precise mode: pmin=n+1 and pmax=100
  // -----------------
  if ( _tgp_mode == NOMAD::TGP_PRECISE ) {
    pmin = n+1;
    pmax = 100;
    return true;
  }

  // TGP user mode: pmin=n+1 and pmax is user decided with constructor 2/2
  // --------------
  if ( _tgp_mode == NOMAD::TGP_USER ) {
    pmin = n+1;
    pmax = _usr_pmax;
    if ( pmax < pmin ) {
      std::ostringstream oss;
      oss << "user pmax value must be > " << n;
      _error_str = oss.str();
      return false;
    }
    return true;
  }

  return true;
}

/*----------------------------------------------*/
/*     get the TGP BTE parameters (private)     */
/*----------------------------------------------*/
bool NOMAD::TGP_Model::get_BTE ( int BTE[3] )
{
  // fast TGP mode:
  // --------------
  // . if p=pmin, then B=2000
  // . if p=pmax, then B=1000
  // . if pmin<p<pmax, a linear relation is used
  // . B is rounded so that 10 divides it
  // . T=3B
  // . E=10 (E divides T-B)
  if ( _tgp_mode == NOMAD::TGP_FAST ) {

    int pmin , pmax;
    get_p_limits ( _n , pmin , pmax );
    
    double a = 1000.0 / (pmin-pmax);
    double b = 2000 - a * pmin;

    BTE[0] = static_cast<int> ( 10.0 * floor ( (a * _p + b) / 10.0 ) );
    BTE[1] = 3*BTE[0];
    BTE[2] = 10;
  }

  // precise TGP mode:
  // -----------------
  else if ( _tgp_mode == NOMAD::TGP_PRECISE ) {
    BTE[0] = 2000;
    BTE[1] = 7000;
    BTE[2] = 2;
  }

  // user mode:
  // ----------
  else {
    BTE[0] = _usr_BTE[0];
    BTE[1] = _usr_BTE[1];
    BTE[2] = _usr_BTE[2];
   
    // check the user BTE parameters:
    if ( BTE[0] <= 0      ||
	 BTE[1] <= 0      ||
	 BTE[2] <= 0      ||
	 BTE[1] <= BTE[0]    ) {
      _error_str = "error with user BTE";
      return false;
    }
  }

  return true;
}

/*--------------------------------------------------------------*/
/*  compute the ranges, check the fixed variables, set the      */
/*  different indexes, and return the number of free variables  */
/*  (private)                                                   */
/*--------------------------------------------------------------*/
int NOMAD::TGP_Model::check_fixed_variables
( const std::list<const NOMAD::Eval_Point *> & X          ,
  bool                                         remove_fv    )
{

  int i;

  if ( !_av_index ) {
    _av_index = new int [_n0];
    for ( i = 0 ; i < _n0 ; ++i )
      _av_index[i] = -1;
  }

  _lb = NOMAD::Point ( _n0 );
  _ub = NOMAD::Point ( _n0 );
  _fv = NOMAD::Point ( _n0 );

  // compute ranges:
  std::list<const NOMAD::Eval_Point *>::const_iterator it , end = X.end();
  for ( it = X.begin() ; it != end ; ++it ) {
    for ( i = 0 ; i < _n0 ; ++i ) {
      if ( !_lb[i].is_defined() || (**it)[i] < _lb[i] )
	_lb[i] = (**it)[i];
      if ( !_ub[i].is_defined() || (**it)[i] > _ub[i] )
	_ub[i] = (**it)[i];
    }
  }

  // compute n (number of free variables), the fixed variables, and the indexes:
  int n = 0;

  std::vector<int> fv_index_tmp;
  for ( i = 0 ; i < _n0 ; ++i ) {
    if ( remove_fv && _lb[i] == _ub[i] ) {
      _fv      [i] = _lb[i];
      _av_index[i] = -1;
    }
    else {
      fv_index_tmp.push_back ( i );
      _av_index[i] = n++;
    }
  }

  // complete the fixed var index construction:
  if ( _fv_index )
    delete [] _fv_index;
  _fv_index = NULL;
  
  if ( n > 0 ) {
    _fv_index = new int[n];
    for ( i = 0 ; i < n ; ++i )
      _fv_index[i] = fv_index_tmp[i];
  }

  return n;
}

/*--------------------------------------------------------*/
/*  filter and sort an interpolation set X (private)      */
/*    . the points are sorted relatively to the distance  */
/*      from the center                                   */
/*    . if the center is not defined (equal to NULL),     */
/*      an alternate center is constructed                */
/*--------------------------------------------------------*/
int NOMAD::TGP_Model::filter_and_sort_X
( const std::vector<const NOMAD::Eval_Point *> & X          ,
  const NOMAD::Point                           * center     ,
  std::list<const NOMAD::Eval_Point *>         & filtered_X   ) const
{
  NOMAD::Point              alt_center;
  const NOMAD::Eval_Point * cur = NULL;
  int                       p0  = X.size() , i;

  // alternate center if center==NULL:
  if ( !center ) {
    int          j;
    NOMAD::Point lb(_n0) , ub(_n0);
    for ( i = 0 ; i < p0 ; ++i ) {  
      cur = X[i];
      if ( test_interpolation_point ( cur ) ) {
	for ( j = 0 ; j < _n0 ; ++j ) {	  
	  if ( !lb[j].is_defined() || (*cur)[j] < lb[j] )
	    lb[j] = (*cur)[j];
	  if ( !ub[j].is_defined() || (*cur)[j] > ub[j] )
	    ub[j] = (*cur)[j];
	}
      }
    }
    alt_center = NOMAD::Point(_n0);
    for ( j = 0 ; j < _n0 ; ++j )
      alt_center[j] = ( lb[j] + ub[j] ) / 2.0;
  }

  // X_tmp is used to sort the points:
  std::multiset<NOMAD::Model_Sorted_Point> tmp_X;

  for ( i = 0 ; i < p0 ; ++i ) {

    cur = X[i];

    // test if the interpolation point is valid for interpolation:
    if ( test_interpolation_point ( cur ) ) {

      NOMAD::Model_Sorted_Point sorted_pt
	( &NOMAD::Cache::get_modifiable_point (*cur) ,
	  (center) ? *center : alt_center );

      tmp_X.insert ( sorted_pt );
    }
  }

  // copy the set X_tmp to filtered_X:
  std::multiset<NOMAD::Model_Sorted_Point>::const_iterator it , end = tmp_X.end();
  for ( it = tmp_X.begin() ; it != end ; ++it )
    filtered_X.push_back ( static_cast<NOMAD::Eval_Point *> ( it->get_point() ) );

  return filtered_X.size();
}

/*-----------------------------------------------------*/
/*  tests to check if an interpolation point is valid  */
/*  (private)                                          */
/*-----------------------------------------------------*/
bool NOMAD::TGP_Model::test_interpolation_point ( const NOMAD::Eval_Point * x ) const
{
  if ( !x || x->size() != _n0 || !x->is_eval_ok() )
    return false;

  int                  m   = _bbot.size();
  const NOMAD::Point & bbo = x->get_bb_outputs();

  if ( bbo.size() != m )
    return false;

  for ( int j = 0 ; j < m ; ++j )
    if ( ( _bbot[j] == NOMAD::OBJ || NOMAD::bbot_is_constraint(_bbot[j])    ) &&
	 ( !bbo[j].is_defined()   || bbo[j].abs() > NOMAD::MODEL_MAX_OUTPUT )    )
      return false;

  // point is valid:
  return true;
}

/*--------------------------------------------*/
/*  compute the models (one for each output)  */
/*--------------------------------------------*/
bool NOMAD::TGP_Model::compute
( std::vector<NOMAD::Eval_Point *> & XX_pts            ,   // IN/OUT
  bool                               compute_Ds2x      ,   // IN
  bool                               compute_improv    ,   // IN
  bool                               pred_outside_bnds   ) // IN
{
  _model_computed = false;

  if ( !_error_str.empty() )
    return false;

  int i , j , index_obj = -1 , n_XX0 = XX_pts.size() , m = _bbot.size();

  // check bbot: there must be exactly one objective:
  for ( i = 0 ; i < m ; ++i ) {
    if ( _bbot[i] == NOMAD::OBJ ) {    
      if ( index_obj < 0 )
	index_obj = i;
      else {
	_error_str = "more than one objective";
	return false;
      }
    }
  }
  if ( index_obj < 0 ) {
    _error_str = "no objective";
    return false;
  }

  // check n_XX0:
  if ( n_XX0 == 0 ) {
    _error_str = "no user-provided prediction point";
    return false;
  }

  // reset XX_pts outputs:
  for ( i = 0 ; i < n_XX0 ; ++i ) {
    for ( j = 0 ; j < m ; ++j )
      XX_pts[i]->set_bb_output ( j , NOMAD::Double() );
    XX_pts[i]->set_eval_status ( NOMAD::EVAL_FAIL );
  }
    
  // create the XX matrix (prediction points):
  std::vector<NOMAD::Eval_Point *> XX_filtered_pts;
  NOMAD::Eval_Point              * cur2;

  // the list of points has to be filtered:
  NOMAD::Double tmp;
  bool          chk;

  for ( i = 0 ; i < n_XX0 ; ++i ) {
    if ( XX_pts[i]->size() == _n0 ) {
      cur2 = XX_pts[i];
      chk  = true;      
      for ( j = 0 ; j < _n0 ; ++j ) {
	tmp = (*cur2)[j];
	if ( !pred_outside_bnds && ( tmp < _lb[j] || tmp > _ub[j] ) ) {
	  chk = false;
	  break;
	}
      }
      if ( chk )
	XX_filtered_pts.push_back ( cur2 );
    }
  }

  if ( _XX ) {
    for ( i = 0 ; i < _n_XX ; ++i )
      delete [] _XX[i];
    delete [] _XX;
  }

  _n_XX = XX_filtered_pts.size();

  if ( _n_XX == 0 ) {
    _error_str = "no prediction point after filtering";
    return false;
  }

  _XX = new double * [_n_XX];
  for ( i = 0 ; i < _n_XX ; ++i ) {
    _XX[i] = new double[_n];
    for ( j = 0 ; j < _n ; ++j )
      _XX[i][j] = (*XX_filtered_pts[i])[_fv_index[j]].value();
  }

  // Xsplit: X+XX: size = nsplit x n:
  int       nsplit = _p + _n_XX;
  double ** Xsplit = new double * [nsplit];
  
  for ( i = 0 ; i < _p ; ++i ) {
    Xsplit[i] = new double [_n];
    for ( j = 0 ; j < _n ; ++j )
      Xsplit[i][j] = _X[i][j];  
  }

  for ( i = _p ; i < nsplit ; ++i ) {
    Xsplit[i] = new double [_n];
    for ( j = 0 ; j < _n ; ++j )
      Xsplit[i][j] = _XX[i-_p][j]; 
  }

  // get the rectangle:
  if ( _tgp_rect )
    delete_matrix ( _tgp_rect );
  _tgp_rect = get_data_rect ( Xsplit , nsplit , _n );

  // TGP parameters:
  Params   tgp_params ( _n );
  double * dparams = NOMAD::TGP_Model::get_TGP_dparams ( _n );
  tgp_params.read_double ( dparams );
  delete [] dparams;

  int BTE[3];
  if ( !get_BTE ( BTE ) ) {
    for ( i = 0 ; i < nsplit ; ++i )
      delete [] Xsplit[i];
    delete [] Xsplit;
    return false;
  }
  
  // display BTE:
#ifdef TGP_DEBUG
  _out << std::endl
       << "BTE={" << BTE[0] << ", " << BTE[1] << ", " << BTE[2] << "}"
       << std::endl;
#endif

  // compute the individual TGP models (one for each output):
  double * ZZ = new double [_n_XX];

  // Ds2x, expected reduction in predictive variance:
  if ( _Ds2x == NULL ) {
    _Ds2x = new double * [m];
    for ( i = 0 ; i < m ; ++i )
      _Ds2x[i] = NULL;
  }
  else {
    for ( i = 0 ; i < m ; ++i )
      if ( _Ds2x[i] ) {
	delete [] _Ds2x[i];
	_Ds2x[i] = NULL;
      }
  }

  // improv, expected improvement of the objective (ranks):
  if ( _improv ) {
    delete [] _improv;
    _improv = NULL;
  }
  if ( compute_improv )
    _improv = new int [_n_XX];

  for ( i = 0 ; i < m ; ++i ) {

    if ( _tgp_models[i] ) {

      _Ds2x[i] = ( compute_Ds2x ) ? new double [_n_XX] : NULL;

      _tgp_models[i]->compute ( _X                              ,
				_XX                             ,
				Xsplit                          ,
				_n                              ,
				_n_XX                           ,
				nsplit                          ,
				&tgp_params                     ,
				_tgp_rect                       ,
				BTE                             ,
				_tgp_linburn                    ,
				_tgp_verb                       ,
				ZZ                              ,
				_Ds2x[i]                        ,
				(i==index_obj) ? _improv : NULL   );

      // set XX_pts outputs #i:
      for ( j = 0 ; j < _n_XX ; ++j ) {
	XX_filtered_pts[j]->set_bb_output   ( i , ZZ[j]      );
	XX_filtered_pts[j]->set_eval_status ( NOMAD::EVAL_OK );
      }

      // check if TGP has been interrupted:
      if ( NOMAD::TGP_Output_Model::get_force_quit() ) {
	_error_str = "TGP has been interrupted with ctrl-c";
	break;
      }
    }
  }

  // clear memory:
  for ( i = 0 ; i < nsplit ; ++i )
    delete [] Xsplit[i];
  delete [] Xsplit;
  delete [] ZZ;

  _model_computed = _error_str.empty();

  return _model_computed;
}

/*--------------------------------------------*/
/*           prediction at one point          */
/*         (x can be of size _n or _n0)       */
/*--------------------------------------------*/
bool NOMAD::TGP_Model::predict ( NOMAD::Eval_Point & x , bool pred_outside_bnds )
{
  if ( !_error_str.empty() )
    return false;

  if ( !_model_computed ) {
    _error_str = "NOMAD::TGP_Model::compute() has not been called";
    return false;
  }

  int i , i0 , ix , m = x.get_m() , nx = x.size();

  // reset point outputs:
  x.set_eval_status ( NOMAD::EVAL_FAIL );
  for ( i = 0 ; i < m ; ++i )
    x.set_bb_output ( i , NOMAD::Double() );
  
  // check dimensions:
  if ( m != static_cast<int>(_bbot.size()) ||
       ( nx != _n0 && nx != _n ) ) {
    _error_str = "predict error: bad x dimensions";
    return false;
  }

  double ZZ , * XX = new double[_n];

  // set the coordinates and check the bounds:
  for ( i = 0 ; i < _n ; ++i ) {
   
    ix = ( nx == _n0 ) ? _fv_index[i] : i;

    if ( !pred_outside_bnds ) {
      i0 = _fv_index[i];
      if ( x[ix] < _lb[i0] || x[ix] > _ub[i0] ) {
	delete [] XX;
	return false; // this is not an error
      }
    }

    XX[i] = x[ix].value();
  }

  // predictions (one for each output):
  for ( i = 0 ; i < m ; ++i )
    if ( _tgp_models[i] ) {
      if ( !_tgp_models[i]->predict ( XX        ,
				      _n        ,
				      ZZ        ,
				      _tgp_rect   ) ) {
	std::ostringstream oss;
	oss << "predict error: problem with model #" << i;
	_error_str = oss.str();
	break;
      }
      x.set_bb_output ( i , ZZ );
    }

  delete [] XX;

  if ( !_error_str.empty() ) {
    x.set_eval_status ( NOMAD::EVAL_FAIL );
    return false;
  }

  x.set_eval_status ( NOMAD::EVAL_OK );
  return true; 
}

/*-----------------------------------------------------------------*/
/*  this function checks if the p x n matrix X is of full rank by  */
/*  applying the Cholesky decomposition to the sym. def. pos.      */
/*  nxn matrix X'X                                                 */
/*  (static, private)                                              */
/*-----------------------------------------------------------------*/
bool NOMAD::TGP_Model::check_full_rank ( double ** X , int p , int n ) {

  int i , j , k , ki , kii , kij , kj , kji , nn12 = n*(n+1)/2;

  // create XTX (X'X):
  double * XTX = new double [nn12];
  for ( i = 0 ; i < n ; ++i ) {
    ki = i*(i+1)/2;
    for ( j = 0 ; j <= i ; ++j ) {
      kij = ki + j;
      XTX[kij] = 0.0;
      for ( k = 0 ; k < p ; ++k )
	XTX[kij] += X[k][i] * X[k][j];
    }
  }

  // create chol:
  double * chol = new double [nn12] , tmp1 , tmp2 , eps = 1e-10;
  
  // Choleski decomposition:
  for ( i = 0 ; i < n ; ++i ) {

    ki   = i*(i+1)/2;
    kii  = ki+i;
    tmp1 = XTX[kii];
    for ( k = 0 ; k < i ; ++k )
      tmp1 -= pow(chol[ki+k],2.0);
    if ( fabs ( tmp1 ) <= eps ) {
      delete [] XTX;
      delete [] chol;
      return false;
    }

    if ( i == n-1 )
      break;

    tmp1 = sqrt(tmp1);
    chol[kii] = tmp1;
    
    for ( j = i+1 ; j < n ; ++j ) {
      kj   = j*(j+1)/2;
      kji  = kj+i;     
      tmp2 = XTX[kji];
      for ( k = 0 ; k < i ; ++k )
	tmp2 -= chol[ki+k]*chol[kj+k];
      chol[kji] = tmp2/tmp1;
    }
  }
    
  delete [] XTX;
  delete [] chol;

  return true;
}

/*----------------------------------------------------------------*/
/*     compute model h and f values given one blackbox output     */
/*----------------------------------------------------------------*/
void NOMAD::TGP_Model::eval_hf ( const NOMAD::Point  & bbo    ,
				 const NOMAD::Double & h_min  ,
				 NOMAD::hnorm_type     h_norm ,
				 NOMAD::Double       & h      ,
				 NOMAD::Double       & f        ) const
{
  f.clear();
  h = 0.0;

  int m = bbo.size();

  if ( m != static_cast<int>(_bbot.size()) )
    throw NOMAD::Exception ( "TGP_Model.cpp" , __LINE__ ,
	  "TGP_Model::eval_hf() called with an invalid bbo argument" );
  
  NOMAD::Double bboi;

  for ( int i = 0 ; i < m ; ++i ) {
    
    bboi = bbo[i];

    if ( bboi.is_defined() ) {
      
      if ( _bbot[i] == NOMAD::EB || _bbot[i] == NOMAD::PEB_E ) {
	  
	if ( bboi > h_min ) {
	  h.clear();
	  return;
	}
      }
	
      else if ( ( _bbot[i] == NOMAD::FILTER ||
		  _bbot[i] == NOMAD::PB     ||
		  _bbot[i] == NOMAD::PEB_P     ) ) {
	if ( bboi > h_min ) {
	  switch ( h_norm ) {
	  case NOMAD::L1:
	    h += bboi;
	    break;
	  case NOMAD::L2:
	    h += bboi * bboi;
	    break;
	  case NOMAD::LINF:
	    if ( bboi > h )
	      h = bboi;
	    break;
	  }
	}
      }
	
      else if ( _bbot[i] == NOMAD::OBJ )
	f = bboi;
    }
   
  }

  if ( h_norm == NOMAD::L2 )
    h = h.sqrt();
}

/*---------------------------------------*/
/*      compute the TGP dparam array     */
/*      (static, private)                */
/*---------------------------------------*/
double * NOMAD::TGP_Model::get_TGP_dparams ( int n ) {

  double * dparams  = new double [n*(n+3)+41];
  int      i , j , k;

  // tree (p <- c(as.numeric(params$tree))):
  dparams[0] = 0.5;
  dparams[1] = 2.0;
  dparams[2] = n+2;
  if ( dparams[2] < 10 )
    dparams[2] = 10;
  dparams[3] = 1.0;
  dparams[4] = n;
  
  // params$meanfn == "linear" :
  dparams[5] = 0.0;

  // params$bprior == "bflat" :
  dparams[6] = 2.0;

  // p <- c(p, as.numeric(params$beta)) : n+1 zeros
  k = 6;
  for ( i = 0 ; i <= n ; ++i )
    dparams[++k] = 0.0;

  // p <- c(p, as.numeric(params$Wi)) : I_{n+1}
  for ( i = 0 ; i <= n ; ++i )
    for ( j = 0 ; j <= n ; ++j )
      dparams[++k] = (i!=j) ? 0.0 : 1.0;

  // p <- c(p, as.numeric(params$s2tau2)):
  dparams[++k] = 1.0;
  dparams[++k] = 1.0;

  // p <- c(p, as.numeric(params$s2.p)):
  dparams[++k] =  5.0;
  dparams[++k] = 10.0;

  // p <- c(p, as.numeric(params$s2.lam)):
  dparams[++k] =  0.2;
  dparams[++k] = 10.0;

  // p <- c(p, as.numeric(params$tau2.p)):
  dparams[++k] =  5.0;
  dparams[++k] = 10.0;

  // p <- c(p, as.numeric(params$tau2.lam)):
  dparams[++k] = 0.2;
  dparams[++k] = 0.1;

  // params$corr == "expsep" :
  dparams[++k] = 1.0;

  // p <- c(p, as.numeric(params$gd)):
  dparams[++k] = 0.1;
  dparams[++k] = 0.5;

  // p <- c(p, as.numeric(params$nug.p)):
  dparams[++k] = 1.0;
  dparams[++k] = 1.0;
  dparams[++k] = 1.0;
  dparams[++k] = 1.0;

  // if (params$nug.lam[1] == "fixed"), p <- c(p, rep(-1, 4)) :
  dparams[++k] = -1.0;
  dparams[++k] = -1.0;
  dparams[++k] = -1.0;
  dparams[++k] = -1.0;

  // p <- c(p, as.numeric(params$gamma)):
  dparams[++k] = 0.0;
  dparams[++k] = 0.2;
  dparams[++k] = 0.7;

  // p <- c(p, as.numeric(params$d.p)):
  dparams[++k] =  1.0;
  dparams[++k] = 20.0;
  dparams[++k] = 10.0;
  dparams[++k] = 10.0;

  // if (params$d.lam[1] == "fixed"), p <- c(p, rep(-1, 4)):
  dparams[++k] = -1.0;
  dparams[++k] = -1.0;
  dparams[++k] = -1.0;
  dparams[++k] = -1.0;

  return dparams;
}

/*--------------------------------------------*/
/*       display the interpolation set X      */
/*--------------------------------------------*/
void NOMAD::TGP_Model::display_X ( const NOMAD::Display & out           ,
				   int                    display_limit   ) const
{
  if ( _p == 0 || !_X ) {
    out << "no interpolation points" << std::endl;
    return;
  }

  int i , j;
  int m  = _bbot.size();
  int i0 = ( display_limit > 0 ) ? _p - display_limit : 0;
  NOMAD::Point x(_n) , bbo(m);

  out << NOMAD::open_block ( "interpolation points (X)");

  if ( i0 > 0 )
    out << "..." << std::endl;
  else if ( i0 < 0 )
    i0 = 0;

  for ( i = i0 ; i < _p ; ++i ) {

    for ( j = 0 ; j < _n ; ++j )
      x[j] = _X[i][j];

    bbo = NOMAD::Point(m);
    if ( _tgp_models )
      for ( j = 0 ; j < m ; ++j )
	if ( _tgp_models[j] )
	  bbo[j] = (_tgp_models[j]->get_Z())[i];
 
    out << "#";
    out.display_int_w ( i , _p );
    out << " x=(";
    x.display ( out , " " , 15 , -1 );
    out << " ) f(x)=[";
    bbo.display ( out , " " , 15 , -1 );
    out << " ]" << std::endl;
  }
  
  std::ostringstream oss;
  oss << "(size=" << _p << ")";
  out << NOMAD::close_block ( oss.str() ) << std::endl;
}

/*---------------------------------------------------------*/
/*  get the XX points with the largest expected reduction  */
/*  in predictive variance, for each output                */
/*  (no duplicates)                                        */
/*---------------------------------------------------------*/
void NOMAD::TGP_Model::get_Ds2x_points ( std::set<int> & pts_indexes ) const
{
  pts_indexes.clear();
  if ( !_Ds2x || _n_XX == 0 )
    return;

  int i , j , k , m = _bbot.size();

  for ( i = 0 ; i < m ; ++i )
    if ( _Ds2x[i] ) {

      NOMAD::Double max;
      k = -1;
      for ( j = 0 ; j < _n_XX ; ++j )
	if ( !max.is_defined() || _Ds2x[i][j] > max ) {
	  max = _Ds2x[i][j];
	  k   = j;
	}

      if ( k >= 0 )
	pts_indexes.insert ( k );
    }
}

/*----------------------------------------------------------------------------*/
/*  get the XX points with the largest expected improvement of the objective  */
/*----------------------------------------------------------------------------*/
void NOMAD::TGP_Model::get_improv_points ( std::list<int> & pts_indexes ) const
{
  pts_indexes.clear();
  if ( !_improv || _n_XX == 0 )
    return;

  int            j;
  NOMAD::Point * XX_pt;
  std::multiset<NOMAD::Model_Sorted_Point> pts;
  std::multiset<NOMAD::Model_Sorted_Point>::const_iterator it , end;

  // 1. sort:
  for ( j = 0 ; j < _n_XX ; ++j ) {   
    XX_pt = new NOMAD::Point ( 1 );
    (*XX_pt)[0] = j;
    pts.insert ( NOMAD::Model_Sorted_Point ( XX_pt , _improv[j] ) );
  }

  // 2. construct pts_indexes (exclude points with improv >= n_XX):
  end = pts.end();
  for ( it = pts.begin() ; it != end ; ++it ) {
    if ( it->get_dist() < _n_XX )
      pts_indexes.push_back ( static_cast<int> ( (*it->get_point())[0].value() ) );
    delete it->get_point();
  }
}

/*--------------------------------------------------------------*/
/*  display the expected improvement of the objective (improv)  */
/*  (better ranks are displayed first)                          */
/*--------------------------------------------------------------*/
void NOMAD::TGP_Model::display_improv ( const NOMAD::Display & out ,
					int          display_limit   ) const
{
  if ( !_improv || _n_XX == 0 ) {
    out << "improv has not been computed" << std::endl;
    return;
  }

  int            j , k;
  NOMAD::Point * XX_pt;
  std::multiset<NOMAD::Model_Sorted_Point> pts;
  std::multiset<NOMAD::Model_Sorted_Point>::const_iterator it , end;

  // 1. sort:
  for ( j = 0 ; j < _n_XX ; ++j ) {
    
    // construct a NOMAD::Point from matrix _XX:
    XX_pt = new NOMAD::Point ( _fv );
    for ( k = 0 ; k < _n0 ; ++k )
      if ( _av_index[k] >= 0 )
	(*XX_pt)[k] = _XX[j][_av_index[k]];
    
    // insert this point in the sorted list:
    pts.insert ( NOMAD::Model_Sorted_Point ( XX_pt , _improv[j] ) );
  }

  // 2. display:
  end = pts.end();
  for ( j = 0 , it = pts.begin() ; it != end ; ++it , ++j ) {
    
    if ( display_limit <= 0 || j < display_limit ) {
      out << "x=( ";
      it->get_point()->display ( out , " " , 6 , -1 );
      out << " ) improv=" << it->get_dist() << std::endl;
    }

    else if ( j == display_limit )
      out << "..." << std::endl;

    delete it->get_point();
  }
}

/*----------------------------------------------------------------*/
/*  display the expected reduction in predictive variance (Ds2x)  */
/*  (larger values are displayed first)                           */
/*----------------------------------------------------------------*/
void NOMAD::TGP_Model::display_Ds2x ( const NOMAD::Display & out ,
				      int          display_limit   ) const
{
  if ( !_Ds2x || _n_XX == 0 ) {
    out << "matrix Ds2x has not been computed" << std::endl;
    return;
  }
    
  int            i , j , k , m = _bbot.size();
  NOMAD::Point * XX_pt;
  std::multiset<NOMAD::Model_Sorted_Point> pts;
  std::multiset<NOMAD::Model_Sorted_Point>::const_iterator it , end;

  for ( i = 0 ; i < m ; ++i ) {

    if ( m > 1 ) {
      std::ostringstream oss;
      oss << "output #" << i;
      out << NOMAD::open_block ( oss.str() );
    }

    if ( _Ds2x[i] ) {

      // 1. sort:
      for ( j = 0 ; j < _n_XX ; ++j ) {

	// construct a NOMAD::Point from matrix _XX:
	XX_pt = new NOMAD::Point ( _fv );
	for ( k = 0 ; k < _n0 ; ++k )
	  if ( _av_index[k] >= 0 )
	    (*XX_pt)[k] = _XX[j][_av_index[k]];
	
	// insert this point in the sorted list:
	pts.insert ( NOMAD::Model_Sorted_Point ( XX_pt , -_Ds2x[i][j] ) );
      }

      // 2. display:
      end = pts.end();
      for ( j = 0 , it = pts.begin() ; it != end ; ++it , ++j ) {
    
	if ( display_limit <= 0 || j < display_limit ) {
	  out << "x=( ";
	  it->get_point()->display ( out , " " , 6 , -1 );
	  out << " ) Ds2x=" << it->get_dist()*-1.0 << std::endl;
	}

	else if ( j == display_limit )
	  out << "..." << std::endl;

	delete it->get_point();
      }

      pts.clear();
    }
    else
      out << "NULL" << std::endl;

    if ( m > 1 )
      out.close_block();
  }
}

/*-------------------------------------------------------*/
/*  display the error stats for the interpolation set X  */
/*-------------------------------------------------------*/
void NOMAD::TGP_Model::display_X_errors ( const NOMAD::Display & out )
{
  if ( _p == 0 || !_X ) {
    out << "no interpolation points" << std::endl;
    return;
  }

  int               i , j , m = _bbot.size();
  NOMAD::Point      min_err(m) , max_err(m) , avg_err(m,0.0) , sd_err(m,0.0);
  NOMAD::Eval_Point x ( _n , m );
  double         ** err = new double * [_p];

  for ( i = 0 ; i < _p ; ++i ) {

    err[i] = new double[m];

    for ( j = 0 ; j < _n ; ++j )
      x[j] = _X[i][j];

    if ( predict ( x , true ) ) {

      for ( j = 0 ; j < m ; ++j )
	if ( _tgp_models[j] ) {

	  // relative error (in %) for point #i and output #j:
	  err[i][j] = ( x.get_bb_outputs()[j].rel_err((_tgp_models[j]->get_Z())[i])
			* 100.0).value();

	  // out << "f=" << (_tgp_models[j]->get_Z())[i] << " "
	  //     << "m=" << x.get_bb_outputs()[j] << " err=" << err[i][j]
	  //     << std::endl;

	  if ( !min_err[j].is_defined() || err[i][j] < min_err[j].value() )
	    min_err[j] = err[i][j];
	  
	  if ( !max_err[j].is_defined() || err[i][j] > max_err[j].value() )
	    max_err[j] = err[i][j];
	  
	  avg_err[j] += err[i][j];
	}
    }
    else {
      for ( j = 0 ; j <= i ; ++j )
	delete [] err[j];
      delete [] err;
      out << "cannot predict interpolation errors ("
	  << _error_str << ")" << std::endl;
      return;
    }
  }

  for ( j = 0 ; j < m ; ++j )
    if ( _tgp_models[j] ) {

      // compute the median error:
      NOMAD::Double med_err;
      {	
	if ( _p == 1 )
	  med_err = err[0][j];
	else if ( _p == 2 )
	  med_err = ( err[0][j] + err[1][j] ) / 2.0;

	else {
	  std::multiset<double> sorted_errors;
	  for ( i = 0 ; i < _p ; ++i )
	    sorted_errors.insert ( err[i][j] );
	  std::multiset<double>::const_iterator it , end = sorted_errors.end();
	  --end;
	  for ( it = sorted_errors.begin() , i = 0 ; it != end ; ++it , ++i ) {
	    if ( i == (_p+1)/2-1 ) {
	      med_err = *it;
	      if ( _p%2==0 ) {
		++it;
		med_err = ( med_err + *it ) / 2.0;
	      }
	      break;
	    }
	  }
	}
      }

      // compute the mean and the standard deviation:
      avg_err[j] /= _p;
      for ( i = 0 ; i < _p ; ++i )
	sd_err[j] += ( avg_err[j] - err[i][j] ).pow2();
      sd_err[j] = (sd_err[j] / _p).sqrt();

      // display:
      if ( m > 1 ) {
	std::ostringstream oss;
	oss << "output #" << j;
	if ( _tgp_models[j]->is_fixed() )
	  oss << " (fixed)";
	else if ( _tgp_models[j]->is_binary() )
	  oss << " (binary)";
	out << NOMAD::open_block ( oss.str() );
      }

      out << "min   : ";
      min_err[j].display ( out , "%6.2f" );
      out << std::endl << "max   : ";
      max_err[j].display ( out , "%6.2f" );
      out << std::endl << "median: ";
      med_err.display ( out , "%6.2f" );
      out << std::endl << "mean  : ";
      avg_err[j].display ( out , "%6.2f" );
      out << std::endl << "sd    : ";
      sd_err[j].display ( out , "%6.2f" );
      out << std::endl;

      if ( m > 1 )
	out.close_block();
    }

  for ( i = 0 ; i < _p ; ++i )
    delete [] err[i];
  delete [] err;
}

/*--------------------------------------------*/
/*                    display                 */
/*--------------------------------------------*/
void NOMAD::TGP_Model::display ( const NOMAD::Display & out ) const
{
  if ( !_error_str.empty() ) {
    out << "error with model" << std::endl;
    return;
  }

  int i , j;

  // fixed variables:
  out << "fixed_var = [ " << _fv << "]" << std::endl;
  out << "av_index  = ";
  if ( _av_index ) {
    out << "[ ";
    for ( i = 0 ; i < _n0 ; ++i )
      out << _av_index[i] << " ";
    out << "]" << std::endl;
  }
  else
    out << "NULL" << std::endl;
  out << "fv_index  = ";
  if ( _fv_index ) {
    out << "[ ";
    for ( i = 0 ; i < _n ; ++i )
      out << _fv_index[i] << " ";
    out << "]" << std::endl;
  }
  else
    out << "NULL" << std::endl;

  // bounds:
  out << "lb        = [ " << _lb << "]" << std::endl
      << "ub        = [ " << _ub << "]" << std::endl
      << std::endl;

  // display X:
  if ( !_X )
    out << "X = NULL" << std::endl;
  else {
    out << "X = [";
    for ( i = 0 ; i < _p ; ++i ) {
      out << "\t";
      for ( j = 0 ; j < _n ; ++j )
	out << std::setw(15) << _X[i][j] << " ";
      out << ( (i==_p-1) ? "]" : ";" ) << std::endl;
    }
    out << "size(X)=" << _p << "x" << _n << std::endl << std::endl;
  }

  // display XX:
  if ( !_XX )
    out << "XX = NULL" << std::endl;
  else {
    out << "XX = [";
    for ( i = 0 ; i < _n_XX ; ++i ) {
      out << "\t";
      for ( j = 0 ; j < _n ; ++j )
	out << std::setw(15) << _XX[i][j] << " ";
      out << ( (i==_n_XX-1) ? "]" : ";" ) << std::endl;
    }
    out << "size(XX)=" << _n_XX << "x" << _n << std::endl << std::endl;
  }

  // display models:
  out << std::endl;
  if ( _tgp_models ) {
    int m = _bbot.size();
    for ( i = 0 ; i < m ; ++i ) {
      if ( _tgp_models[i] ) {
	std::ostringstream oss;
	oss << "model #" << i;
	out.open_block ( oss.str() );
	_tgp_models[i]->display ( out );
	out.close_block();
      }
      else
	out << "model #" << i << ": NULL" << std::endl;
      out << std::endl;
    }
  }
  else
    out << "no models" << std::endl << std::endl;
}

/*----------------------------------------------------------------*/
/*           access to the width of the interpolation set         */
/*----------------------------------------------------------------*/
#ifdef MODEL_STATS
NOMAD::Double NOMAD::TGP_Model::get_Yw ( void ) const
{
  NOMAD::Double Yw , tmp;
  for ( int i = 0 ; i < _n0 ; ++i ) {
    tmp = _ub[i]-_lb[i];
    if ( !Yw.is_defined() || tmp > Yw )
      Yw = tmp;
  }
  return Yw;
}
#endif

#endif
