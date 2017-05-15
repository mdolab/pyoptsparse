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
  \file   TGP_Model_Search.cpp
  \brief  TGP Model search (implementation)
  \author Sebastien Le Digabel
  \date   2011-02-17
  \see    TGP_Model_Search.hpp
*/

#ifndef USE_TGP

int TGP_MODEL_SEARCH_DUMMY; // avoids that TGP_Model_Search.o has no symbols with ranlib

#else

#include "TGP_Model_Search.hpp"

/*-----------------------------------*/
/*           reset (virtual)         */
/*-----------------------------------*/
void NOMAD::TGP_Model_Search::reset ( void )
{
  if ( _model )
    delete _model;
  _model = NULL;
}

/*--------------------------------------------------------*/
/*  delete a list of points: one version for points, and  */
/*  one version for evaluation points (static, private)   */
/*--------------------------------------------------------*/
void NOMAD::TGP_Model_Search::clear_pts ( std::vector<NOMAD::Point *> & pts )
{
  size_t k , n = pts.size();
  for ( k = 0 ; k < n ; ++k )
    delete pts[k];
  pts.clear();
}

void NOMAD::TGP_Model_Search::clear_pts ( std::vector<NOMAD::Eval_Point *> & pts )
{
  size_t k , n = pts.size();
  for ( k = 0 ; k < n ; ++k )
    delete pts[k];
  pts.clear();
}

/*------------------------------------------------------------------*/
/*                             the search                           */
/*------------------------------------------------------------------*/
/* Search parameters:                                               */
/* ------------------                                               */
/*                                                                  */
/*  . MODEL_SEARCH: flag to activate the model search (MS)          */
/*                  (here its value is NOMAD::TGP_MODEL)            */
/*                                                                  */
/*  . MODEL_SEARCH_OPTIMISTIC: if true, the direction from the      */
/*                             model center to the trial point      */
/*                             is computed and prossibly used       */
/*                             in the speculative search            */
/*                             default=yes                          */
/*                                                                  */
/*  . MODEL_SEARCH_PROJ_TO_MESH: project or not to mesh             */
/*                               default=yes                        */
/*                                                                  */
/*  . MODEL_SEARCH_MAX_TRIAL_PTS: limit on the number of trial      */
/*                                points for one search             */
/*                                default=10                        */
/*                                                                  */
/*  . MODEL_TGP_MODE: TGP mode (FAST or PRECISE)                    */
/*                    default=FAST                                  */
/*                                                                  */
/*------------------------------------------------------------------*/
void NOMAD::TGP_Model_Search::search ( NOMAD::Mads              & mads           ,
				       int                      & nb_search_pts  ,
				       bool                     & stop           ,
				       NOMAD::stop_type         & stop_reason    ,
				       NOMAD::success_type      & success        ,
				       bool                     & count_search   ,
				       const NOMAD::Eval_Point *& new_feas_inc   ,
				       const NOMAD::Eval_Point *& new_infeas_inc   )
{
  new_feas_inc  = new_infeas_inc = NULL;
  nb_search_pts = 0;
  success       = NOMAD::UNSUCCESSFUL;
  count_search  = false;

  _one_search_stats.reset();

  const NOMAD::Display    & out = _p.out();
  NOMAD::dd_type display_degree = out.get_search_dd();
  int               display_lim = 15;
 
  if ( stop ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << "TGP_Model_Search::search(): not performed (stop flag is active)"
	  << std::endl;
    return;
  }

  // active cache (we accept only true function evaluations):
  const NOMAD::Cache & cache = mads.get_cache();
  if ( cache.get_eval_type() != NOMAD::TRUTH ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << "TGP_Model_Search::search(): not performed on surrogates"
	  << std::endl;
    return;
  }

  // check that there is one objective exactly:
  const std::list<int> & index_obj_list = _p.get_index_obj();

  if ( index_obj_list.empty() ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << "TGP_Model_Search::search(): not performed with no objective function"
	  << std::endl;
    return;
  }
  if ( index_obj_list.size() > 1 ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << "TGP_Model_Search::search(): not performed with biobjective"
	  << std::endl;
    return;
  }

  // active barrier:
  const NOMAD::Barrier & barrier = mads.get_true_barrier();

  // get the incumbent:
  const NOMAD::Eval_Point * incumbent = barrier.get_best_feasible();
  if ( !incumbent )
    incumbent = barrier.get_best_infeasible();
  if ( !incumbent ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << "TGP_Model_Search::search(): no incumbent"
	  << std::endl;
    return;
  }

  // get and check the signature, and compute the dimension:
  NOMAD::Signature * signature = incumbent->get_signature();

  if ( !signature ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << "TGP_Model_Search::search(): no signature"
	  << std::endl;
    return;
  }

  int n = signature->get_n();

  if ( n != incumbent->size() ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << "TGP_Model_Search::search(): incompatible signature"
	  << std::endl;
    return;
  }

  // initial displays:
  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    std::ostringstream oss;
    oss << NOMAD::MODEL_SEARCH << " #"
 	<< _all_searches_stats.get_MS_nb_searches();
    out << std::endl << NOMAD::open_block ( oss.str() ) << std::endl;
  }

  // from this point the search is counted:
  count_search = true;
  _one_search_stats.add_MS_nb_searches();

// C.Tribes august 26, 2014 --- not needed
//  // mesh:
//  int          mesh_index = NOMAD::Mesh::get_mesh_index();
      NOMAD::Point delta_m;
    // B.Talgorn march 7, 2014 ---
    //NEWMESH
    //signature->get_mesh().get_delta_m ( delta_m , mesh_index );
    signature->get_mesh()->get_delta(delta_m);

  // initial displays:
  if ( display_degree == NOMAD::FULL_DISPLAY ) {
#ifdef TGP_DEBUG
    out << "seed                  : "
	<< _p.get_seed() << std::endl;
#endif
    out << "number of cache points: "   << cache.size()    << std::endl
	// C.Tribes august 26, 2014  --- not needed
    //   << "mesh index            : "   << mesh_index      << std::endl
	<< "mesh size parameter   : ( " << delta_m << " )" << std::endl
	<< "incumbent             : ( ";
    incumbent->NOMAD::Point::display
      ( out , " " , 2 , NOMAD::Point::get_display_limit() );
    out << " )" << std::endl;
  }

  // construct the model:
  NOMAD::Stats                   & stats = mads.get_stats();
  bool                             compute_Ds2x;
  std::vector<NOMAD::Eval_Point *> XX;
  std::string                      error_str;

  if ( !model_construction ( cache          ,
			     *incumbent     ,
			     delta_m        ,
			     out            ,
			     display_degree ,
			     display_lim    ,
			     stats          ,
			     compute_Ds2x   ,
			     XX             ,
			     stop           ,
			     stop_reason    ,
			     error_str        ) ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << NOMAD::close_block ( "failure: " + error_str )
	  << std::endl; 
    return;
  }

  // trial_pts = oracle_pts + Ds2x_pts + improv_pts
  //
  //  oracle_pts: given by the model optimization
  //  Ds2x_pts  : XX points with the largest expected reduction in predictive variance
  //  improv_pts: XX points with the largest expected improvement for the objective

  int max_pts = _p.get_model_search_max_trial_pts();

  /*-----------------------*/
  /*  oracle points (1/3)  */
  /*-----------------------*/
  std::vector<NOMAD::Point *> oracle_pts;
  if ( !create_oracle_pts ( cache          ,
			    *incumbent     ,
			    delta_m        ,
			    out            ,
			    display_degree ,
			    display_lim    ,
			    XX             ,
			    oracle_pts     ,
			    stop           ,
			    stop_reason      ) && stop ) {
    
    // delete XX and oracle_pts:
    NOMAD::TGP_Model_Search::clear_pts ( XX         );
    NOMAD::TGP_Model_Search::clear_pts ( oracle_pts );

    // quit:
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << NOMAD::close_block ( "algorithm stop" )
	  << std::endl; 
    return;
  }

  /*---------------------*/
  /*  Ds2x points (2/3)  */
  /*---------------------*/
  std::vector<NOMAD::Point *> Ds2x_pts;
  if ( compute_Ds2x )
    create_Ds2x_pts ( XX , out , display_degree , display_lim , Ds2x_pts );

  /*-----------------------*/
  /*  improv points (3/3)  */
  /*-----------------------*/
  std::vector<NOMAD::Point *> improv_pts;
  create_improv_pts ( XX             ,
		      *incumbent     ,
		      max_pts        ,
		      out            ,
		      display_degree ,
		      display_lim    ,
		      improv_pts       );

  // create the complete list of trial points:
  // -----------------------------------------
  std::vector<NOMAD::Point *> trial_pts;
  create_trial_pts ( oracle_pts     ,
		     Ds2x_pts       ,
		     improv_pts     ,
		     *incumbent     ,
		     max_pts        ,
		     out            ,
		     display_degree ,
		     trial_pts        );

  // evaluator control:
  NOMAD::Evaluator_Control & ev_control = mads.get_evaluator_control();

  // add the trial points to the evaluator control for evaluation:
  int i , nop = trial_pts.size();
  for ( i = 0 ; i < nop ; ++i )
    register_point ( *trial_pts[i]  ,
		     *signature     ,
		     *incumbent     ,
		     // C.Tribes august 26, 2014 --- not needed
             // mesh_index     ,
		     display_degree ,
		     ev_control       );

  // display the evaluator control list of points:
  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    out << std::endl << NOMAD::open_block ( "list of trial points" );
    const std::set<NOMAD::Priority_Eval_Point> & lop = ev_control.get_eval_lop();
    std::set<NOMAD::Priority_Eval_Point>::const_iterator it , end = lop.end();
    nop = lop.size();
    for ( it = lop.begin() , i = 0 ; it != end ; ++it , ++i ) {
      out << "#";
      out.display_int_w ( i , nop );
      out << " x=( ";
      it->get_point()->NOMAD::Point::display ( out , " " , 15 , -1 );
      out << " )" << std::endl;
    }
    out.close_block();
  }

  // delete XX and the trial points
  // (do not delete Ds2x_pts and improv_pts because
  //  they are XX points, contrary to oracle_pts):
  NOMAD::TGP_Model_Search::clear_pts ( XX           );
  NOMAD::TGP_Model_Search::clear_pts ( oracle_pts   );

  nb_search_pts = ev_control.get_nb_eval_points();

  if ( nb_search_pts == 0 ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << std::endl << "no trial point" << std::endl;
  }

  else {

    _one_search_stats.update_MS_max_search_pts ( nb_search_pts );

    // evaluate the trial points:
    // --------------------------
    int bbe        = stats.get_bb_eval();
    int sgte_eval  = stats.get_sgte_eval ();
    int cache_hits = stats.get_cache_hits();

    new_feas_inc = new_infeas_inc = NULL;

    ev_control.disable_model_eval_sort();

    std::list<const NOMAD::Eval_Point *> * evaluated_pts = NULL;
    if ( display_degree == NOMAD::FULL_DISPLAY )
      evaluated_pts = new std::list<const NOMAD::Eval_Point *>;

    ev_control.eval_list_of_points ( _type                   ,
				     mads.get_true_barrier() ,
				     mads.get_sgte_barrier() ,
				     mads.get_pareto_front() ,
				     stop                    ,
				     stop_reason             ,
				     new_feas_inc            ,
				     new_infeas_inc          ,
				     success                 ,
				     evaluated_pts             );

    // display the prediction error for the evaluated points:
    if ( display_degree == NOMAD::FULL_DISPLAY ) {
      display_eval_pred_errors ( *evaluated_pts , out );
      delete evaluated_pts;
    }

    ev_control.enable_model_eval_sort();

    // update stats:
    _one_search_stats.add_MS_bb_eval    ( stats.get_bb_eval   () - bbe        );
    _one_search_stats.add_MS_sgte_eval  ( stats.get_sgte_eval () - sgte_eval  );
    _one_search_stats.add_MS_cache_hits ( stats.get_cache_hits() - cache_hits );

    if ( success == NOMAD::FULL_SUCCESS )
      _one_search_stats.add_MS_success();

    _one_search_stats.add_MS_pts ( nb_search_pts );
  }

  // update stats objects:
  stats.update_model_stats   ( _one_search_stats );
  _all_searches_stats.update ( _one_search_stats );

  // final display:
  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    std::ostringstream oss;
    oss << "end of " << NOMAD::MODEL_SEARCH << " (" << success << ")";
    out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
  }
}

/*---------------------------------------------------------------*/
/*         create XX a list of prediction points (private)       */
/*---------------------------------------------------------------*/
void NOMAD::TGP_Model_Search::set_XX
( const NOMAD::Cache               & cache     ,
  int                                n         ,
  int                                m         ,
  const NOMAD::Point               & incumbent ,
  const NOMAD::Point               & delta_m   ,
  std::vector<NOMAD::Eval_Point *> & XX          ) const
{
  // . we begin with 1999 points, filter them to remove points
  //     that appear more than once or in the cache.
  // . then we prune the list to 499 points.
  // . we do not begin directly with 500 points because it
  //    gives more flexibility with the projection.
  // . we terminate by adding the incumbent.
  // --> we obtain a list of at most 500 prediction points.

  int  j , i = 0 , n_XX = 1999; // 2000-1
  bool proj_to_mesh = _p.get_model_search_proj_to_mesh();
  bool remove_pt;

  // XX is made of n_XX-1 LH points inside the model bounds
  // (and the incumbent will be added at the end):
  NOMAD::LH_Search::LH_points ( n, m, n_XX, _model->get_lb(), _model->get_ub(), XX );

  while ( i < n_XX ) {

    remove_pt = false;   

    // project to the mesh:
    if ( proj_to_mesh )
      XX[i]->project_to_mesh ( incumbent , delta_m , _p.get_lb() , _p.get_ub() );

    // remove if point is in cache:
    if ( cache.find ( *XX[i] ) )
      remove_pt = true;

    // check if this point is already in XX
    // (may occur only if the point has been projected):
    if ( proj_to_mesh && !remove_pt ) {

      if ( incumbent == (*XX[i]) )
	remove_pt = true;
      
      else
	for ( j = 0 ; j < i ; ++j )
	  if ( XX[j]->NOMAD::Point::operator == ( *XX[i] ) ) {
	    remove_pt = true;
	    break;
	  }
    }

    // remove the point:
    if ( remove_pt ) {
      delete XX[i];
      --n_XX;
      if ( i != n_XX )
	XX[i] = XX[n_XX];
      XX.resize ( n_XX );
    }
    else
      ++i;
  }

  // reduce to 500-1 points (we eliminate randomly):
  while ( n_XX >= 500 ) {

    i = rand()%n_XX;
 
    delete XX[i];
    --n_XX;
    if ( i != n_XX )
      XX[i] = XX[n_XX];
    XX.resize(n_XX);
  }

  // add the incumbent as the last point of XX:
  XX.push_back ( new NOMAD::Eval_Point ( n , m ) );
  XX[n_XX]->NOMAD::Point::operator = ( incumbent );
}

/*---------------------------------------------------------------------*/
/*  create the complete list of trial points (oracle + Ds2x + improv)  */
/*  (private)                                                          */
/*---------------------------------------------------------------------*/
void NOMAD::TGP_Model_Search::create_trial_pts
( const std::vector<NOMAD::Point *> & oracle_pts     ,         // IN
  const std::vector<NOMAD::Point *> & Ds2x_pts       ,         // IN
  const std::vector<NOMAD::Point *> & improv_pts     ,         // IN
  const NOMAD::Point                & incumbent      ,         // IN
  int                                 max_pts        ,         // IN
  const NOMAD::Display              & out            ,         // IN
  NOMAD::dd_type                      display_degree ,         // IN
  std::vector<NOMAD::Point *>       & trial_pts        ) const // OUT
{
  bool   found;
  size_t i , j ,
    n1 = oracle_pts.size() ,
    n2 = Ds2x_pts.size  () ,
    n3 = improv_pts.size();

  std::vector<NOMAD::Point *> l2 , l3;

  // 1. remove duplicates:
  // ---------------------
  //
  // . oracle_pts are not XX points
  // . Ds2x_pts and improv_pts are XX points
  // . there is no duplicate in each list separately
  //
  for ( i = 0 ; i < n2 ; ++i ) {
    found = false;

    if ( *Ds2x_pts[i] == incumbent )
      found = true;
    else {
      for ( j = 0 ; j < n1 ; ++j )
	if ( *Ds2x_pts[i] == *oracle_pts[j] ) {
	  found = true;
	  break;
	}
    }
    if ( !found )
      l2.push_back ( Ds2x_pts[i] );
  }

  n2 = l2.size();

  for ( i = 0 ; i < n3 ; ++i ) {
    found = false;
    if ( *improv_pts[i] == incumbent )
      found = true;
    else {
      for ( j = 0 ; j < n1 ; ++j )
	if ( *improv_pts[i] == *oracle_pts[j] ) {
	  found = true;
	  break;
	}
    }
    if ( !found ) {
      for ( j = 0 ; j < n2 ; ++j )
	if ( improv_pts[i] == l2[j] ) {
	  found = true;
	  break;
	}
    }
    if ( !found )
      l3.push_back ( improv_pts[i] );
  }
    
  n3 = l3.size();

  // 2. construct the list of trial points:
  // ------------------------------------
  trial_pts.clear();

  int nb_pts = static_cast<int> ( n1 + n2 + n3 );

  // no need to reduce the number of trial points:
  if ( max_pts <= 0 || nb_pts <= max_pts ) {


    // 1. oracle points:
    for ( i = 0 ; i < n1 ; ++i )
      trial_pts.push_back ( oracle_pts[i] );
    
    // 2, improv points:
    for ( i = 0 ; i < n3 ; ++i )
      trial_pts.push_back ( l3[i] );    

    // 3. Ds2x points:
    for ( i = 0 ; i < n2 ; ++i )
      trial_pts.push_back ( l2[i] );
    
  }

  // reduce the list to max_pts points:
  else {

    nb_pts = 0;

    size_t i1 = 0 , i2 = 0 , i3 = 0;

    while ( true ) {

      // one point from the oracle points:
      if ( i1 < n1 ) {
	trial_pts.push_back ( oracle_pts[i1++] );
	++nb_pts;
	if ( nb_pts == max_pts )
	  break;
      }

      // two from the improv points:
      for ( i = 0 ; i < 2 ; ++i ) {
	if ( i3 < n3 ) {
	  trial_pts.push_back ( l3[i3++] );
	  ++nb_pts;
	  if ( nb_pts == max_pts )
	    break;
	}
      }
      if ( nb_pts == max_pts )
	break;

      // one from the Ds2x points:
      if ( i2 < n2 ) {
	trial_pts.push_back ( l2[i2++] );
	++nb_pts;
	if ( nb_pts == max_pts )
	  break;
      }
    }
  }

  // 3. display the list of trial points:
  // ------------------------------------
  // if ( display_degree == NOMAD::FULL_DISPLAY ) {   
  //   out.open_block ( "list of trial points (debug)" );
  //   n1 = trial_pts.size();   
  //   for ( i = 0 ; i < n1 ; ++i ) {
  //     out << "#";
  //     out.display_int_w ( i , n1 );
  //     out << " x=( " << *trial_pts[i] << " )" << std::endl;
  //   }
  //   out.close_block();
  // }
}

/*---------------------------------------------------------*/
/*       create the list of improv points (private)        */
/*---------------------------------------------------------*/
void NOMAD::TGP_Model_Search::create_improv_pts
( const std::vector<NOMAD::Eval_Point *> & XX             ,         // IN
  const NOMAD::Point                     & incumbent      ,         // IN
  int                                      max_pts        ,         // IN
  const NOMAD::Display                   & out            ,         // IN
  NOMAD::dd_type                           display_degree ,         // IN
  int                                      display_lim    ,         // IN
  std::vector<NOMAD::Point *>            & improv_pts       ) const // OUT
{
  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << NOMAD::open_block ( "improv points construction" )
	<< std::endl;

  // reset improv points:
  improv_pts.clear();

  // display improv points directly from the NOMAD::TGP_Model object:
#ifdef TGP_DEBUG
  out << NOMAD::open_block ( "expected improvement of the objective (improv)" );
  _model->display_improv ( out , display_lim );
  out << NOMAD::close_block() << std::endl;
#endif

  std::list<int> pts_indexes;
  _model->get_improv_points ( pts_indexes );

  if ( pts_indexes.empty() ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << std::endl
	  << NOMAD::close_block ( "no improv candidate" )
	  << std::endl;
    return;
  }

  std::list<int>::const_iterator it , end = pts_indexes.end();
  
  // with constraints, the list is re-sorted in order to put the
  // predicred feasible points first:
  if ( _p.has_constraints() ) {
  
    std::list<int> feas_pts , infeas_pts;

    NOMAD::Double         h , f;
    const NOMAD::Double & h_min = _p.get_h_min();

    for ( it = pts_indexes.begin() ; it != end ; ++it ) {
      if ( predict ( *XX[*it] , h , f ) && h <= h_min )
	feas_pts.push_back ( *it );
      else
	infeas_pts.push_back ( *it );
    }

    pts_indexes.clear();
    
    end = feas_pts.end();
    for ( it = feas_pts.begin() ; it != end ; ++it )
      pts_indexes.push_back ( *it );
    
    end = infeas_pts.end();
    for ( it = infeas_pts.begin() ; it != end ; ++it )
      pts_indexes.push_back ( *it );

    end = pts_indexes.end();
  }

  // compute max_index just for the display:
  int i , j , max_index = -1 , ni = -1;

  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    ni = pts_indexes.size();
    if ( max_pts > 0 && ni > max_pts )
      ni = max_pts;

    for ( it = pts_indexes.begin() ; it != end ; ++it ) {
      if ( *it > max_index )
	max_index = *it;
    }
  }

  // add the points to improv_pts:
  bool rejected;
  i = j = 0;
  for ( it = pts_indexes.begin() ; it != end ; ++it ) {
    
    // we check the max number of points:
    rejected = ( max_pts > 0 && max_pts == i );

    // we reject the point if it is the incumbent:
    if ( !rejected && incumbent == *XX[*it] )
      rejected = true;

    // we add the point:
    if ( !rejected ) {
      improv_pts.push_back ( XX[*it] );
      ++i;
    }

    // display:
    if ( display_degree == NOMAD::FULL_DISPLAY ) {
      if ( display_lim <= 0 || j < display_lim ) {
	if ( rejected )
	  out << "rejected candidate ";
	else {
	  out << "improv candidate #";     
	  out.display_int_w ( i-1 , ni );
	}
	out << " (XX point #";
	out.display_int_w ( *it , max_index );
	out << "): x=( ";
	XX[*it]->NOMAD::Point::display ( out , " " , 6 , -1 );
	out << " )";
	if ( rejected ) {
	  if ( max_pts > 0 && max_pts == i )
	    out << " (max number of points)";
	  else
	    out << " (candidate==incumbent)";
	}
	else
	  out << " improv=" << _model->get_improv(*it);
	out << std::endl;
      }
      if ( display_lim > 0 && j == display_lim )
	out << "..." << std::endl;
      ++j;
    }

    // if no display, stop the loop if there is already too many points:
    else if ( max_pts > 0 && max_pts == i )
      break;
  }

  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << std::endl
	<< NOMAD::close_block ( "end of improv points construction" )
	<< std::endl;
}

/*------------------------------------------------------------*/
/*  create the list of Ds2x points, the points that maximize  */
/*  the expected reduction in predictive variance             */
/*  (private)                                                 */
/*------------------------------------------------------------*/
void NOMAD::TGP_Model_Search::create_Ds2x_pts
( const std::vector<NOMAD::Eval_Point *> & XX             ,         // IN
  const NOMAD::Display                   & out            ,         // IN
  NOMAD::dd_type                           display_degree ,         // IN
  int                                      display_lim    ,         // IN
  std::vector<NOMAD::Point *>            & Ds2x_pts         ) const // OUT
{
  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << NOMAD::open_block ( "Ds2x points construction" )
	<< std::endl;

  // reset Ds2x points:
  Ds2x_pts.clear();

  // display Ds2x points directly from the NOMAD::TGP_Model object:
#ifdef TGP_DEBUG
  out << NOMAD::open_block ( "expected reduction in predictive variance (Ds2x)" );
  _model->display_Ds2x ( out , display_lim );
  out << NOMAD::close_block() << std::endl;
#endif
    
  std::set<int> pts_indexes;
  _model->get_Ds2x_points ( pts_indexes );

  if ( pts_indexes.empty() ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << std::endl
	  << NOMAD::close_block ( "no Ds2x candidate" )
	  << std::endl;
    return;
  }

  std::set<int>::const_iterator it , end = pts_indexes.end();
  int i , max_index = -1 , m = _p.get_bb_nb_outputs() , n_XX = XX.size();

  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    for ( it = pts_indexes.begin() ; it != end ; ++it ) {
      if ( *it > max_index )
	max_index = *it;
    }
  }

  for ( it = pts_indexes.begin() , i = 0 ; it != end ; ++it , ++i ) {

    Ds2x_pts.push_back ( XX[*it] );

    if ( display_degree == NOMAD::FULL_DISPLAY ) {
      out << "Ds2x candidate #";     
      out.display_int_w ( i , m );
      out << " (XX point #";
      out.display_int_w ( *it , max_index );
      out << "): x=( ";
      XX[*it]->NOMAD::Point::display ( out , " " , 6 , -1 );
      out << " )";
      if ( *it == n_XX - 1 )
	out << " (rejected: candidate==incumbent)";
      out << std::endl;
    }
  }

  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << std::endl
	<< NOMAD::close_block ( "end of Ds2x points construction" )
	<< std::endl;
}

/*-------------------------------------------------------------------*/
/*  create a list of oracle points, given by the model optimization  */
/*  (private)                                                        */
/*-------------------------------------------------------------------*/
bool NOMAD::TGP_Model_Search::create_oracle_pts
( const NOMAD::Cache                     & cache          ,
  const NOMAD::Point                     & incumbent      ,
  const NOMAD::Point                     & delta_m        ,
  const NOMAD::Display                   & out            ,
  NOMAD::dd_type                           display_degree ,
  int                                      display_lim    ,
  const std::vector<NOMAD::Eval_Point *> & XX             ,
  std::vector<NOMAD::Point *>            & oracle_pts     ,
  bool                                   & stop           ,
  NOMAD::stop_type                       & stop_reason      )
{
  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << NOMAD::open_block ( "oracle points construction" )
	<< std::endl;

  // reset oracle points:
  NOMAD::TGP_Model_Search::clear_pts ( oracle_pts );

  int i , n_XX = XX.size();

  // starting points selection:
  const NOMAD::Eval_Point * x0s[3];
  x0s[0] = x0s[1] = x0s[2] = NULL;

  // open display block for model predictions:
#ifdef TGP_DEBUG
  out << NOMAD::open_block ( "TGP predictions (XX+ZZ)");
  int i0 = ( display_lim > 0 ) ? n_XX - display_lim : 0;
  if ( i0 > 0 )
    out << "..." << std::endl;
#endif

  NOMAD::Double         f_model , h_model;
  const NOMAD::Double & h_min  = _p.get_h_min();
  NOMAD::hnorm_type     h_norm = _p.get_h_norm();

  for ( i = 0 ; i < n_XX ; ++i ) {

    // compute model h and f values:
    _model->eval_hf ( XX[i]->get_bb_outputs() ,
		      h_min                   ,
		      h_norm                  ,
		      h_model                 ,
		      f_model                   );    

    if ( h_model.is_defined() && f_model.is_defined() ) {

      XX[i]->set_f ( f_model );
      XX[i]->set_h ( h_model );

      // feasible point:
      if ( XX[i]->is_feasible ( h_min ) ) {
	if ( !x0s[0] || f_model < x0s[0]->get_f() )
	  x0s[0] = XX[i];
      }

      // infeasible point:
      else {
	if ( !x0s[1] || h_model < x0s[1]->get_h() )
	  x0s[1] = XX[i];
	
	if ( !x0s[2] || f_model < x0s[2]->get_f() )
	  x0s[2] = XX[i];
      }
    }

    // display model predictions:
#ifdef TGP_DEBUG
    if ( i >= i0 ) {
      out << "#";
      out.display_int_w ( i , n_XX );
      out << " x=(";
      XX[i]->NOMAD::Point::display ( out , " " , 15 , -1 );
      out << " ) m(x)=[";
      XX[i]->get_bb_outputs().display ( out , " " , 15 , -1 );
      out << " ]";
	      
      if ( h_model.is_defined() && f_model.is_defined() )
	out << " hm=" << h_model << " fm=" << f_model;
      else
	out << " no model value";
      out << std::endl;
    }
#endif	    
      
#ifdef MODEL_STATS
    if ( XX[i] && f_model.is_defined() && h_model.is_defined() ) {
      XX[i]->set_mod_use ( 1                 ); // 1 for model search
      XX[i]->set_Yw      ( _model->get_Yw () );
      XX[i]->set_nY      ( p                 );
      XX[i]->set_mh      ( h_model           );
      XX[i]->set_mf      ( f_model           );
    }
#endif
  }

#ifdef TGP_DEBUG
  // close display block for model predictions:
  {
    std::ostringstream oss;
    oss << "(size=" << n_XX << ")";
    out << NOMAD::close_block ( oss.str() ) << std::endl;
  }

  // compute and display prediction errors:
  out << NOMAD::open_block ( "prediction relative errors on X(%)" );
  _model->display_X_errors ( out );
  out << NOMAD::close_block() << std::endl;

#endif

  if ( !x0s[0] && !x0s[1] && !x0s[2] ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << std::endl
	  << NOMAD::close_block ( "oracle points error: no model starting point" )
	  << std::endl;
    return false;
  }

  // display starting points:
  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    out << std::endl
 	<< NOMAD::open_block ( "model starting points" );
    
    for ( i = 0 ; i < 3 ; ++i ) {

      out << "#" << i << ": ";
      if ( !x0s[i] )
	out << "NULL" << std::endl;
      else {
	
	out << " x=(";
	x0s[i]->NOMAD::Point::display ( out , " " , 15 , -1 );
	out << " ) m(x)=[";
	x0s[i]->get_bb_outputs().display ( out , " " , 15 , -1 );
	out << " ]"
	    << " hm=" << std::setw(15) << x0s[i]->get_h()
	    << " fm=" << std::setw(15) << x0s[i]->get_f()
	    << std::endl;
      }
    }
    out << NOMAD::close_block() << std::endl;
  }

  // optimize model:
  // ---------------
  NOMAD::Point * xf = NULL , * xi = NULL;
  NOMAD::Clock clock;

  bool optimization_ok = optimize_model ( x0s            ,
					  out            ,
					  display_degree ,
					  xf             ,
					  xi             ,
					  stop           ,
					  stop_reason      );

  _one_search_stats.add_optimization_time ( clock.get_CPU_time() );
 
  if ( stop || !optimization_ok || ( xf == NULL && xi == NULL ) ) {
    std::string error_str;
    if ( xf == NULL && xi == NULL )
      error_str = "no model optimization solution";
    else {
      if ( xf ) delete xf;
      if ( xi ) delete xi;
      error_str = ( stop ) ? "algorithm stop" : "model optimization error";
    }

    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << std::endl
	  << NOMAD::close_block ( "oracle points error: " + error_str )
	  << std::endl;
    return false;
  }
 
  // project and check xf and xi:
  if ( xf && !check_oracle_point ( cache          ,
				   incumbent      ,
				   delta_m        ,
				   out            ,
				   display_degree ,
				   *xf              ) ) {
    delete xf;
    xf = NULL;
  }

  if ( xi && !check_oracle_point ( cache          ,
				   incumbent      ,
				   delta_m        ,
				   out            ,
				   display_degree ,
				   *xi              ) ) {
    delete xi;
    xi = NULL;
  }

  // add xf and xi in the list of oracle points:
  if ( xf )
    oracle_pts.push_back ( xf );

  if ( xi ) {

    // check that xi != xf:
    if ( xf && *xf == *xi ) {
      delete xi;
      xi = NULL;
    }

    if ( xi )
      oracle_pts.push_back ( xi );
  }

  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << std::endl
	<< NOMAD::close_block ( "end of oracle points construction" )
	<< std::endl;

  return true;
}

/*------------------------------------------------------*/
/*  project and accept or reject an oracle trial point  */
/*  (private)                                           */
/*------------------------------------------------------*/
bool NOMAD::TGP_Model_Search::check_oracle_point
( const NOMAD::Cache   & cache          ,
  const NOMAD::Point   & incumbent      ,
  const NOMAD::Point   & delta_m        ,
  const NOMAD::Display & out            ,
  NOMAD::dd_type         display_degree ,
  NOMAD::Point         & x                )
{
  bool proj_to_mesh = _p.get_model_search_proj_to_mesh();

  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    out << std::endl << "oracle candidate";
    if ( proj_to_mesh )
      out << " (before projection)";
    out << ": ( " << x << " )" << std::endl;
  }

  // projection to mesh:
  if ( proj_to_mesh ) {
    x.project_to_mesh ( incumbent , delta_m , _p.get_lb() , _p.get_ub() );
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << "oracle candidate (after projection) : ( "
	  << x << " )" << std::endl;
  }

  // compare x and incumbent coordinates:
  if ( x == incumbent ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << "oracle candidate rejected (candidate==incumbent)" << std::endl;
    return false;
  }
  
  // two evaluations points are created in order to:
  //   1. check if the candidate is in cache
  //   2. have a prediction at x and at the incumbent:
  int n = x.size() , m = _p.get_bb_nb_outputs();

  NOMAD::Eval_Point * tk = new NOMAD::Eval_Point ( n , m ); // trial point
  tk->Point::operator = ( x );

  // check if the point is in cache:
  if ( cache.find ( *tk ) ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << "oracle candidate rejected (found in cache)" << std::endl;
    delete tk;
    return false;
  }

  NOMAD::Eval_Point * ic = new NOMAD::Eval_Point ( n , m ); // incumbent copy
  ic->Point::operator = ( incumbent );

  // model predictions (in order to accept or reject the trial point):
  bool pred_error = !_model->predict ( *ic , true ) || // pred_outside_bnds = true
                    !_model->predict ( *tk , true );   // pred_outside_bnds = true

  const NOMAD::Double & h_min  = _p.get_h_min();
  NOMAD::hnorm_type     h_norm = _p.get_h_norm();

  NOMAD::Double h0 , f0; // model values of f and h at the center
  NOMAD::Double h1 , f1; // model values of f and h at the trial point

  if ( !pred_error )
    _model->eval_hf ( tk->get_bb_outputs() , h_min , h_norm , h1 , f1 );
 
  delete tk;

  if ( pred_error || !h1.is_defined() || !f1.is_defined() ) {

    if ( display_degree == NOMAD::FULL_DISPLAY ) {
      if ( pred_error )
 	out << "prediction error: oracle candidate rejected";
      else
 	out << "no model value (EB constraint violated?): oracle candidate rejected";
      out << std::endl;
    }

    delete ic;

    return false;
  }

  // accept or reject the trial point:
  bool accept_point = false;
  _model->eval_hf ( ic->get_bb_outputs(), h_min, h_norm, h0, f0 );
  if ( !h0.is_defined() || !f0.is_defined() )
    accept_point = true;
 
  delete ic;

  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << std::endl
 	<< "incumbent prediction  : h0=" << h0 << " f0=" << f0 << std::endl
 	<< "trial point prediction: h1=" << h1 << " f1=" << f1 << std::endl;

  if ( !accept_point )
    accept_point = (f1 < f0) || (h1 < h0);

  if ( !accept_point ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << "oracle candidate rejected" << std::endl;
    _one_search_stats.add_MS_rejected();
    return false;
  }

  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << "oracle candidate accepted" << std::endl;

  return true;
}

/*--------------------------------------------------------*/
/*  insert a trial point in the evaluator control object  */
/*  (private)                                             */
/*--------------------------------------------------------*/
void NOMAD::TGP_Model_Search::register_point
( NOMAD::Point               x              ,
  NOMAD::Signature         & signature      ,
  const NOMAD::Point       & incumbent      ,
  // C.Tribes august 26, 2014 --- mesh_index not needed
  // int                        mesh_index     ,
  NOMAD::dd_type             display_degree ,
  NOMAD::Evaluator_Control & ev_control       ) const
{
  int n = x.size();

  NOMAD::Eval_Point * tk = new NOMAD::Eval_Point ( n , _p.get_bb_nb_outputs() );

  // if the search is optimistic, a direction is computed (this
  // will be used in case of success in the speculative search):
  if ( _p.get_model_search_optimistic() ) {
    NOMAD::Direction dir ( n , 0.0 , NOMAD::MODEL_SEARCH_DIR );
    dir.Point::operator = ( x - incumbent );
    tk->set_direction  ( &dir );
  }

  tk->set_signature  ( &signature  );
    // B.Talgorn march 7, 2014 ---
    //NEWMESH : mesh_index is not needed anymore. It is included in the signature.
    //tk->set_mesh_index ( &mesh_index );
    tk->Point::operator = ( x );

  // add the new point to the list of search trial points:
  ev_control.add_eval_point ( tk                      ,
			      display_degree          ,
			      _p.get_snap_to_bounds() ,
			      NOMAD::Double()         ,
			      NOMAD::Double()         ,
			      NOMAD::Double()         ,
			      NOMAD::Double()           );
#ifdef MODEL_STATS
  if ( tk ) {
    NOMAD::Double h1 , f1;
    if ( _model->predict ( *tk , true ) ) // pred_outside_bnds = true
      _model->eval_hf ( tk->get_bb_outputs() ,
			_p.get_h_min()       ,
			_p.get_h_norm()      ,
			h1                   ,
			f1                     );
    if ( h1.is_defined() && f1.is_defined() ) {
      tk->set_mod_use ( 1                ); // 1 for model search
      tk->set_Yw      ( _model->get_Yw() );
      tk->set_nY      ( model.get_p ()   );
      tk->set_mh      ( h1               );
      tk->set_mf      ( f1               );
    }
  }
#endif
}

/*---------------------------------------------------------------*/
/*                    model construction (private)               */
/*---------------------------------------------------------------*/
bool NOMAD::TGP_Model_Search::model_construction
( const NOMAD::Cache               & cache          ,
  const NOMAD::Point               & incumbent      ,
  const NOMAD::Point               & delta_m        ,
  const NOMAD::Display             & out            ,
  NOMAD::dd_type                     display_degree ,
  int                                display_lim    ,
  NOMAD::Stats                     & stats          ,
  bool                             & compute_Ds2x   ,
  std::vector<NOMAD::Eval_Point *> & XX             ,
  bool                             & stop           ,
  NOMAD::stop_type                 & stop_reason    ,
  std::string                      & error_str        )
{
  int i , n_XX , n = incumbent.size();

  compute_Ds2x = false;

  // reset XX:
  {
    n_XX = XX.size();
    for ( i = 0 ; i < n_XX ; ++i )
      delete XX[i];
    XX.clear();
    n_XX = 0;
  }

  error_str.clear();

  if ( stop )
    return false;

  NOMAD::Clock clock;

  // TGP model creation:
  if ( _model )
    delete _model;

  NOMAD::TGP_mode_type tgp_mode = _p.get_model_tgp_mode();

  _model = new NOMAD::TGP_Model ( n                       ,
				  _p.get_bb_output_type() ,
				  out                     ,
				  tgp_mode                  );

  // construct interpolation set (X):
  // --------------------------------
  if ( !_model->set_X ( cache          ,
			&incumbent     ,
			_p.get_seed()  ,
			true             ) ) { // remove_fv = true

    if ( _model->get_nep_flag() )
      _one_search_stats.add_not_enough_pts();

    stats.update_model_stats   ( _one_search_stats );
    _all_searches_stats.update ( _one_search_stats );

    error_str = _model->get_error_str();

    delete _model;
    _model = NULL;

    return false;
  }

  // create the list of prediction points (XX)
  // (they are used to get starting points
  //  and to get IMPROV and DS2X stats):
  // -----------------------------------
  set_XX ( cache                          ,
	   n                              ,
	   _p.get_bb_output_type().size() ,
	   incumbent                      ,
	   delta_m                        ,
	   XX                               );

  n_XX = XX.size();

  // display sets X and XX:
  // ----------------------
#ifdef TGP_DEBUG

  // X:
  _model->display_X ( out , display_lim );

  // XX (only the 10 last points are displayed):
  out << NOMAD::open_block ( "prediction points (XX)");
  int i0 = ( display_lim > 0 ) ? n_XX - display_lim : 0;
  if ( i0 > 0 )
    out << "..." << std::endl;
  else if ( i0 < 0 )
    i0 = 0;
  for ( i = i0 ; i < n_XX ; ++i ) {
    out << "#";
    out.display_int_w ( i , n_XX );
    out << " x=(";
    XX[i]->NOMAD::Point::display ( out , " " , 15 , -1 );
    out << " )" << std::endl;
  }
  std::ostringstream oss;
  oss << "(size=" << n_XX << ")";
  out << NOMAD::close_block ( oss.str() ) << std::endl;
#endif 

  // model construction:
  // -------------------
  int p = _model->get_p(); // number of points in X

  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    out << "TGP model construction (p="
	<< p << ") ...";
    out.flush();
  }

  // decide compute_Ds2x flag:
  compute_Ds2x = true;
  if ( tgp_mode == NOMAD::TGP_PRECISE && n_XX > 100 )
    compute_Ds2x = false;

  if ( !_model->compute ( XX           ,
			  compute_Ds2x ,
			  true         ,       // compute_improv    = true
			  false          ) ) { // pred_outside_bnds = false
    
    _one_search_stats.add_construction_error();

    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << "... error" << std::endl;

    // delete XX points:
    for ( i = 0 ; i < n_XX ; ++i )
      delete XX[i];
    XX.clear();
    n_XX = 0;

    stats.update_model_stats   ( _one_search_stats );
    _all_searches_stats.update ( _one_search_stats );

    error_str = _model->get_error_str();

    delete _model;
    _model = NULL;

    // check if ctrl-c has been pressed:
    if ( NOMAD::TGP_Output_Model::get_force_quit() ) {
      stop        = true;
      stop_reason = NOMAD::CTRL_C;
    }

    return false;
  }
  
  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << "... OK" << std::endl << std::endl;

  // update model stats:
  _one_search_stats.add_construction_time ( clock.get_CPU_time() );
  _one_search_stats.update_nY             ( p                    );
  _one_search_stats.add_nb_truth();
  _one_search_stats.add_nb_TGP();

  return true;
}

/*---------------------------------------------------------------*/
/*                    optimize a model (private)                 */
/*---------------------------------------------------------------*/
bool NOMAD::TGP_Model_Search::optimize_model
( const NOMAD::Eval_Point * x0s[3]         ,
  const NOMAD::Display    & out            ,
  NOMAD::dd_type            display_degree ,
  NOMAD::Point           *& xf             ,
  NOMAD::Point           *& xi             ,
  bool                    & stop           ,
  NOMAD::stop_type        & stop_reason      )
{
  // reset xf and xi:
  if ( xf ) delete xf; xf = NULL;
  if ( xi ) delete xi; xi = NULL;

  if ( !_model ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out <<  std::endl << "model optimization error (no model)"
	  << std::endl;
    return false;
  }

  std::string error_str;
  bool        error = false;
  int         i , n = _model->get_n0();

  // model bounds:
  const NOMAD::Point & lb = _model->get_lb();
  const NOMAD::Point & ub = _model->get_ub();

  // initial displays:
  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << std::endl << NOMAD::open_block ( "model optimization" );

  // parameters creation:
  NOMAD::Parameters model_param ( out );
  
    // C.Tribes august 26, 2014 --- Seed is provided once when reading parameter file or NOMAD::Parameters::set_SEED function
	// random seed:
	// model_param.set_SEED ( _p.get_seed() + 10*_all_searches_stats.get_MS_nb_searches() );

  // number of variables:
  model_param.set_DIMENSION ( n );
  
  // blackbox outputs:
  model_param.set_BB_OUTPUT_TYPE ( _p.get_bb_output_type() );

	// C.Tribes august 26, 2014 --- all variables are treated as continuous (integer and binary. Categoricals disable models anyway). Projection on mesh is done later.
	// C.Tribes august 26, --- this change prevents an exception when fixed variables are detected
	// model_param.set_BB_INPUT_TYPE ( _p.get_bb_input_type() );
	// blackbox inputs:
	// Default: all variables are treated as continuous

  // barrier parameters:
  model_param.set_H_MIN  ( _p.get_h_min () );
  model_param.set_H_NORM ( _p.get_h_norm() );

  // starting points:
  for ( i = 0 ; i < 3 ; ++i )
    if ( x0s[i] )
      model_param.set_X0 ( *x0s[i] );

  // fixed variables:
  for ( i = 0 ; i < n ; ++i )
    if ( lb[i] == ub[i] || _p.variable_is_fixed(i) )
      model_param.set_FIXED_VARIABLE(i);

  // no model search and no model ordering:
  model_param.set_MODEL_SEARCH    ( false );
  model_param.set_MODEL_EVAL_SORT ( false );

  // display:
  model_param.set_DISPLAY_DEGREE ( NOMAD::NO_DISPLAY );
  if ( display_degree == NOMAD::FULL_DISPLAY ) {

    model_param.set_DISPLAY_DEGREE ( NOMAD::NORMAL_DISPLAY );
    // model_param.set_DISPLAY_DEGREE ( NOMAD::FULL_DISPLAY );

    if ( n <= 5 )
      model_param.set_DISPLAY_STATS ( "bbe mesh_index ( %14.12gsol ) %14.12gobj" );
    else if ( n <= 10 )
      model_param.set_DISPLAY_STATS ( "bbe mesh_index ( sol ) obj" );
    else
      model_param.set_DISPLAY_STATS ( "bbe obj" );
  }
// C.Tribes august 26, 2014 --- modif for OrthogonalMesh
//  // mesh:
//  int mesh_index       = NOMAD::Mesh::get_mesh_index();
//  int min_mesh_index   = NOMAD::Mesh::get_min_mesh_index();
//  int max_mesh_index   = NOMAD::Mesh::get_max_mesh_index();
//  int max_halton_index = NOMAD::Mesh::get_max_halton_index();
//
//  NOMAD::Mesh::init ( 4.0 , 1 , -1 , 0 );
	// mesh: use isotropic mesh
	model_param.set_ANISOTROPIC_MESH ( false );
	model_param.set_MESH_UPDATE_BASIS ( 4.0 );
	model_param.set_MESH_COARSENING_EXPONENT ( 1 );
	model_param.set_MESH_REFINING_EXPONENT ( -1 );
	model_param.set_INITIAL_MESH_INDEX ( 0 );
	

  // searches:
  // model_param.set_LH_SEARCH ( 1000 , 100 );
  // model_param.set_OPPORTUNISTIC_LH ( true );
  // model_param.set_VNS_SEARCH ( true );

  // maximum number of evaluations (2000 or 10000):
  model_param.set_MAX_BB_EVAL
    ( ( _p.get_model_tgp_mode() == NOMAD::TGP_PRECISE ) ? 10000 : 2000 );

  // min mesh size:
  // model_param.set_MAX_MESH_INDEX ( 30 );
  // model_param.set_MIN_MESH_SIZE ( NOMAD::Double ( 1e-8 ) , false );

  model_param.set_SNAP_TO_BOUNDS ( true );
  // model_param.set_SNAP_TO_BOUNDS ( false );

  // disable user calls:
  model_param.set_USER_CALLS_ENABLED ( false );

  // set flags:
  bool flag_check_bimads , flag_reset_mesh , flag_reset_barriers , flag_p1_active;
  NOMAD::Mads::get_flags ( flag_check_bimads   ,
 			   flag_reset_mesh     ,
 			   flag_reset_barriers ,
 			   flag_p1_active        );

  NOMAD::Mads::set_flag_check_bimads   ( true  );
  NOMAD::Mads::set_flag_reset_mesh     ( true  );
  NOMAD::Mads::set_flag_reset_barriers ( true  );
  NOMAD::Mads::set_flag_p1_active      ( false );

  // bounds:
  model_param.set_LOWER_BOUND ( lb );
  model_param.set_UPPER_BOUND ( ub );

  try {

    // parameters validation:
    model_param.check();

    // out << "TGP PARAMETERS:" << std::endl << model_param << std::endl;

    // model evaluator creation:
    NOMAD::TGP_Model_Evaluator ev ( model_param , *_model );

    // algorithm creation and execution:
    NOMAD::Mads    mads ( model_param , &ev );
    NOMAD::stop_type st = mads.run();
      
    // check the stopping criterion:
    if ( st == NOMAD::CTRL_C || st == NOMAD::MAX_CACHE_MEMORY_REACHED ) {
      std::ostringstream oss;
      oss << "model optimization: " << st;
      error_str   = oss.str();
      error       = true;
      stop        = true;
      stop_reason = st;
    }
    
    else if ( st == NOMAD::MAX_BB_EVAL_REACHED )
      _one_search_stats.add_MS_max_bbe();
   
    // update the stats on the number of model evaluations:
    _one_search_stats.update_MS_model_opt ( mads.get_stats().get_bb_eval() );

    // get the solution(s):
    const NOMAD::Eval_Point * best_feas   = mads.get_best_feasible  ();
    const NOMAD::Eval_Point * best_infeas = mads.get_best_infeasible();

    if ( best_feas )
      xf = new NOMAD::Point ( *best_feas );
    if ( best_infeas )
      xi = new NOMAD::Point ( *best_infeas );
	
    if ( !xf && !xi ) {
      error     = true;
      error_str = "optimization error: no solution";
    }
  }
  catch ( std::exception & e ) {
    error     = true;
    error_str = std::string ( "optimization error: " ) + e.what();
  }

  // reset flags:
  NOMAD::Mads::set_flag_check_bimads   ( flag_check_bimads   );
  NOMAD::Mads::set_flag_reset_mesh     ( flag_reset_mesh     );
  NOMAD::Mads::set_flag_reset_barriers ( flag_reset_barriers );
  NOMAD::Mads::set_flag_p1_active      ( flag_p1_active      );

// C.Tribes august 26, 2014 --- not needed with orthogonalMesh
//  // reset mesh to what it was before:
//  NOMAD::Mesh::init ( _p.get_mesh_update_basis().value() ,
// 		      _p.get_mesh_coarsening_exponent()  , 
// 		      _p.get_mesh_refining_exponent()    ,
// 		      _p.get_initial_mesh_index()          );
//
//  NOMAD::Mesh::set_max_halton_index ( max_halton_index );
//
//  NOMAD::Mesh::set_mesh_index ( min_mesh_index );
//  NOMAD::Mesh::set_mesh_index ( max_mesh_index );
//  NOMAD::Mesh::set_mesh_index ( mesh_index );

  // close display block:
  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    if ( error )
      out.close_block ( error_str );
    else
      out.close_block();
  }

  return !error;
}

/*---------------------------------------------------------*/
/*  display the prediction error for the evaluated points  */
/*  (private)                                              */
/*---------------------------------------------------------*/
bool NOMAD::TGP_Model_Search::predict ( const NOMAD::Point & x ,
					NOMAD::Double      & h ,
					NOMAD::Double      & f   ) const
{
  h.clear();
  f.clear();

  if ( !_model )
    return false;

  NOMAD::Eval_Point y ( x.size() , _p.get_bb_nb_outputs() );
  y.NOMAD::Point::operator = ( x );

  if ( _model->predict ( y , true ) ) // pred_outside_bnds = true
    _model->eval_hf ( y.get_bb_outputs() , _p.get_h_min() , _p.get_h_norm() , h , f );

  return ( h.is_defined() && f.is_defined() );
}

/*---------------------------------------------------------*/
/*  display the prediction error for the evaluated points  */
/*  (private)                                              */
/*---------------------------------------------------------*/
void NOMAD::TGP_Model_Search::display_eval_pred_errors
( const std::list<const NOMAD::Eval_Point *> & evaluated_pts ,
  const NOMAD::Display                       & out             )
{
  if ( !_model )
    return;

  int                   i ,  j = 0;
  int                   nb_pts = evaluated_pts.size()   ,
                        n      = _model->get_n0()       ,
                        m      = _p.get_bb_nb_outputs();
  const NOMAD::Double & h_min  = _p.get_h_min();
  NOMAD::hnorm_type     h_norm = _p.get_h_norm();
  NOMAD::Double         h , f , hm , fm , err;
  NOMAD::Eval_Point     x ( n , m );

  out << std::endl << NOMAD::open_block ( "evaluated points" );
  std::list<const NOMAD::Eval_Point *>::const_iterator it , end = evaluated_pts.end();
  for ( it = evaluated_pts.begin() ; it != end ; ++it ) {

    if ( !(*it) )
      continue;

    h = (*it)->get_h();
    f = (*it)->get_f();

    for ( i = 0 ; i < n ; ++i )
      x[i] = (**it)[i];

    _model->predict ( x , true ); // pred_outside_bnds = true
    _model->eval_hf ( x.get_bb_outputs(), h_min , h_norm , hm , fm );

    out << "#";
    out.display_int_w ( j++ , nb_pts );

    out << " x=(";
    if ( n < 5 )
      (*it)->NOMAD::Point::display ( out , " " , 6 , -1 );
    else
      (*it)->NOMAD::Point::display ( out , " " , -1 , 20 );
    out << ")";
    if ( _p.has_constraints() ) {
      err = (h.is_defined() && hm.is_defined()) ? h.rel_err(hm) * 100.0 : 100.0;
      out << " [h=";
      h.display ( out , "%11.3g" );
      out << " hm=";
      hm.display ( out , "%11.3g" );
      out << " err_h=";
      err.display ( out , "%.2f" );
      out << "%]";
    }
    err = (f.is_defined() && fm.is_defined()) ? f.rel_err(fm) * 100.0 : 100.0;
    out << " [f=";
    f.display ( out , "%11.3g" );
    out << " fm=";
    fm.display ( out , "%11.3g" );
    out << " err_f=";
    err.display ( out , "%.2f" );
    out	<< "%]" << std::endl;
  }

  out.close_block();
}

#endif
