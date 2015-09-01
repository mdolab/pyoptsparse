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
  \file   Cache_Search.cpp
  \brief  NOMAD::Search subclass for the cache search (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-08
  \see    Cache_Search.hpp
*/
#include "Cache_Search.hpp"

/*---------------------------------------------------------*/
/*                       the search                        */
/*---------------------------------------------------------*/
void NOMAD::Cache_Search::search ( NOMAD::Mads              & mads           ,
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

  NOMAD::Evaluator_Control & ev_control    = mads.get_evaluator_control();
  const Cache              & cache         = mads.get_cache();
  int                        nb_extern_pts = cache.get_nb_extern_points();

  // do not perform the search if the number of extern points did not change:
  if ( stop || nb_extern_pts == 0 || nb_extern_pts == _last_search_nb_extern_pts )
    return;

  count_search = true;

  // initial display:
  const NOMAD::Display    & out = _p.out();
  NOMAD::dd_type display_degree = out.get_search_dd();
  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    std::ostringstream oss;
    oss << NOMAD::CACHE_SEARCH;
    out << std::endl << NOMAD::open_block ( oss.str() ) << std::endl;
  }

  std::list<const NOMAD::Eval_Point*> list_of_extern_pts;
  const NOMAD::Eval_Point * extern_pt = cache.get_and_remove_extern_point();
  NOMAD::Eval_Point       * pt        = NULL;
  int                       n         = _p.get_dimension();

  // add the extern points to the list of points to be evaluated:
  while ( extern_pt ) {

    list_of_extern_pts.push_front ( extern_pt );

    pt = &NOMAD::Cache::get_modifiable_point ( *extern_pt );
    
    if ( extern_pt->get_signature() )
      pt->set_signature ( extern_pt->get_signature() );
    else if ( extern_pt->size() == n )
      pt->set_signature ( _p.get_signature() );
    
    if ( pt->get_signature() )
      ev_control.add_eval_point ( pt              ,
				  display_degree  ,
				  false           ,
				  NOMAD::Double() ,
				  NOMAD::Double() ,
				  NOMAD::Double() ,
				  NOMAD::Double()   );
    else {
      if ( display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY )
	out << std::endl << "Warning (Cache_Search.cpp, " << __LINE__
	    << "): could not use the point " << *pt
	    << "(no signature)" << std::endl;
    }

    extern_pt = cache.get_and_remove_extern_point();
  }

  nb_search_pts = ev_control.get_nb_eval_points();

  // display all the search points:
  if ( display_degree == NOMAD::FULL_DISPLAY )
    ev_control.display_eval_lop ( NOMAD::CACHE_SEARCH );

  // eval_list_of_points:
  new_feas_inc = new_infeas_inc = NULL;
  ev_control.eval_list_of_points ( NOMAD::CACHE_SEARCH     ,
				   mads.get_true_barrier() ,
				   mads.get_sgte_barrier() ,
				   mads.get_pareto_front() ,
				   stop                    ,
				   stop_reason             ,
				   new_feas_inc            ,
				   new_infeas_inc          ,
				   success                   );

  // the method cache.get_and_remove_extern_point() removes the first extern
  // point from the cache. If the search is opportunistic and if there are
  // extern points that have not been treated, then they must be put back
  // in the cache list of extern points:
  {
    std::list<const NOMAD::Eval_Point*>::const_iterator
      it , end = list_of_extern_pts.end();
    for ( it = list_of_extern_pts.begin() ; it != end ; ++it )
      if ( !(*it)->get_current_run() )
	cache.insert_extern_point ( **it );
  }

  // update _last_search_nb_extern_pts:
  _last_search_nb_extern_pts = cache.get_nb_extern_points();

  // final display:
  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    std::ostringstream oss;
    oss << "end of cache search (" << success << ")";
    out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
  }
}
