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
 \file   Phase_One_Search.cpp
 \brief  NOMAD::Search subclass for the phase one (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-09
 \see    Phase_One_Search.hpp
 */
#include "Phase_One_Search.hpp"

/*-------------------------------------------------------------*/
/*                         phase one search                    */
/*                 (try to satisfy EB constraints)             */
/*-------------------------------------------------------------*/
void NOMAD::Phase_One_Search::search ( NOMAD::Mads              & mads           ,
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
    stop          = false;
    count_search  = true;
    
    const NOMAD::Display    & out = _p.out();
    NOMAD::dd_type display_degree = out.get_search_dd();
    
    // initial display:
    if ( display_degree == NOMAD::FULL_DISPLAY )
    {
        std::ostringstream oss;
        oss << NOMAD::P1_SEARCH;
        out << std::endl << NOMAD::open_block ( oss.str() ) << std::endl;
    }
    
    // stats:
    NOMAD::Stats & stats = mads.get_stats();
    
    // counters:
    int old_bbe = stats.get_bb_eval();
    int old_it  = stats.get_iterations();
    
    // Evaluator_Control:
    NOMAD::Evaluator_Control & ev_control = mads.get_evaluator_control();
    
    // save and modify parameters:
    std::string old_display_degree;
    _p.out().get_display_degree ( old_display_degree );
    const std::list<std::string>         old_ds = _p.get_display_stats();
    NOMAD::Double                     old_VNS_trigger = _p.get_VNS_trigger();
    const std::string             old_stats_file_name = _p.get_stats_file_name();
    const std::string                    old_sol_file = _p.get_solution_file();
    const std::list<std::string>       old_stats_file = _p.get_stats_file();
    const NOMAD::Point                   old_f_target = _p.get_f_target();
    NOMAD::Double                             old_lct = _p.get_L_curve_target();
    bool                                      old_sif = _p.get_stop_if_feasible();
    const std::vector<NOMAD::bb_output_type> old_bbot = _p.get_bb_output_type();
    std::vector<NOMAD::bb_output_type>        p1_bbot = old_bbot;
    
    
    if ( display_degree == NOMAD::NORMAL_DISPLAY)  // Normal display -> minimal display for Phase one
        _p.set_DISPLAY_DEGREE( NOMAD::MINIMAL_DISPLAY);
    else if (display_degree == NOMAD::FULL_DISPLAY)
        _p.set_DISPLAY_DEGREE( NOMAD::FULL_DISPLAY);// Full display -> full display for Phase one
    
    
    
    int  m   = static_cast<int> ( old_bbot.size() );
    int  cnt = 0;
    for ( int i = 0 ; i < m ; ++i )
    {
        if ( old_bbot[i] == NOMAD::EB )
        {
            p1_bbot[i] = NOMAD::OBJ;
            ++cnt;
        }
        else if ( old_bbot[i] == NOMAD::OBJ )
            p1_bbot[i] = NOMAD::UNDEFINED_BBO;
    }
    
    if ( cnt == 0 ) {
        stop        = true;
        stop_reason = NOMAD::P1_FAIL;
        return;
    }
    
    _p.set_F_TARGET         ( NOMAD::Point ( cnt , 0.0 ) );
    _p.set_L_CURVE_TARGET   ( NOMAD::Double()            );
    _p.set_STOP_IF_FEASIBLE ( false                      );
    _p.set_VNS_SEARCH       ( false                      );
    _p.set_BB_OUTPUT_TYPE   ( p1_bbot                    );
    _p.set_SOLUTION_FILE    ( ""                         );
    _p.reset_stats_file();
    
    // DISPLAY_STATS  and STATS_FILE
    {
        std::list<std::string>                 ds    = old_ds;
        std::list<std::string>                 sf    = old_stats_file;
        ds.push_back ( " (PhaseOne)" );
        _p.set_DISPLAY_STATS ( ds );
        sf.push_back ( " (PhaseOne)" );
        _p.set_STATS_FILE ( old_stats_file_name , sf );
        
    }
    
    _p.check ( false ,    // remove_history_file  = false
              true ,    // remove_solution_file = true
              false    ); // remove_stats_file    = false
    
    // modify evaluator:
    NOMAD::Evaluator           * old_ev = ev_control.get_evaluator();
    NOMAD::Phase_One_Evaluator * p1ev   = new NOMAD::Phase_One_Evaluator ( _p , *old_ev );
    ev_control.set_evaluator ( p1ev );
    
    // disable the Pareto front:
    NOMAD::Pareto_Front * old_pareto_front = mads.get_pareto_front();
    mads.set_pareto_front ( NULL );
    
    int old_eval = stats.get_eval();
    
    // run MADS with modified parameters:
    // ----------------------------------
    
    // C. Tribes  march 2014 ---- these flags are mads static and must be put back to their original value after running mads (see below)
    // get flags:
    bool flag_check_bimads , flag_reset_mesh , flag_reset_barriers , flag_p1_active;
    NOMAD::Mads::get_flags ( flag_check_bimads   ,
                            flag_reset_mesh     ,
                            flag_reset_barriers ,
                            flag_p1_active        );
    
    
    // set flags:
    NOMAD::Mads::set_flag_check_bimads   ( false );
    NOMAD::Mads::set_flag_reset_mesh     ( false );
    NOMAD::Mads::set_flag_p1_active      ( true  );
    NOMAD::Mads::set_flag_reset_barriers ( true  );
	   
    // run:
    stop_reason = mads.run();
    
    // reset stopping condition:
    if ( stop_reason == NOMAD::F_TARGET_REACHED )
    {
        
        // stop if feasible:
        if ( old_sif )
        {
            stop        = true;
            stop_reason = NOMAD::FEAS_REACHED;
        }
        
        // continue:
        else {
            stop        = false;
            stop_reason = NOMAD::NO_STOP;
        }
    }
    else
        stop = true;
    
    // reset flags to there previous state :
    // C. Tribes  march 2014 ---- these flags are mads static and must be put back to their original value after running mads
    NOMAD::Mads::set_flag_check_bimads ( flag_check_bimads  );
    NOMAD::Mads::set_flag_reset_mesh   ( flag_reset_mesh  );
    NOMAD::Mads::set_flag_p1_active    ( flag_p1_active );
    NOMAD::Mads::set_flag_reset_barriers    ( flag_reset_barriers );
    
    //	NOMAD::Mads::set_flag_check_bimads ( true  );
    //	NOMAD::Mads::set_flag_reset_mesh   ( true  );
    //	NOMAD::Mads::set_flag_p1_active    ( false );
    
    // number of search points:
    nb_search_pts = stats.get_eval() - old_eval;
    
    // restore evaluator:
    ev_control.set_evaluator ( old_ev );
    delete p1ev;
    
    // restore the Pareto front:
    mads.set_pareto_front ( old_pareto_front );
    
    // restore parameters:
    _p.set_VNS_SEARCH     ( old_VNS_trigger    );
    _p.set_F_TARGET       ( old_f_target       );
    _p.set_L_CURVE_TARGET ( old_lct            );
    _p.set_BB_OUTPUT_TYPE ( old_bbot           );
    _p.set_DISPLAY_DEGREE ( old_display_degree );
    _p.set_SOLUTION_FILE  ( old_sol_file       );
    _p.reset_stats_file();
    _p.set_STATS_FILE     ( old_stats_file_name , old_stats_file );
    _p.set_DISPLAY_STATS (old_ds);
    
    _p.check ( false ,    // remove_history_file  = false
              true  ,    // remove_solution_file = true
              false    ); // remove_stats_file    = true
    
    
    // counters:
    stats.add_p1_iterations ( stats.get_iterations() - old_it  );
    stats.add_p1_bbe        ( stats.get_bb_eval   () - old_bbe );
    
    // for the update of new_feas_inc and new_infeas_inc (1/2):
    const NOMAD::Barrier    & active_barrier           = mads.get_active_barrier();
    const NOMAD::Eval_Point * old_feasible_incumbent   = NULL;
    const NOMAD::Eval_Point * old_infeasible_incumbent = NULL;
    old_feasible_incumbent   = active_barrier.get_best_feasible();
    old_infeasible_incumbent = active_barrier.get_best_infeasible();
    
    // update the barriers and compute the true values
    // of f and h for all evaluated points:
    NOMAD::Barrier & true_barrier = mads.get_true_barrier();
    NOMAD::Barrier & sgte_barrier = mads.get_sgte_barrier();
    
    true_barrier.reset();
    sgte_barrier.reset();
    
    // scan the active cache:
    const NOMAD::Cache & active_cache = mads.get_cache();
    if ( active_cache.empty() )
    {
        stop        = true;
        stop_reason = NOMAD::P1_FAIL;
        return;
    }
    
    const NOMAD::Eval_Point * cur = active_cache.begin();
    while ( cur )
    {
        
        if ( cur->is_eval_ok() && cur->get_signature() )
        {
            
            NOMAD::Eval_Point * modifiable_x = &NOMAD::Cache::get_modifiable_point ( *cur );
            
            modifiable_x->set_direction          ( NULL                              );
            modifiable_x->set_poll_center_type   ( NOMAD::UNDEFINED_POLL_CENTER_TYPE );
            modifiable_x->set_user_eval_priority ( NOMAD::Double()                   );
            modifiable_x->set_rand_eval_priority ( NOMAD::Double()                   );
            
            old_ev->compute_f ( *modifiable_x );
            old_ev->compute_h ( *modifiable_x );
            
            // insertion in barrier:
            (( cur->get_eval_type() == NOMAD::TRUTH ) ? true_barrier : sgte_barrier).insert (*cur);
        }
        
        cur = active_cache.next();
    }
    
    success = active_barrier.get_success();
    
    if ( !stop && success == NOMAD::UNSUCCESSFUL )
    {
        stop        = true;
        stop_reason = NOMAD::P1_FAIL;
        return;
    }
    
    true_barrier.update_and_reset_success();
    sgte_barrier.update_and_reset_success();
    
    // update of new_feas_inc and new_infeas_inc (2/2):
    const NOMAD::Eval_Point * bf = active_barrier.get_best_feasible  ();
    const NOMAD::Eval_Point * bi = active_barrier.get_best_infeasible();
    if ( bf && bf != old_feasible_incumbent )
        new_feas_inc = bf;
    if ( bi && bi != old_infeasible_incumbent )
        new_infeas_inc = bi;
    
    // final displays:
    if ( bf && bf->get_current_run() )
    {
        
        // solution file:
        ev_control.write_solution_file ( *bf );
        
        // stats_file:
        const std::string & stats_file_name = _p.get_stats_file_name();
        if ( !stats_file_name.empty() && display_degree > NOMAD::NO_DISPLAY) 
            ev_control.stats_file ( stats_file_name , bf , true , NULL );
        
        // display_stats
        if (display_degree > NOMAD::NO_DISPLAY)
            ev_control.display_stats(false, out, old_ds , bf , true , NULL );
    }
    
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out << std::endl << NOMAD::close_block ( "end of phase one" );
}
