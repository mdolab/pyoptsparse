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
 \file   VNS_Search.cpp
 \brief  VNS search (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-12
 \see    VNS_Search.hpp
 */
#include "VNS_Search.hpp"

/*---------------------------------------------------------*/
/*                           reset                         */
/*---------------------------------------------------------*/
void NOMAD::VNS_Search::reset ( void )
{
    _k = _k_max   = 1;
    _old_x        = NULL;
    
}

/*---------------------------------------------------------*/
/*                       the search                        */
/*       VNS: x --[shaking(k)]--> x' --[descent]--> x"     */
/*---------------------------------------------------------*/
void NOMAD::VNS_Search::search ( NOMAD::Mads              & mads           ,
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
    count_search  = !stop;
    
    if ( stop )
        return;
    
    // initial display:
    const NOMAD::Display    & out = _p.out();
    NOMAD::dd_type display_degree = out.get_search_dd();
    if ( display_degree == NOMAD::FULL_DISPLAY )
    {
        std::ostringstream oss;
        oss << NOMAD::VNS_SEARCH;
        out << std::endl << NOMAD::open_block ( oss.str() ) << std::endl;
    }
    
    bool opt_only_sgte = _p.get_opt_only_sgte();
    
    // the barriers:
    NOMAD::Barrier       & true_barrier   = mads.get_true_barrier();
    NOMAD::Barrier       & sgte_barrier   = mads.get_sgte_barrier();
    const NOMAD::Barrier & active_barrier = mads.get_active_barrier();
    
    // point x:
    NOMAD::Double             best_f;
    bool                      x_feas = true;
    const NOMAD::Eval_Point * x      = active_barrier.get_best_feasible();
    if ( x )
        best_f = x->get_f();
    else
    {
        x      = active_barrier.get_best_infeasible();
        x_feas = false;
    }
    
    if ( !x )
    {
        if ( display_degree == NOMAD::FULL_DISPLAY )
            out.close_block ( "end of VNS search (no incumbent)" );
        return;
    }
    
    // update _k and _old_x:
    if ( x == _old_x )
    {
        ++_k;
        if ( _k > _k_max )
            _k_max = _k;
    }
    else
        _k = 1;
    
    _old_x = x;
    
    // get the signature:
    NOMAD::Signature * signature = x->get_signature();
    if ( !signature )
    {
        if ( display_degree == NOMAD::FULL_DISPLAY )
            out.close_block ( "end of VNS search (no signature)" );
        return;
    }
    
    int n = signature->get_n();
    if ( n != x->size() )
    {
        if ( display_degree == NOMAD::FULL_DISPLAY )
            out.close_block ( "end of VNS search (incompatible signature)" );
        return;
    }
    
    // shaking: get ONE direction from the signature:
    NOMAD::Direction dir;
    
    signature->get_one_direction ( dir , _k - 1);
    
    
    // shaking: construct x':
    NOMAD::Point xp = *x + dir;
    
    // shaking: the perturbation is tried twice with dir and -dir
    //          (in case x == x + dir after snapping)
    for ( int nbt = 0 ; nbt < 2 ; ++nbt )
    {
        
        // treat xp: periodic variables or bounds:
        if ( _p.has_periodic_variables() )
        {
            NOMAD::Direction * tmp_dir = NULL;
            signature->treat_periodic_variables ( xp , NULL , tmp_dir );
        }
        else
            signature->snap_to_bounds ( xp , NULL );
        
        if ( xp == *x )
        {
            
            // no third try: the search fails
            if ( nbt == 1 )
            {
                if ( display_degree == NOMAD::FULL_DISPLAY )
                    out.close_block ( "end of VNS search (shaking failed)" );
                return;
            }
            
            // 2nd try (-dir instead of dir):
            xp = *x - dir;
        }
    }
    
    // Current mesh indices
    const NOMAD::Point         old_mesh_indices = signature->get_mesh()->get_mesh_indices ( );
    const NOMAD::Point            old_delta_min = signature->get_mesh()->get_min_mesh_size();
    
    
    // stats:
    NOMAD::Stats & stats = mads.get_stats();
    
    // current number of blackbox evaluations:
    int  bbe             = stats.get_bb_eval();
    int  blk_eva		 = stats.get_block_eval();
    int  sgte_eval       = stats.get_sgte_eval();
    int  mads_iterations = stats.get_iterations();
    bool has_sgte        = _p.has_sgte();
    
    // displays:
    if ( display_degree == NOMAD::FULL_DISPLAY )
    {
        out << "        it = " << mads_iterations << std::endl
        << "       bbe = " << bbe << std::endl
        << "   blk_eva = " << blk_eva << std::endl;
        if ( has_sgte )
            out << " sgte_eval = " << sgte_eval << std::endl;
        out << "mesh_indices = ( " << old_mesh_indices << " ) " << std::endl
        << "         k = " << _k << std::endl
        << "      kmax = " << _k_max << std::endl
        << "         x = ( ";
        x->Point::display ( out , " " , 5 , _p.get_point_display_limit() );
        out << " ) f=" << x->get_f() << " h=" << x->get_h() << std::endl
        << "       dir = ( ";
        dir.Point::display ( out , " " , 5 , _p.get_point_display_limit() );
        out << " ) |dir|=";
        NOMAD::Double norm = dir.norm();
        out << norm << std::endl;
        out << "        xp = ( ";
        xp.display ( out , " " , 5 , _p.get_point_display_limit() );
        out << " )" << std::endl << std::endl;
        out << "bb_eval (before+VNS only) objective_value"
        << std::endl << std::endl;
    }
    
    // save parameters that are going to be modified:
    // ----------------------------------------------
    std::string old_display_degree;
    _p.out().get_display_degree ( old_display_degree );
    
    NOMAD::model_params_type old_mp;
    _p.get_model_parameters ( old_mp );
    
    bool                                old_ses = _p.get_sgte_eval_sort();
    bool                                old_sif = _p.get_stop_if_feasible();
    int                            old_max_time = _p.get_max_time();
    int                             old_max_bbe = _p.get_max_bb_eval();
    int                            old_max_eval = _p.get_max_eval();
    int                       old_max_sgte_eval = _p.get_max_sgte_eval();
    int                              old_max_it = _p.get_max_iterations();
    int                             old_max_cfi = _p.get_max_consecutive_failed_iterations();
    int                               old_LH_p0 = _p.get_LH_search_p0();
    int                               old_LH_pi = _p.get_LH_search_pi();
    bool                             old_opp_LH = _p.get_opportunistic_LH();
    bool                                 old_CS = _p.get_cache_search();
    bool                             old_opp_CS = _p.get_opportunistic_cache_search();
    int                             old_max_sbe = _p.get_max_sim_bb_eval();
    NOMAD::Double                       old_sst = _p.get_stat_sum_target();
    NOMAD::Double                       old_lct = _p.get_L_curve_target();
    NOMAD::Double                   old_trigger = _p.get_VNS_trigger();
    NOMAD::Point                         old_ft = _p.get_f_target();
    const std::list<std::string>         old_ds = _p.get_display_stats();
    const std::list<std::string> old_stats_file = _p.get_stats_file();
    const std::string       old_stats_file_name = _p.get_stats_file_name();
    const std::string              old_sol_file = _p.get_solution_file();
    const std::string              old_his_file = _p.get_history_file();
    bool                                old_uce = _p.get_user_calls_enabled();
    bool                                old_epe = _p.get_extended_poll_enabled();
    const std::vector<NOMAD::bb_output_type> old_bbot = _p.get_bb_output_type();
    
    // save list of starting points:
    std::string x0_cache_file = _p.get_x0_cache_file();
    std::vector<NOMAD::Point *> x0s;
    {
        const std::vector<NOMAD::Point *> & x0s_tmp = _p.get_x0s();
        size_t nx0 = x0s_tmp.size() , k;
        for ( k = 0 ; k < nx0 ; ++k )
            x0s.push_back ( new Point ( *x0s_tmp[k] ) );
    }
    
    // modify parameters:
    // ------------------
    _p.set_DISPLAY_DEGREE(NOMAD::NO_DISPLAY);
    
    _p.set_SOLUTION_FILE  ( ""       );
    _p.set_LH_SEARCH      ( 0 , 0    );
    _p.set_VNS_SEARCH     ( false    );
    _p.set_CACHE_SEARCH   ( false    );
    _p.set_MAX_ITERATIONS ( -1       );
    _p.set_MAX_CONSECUTIVE_FAILED_ITERATIONS ( -1 );
    _p.reset_X0();
    _p.reset_stats_file();
    
    
    if ( has_sgte )
    {
        _p.set_OPT_ONLY_SGTE   ( true            );
        _p.set_MODEL_SEARCH    ( NOMAD::NO_MODEL );
        _p.set_MODEL_EVAL_SORT ( NOMAD::NO_MODEL );
    }
    
    _p.set_USER_CALLS_ENABLED    ( false );
    _p.set_EXTENDED_POLL_ENABLED ( false );
    
    // DISPLAY_STATS:
    {
        if ( has_sgte )
            _p.set_DISPLAY_STATS ( NOMAD::itos(sgte_eval) + "+SGTE OBJ (VNS--surrogate)" );
        else
        {
            std::list<std::string>                 ds    = old_ds;
            std::list<std::string>::iterator       it    = ds.begin();
            std::list<std::string>::const_iterator end   = ds.end();
            std::string                            s_bbe = NOMAD::itos(bbe) + "+";
            std::string                            s_blk = NOMAD::itos(blk_eva) + "+";
            while ( it != end )
            {
                if ( *it == "BBE" )
                    ds.insert ( it , s_bbe );
                if ( *it == "BLK_EVA" )
                    ds.insert ( it , s_blk );
                ++it;
            }
            ds.push_back ( " (VNS)" );
            _p.set_DISPLAY_STATS ( ds );
        }
    }
    
    // STATS_FILE:
    {
        std::list<std::string>                 sf    = old_stats_file;
        std::list<std::string>::iterator       it    = sf.begin();
        std::list<std::string>::const_iterator end   = sf.end();
        std::string                            s_bbe = NOMAD::itos(bbe) + "+";
        std::string                            s_blk	= NOMAD::itos(blk_eva) + "+";
        while ( it != end )
        {
            if ( *it == "BBE" )
                sf.insert ( it , s_bbe );
            if ( *it == "BLK_EVA" )
                sf.insert ( it , s_blk );
            ++it;
        }
        sf.push_back ( " (VNS)" );
        
        _p.set_STATS_FILE( old_stats_file_name, sf);
    }
    
    
    // Mesh size at current mads iterate can be used as termination criterion for vns search.
    // Mesh indices are reinitialized during p.check()
    NOMAD::Point delta=signature->get_mesh()->get_delta ( );
    signature->get_mesh()->set_min_mesh_sizes( delta );
    
    
    // X0:
    _p.set_EXTERN_SIGNATURE ( signature );
    _p.set_X0 ( xp );
    
    // MAX_BB_EVAL:
    if ( old_max_bbe < 0 )
        _p.set_MAX_BB_EVAL ( 100 * n );
    else
        _p.set_MAX_BB_EVAL ( old_max_bbe - bbe );
    
    // MAX_SGTE_EVAL:
    if ( old_max_sgte_eval > 0 )
        _p.set_MAX_SGTE_EVAL ( old_max_sgte_eval - sgte_eval );
    
    // MAX_EVAL:
    if ( old_max_eval > 0 )
        _p.set_MAX_EVAL ( old_max_eval - stats.get_eval() );
    
    // MAX_SIM_BB_EVAL:
    if ( old_max_sbe > 0 )
        _p.set_MAX_SIM_BB_EVAL ( old_max_sbe - stats.get_sim_bb_eval() );
    
    // STAT_SUM_TARGET:
    if ( old_sst.is_defined() )
        _p.set_STAT_SUM_TARGET ( old_sst - stats.get_stat_sum() );
    
    // MAX_TIME:
    if ( old_max_time > 0 )
        _p.set_MAX_TIME ( old_max_time - stats.get_real_time() );
    
    // L_CURVE_TARGET:
    if ( !has_sgte )
        _p.set_L_CURVE_TARGET ( best_f );
    
    // F_TARGET and STOP_IF_FEASIBLE:
    if ( has_sgte )
    {
        _p.reset_f_target();
        _p.set_STOP_IF_FEASIBLE ( false );
    }
    
    // check the parameters:
    _p.check ( false ,    // remove_history_file  = false
              false ,    // remove_solution_file = false
              false   ); // remove_stats_file    = false
    
    // Evaluator_Control:
    NOMAD::Evaluator_Control & ev_control = mads.get_evaluator_control();
    
    // descent: run MADS:
    // ------------------
    NOMAD::Mads VNS_mads ( _p                           ,
                          ev_control.get_evaluator  () ,
                          mads.get_extended_poll    () ,
                          &ev_control.get_cache     () ,
                          &ev_control.get_sgte_cache()   );
    NOMAD::Mads::set_flag_reset_mesh     ( false );
    NOMAD::Mads::set_flag_reset_barriers ( true  );
    
    NOMAD::stop_type st = VNS_mads.run();
    
    NOMAD::Mads::set_flag_reset_mesh ( true );
    
    // update stats:
    {
        const NOMAD::Stats & VNS_stats = VNS_mads.get_stats();
        stats.update            ( VNS_stats , true ); // for_search = true
        stats.add_VNS_bb_eval   ( VNS_stats.get_bb_eval  () );
        stats.add_VNS_sgte_eval ( VNS_stats.get_sgte_eval() );
    }
    
    // check MADS stopping criteria:
    if ( st == NOMAD::CTRL_C                   ||
        st == NOMAD::ERROR                    ||
        st == NOMAD::UNKNOWN_STOP_REASON      ||
        st == NOMAD::FEAS_REACHED             ||
        st == NOMAD::MAX_CACHE_MEMORY_REACHED ||
        st == NOMAD::STAT_SUM_TARGET_REACHED  ||
        st == NOMAD::MAX_SGTE_EVAL_REACHED    ||
        st == NOMAD::F_TARGET_REACHED         ||
        st == NOMAD::MAX_SIM_BB_EVAL_REACHED  ||
        st == NOMAD::MAX_TIME_REACHED         ||
        (st == NOMAD::MAX_BB_EVAL_REACHED && old_max_bbe > 0 ) )
    {
        stop_reason = st;
        stop        = true;
    }
    
    // Pareto front:
    NOMAD::Pareto_Front * pareto_front = mads.get_pareto_front();
    
    // restore starting points:
    {
        _p.reset_X0();
        size_t nx0 = x0s.size();
        
        if ( nx0 > 0 )
        {
            
            for ( size_t k = 0 ; k < nx0 ; ++k )
            {
                _p.set_X0 ( *x0s[k] );
                delete x0s[k];
            }
        }
        else if ( !x0_cache_file.empty() )
            _p.set_X0 ( x0_cache_file );
    }
    
    // restore other saved parameters:
    _p.set_model_parameters           ( old_mp                               );
    _p.set_USER_CALLS_ENABLED         ( old_uce                              );
    _p.set_EXTENDED_POLL_ENABLED      ( old_epe                              );
    _p.set_VNS_SEARCH                 ( old_trigger                          );
    _p.set_F_TARGET                   ( old_ft                               );
    _p.set_STOP_IF_FEASIBLE           ( old_sif                              );
    _p.set_L_CURVE_TARGET             ( old_lct                              );
    _p.set_DISPLAY_DEGREE             ( old_display_degree                   );
    _p.set_DISPLAY_STATS              ( old_ds                               );
    _p.set_STATS_FILE                 ( old_stats_file_name , old_stats_file );
    _p.set_SOLUTION_FILE              ( old_sol_file                         );
    _p.set_MAX_BB_EVAL                ( old_max_bbe                          );
    _p.set_MAX_EVAL                   ( old_max_eval                         );
    _p.set_MAX_SGTE_EVAL              ( old_max_sgte_eval                    );
    _p.set_MAX_ITERATIONS             ( old_max_it                           );
    _p.set_MAX_CONSECUTIVE_FAILED_ITERATIONS ( old_max_cfi                   );
    _p.set_STAT_SUM_TARGET            ( old_sst                              );
    _p.set_LH_SEARCH                  ( old_LH_p0 , old_LH_pi                );
    _p.set_OPPORTUNISTIC_LH           ( old_opp_LH                           );
    _p.set_CACHE_SEARCH               ( old_CS                               );
    _p.set_OPPORTUNISTIC_CACHE_SEARCH ( old_opp_CS                           );
    _p.set_MAX_SIM_BB_EVAL            ( old_max_sbe                          );
    _p.set_MAX_TIME                   ( old_max_time                         );
    _p.set_SGTE_EVAL_SORT             ( old_ses                              );
    _p.set_OPT_ONLY_SGTE              ( opt_only_sgte                        );
    _p.set_BB_OUTPUT_TYPE			  ( old_bbot							 );
    
    
    _p.check ( false ,    // remove_history_file  = false
              false ,    // remove_solution_file = false
              false   ); // remove_stats_file    = false
    
    // restore min mesh sizes and mesh indices
    signature->get_mesh()->set_min_mesh_sizes( old_delta_min );
    signature->get_mesh()->set_mesh_indices( old_mesh_indices ); // Needed because mesh indices reinitialized during p.check()
    
    
    // surrogate evaluations: perform only one true evaluation:
    if ( has_sgte && !opt_only_sgte )
    {
        
        if ( !stop )
        {
            
            // remember old best surrogates incumbents:
            const NOMAD::Eval_Point * old_sgte_bf = sgte_barrier.get_best_feasible  ();
            const NOMAD::Eval_Point * old_sgte_bi = sgte_barrier.get_best_infeasible();
            
            // update the surrogate barrier
            // (no need to invoke Evaluator_Control::process_barrier_points() here
            //  since only surrogate evaluations have been made):
            sgte_barrier.insert ( VNS_mads.get_sgte_barrier() );
            NOMAD::success_type sgte_succ = sgte_barrier.get_success();
            sgte_barrier.update_and_reset_success();
            
            // we generate only a true trial point if the
            // surrogates improved the surrogate barrier:
            if ( sgte_succ != NOMAD::UNSUCCESSFUL )
            {
                
                // choose the best surrogate point(s) where to evaluate the true function:
                const NOMAD::Eval_Point * sgte_bf = sgte_barrier.get_best_feasible  ();
                const NOMAD::Eval_Point * sgte_bi = sgte_barrier.get_best_infeasible();
                
                std::list<const NOMAD::Eval_Point *> candidates;
                
                if ( sgte_bf && ( !x_feas || sgte_bf != old_sgte_bf ) )
                    candidates.push_back ( sgte_bf );
                
                if ( sgte_bi && sgte_bi != old_sgte_bi )
                    candidates.push_back ( sgte_bi );
                
                // generate the new trial points:
                NOMAD::Eval_Point * sk;
                std::list<const NOMAD::Eval_Point *>::const_iterator
                it , end = candidates.end();
                for ( it = candidates.begin() ; it != end ; ++it )
                {
                    
                    // display:
                    if ( display_degree == NOMAD::FULL_DISPLAY )
                        out << std::endl << "VNS surrogate candidate: "
                        << **it << std::endl;
                    
                    sk = new NOMAD::Eval_Point;
                    sk->set ( n , _p.get_bb_nb_outputs() );
                    sk->set_signature  ( signature   );
                    sk->Point::operator = ( **it );
                    
                    // add the new point to the list of search trial points:
                    ev_control.add_eval_point ( sk                      ,
                                               display_degree          ,
                                               _p.get_snap_to_bounds() ,
                                               NOMAD::Double()         ,
                                               NOMAD::Double()         ,
                                               NOMAD::Double()         ,
                                               NOMAD::Double()           );
                }
                
                // eval_list_of_points:
                // --------------------
                success = NOMAD::UNSUCCESSFUL;
                new_feas_inc = new_infeas_inc = NULL;
                
                ev_control.eval_list_of_points ( _type          ,
                                                true_barrier   ,
                                                sgte_barrier   ,
                                                pareto_front   ,
                                                stop           ,
                                                stop_reason    ,
                                                new_feas_inc   ,
                                                new_infeas_inc ,
                                                success          );
                
                // number of search points (0 or 1 or 2):
                nb_search_pts = static_cast<int> ( candidates.size() );
            }
        }
    }
    
    // true evaluations (or surrogate evaluations if opt_only_sgte==true):
    else 
    {
        
        // for the update of new_feas_inc and new_infeas_inc (1/2):
        const NOMAD::Eval_Point * old_feasible_incumbent   =
        active_barrier.get_best_feasible();
        const NOMAD::Eval_Point * old_infeasible_incumbent =
        active_barrier.get_best_infeasible();
        
        // update barriers and process VNS search points:
        NOMAD::success_type sgte_succ
        = ev_control.process_barrier_points ( sgte_barrier                ,
                                             VNS_mads.get_sgte_barrier() ,
                                             pareto_front                ,
                                             display_degree              ,
                                             NOMAD::VNS_SEARCH             );
        NOMAD::success_type true_succ
        = ev_control.process_barrier_points ( true_barrier                ,
                                             VNS_mads.get_true_barrier() ,
                                             pareto_front                ,
                                             display_degree              ,
                                             NOMAD::VNS_SEARCH             );
        
        // update of new_feas_inc and new_infeas_inc (2/2):
        const NOMAD::Eval_Point * bf = active_barrier.get_best_feasible  ();
        const NOMAD::Eval_Point * bi = active_barrier.get_best_infeasible();
        if ( bf && bf != old_feasible_incumbent )
            new_feas_inc = bf;
        
        if ( bi && bi != old_infeasible_incumbent )
        {
            new_infeas_inc = bi;
            // check the PEB constraints: if we have a new best infeasible
            // incumbent from another infeasible incumbent
            // ( active_barrier.check_PEB_constraints() ):
            if ( _p.get_barrier_type() == NOMAD::PEB_P )
                ( ( _p.get_opt_only_sgte() ) ? sgte_barrier : true_barrier ).check_PEB_constraints ( *new_infeas_inc , display_degree==NOMAD::FULL_DISPLAY );
        }		
        
        // number of search points and success:
        if ( opt_only_sgte )
        {
            nb_search_pts = VNS_mads.get_stats().get_sgte_eval();
            success       = sgte_succ;
        }
        else 
        {
            nb_search_pts = VNS_mads.get_stats().get_eval();
            success       = true_succ;
        }
        
        // solution file:
        if ( bf )
            ev_control.write_solution_file ( *bf );
        
    }
    
    // final display:
    if ( display_degree == NOMAD::FULL_DISPLAY )
    {
        std::ostringstream oss;
        oss << "end of VNS search (" << success << ")";
        out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
    }
}
