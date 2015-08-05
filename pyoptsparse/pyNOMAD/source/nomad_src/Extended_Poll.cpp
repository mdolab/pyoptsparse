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
 \file   Extended_Poll.cpp
 \brief  Extended poll for categorical variables (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-14
 \see    Extended_Poll.hpp
 */
#include "Extended_Poll.hpp"

/*----------------------------------------------------------------*/
/*                            destructor                          */
/*----------------------------------------------------------------*/
NOMAD::Extended_Poll::~Extended_Poll ( void )
{
    std::set<NOMAD::Signature_Element>::const_iterator it , end = _signatures.end();
    for ( it = _signatures.begin() ; it != end ; ++it )
        delete (*it).get_signature();
    poll_reset();
}

/*----------------------------------------------------------------*/
/*                               reset                            */
/*----------------------------------------------------------------*/
void NOMAD::Extended_Poll::reset ( void )
{
    // successful directions:
    std::set<NOMAD::Signature_Element>::const_iterator it , end = _signatures.end();
    for ( it = _signatures.begin() ; it != end ; ++it )
    {
        (*it).get_signature()->reset_feas_success_dir();
        (*it).get_signature()->reset_infeas_success_dir();
    }
    
    // poll_reset:
    poll_reset();
}

/*----------------------------------------------------------------*/
/*         poll reset: before the extended poll is launched       */
/*----------------------------------------------------------------*/
void NOMAD::Extended_Poll::poll_reset ( void )
{
    _poll_signatures.clear();
    std::list<NOMAD::Eval_Point *>::const_iterator end = _extended_points.end();
    for ( std::list<NOMAD::Eval_Point *>::iterator  it = _extended_points.begin() ;
         it != end ; ++it )
        if ( !(*it)->is_in_cache() )
            delete *it;
    _extended_points.clear();
}

/*----------------------------------------------------------------*/
/*  get, check and register the extended point and its signature  */
/*  created by the user in construct_extended_points()            */
/*----------------------------------------------------------------*/
void NOMAD::Extended_Poll::add_extended_poll_point ( NOMAD::Point     & ep ,
                                                    NOMAD::Signature & s    )
{
    // new signature:
    // --------------
    NOMAD::Signature * new_s = new NOMAD::Signature (s);
    
    // . 's' can be standard, but its copy 'new_s' is not
    
    // . a standard signature will never be inserted into _signatures,
    //   since standard signatures are handled and deleted by the Parameters class
    
    {
        // signature already registered ?
        NOMAD::Signature_Element se ( new_s );
        std::set<NOMAD::Signature_Element>::const_iterator it = _signatures.find ( se );
        
        // signature already registered:
        if ( it != _signatures.end() )
        {
            
            // success directions eventually included in new_s are not considered
            // since new_s is the copy of s, which is user provided
            
            delete new_s;
            new_s = it->get_signature();
        }
        
        // new signature to register:
        else
            _signatures.insert ( se );
        
        _poll_signatures.insert ( NOMAD::Signature_Element ( new_s ) );
    }
    
    // new eval point:
    // ---------------
    NOMAD::Eval_Point * pt  = new NOMAD::Eval_Point;
    pt->set				 ( ep , _p.get_bb_nb_outputs() );
    pt->set_signature	 ( new_s                       );
    
    for ( int i = 0 ; i < pt->get_n() ; ++i )
    {
        if (  (pt->get_signature()->get_input_type())[i] != NOMAD::CONTINUOUS && ! (*pt)[i].is_integer() )
            
            throw NOMAD::Exception ( "Extended_Poll.cpp" , __LINE__ ,
                                    "NOMAD::Extended_Poll::add_extended_points(): the categorical variables of the added point must be an integer." );
    }
    
    
    _extended_points.push_back ( pt );
}

/*----------------------------------------------------------------*/
/*                 check the extended poll trigger                */
/*----------------------------------------------------------------*/
/*  . return true if the extended poll has to be performed        */
/*  . private method                                              */
/*----------------------------------------------------------------*/
bool NOMAD::Extended_Poll::check_trigger ( const NOMAD::Eval_Point * old_bf ,
                                          const NOMAD::Eval_Point * old_bi ,
                                          const NOMAD::Eval_Point * y        ) const
{
    if ( !y->is_in_cache()        ||
        !y->is_eval_ok()         ||
        !y->get_f().is_defined() ||
        !y->get_h().is_defined()    )
        return false;
    
    // y is feasible:
    // --------------
    if ( y->is_feasible ( _p.get_h_min() ) )
    {
        
        if ( !old_bf )
            return true;
        
        return check_trigger_on_f ( old_bf->get_f() , y->get_f() );
    }
    
    // y is infeasible:
    // ----------------
    if ( !old_bf && !old_bi )
        return true;
    
    if ( !old_bf )
        return ( y->get_h() < old_bi->get_h() );
    
    if ( !old_bi )
        return check_trigger_on_f ( old_bf->get_f() , y->get_f() );
    
    if ( y->get_h() >= old_bi->get_h() )
        return false;
    
    // y is infeasible, and old best feasible and best infeasible solutions are
    // available: the extended poll will be performed if the y point in the
    // (h,f) space is below the line joining [ h(old_bf) , f(old_bf)+trigger ]
    // to [ h(old_bi) , f(old_bi)+trigger ] :
    const NOMAD::Double & hA  = old_bf->get_h();
    NOMAD::Double         fA  = old_bf->get_f();
    
    const NOMAD::Double & hB  = old_bi->get_h();
    NOMAD::Double         fB  = old_bi->get_f();
    
    const NOMAD::Double & hy  = y->get_h();
    const NOMAD::Double & fy  = y->get_f();
    
    const NOMAD::Double & ept = _p.get_extended_poll_trigger();
    
    if ( _p.get_relative_ept() && fA != 0.0 && fB != 0.0 && fy != 0.0 )
    {
        fA = fA + fA.abs() * ept;
        fB = fB + fB.abs() * ept;
    }
    else
    {
        fA = fA + ept;
        fB = fB + ept;
    }
    
    // line joining [h(A),f(A)] to [h(B),f(B)]: f=a*h+b :
    NOMAD::Double a = (fA-fB) / (hA-hB);
    NOMAD::Double b = fA - a * hA;
    
    return fy < a*hy + b;
}

/*-------------------------------------------------------------------*/
/*  check only the f values for the extended poll trigger (private)  */
/*-------------------------------------------------------------------*/
bool NOMAD::Extended_Poll::check_trigger_on_f  ( const NOMAD::Double & old_f ,
                                                const NOMAD::Double & new_f   ) const
{
    if ( new_f <= old_f )
        return true;
    
    // relative comparison (both values are != 0):
    if ( _p.get_relative_ept() && old_f != 0.0 && new_f != 0.0 )
        return ( new_f < old_f + old_f.abs() * _p.get_extended_poll_trigger() );
    
    // absolute comparison:
    return ( new_f < old_f + _p.get_extended_poll_trigger() );
}

/*----------------------------------------------------------------*/
/*          descent from the extended poll center (private)       */
/*----------------------------------------------------------------*/
void NOMAD::Extended_Poll::descent ( const NOMAD::Eval_Point  * y               ,
                                    Mads                     & mads            ,
                                    int                      & nb_ext_poll_pts ,
                                    bool                     & stop            ,
                                    NOMAD::stop_type         & stop_reason     ,
                                    NOMAD::success_type      & success         ,
                                    const NOMAD::Eval_Point *& new_feas_inc    ,
                                    const NOMAD::Eval_Point *& new_infeas_inc    )
{
    
    
    NOMAD::Stats         & stats          = mads.get_stats();
    bool                   has_sgte       = _p.has_sgte();
    bool                   opt_only_sgte  = _p.get_opt_only_sgte();
    const NOMAD::Display & out            = _p.out();
    NOMAD::dd_type         display_degree = out.get_poll_dd();
    
    // get the signature:
    NOMAD::Signature * signature = y->get_signature();
    
    
    // displays:
    if ( display_degree == NOMAD::FULL_DISPLAY )
    {
        
        std::ostringstream oss;
        oss << NOMAD::EXTENDED_POLL << " descent";
        if ( has_sgte )
            oss << " (on surrogates)";
        
        out << std::endl << NOMAD::open_block ( oss.str() ) << std::endl;
        
        out << "       iteration = " << stats.get_iterations() << std::endl
        << " blackbox eval.  = " << stats.get_bb_eval() << std::endl;
        if ( has_sgte )
            out << "      sgte eval. = " << stats.get_sgte_eval() << std::endl;
        out << "     mesh indices = (" << signature->get_mesh()->get_mesh_indices() << " )" << std::endl
        << "ext. poll center = ( ";
        y->Point::display ( out , " " , 5 , _p.get_point_display_limit() );
        out << " ) f=" << y->get_f() << " h=" << y->get_h() << std::endl << std::endl;
    }
    
    
    
    // create the descent parameters:
    NOMAD::Parameters descent_p ( signature , _p.out() );
    set_descent_parameters ( y , stats , descent_p );
    
    
    // mesh indices before the descent     --> restore after descent
    NOMAD::Point           mesh_indices     = signature->get_mesh()->get_mesh_indices();
    
    // limit mesh index before the descent --> restore after descent
    int limit_mesh_index = signature->get_mesh()->get_limit_mesh_index();
    
    
    // Set mesh indices to 0
    NOMAD::Point delta( signature->get_n(), 0 );
    descent_p.get_signature()->get_mesh()->set_mesh_indices( delta );
    
    // Use best_feasible or best_infeasible limit mesh index as a termination criterion for descent
    const NOMAD::Eval_Point * old_bf = mads.get_best_feasible();
    const NOMAD::Eval_Point * old_bi = mads.get_best_infeasible();
    int l1=0,l2=0;
    if ( old_bf )
        l1 = static_cast<int>((old_bf->get_signature()->get_mesh()->get_min_mesh_indices())[0].value());  // index same for all variables when using categorical variable
    else if ( old_bi )
        l2 = static_cast<int>((old_bi->get_signature()->get_mesh()->get_min_mesh_indices())[0].value());
    descent_p.get_signature()->get_mesh()->set_limit_mesh_index( min( l1, l2) );
    
    
    
    // Evaluator_Control object:
    NOMAD::Evaluator_Control & ev_control = mads.get_evaluator_control();
    
    // descent: run MADS:
    // ------------------
    NOMAD::Mads EP_mads ( descent_p                    ,
                         ev_control.get_evaluator  () ,
                         NULL                         ,
                         &ev_control.get_cache     () ,
                         &ev_control.get_sgte_cache()   );
#ifdef DEBUG
    out << std::endl << NOMAD::open_block ( "MADS run (ext. poll)" ) << std::endl;
#endif
    
    NOMAD::Mads::set_flag_reset_barriers ( true  );
    NOMAD::Mads::set_flag_reset_mesh     ( false );
    
    NOMAD::stop_type st = EP_mads.run();
    
    NOMAD::Mads::set_flag_reset_mesh ( true );
    
#ifdef DEBUG
    out << std::endl << NOMAD::close_block ( "end of run (ext. poll)" ) << std::endl;
#endif
    
    // Restore mesh indices and limit mesh index (termination criterion)
    signature->get_mesh()->set_mesh_indices( mesh_indices );
    signature->get_mesh()->set_limit_mesh_index( limit_mesh_index );
    
    
    // update stats:
    const NOMAD::Stats & EP_stats = EP_mads.get_stats();
    stats.update               ( EP_stats , true ); // for_search = true
    stats.add_ext_poll_bb_eval ( EP_stats.get_bb_eval() );
    stats.add_ext_poll_descent ();
    
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
        st == NOMAD::MAX_BB_EVAL_REACHED         )
    {
        stop_reason = st;
        stop        = true;
    }
    
    // Pareto front:
    NOMAD::Pareto_Front * pareto_front = mads.get_pareto_front();
    
    // the barriers:
    NOMAD::Barrier & true_barrier = mads.get_true_barrier();
    NOMAD::Barrier & sgte_barrier = mads.get_sgte_barrier();
    
    // surrogate evaluations: perform at most one true evaluation:
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
            sgte_barrier.insert ( EP_mads.get_sgte_barrier() );
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
                
                if ( sgte_bf && ( !y->is_feasible(_p.get_h_min()) || sgte_bf != old_sgte_bf ) )
                    candidates.push_back ( sgte_bf );
                
                if ( sgte_bi && sgte_bi != old_sgte_bi )
                    candidates.push_back ( sgte_bi );
                
                // generate the new trial points:
                NOMAD::Eval_Point * sk;
                std::list<const NOMAD::Eval_Point *>::const_iterator it , end = candidates.end();
                for ( it = candidates.begin() ; it != end ; ++it )
                {
                    
                    // display:
                    if ( display_degree == NOMAD::FULL_DISPLAY )
                        out << std::endl << "ext. poll surrogate candidate: " << **it << std::endl;
                    
                    sk = new NOMAD::Eval_Point;
                    sk->set ( signature->get_n() , _p.get_bb_nb_outputs() );
                    sk->set_signature  ( signature   );
                    sk->Point::operator = ( **it );
                    
                    // add the new point to the list of trial points:
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
                
                ev_control.eval_list_of_points ( NOMAD::EXTENDED_POLL ,
                                                true_barrier         ,
                                                sgte_barrier         ,
                                                pareto_front         ,
                                                stop                 ,
                                                stop_reason          ,
                                                new_feas_inc         ,
                                                new_infeas_inc       ,
                                                success                );
                
                // number of search points (0 or 1 or 2):
                nb_ext_poll_pts += static_cast<int> ( candidates.size() );
            }
        }
    }
    
    // true evaluations (or surrogate evaluations if opt_only_sgte==true):
    else
    {
        
        const NOMAD::Barrier & active_barrier = mads.get_active_barrier();
        
        // for the update of new_feas_inc and new_infeas_inc (1/2):
        const NOMAD::Eval_Point * old_feasible_incumbent
        = active_barrier.get_best_feasible();
        const NOMAD::Eval_Point * old_infeasible_incumbent
        = active_barrier.get_best_infeasible();
        
        // update barriers and process extended poll points:
        NOMAD::success_type sgte_succ
        = ev_control.process_barrier_points ( sgte_barrier               ,
                                             EP_mads.get_sgte_barrier() ,
                                             pareto_front               ,
                                             display_degree             ,
                                             NOMAD::EXTENDED_POLL         );
        
        NOMAD::success_type true_succ
        = ev_control.process_barrier_points ( true_barrier               ,
                                             EP_mads.get_true_barrier() ,
                                             pareto_front               ,
                                             display_degree             ,
                                             NOMAD::EXTENDED_POLL         );
        
        // update of new_feas_inc and new_infeas_inc (2/2):
        const NOMAD::Eval_Point * bf = active_barrier.get_best_feasible  ();
        const NOMAD::Eval_Point * bi = active_barrier.get_best_infeasible();
        if ( bf && bf != old_feasible_incumbent )
            new_feas_inc = bf;
        if ( bi && bi != old_infeasible_incumbent )
            new_infeas_inc = bi;
        
        // number of extended poll points and success:
        if ( opt_only_sgte )
        {
            nb_ext_poll_pts += EP_mads.get_stats().get_sgte_eval();
            success          = sgte_succ;
        }
        else
        {
            nb_ext_poll_pts += EP_mads.get_stats().get_eval();
            success          = true_succ;
        }
    }
    
    // final display:
    if ( display_degree == NOMAD::FULL_DISPLAY )
    {
        std::ostringstream oss;
        oss << "end of " << NOMAD::EXTENDED_POLL << " descent (" << success << ")";
        out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
    }
}

/*----------------------------------------------------------------*/
/*               create the descent parameters (private)          */
/*----------------------------------------------------------------*/
void NOMAD::Extended_Poll::set_descent_parameters
( const NOMAD::Eval_Point * y         ,
 const NOMAD::Stats      & stats     ,
 NOMAD::Parameters       & descent_p   ) const
{
    
    // extended poll center signature
    // (will be the temporary standard signature):
    NOMAD::Signature * epc_signature = y->get_signature();
    if ( !epc_signature )
        throw NOMAD::Exception ( "Extended_Poll.cpp" , __LINE__ ,
                                "Extended_Poll::set_descent_parameters(): ext. poll center has no signature" );
    
    // we set all the parameters:
    descent_p.set_DIMENSION         ( epc_signature->get_n()                  );
    descent_p.set_BB_INPUT_TYPE     ( epc_signature->get_input_types()        );
    descent_p.set_LOWER_BOUND       ( epc_signature->get_lb()                 );
    descent_p.set_UPPER_BOUND       ( epc_signature->get_ub()                 );
    descent_p.set_FIXED_VARIABLE    ( epc_signature->get_fixed_variables()    );
    descent_p.set_PERIODIC_VARIABLE ( epc_signature->get_periodic_variables() );
    descent_p.set_VARIABLE_GROUP    ( epc_signature->get_var_groups()         );
    descent_p.set_BB_OUTPUT_TYPE    ( _p.get_bb_output_type() );
    descent_p.set_DIRECTION_TYPE ( _p.get_direction_types() );
    descent_p.set_SEC_POLL_DIR_TYPE ( _p.get_sec_poll_dir_types() );
    
    
    {
        const std::list<std::string> & bb_exe = _p.get_bb_exe();
        descent_p.set_BB_EXE ( bb_exe );
        std::list<std::string>::const_iterator it , end = bb_exe.end();
        for ( it = bb_exe.begin() ; it != end ; ++it )
            descent_p.set_SGTE_EXE ( *it , _p.get_sgte_exe ( *it ) );
    }
    
    descent_p.set_PROBLEM_DIR    ( _p.get_problem_dir()    );
    descent_p.set_TMP_DIR        ( _p.get_tmp_dir()        );
    descent_p.set_SGTE_COST      ( _p.get_sgte_cost()      );
    descent_p.set_SGTE_EVAL_SORT ( _p.get_sgte_eval_sort() );
    descent_p.set_X0             ( *y                      );
    
    bool has_sgte = _p.has_sgte();
    if ( has_sgte )
    {
        descent_p.reset_f_target();
        descent_p.set_HAS_SGTE         ( true            );
        descent_p.set_OPT_ONLY_SGTE    ( true            );
        descent_p.set_STOP_IF_FEASIBLE ( false           );
        descent_p.set_MODEL_SEARCH     ( false );
        descent_p.set_MODEL_EVAL_SORT  ( NOMAD::NO_MODEL );
    }
    else
    {
        descent_p.set_F_TARGET         ( _p.get_f_target()         );
        descent_p.set_STOP_IF_FEASIBLE ( _p.get_stop_if_feasible() );
        descent_p.set_MODEL_EVAL_SORT ( _p.get_model_eval_sort());
        descent_p.set_MODEL_SEARCH (_p.has_model_search());
        
        
    }
    
    descent_p.set_LH_SEARCH   ( 0 , 0 );
    
    int bbe       = stats.get_bb_eval();
    int sgte_eval = stats.get_sgte_eval();
    
    {
        int p_max_bbe = _p.get_max_bb_eval();
        if ( p_max_bbe > 0 )
            descent_p.set_MAX_BB_EVAL ( p_max_bbe - bbe );
        
        int p_max_sgte_eval = _p.get_max_sgte_eval();
        if ( p_max_sgte_eval > 0 )
            descent_p.set_MAX_SGTE_EVAL ( p_max_sgte_eval - sgte_eval );
        
        int p_max_eval = _p.get_max_eval();
        if ( p_max_eval > 0 )
            descent_p.set_MAX_EVAL ( p_max_eval - stats.get_eval() );
        
        int p_max_sbe = _p.get_max_sim_bb_eval();
        if ( p_max_sbe > 0 )
            descent_p.set_MAX_SIM_BB_EVAL ( p_max_sbe - stats.get_sim_bb_eval() );
        
        int p_max_time = _p.get_max_time();
        if ( p_max_time > 0 )
            descent_p.set_MAX_TIME ( p_max_time - stats.get_real_time() );
        
        NOMAD::Double p_sst = _p.get_stat_sum_target();
        if ( p_sst.is_defined() )
            descent_p.set_STAT_SUM_TARGET ( p_sst - stats.get_stat_sum() );
    }
    
    descent_p.set_OPPORTUNISTIC_EVAL    ( _p.get_opportunistic_eval()    );
    descent_p.set_BB_INPUT_INCLUDE_SEED ( _p.get_bb_input_include_seed() );
    descent_p.set_BB_INPUT_INCLUDE_TAG  ( _p.get_bb_input_include_tag()  );
    descent_p.set_BB_REDIRECTION        ( _p.get_bb_redirection()        );
    
    descent_p.set_EXTENDED_POLL_ENABLED ( false );
    descent_p.set_USER_CALLS_ENABLED    ( false );
    
    descent_p.set_H_MAX_0                      ( _p.get_h_max_0()                      );
    descent_p.set_H_MIN                        ( _p.get_h_min()                        );
    descent_p.set_H_NORM                       ( _p.get_h_norm()                       );
    descent_p.set_RHO                          ( _p.get_rho()                          );
    descent_p.set_SNAP_TO_BOUNDS               ( _p.get_snap_to_bounds()               );
    descent_p.set_MAX_CACHE_MEMORY             ( _p.get_max_cache_memory()             );
    descent_p.set_SPECULATIVE_SEARCH           ( _p.get_speculative_search()           );
    descent_p.set_OPPORTUNISTIC_LUCKY_EVAL     ( _p.get_opportunistic_lucky_eval()     );
    descent_p.set_OPPORTUNISTIC_MIN_EVAL       ( _p.get_opportunistic_min_eval()       );
    descent_p.set_OPPORTUNISTIC_MIN_F_IMPRVMT  ( _p.get_opportunistic_min_f_imprvmt()  );
    descent_p.set_OPPORTUNISTIC_MIN_NB_SUCCESS ( _p.get_opportunistic_min_nb_success() );
    
    if (_p.eval_points_as_block())
        descent_p.set_BB_MAX_BLOCK_SIZE(	_p.get_bb_max_block_size()		);
    
    
    descent_p.set_CACHE_FILE        ( _p.get_cache_file()        );
    descent_p.set_SGTE_CACHE_FILE   ( _p.get_sgte_cache_file()   );
    descent_p.set_CACHE_SAVE_PERIOD ( _p.get_cache_save_period() );
    
    descent_p.set_ADD_SEED_TO_FILE_NAMES ( _p.get_add_seed_to_file_names() );
    
    descent_p.set_DISPLAY_ALL_EVAL(_p.get_display_all_eval());
    if ( _p.out().get_poll_dd() == NOMAD::FULL_DISPLAY )
        descent_p.set_DISPLAY_DEGREE ( NOMAD::NORMAL_DISPLAY );
    else if (_p.out().get_poll_dd() == NOMAD::NORMAL_DISPLAY )
        descent_p.set_DISPLAY_DEGREE ( NOMAD::MINIMAL_DISPLAY );
    else
        descent_p.set_DISPLAY_DEGREE ( _p.out().get_poll_dd());
    
    // Stats style modified
    if ( has_sgte )
        descent_p.set_DISPLAY_STATS ( NOMAD::itos(sgte_eval) + "+SGTE OBJ (ExtendedPoll---surrogate)" );
    else
    {
        std::list<std::string> ds = _p.get_display_stats();
        std::list<std::string>::iterator       it    = ds.begin();
        std::list<std::string>::const_iterator end   = ds.end();
        std::string  s_bbe = NOMAD::itos(bbe) + "+";
        while ( it != end )
        {
            if (*it == "BBE")
                ds.insert ( it , s_bbe );
            ++it;
        }
        ds.push_back ( " (ExtendedPoll)" );
        descent_p.set_DISPLAY_STATS ( ds );
    }
    
    // STATS_FILE:
    if ( has_sgte )
        descent_p.set_STATS_FILE ( _p.get_stats_file_name() , NOMAD::itos(sgte_eval) + "+SGTE OBJ (ExtendedPoll---surrogate)" );
    else
    {
        std::list<std::string>                 sf    = _p.get_stats_file();
        std::list<std::string>::iterator       it    = sf.begin();
        std::list<std::string>::const_iterator end   = sf.end();
        std::string                            s_bbe = NOMAD::itos(bbe) + "+";
        while ( it != end )
        {
            if ( *it == "BBE" )
                sf.insert ( it , s_bbe );
            ++it;
        }
        sf.push_back ( " (ExtendedPoll)" );
        descent_p.set_STATS_FILE ( _p.get_stats_file_name() , sf );
    }
    
    // Mesh:
    {
        const OrthogonalMesh * mesh = epc_signature->get_mesh();
        descent_p.set_MIN_MESH_SIZE ( mesh->get_min_mesh_size() );
        descent_p.set_MIN_POLL_SIZE ( mesh->get_min_poll_size() );
        descent_p.set_INITIAL_POLL_SIZE ( mesh->get_initial_poll_size() , false );
        
    }
    
    // check the parameters:
    try
    {
        descent_p.check ( false ,    // remove_history_file  = false
                         false ,    // remove_solution_file = false
                         false   ); // remove_stats_file    = false
    }
    catch ( NOMAD::Exception & e )
    {
        std::ostringstream err;
        err << "-- " << e.what();
        throw NOMAD::Exception ( "Extended_Poll.cpp" , __LINE__ , err.str() );
    }
    
}

/*----------------------------------------------------------------*/
/*            evaluation of an extended poll point (private)      */
/*----------------------------------------------------------------*/
const NOMAD::Eval_Point * NOMAD::Extended_Poll::eval_epp
( NOMAD::Eval_Point        * y              ,
 Mads                     & mads           ,
 bool                     & stop           ,
 NOMAD::stop_type         & stop_reason    ,
 NOMAD::success_type      & success        ,
 const NOMAD::Eval_Point *& new_feas_inc   ,
 const NOMAD::Eval_Point *& new_infeas_inc ) const
{
    NOMAD::Evaluator_Control & ev_control     = mads.get_evaluator_control();
    const NOMAD::Display     & out            = _p.out();
    NOMAD::dd_type             display_degree = out.get_poll_dd();
    
    // initial display:
    if ( display_degree == NOMAD::FULL_DISPLAY )
    {
        out << std::endl << NOMAD::open_block ( "extended poll point eval" ) << std::endl
        << "extended poll point = ( ";
        y->Point::display ( out , " " , 5 , _p.get_point_display_limit() );
        out << " )" << std::endl;
    }
    
    // add the eval point to the evaluator control:
    ev_control.add_eval_point ( y                       ,
                               display_degree          ,
                               _p.get_snap_to_bounds() ,
                               NOMAD::Double()         ,
                               NOMAD::Double()         ,
                               NOMAD::Double()         ,
                               NOMAD::Double()           );
    
    // get the stats:
    NOMAD::Stats & stats   = mads.get_stats();
    int            old_bbe = stats.get_bb_eval();
    
    // eval list of points:
    new_feas_inc = new_infeas_inc = NULL;
    std::list<const NOMAD::Eval_Point *> evaluated_pts;
    
    ev_control.eval_list_of_points ( NOMAD::EXTENDED_POLL    ,
                                    mads.get_true_barrier() ,
                                    mads.get_sgte_barrier() ,
                                    mads.get_pareto_front() ,
                                    stop                    ,
                                    stop_reason             ,
                                    new_feas_inc            ,
                                    new_infeas_inc          ,
                                    success                 ,
                                    &evaluated_pts            );
    
    // update the number of extended poll blackbox evaluations:
    stats.add_ext_poll_bb_eval ( stats.get_bb_eval() - old_bbe );
    
    // final display:
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out << std::endl << NOMAD::close_block() << std::endl;
    
    // return the evaluated point:
    return ( evaluated_pts.size() != 1 ) ? NULL : *evaluated_pts.begin();
}

/*-----------------------------------------*/
/*         run the extended poll           */
/*-----------------------------------------*/
void NOMAD::Extended_Poll::run ( Mads                     & mads            ,
                                int                      & nb_ext_poll_pts ,
                                bool                     & stop            ,
                                NOMAD::stop_type         & stop_reason     ,
                                NOMAD::success_type      & success         ,
                                const NOMAD::Eval_Point *& new_feas_inc    ,
                                const NOMAD::Eval_Point *& new_infeas_inc   )
{
    nb_ext_poll_pts = 0;
    success         = NOMAD::UNSUCCESSFUL;
    new_feas_inc    = new_infeas_inc = NULL;
    
    if ( stop || _extended_points.empty() )
        return;
    
    const NOMAD::Display & out            = _p.out();
    NOMAD::dd_type         display_degree = out.get_poll_dd();
    NOMAD::Eval_Point    * cur;
    
    // phase 1: evaluate the extended poll points in order to sort them
    // -------- (based on the surrogates or on the true function):
    if ( _extended_points.size() > 1 )
    {
        
        bool has_sgte           = _p.has_sgte();
        bool old_sgte_eval_sort = _p.get_sgte_eval_sort();
        
        // phase 1 initial display:
        if ( display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss;
            oss << "extended poll pts sorting";
            if ( has_sgte )
                oss << " (on surrogates)";
            out << std::endl << NOMAD::open_block ( oss.str() ) << std::endl;
        }
        
        if ( has_sgte ) {
            _p.set_SGTE_EVAL_SORT ( false ); // this ensures that only surrogate
            _p.force_check_flag();           // evaluations are performed
        }
        
        NOMAD::Evaluator_Control & ev_control = mads.get_evaluator_control();
        
        // loop on the extended poll points:
        std::list<NOMAD::Eval_Point *>::const_iterator end = _extended_points.end();
        for ( std::list<NOMAD::Eval_Point *>::iterator it = _extended_points.begin() ;
             it != end ; ++it )
        {
            cur = *it;
            
            if ( has_sgte )
                cur->set_eval_type ( NOMAD::SGTE );
            
            ev_control.add_eval_point ( cur                     ,
                                       display_degree          ,
                                       _p.get_snap_to_bounds() ,
                                       NOMAD::Double()         ,
                                       NOMAD::Double()         ,
                                       NOMAD::Double()         ,
                                       NOMAD::Double()           );
        }
        
        _extended_points.clear();
        
        // get the stats:
        NOMAD::Stats & stats   = mads.get_stats();
        int            old_bbe = stats.get_bb_eval();
        
        // number of eval points:
        nb_ext_poll_pts = ev_control.get_nb_eval_points();
        
        // eval list of points:
        new_feas_inc = new_infeas_inc = NULL;
        std::list<const NOMAD::Eval_Point *> evaluated_pts;
        
        ev_control.eval_list_of_points ( NOMAD::EXTENDED_POLL    ,
                                        mads.get_true_barrier() ,
                                        mads.get_sgte_barrier() ,
                                        mads.get_pareto_front() ,
                                        stop                    ,
                                        stop_reason             ,
                                        new_feas_inc            ,
                                        new_infeas_inc          ,
                                        success                 ,
                                        &evaluated_pts            );
        if ( has_sgte )
        {
            if ( !_p.get_opt_only_sgte() ) {
                success      = NOMAD::UNSUCCESSFUL;
                new_feas_inc = new_infeas_inc = NULL;
            }
            _p.set_SGTE_EVAL_SORT ( old_sgte_eval_sort );
            _p.force_check_flag();
        }
        
        // update the number of extended poll blackbox evaluations:
        stats.add_ext_poll_bb_eval ( stats.get_bb_eval() - old_bbe );
        
        // sort the evaluated extended poll points:
        if ( success != NOMAD::FULL_SUCCESS )
            sort_epp ( evaluated_pts );
        
        // the extended poll is terminated in case of success:
        if ( stop || success == NOMAD::FULL_SUCCESS || new_feas_inc || new_infeas_inc )
        {
            if ( display_degree == NOMAD::FULL_DISPLAY )
                out << std::endl
                << NOMAD::close_block ( "end of ext poll pts sorting (success)" )
                << std::endl;
            return;
        }
        
        // phase 1 final display:
        if ( display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss;
            oss << "end of ext poll pts sorting";
            if ( has_sgte )
                oss << " (on surrogates)";
            out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
        }
        
    } // end of phase 1
    // --------------
    
    /*--------------------------------------------------------------*/
    
    // phase 2: execute the extended poll for each extended point:
    // --------
    {
        const NOMAD::Eval_Point * old_bf = mads.get_best_feasible();
        const NOMAD::Eval_Point * old_bi = mads.get_best_infeasible();
        const NOMAD::Eval_Point * y;
        
        while ( !_extended_points.empty() )
        {
            
            cur = *_extended_points.begin();
            
            // the point has already been evaluated during
            // the extended poll points sorting:
            if ( cur->is_in_cache() && cur->get_eval_type() == NOMAD::TRUTH )
                y = cur;
            
            // the point has to be evaluated:
            else
            {
                
                y = eval_epp ( cur            ,
                              mads           ,
                              stop           ,
                              stop_reason    ,
                              success        ,
                              new_feas_inc   ,
                              new_infeas_inc   );
                
                ++nb_ext_poll_pts;
                
                // the extended poll is terminated in case of success:
                if ( !y                             ||
                    stop                           ||
                    success == NOMAD::FULL_SUCCESS ||
                    new_feas_inc                   ||
                    new_infeas_inc                    )
                    break;
            }
            
            _extended_points.pop_front();
            
            // perform the extended poll descent ?
            if ( display_degree == NOMAD::FULL_DISPLAY )
            {
                out << std::endl
                << "extended poll center: ( ";
                y->Point::display ( out , " " , 5 , _p.get_point_display_limit() );
                out << " )" << std::endl << std::endl
                << "perform extended poll descent ...";
            }
            if ( check_trigger ( old_bf , old_bi , y ) )
            {
                if ( display_degree == NOMAD::FULL_DISPLAY )
                    out << "... yes" << std::endl;
                descent ( y               ,
                         mads            ,
                         nb_ext_poll_pts ,
                         stop            ,
                         stop_reason     ,
                         success         ,
                         new_feas_inc    ,
                         new_infeas_inc    );
                
                // the extended poll is terminated in case of success:
                if ( stop || success == NOMAD::FULL_SUCCESS || new_feas_inc || new_infeas_inc )
                    break;
            }
            else if ( display_degree == NOMAD::FULL_DISPLAY )
                out << "... no" << std::endl;
        }
    }  // end of phase 2
    // --------------
    
    // clean the extended points that have not been considered:
    std::list<NOMAD::Eval_Point *>::const_iterator end = _extended_points.end();
    for ( std::list<NOMAD::Eval_Point *>::iterator it = _extended_points.begin() ;
         it != end ; ++it )
        if ( !(*it)->is_in_cache() )
            delete *it;
    _extended_points.clear();
}

/*----------------------------------------------------------------*/
/*  sort the extended poll points after they have been evaluated  */
/*  (private)                                                     */
/*----------------------------------------------------------------*/
void NOMAD::Extended_Poll::sort_epp
( const std::list<const NOMAD::Eval_Point *> & evaluated_pts )
{
    const NOMAD::Display                 & out            = _p.out();
    NOMAD::dd_type                         display_degree = out.get_poll_dd();
    const NOMAD::Double                  & h_min          = _p.get_h_min();
    std::set<NOMAD::Priority_Eval_Point>   sorted_pts;
    
    // 1. loop on the evaluated points:
    std::list<const NOMAD::Eval_Point *>::const_iterator it1 , end1 = evaluated_pts.end();
    for ( it1 = evaluated_pts.begin() ; it1 != end1 ; ++it1 )
    {
        
        // creation of a Priority_Eval_Point:
        NOMAD::Priority_Eval_Point pep ( *it1 , h_min );
        
        // surrogate values for f and h:
        if ( (*it1)->get_eval_type() == NOMAD::SGTE )
        {
            pep.set_f_sgte ( (*it1)->get_f() );
            pep.set_h_sgte ( (*it1)->get_h() );
        }
        
        // insertion in the sorted list of points:
        sorted_pts.insert ( pep );
    }
    
    // 2. loop on the sorted points:
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out << std::endl << NOMAD::open_block ( "sorted ext poll pts" ) << std::endl;
    
    const NOMAD::Eval_Point * cur;
    NOMAD::Eval_Point       * y;
    int i = 0 , nb_pts = static_cast<int> ( sorted_pts.size() );
    std::set<NOMAD::Priority_Eval_Point>::const_iterator
    it2 , end2 = sorted_pts.end();
    
    for ( it2 = sorted_pts.begin() ; it2 != end2 ; ++it2 )
    {
        
        // we copy y=cur an create a new Eval_Point because cur can be a surrogate
        // point already in the surrogate cache
        
        cur = (*it2).get_point();
        
        y = new NOMAD::Eval_Point;
        y->set ( cur->size() , _p.get_bb_nb_outputs() );
        y->set_signature  ( cur->get_signature () );
        y->set_direction  ( cur->get_direction () );
        y->Point::operator = ( *cur );
        
        // display:
        if ( display_degree == NOMAD::FULL_DISPLAY ) 
        {
            out << "point #";
            out.display_int_w ( ++i , nb_pts );
            out << "/" << nb_pts << " : ( ";
            y->Point::display ( out , " " , 5 , _p.get_point_display_limit() );
            out << " )" << std::endl;
        }
        
        // insertion in _extended_points:
        _extended_points.push_back ( y );
    }
    
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out << std::endl << NOMAD::close_block() << std::endl;
}

/*--------------------------------------------------------------------*/
/*  set the neighbors executable name for the default implementation  */
/*--------------------------------------------------------------------*/
bool NOMAD::Extended_Poll::set_neighbors_exe ( std::string & error_str )
{
    error_str.clear();    
    
    _neighbors_exe = _p.get_neighbors_exe();
    
    if ( _neighbors_exe.empty() ) {
        error_str = "categorical variables: parameter NEIGHBORS_EXE is undefined";
        return false;
    }
    
    const std::string problem_dir = _p.get_problem_dir();
    
    std::list<std::string> neighbors_exe_words;
    NOMAD::get_words ( _neighbors_exe , neighbors_exe_words );
    
    // _neighbors_exe is composed of several words (it is a command):
    if ( neighbors_exe_words.size() > 1 ) 
    {
        
        _neighbors_exe.clear();
        
        std::list<std::string>::const_iterator it  = neighbors_exe_words.begin() ,
        end = neighbors_exe_words.end();
        while (true) {
            
            if ( (*it)[0] != '$' ) {
                _neighbors_exe += "\"" + problem_dir;
                _neighbors_exe += *it + "\"";
            }
            else
                _neighbors_exe += it->substr ( 1 , it->size()-1 );
            
            ++it;
            
            if ( it == end )
                break;
            
            _neighbors_exe += " ";
        }
    }
    
    // _neighbors_exe is just composed of one name (it is an executable):
    else
    {
        
        if ( _neighbors_exe[0] != '$' )
            _neighbors_exe = problem_dir + _neighbors_exe;
        else
            _neighbors_exe = _neighbors_exe.substr ( 1 , _neighbors_exe.size()-1 );
        
        if ( !NOMAD::check_exe_file ( _neighbors_exe ) )
        {
            error_str =   "categorical variables: \'" + _neighbors_exe
            + "\' is not a valid executable file";
            return false;
        }
        
        if ( _neighbors_exe[0] != '$' )
            _neighbors_exe = "\"" + _neighbors_exe + "\"";
    }
    
    return true;
}

/*----------------------------------------------------------------------*/
/*  construct the extended poll points: this is the default version of  */
/*  this virtual function: it calls the executable defined by the       */
/*  NEIGHBORS_EXE parameter                                             */
/*----------------------------------------------------------------------*/
void NOMAD::Extended_Poll::construct_extended_points ( const NOMAD::Eval_Point & xk )
{
    if ( _neighbors_exe.empty() )
        throw NOMAD::Exception ( "Extended_Poll.cpp" , __LINE__ ,
                                "NOMAD::Extended_Poll::construct_extended_points(): no NEIGHBORS_EXE executable (batch mode) or no subclass implementation of the method (library mode)" );
    
    if ( !xk.is_complete() )
        throw NOMAD::Exception ( "Extended_Poll.cpp" , __LINE__ ,
                                "NOMAD::Extended_Poll::construct_extended_points(): bad extended poll center");
    
    NOMAD::Signature * signature = _p.get_signature();
    if ( !signature )
        throw NOMAD::Exception ( "Extended_Poll.cpp" , __LINE__ ,
                                "NOMAD::Extended_Poll::construct_extended_points(): no signature" );
    
    std::string tmp_dir = _p.get_tmp_dir();
    
    std::ostringstream oss;
    oss << "." << _p.get_seed() << "." << xk.get_tag() << ".neighbors.";
    const std::string & sint = oss.str();
    
    // input file writing:
    // -------------------
    std::string input_file_name =
    tmp_dir + NOMAD::BLACKBOX_INPUT_FILE_PREFIX 
    + sint + NOMAD::BLACKBOX_INPUT_FILE_EXT;
    
    std::string output_file_name =
    tmp_dir + NOMAD::BLACKBOX_OUTPUT_FILE_PREFIX
    + sint + NOMAD::BLACKBOX_OUTPUT_FILE_EXT;
    
    std::ofstream fout ( input_file_name.c_str() );
    if ( fout.fail() )
    {
        remove ( input_file_name.c_str () );
        std::string err = "could not open file neighbors input file " + input_file_name;
        throw NOMAD::Exception ( "Extended_Poll.cpp" , __LINE__ , err );
    }
    
    fout.setf ( std::ios::fixed );
    fout.precision ( NOMAD::DISPLAY_PRECISION_BB );
    xk.Point::display ( fout , " " , -1 , -1 );
    fout << std::endl;
    
    fout.close();
    
    if ( fout.fail() ) 
    {
        remove ( input_file_name.c_str () );
        std::string err = "could not write file neighbors input file " + input_file_name;
        throw NOMAD::Exception ( "Extended_Poll.cpp" , __LINE__ , err );
    }
    
    // system call to get the neighbors:
    // ---------------------------------
    
    std::string cmd = _neighbors_exe + " " + input_file_name + " > " + output_file_name;
    
#ifdef DEBUG
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank);
    _p.out() << "command(rank=" << rank
    << ") = \'" << cmd << "\'" << std::endl;
#else
    _p.out() << "command=\'" << cmd << "\'" << std::endl;
#endif
#endif
    
    // the call:
    if ( ( system ( cmd.c_str() ) ) != 0 )
    {
        remove ( input_file_name.c_str () );
        remove ( output_file_name.c_str() );
        std::string err = "error with command " + cmd;
        throw NOMAD::Exception ( "Extended_Poll.cpp" , __LINE__ , err );
    }
    
    // reading of the output file:
    // ---------------------------
    
    std::ifstream fin ( output_file_name.c_str() );
    
    if ( fin.fail() ) 
    {
        remove ( input_file_name.c_str () );
        remove ( output_file_name.c_str() );
        std::string err = "could not open neighbors output file " + output_file_name;
        throw NOMAD::Exception ( "Extended_Poll.cpp" , __LINE__ , err );
    }
    
    int n = xk.size();
    
    while ( true )
    {
        NOMAD::Point y(n);
        try 
        {      
            fin >> y;
        }
        catch ( NOMAD::Point::Bad_Input & )
        {
            if ( y.is_defined() ) {
                remove ( input_file_name.c_str () );
                remove ( output_file_name.c_str() );
                std::string err = "error with neighbor in file " + output_file_name;
                throw NOMAD::Exception ( "Extended_Poll.cpp" , __LINE__ , err );
            }
            break;
        }
        
        if ( !y.is_complete() )
        {
            remove ( input_file_name.c_str () );
            remove ( output_file_name.c_str() );
            std::string err = "error with neighbor in file " + output_file_name;
            throw NOMAD::Exception ( "Extended_Poll.cpp" , __LINE__ , err );
        }
        
        add_extended_poll_point ( y , *signature );
    }
    
    fin.close();
    
    // delete the input and output files:
    // ----------------------------------
    remove ( input_file_name.c_str () );
    remove ( output_file_name.c_str() );
}
