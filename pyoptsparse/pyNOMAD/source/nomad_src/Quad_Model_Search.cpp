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
 \file   Quad_Model_Search.cpp
 \brief  Quadratic Model search (implementation)
 \author Sebastien Le Digabel
 \date   2010-08-30
 \see    Quad_Model_Search.hpp
 */
#include "Quad_Model_Search.hpp"

/*----------------------------------------------------------------*/
/*                            the search                          */
/*----------------------------------------------------------------*/
/* Search parameters:                                             */
/* ------------------                                             */
/*                                                                */
/*  . MODEL_SEARCH: flag to activate the model search (MS)        */
/*                  (here its value is NOMAD::QUADRATIC_MODEL)    */
/*                                                                */
/*  . MODEL_SEARCH_OPTIMISTIC: if true, the direction from the    */
/*                             model center to the trial point    */
/*                             is computed and prossibly used     */
/*                             in the speculative search          */
/*                             default=yes                        */
/*                                                                */
/*  . MODEL_SEARCH_PROJ_TO_MESH: project or not to mesh           */
/*                               default=yes                      */
/*                                                                */
/*  . MODEL_SEARCH_MAX_TRIAL_PTS: limit on the number of trial    */
/*                                points for one search (in       */
/*                                {1,2,3,4} and with default=4    */
/*                                for quadratic models)           */
/*                                                                */
/*  . MODEL_RADIUS_FACTOR (r): Y points are in B(xk,r.Delta^p_k)  */
/*                             default=2.0                        */
/*                                                                */
/*  . MODEL_MAX_Y_SIZE: limit on the size of Y; if equal to       */
/*                      (n+1)(n+2)/2, regression is never used    */
/*                      default=500                               */
/*                                                                */
/*  . MODEL_EVAL_SORT: if true, all evaluation points are sorted  */
/*                     according to model values; for the model   */
/*                     search, this is used to order the up to 4  */
/*                     trial points.                              */
/*                     default=yes                                */
/*                                                                */
/*  . MODEL_USE_WP: if true, well-poisedness strategy is applied  */
/*                  for regression models                         */
/*                                                                */
/*  . construct model centered around best_feas and best_infeas   */
/*    or just around best_feas;                                   */
/*    default=around them both, not modifiable                    */
/*                                                                */
/* SVD decomposition parameters:                                  */
/* -----------------------------                                  */
/*  . SVD_EPS    : epsilon; default=1e-13                         */
/*  . SVD_MAX_MPN: max matrix size: default: m+n <= 1500          */
/*                                                                */
/*----------------------------------------------------------------*/
void NOMAD::Quad_Model_Search::search ( NOMAD::Mads              & mads           ,
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
    
    if ( stop )
    {
        if ( display_degree == NOMAD::FULL_DISPLAY )
            out << "Quad_Model_Search::search(): not performed (stop flag is active)"
            << std::endl;
        return;
    }
    
    // black-box output types:
    const std::vector<NOMAD::bb_output_type> & bbot           = _p.get_bb_output_type();
    const std::list<int>                     & index_obj_list = _p.get_index_obj();
    
    if ( index_obj_list.empty() )
    {
        if ( display_degree == NOMAD::FULL_DISPLAY )
            out << "Quad_Model_Search::search(): not performed with no objective function"
            << std::endl;
        return;
    }
    
    
    // initial displays:
    if ( display_degree == NOMAD::FULL_DISPLAY )
    {
        std::ostringstream oss;
        oss << NOMAD::MODEL_SEARCH << " #"
        << _all_searches_stats.get_MS_nb_searches();
        out << std::endl << NOMAD::open_block ( oss.str() ) << std::endl;
    }
    
    // surrogate or truth model evaluations:
    NOMAD::eval_type ev_type =
    ( _p.get_opt_only_sgte() ) ? NOMAD::SGTE : NOMAD::TRUTH;
    
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out << "model construction for " << ev_type << std::endl;
    
    // active cache:
    const NOMAD::Cache & cache = mads.get_cache();
    
    // active barrier:
    const NOMAD::Barrier & barrier =
    ( ev_type == NOMAD::SGTE ) ? mads.get_sgte_barrier() : mads.get_true_barrier();
    
    // current incumbents: xk[0]=best_feas and xk[1]=best_infeas:
    const NOMAD::Eval_Point * xk[2];
    xk[0] = barrier.get_best_feasible  ();
    xk[1] = barrier.get_best_infeasible();
    
    if ( !xk[0] && !xk[1] )
    {
        if ( display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss;
            oss << "end of " << NOMAD::MODEL_SEARCH << " (no incumbent)";
            out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
        }
        return;
    }
    
    // from this point the search is counted:
    count_search = true;
    _one_search_stats.add_MS_nb_searches();
    
    // display the number of cache points:
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out << std::endl << "number of points in cache: "
        << cache.size() << std::endl;
    
    // stats:
    NOMAD::Stats & stats = mads.get_stats();
    
    // number of interpolation points:
    int nY[2];
    nY[0] = nY[1] = -1;
    
    int min_Y_size = _p.get_model_quad_min_Y_size();
    int max_Y_size = _p.get_model_quad_max_Y_size();
    
    // use or not well-poisedness:
    bool use_WP = _p.get_model_quad_use_WP();
    
    // flag to detect model errors:
    bool                       model_ok   = false;
    NOMAD::Evaluator_Control & ev_control = mads.get_evaluator_control();
    
    // main loop on the two incumbents (feasible and infeasible):
    // ---------
    for ( int i_inc = 0 ; i_inc < 2 ; ++i_inc )
    {
        
        if ( xk[i_inc] )
        {
            
            // display the model center:
            if ( display_degree == NOMAD::FULL_DISPLAY )
            {
                out << std::endl << "model center";
                if ( xk[0] && xk[1] )
                    out << " (" << i_inc+1 << "/2)";
                out << ": " << *xk[i_inc] << std::endl;
            }
            
            // get and check the signature:
            NOMAD::Signature * signature = xk[i_inc]->get_signature();
            if ( !signature )
            {
                if ( display_degree == NOMAD::FULL_DISPLAY ) {
                    std::ostringstream oss;
                    oss << "end of " << NOMAD::MODEL_SEARCH << " (no signature)";
                    out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
                }
                stats.update_model_stats   ( _one_search_stats );
                _all_searches_stats.update ( _one_search_stats );
                return;
            }
            
            // current mesh index:
            NOMAD::Point mesh_indices = signature->get_mesh()->get_mesh_indices();
            
            
            int n = signature->get_n();
            if ( n != xk[i_inc]->size() )
            {
                if ( display_degree == NOMAD::FULL_DISPLAY )
                {
                    std::ostringstream oss;
                    oss << "end of " << NOMAD::MODEL_SEARCH << " (incompatible signature)";
                    out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
                }
                stats.update_model_stats   ( _one_search_stats );
                _all_searches_stats.update ( _one_search_stats );
                return;
            }
            
            // compute the interpolation radius: points in Y must be at
            // a max distance of ms_radius_factor times Delta^k:
            NOMAD::Point Delta , delta;
            signature->get_mesh()->get_Delta ( Delta );
            signature->get_mesh()->get_delta ( delta );
            
            
            NOMAD::Point interpolation_radius = Delta;
            interpolation_radius *= _p.get_model_quad_radius_factor();
            
            if ( display_degree == NOMAD::FULL_DISPLAY )
                out << "mesh indices         : ("   << mesh_indices  << ")"    << std::endl
                << "mesh size parameter  : ( " << delta << " )" << std::endl
                << "poll size parameter  : ( " << Delta << " )" << std::endl
                << "interpolation radius : ( " << interpolation_radius
                << " )" << std::endl;
            
            // creation of the model:
            NOMAD::Quad_Model model ( out , bbot , cache , *signature );
            
            
            
            NOMAD::Clock      clock;
            
            // construct interpolation set Y:
            // ------------------------------
            model.construct_Y ( *xk[i_inc]           ,
                               interpolation_radius ,
                               max_Y_size             );
            
            nY[i_inc] = model.get_nY();
            
            if ( display_degree == NOMAD::FULL_DISPLAY )
            {
                out << "number of points in Y: " << nY[i_inc]
                << " (p=" << nY[i_inc]-1;
                if ( nY[i_inc] < 2 )
                    out << ", not enough";
                out << ")" << std::endl;
            }
            
            if ( nY[i_inc] < 2 )
                _one_search_stats.add_not_enough_pts();
            else
            {
                
#ifdef DEBUG
                out << std::endl;
                model.display_Y ( out , "unscaled interpolation set Y" );
#endif
                
                // define scaling:
                // ---------------
                // The min box around the interpolation set Y
                //   is scaled to [-r;r] with r=MODEL_RADIUS_FACTOR.
                model.define_scaling ( _p.get_model_quad_radius_factor() );
                
#ifdef DEBUG
                out << std::endl;
                model.display_Y ( out , "scaled interpolation set Ys" );
#endif
                
                // error check:
                if ( model.get_error_flag() )
                    _one_search_stats.add_construction_error();
                
                // not enough points (this is not counted as an error):
                else if ( nY[i_inc] < 2                                        ||
                         ( min_Y_size < 0 && nY[i_inc] <= model.get_nfree() )    )
                    _one_search_stats.add_not_enough_pts();
                
                // no error and enough points in Y:
                else
                {
                    
                    // construct model:
                    // ----------------
                    model.construct ( use_WP , NOMAD::SVD_EPS , NOMAD::SVD_MAX_MPN , max_Y_size );
                    _one_search_stats.add_construction_time ( clock.get_CPU_time() );
                    _one_search_stats.update_nY ( model.get_nY() );
                    
                    // display model characteristics:
#ifdef DEBUG
                    out << std::endl;
                    model.display_model_coeffs ( out );
                    out << std::endl;
                    model.display_Y_error ( out );
#endif
                    // count model:
                    if ( ev_type == NOMAD::TRUTH )
                        _one_search_stats.add_nb_truth();
                    else
                        _one_search_stats.add_nb_sgte();
                    
                    switch ( model.get_interpolation_type() )
                    {
                        case NOMAD::MFN:
                            _one_search_stats.add_nb_MFN();
                            break;
                        case NOMAD::WP_REGRESSION:
                            _one_search_stats.add_nb_WP_regression();
                            break;
                        case NOMAD::REGRESSION:
                            _one_search_stats.add_nb_regression();
                            break;
                        default:
                            break;
                    }
                    
                    // check model error flag:
                    const NOMAD::Double & cond = model.get_cond();
                    if ( model.get_error_flag()     ||
                        !cond.is_defined()         ||
                        cond > NOMAD::SVD_MAX_COND    )
                    {
                        if ( model.get_error_flag() )
                            _one_search_stats.add_construction_error();
                        else
                            _one_search_stats.add_bad_cond();
                    }
                    else
                    {
                        
                        model_ok = true;
                        
                        // optimize model:
                        // ---------------
                        NOMAD::Point xf , xi;
                        
                        clock.reset();
                        
                        bool optimization_ok = optimize_model ( model          ,
                                                               xk             ,
                                                               i_inc          ,
                                                               display_degree ,
                                                               out            ,
                                                               xf             ,
                                                               xi             ,
                                                               stop           ,
                                                               stop_reason      );
                        
                        _one_search_stats.add_optimization_time ( clock.get_CPU_time() );
                        
                        if ( optimization_ok )
                        {
                            
                            // get solution(s), project to mesh (+round for integers), and create trial points:
                            // ----------------------------------------------------------
                            if ( xf.is_defined() )
                                create_trial_point ( ev_control     ,
                                                    xf             ,
                                                    model          ,
                                                    *signature     ,
                                                    mesh_indices   ,
                                                    delta          ,
                                                    display_degree ,
                                                    out              );
                            
                            if ( xi.is_defined() )
                                create_trial_point ( ev_control     ,
                                                    xi             ,
                                                    model          ,
                                                    *signature     ,
                                                    mesh_indices   ,
                                                    delta          ,
                                                    display_degree ,
                                                    out              );
                        }
                        else
                            _one_search_stats.add_MS_opt_error();
                    }
                }
            }
        }
    } // end of main loop
    
    // check the number of times that not enough points could be considered:
    if ( nY[0] <= 1 && nY[1] <= 1 ) {
        if ( display_degree == NOMAD::FULL_DISPLAY ) {
            std::ostringstream oss;
            oss << "end of " << NOMAD::MODEL_SEARCH
            << " (not enough points)";
            out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
        }
        stats.update_model_stats   ( _one_search_stats );
        _all_searches_stats.update ( _one_search_stats );
        return;
    }
    
    // check if no model has been computed:
    if ( !model_ok )
    {
        if ( display_degree == NOMAD::FULL_DISPLAY ) {
            std::ostringstream oss;
            oss << "end of " << NOMAD::MODEL_SEARCH
            << " (model computation or optimization error)";
            out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
        }
        stats.update_model_stats   ( _one_search_stats );
        _all_searches_stats.update ( _one_search_stats );
        return;
    }
    
    nb_search_pts = ev_control.get_nb_eval_points();
    
    // reduce the list of trial points from a maximum
    // of 4 points to MODEL_SEARCH_MAX_TRIAL_PTS points:
    int max_trial_pts = _p.get_model_search_max_trial_pts();
    if ( max_trial_pts > 4 )
        max_trial_pts = 4;
    if ( nb_search_pts > max_trial_pts )
    {
        ev_control.reduce_eval_lop ( max_trial_pts );
        nb_search_pts = ev_control.get_nb_eval_points();
        if ( display_degree == NOMAD::FULL_DISPLAY )
        {
            out << "the list of trial points is reduced to "
            << nb_search_pts << " point";
            if ( nb_search_pts > 1 )
                out << "s";
            out << std::endl;
        }
    }
    
    _one_search_stats.update_MS_max_search_pts ( nb_search_pts );
    
    // evaluate the trial points:
    // --------------------------
    int bbe        = stats.get_bb_eval();
    int sgte_eval  = stats.get_sgte_eval ();
    int cache_hits = stats.get_cache_hits();
    
    new_feas_inc = new_infeas_inc = NULL;
    
    ev_control.disable_model_eval_sort();
    
    ev_control.eval_list_of_points ( _type                   ,
                                    mads.get_true_barrier() ,
                                    mads.get_sgte_barrier() ,
                                    mads.get_pareto_front() ,
                                    stop                    ,
                                    stop_reason             ,
                                    new_feas_inc            ,
                                    new_infeas_inc          ,
                                    success                   );
    
    ev_control.enable_model_eval_sort();
    
    // update stats:
    _one_search_stats.add_MS_bb_eval    ( stats.get_bb_eval   () - bbe        );
    _one_search_stats.add_MS_sgte_eval  ( stats.get_sgte_eval () - sgte_eval  );
    _one_search_stats.add_MS_cache_hits ( stats.get_cache_hits() - cache_hits );
    
    if ( success == NOMAD::FULL_SUCCESS )
        _one_search_stats.add_MS_success();
    
    _one_search_stats.add_MS_pts ( nb_search_pts );
    
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
/*        project to mesh and create a trial point (private)     */
/*---------------------------------------------------------------*/
void NOMAD::Quad_Model_Search::create_trial_point
( NOMAD::Evaluator_Control & ev_control     ,
 NOMAD::Point               x              ,
 const NOMAD::Quad_Model  & model          ,
 NOMAD::Signature         & signature      ,
 const NOMAD::Point       & mesh_indices   ,
 const NOMAD::Point       & delta        ,
 NOMAD::dd_type             display_degree ,
 const NOMAD::Display     & out              )
{
    
    bool proj_to_mesh = _p.get_model_search_proj_to_mesh();
    
    if ( display_degree == NOMAD::FULL_DISPLAY ) {
        out << "candidate";
        if ( proj_to_mesh )
            out << " (before projection)";
        out << ": ( " << x << " )" << std::endl;
    }
    
    // model center:
    NOMAD::Point center = model.get_center();
    
    // model search point:
    int n = x.size();
    
    // projection to mesh:
    if ( proj_to_mesh )
    {
        x.project_to_mesh ( center , delta , _p.get_lb() , _p.get_ub() );
        if ( display_degree == NOMAD::FULL_DISPLAY )
            out << "candidate (after projection) : ( "
            << x << " )" << std::endl;
    }
    
    
    // Round for integer and binary variables:
    bool has_integer=false;
    bool has_binary=false;
    for (int i=0;i<n;i++)
    {
        if ( _p.get_bb_input_type()[i] == NOMAD::INTEGER )
        {
            has_integer=true;
            if ( x[i] >= 0.0 )
                x[i] = x[i].NOMAD::Double::ceil();
            else
                x[i] = x[i].NOMAD::Double::floor();
        }
        // binary variables:
        else if ( _p.get_bb_input_type()[i] == NOMAD::BINARY )
        {
            has_binary=true;
            if ( x[i]!= 0.0 )
                x[i] = 1.0;
        }
    }
    if ( has_integer && display_degree == NOMAD::FULL_DISPLAY )
        out << "candidate (after rounding integer) : ( "
        << x << " )" << std::endl;
    
    if ( has_binary && display_degree == NOMAD::FULL_DISPLAY )
        out << "candidate (after rounding binary) : ( "
        << x << " )" << std::endl;
    
    
    // compare x and center:
    if ( x == center )
    {
        if ( display_degree == NOMAD::FULL_DISPLAY )
            out << "candidate rejected (candidate==model center)" << std::endl;
        return;
    }
    
    
    
    NOMAD::Eval_Point * tk = new NOMAD::Eval_Point;
    
    // if the search is optimistic, a direction is computed (this
    // will be used in case of success in the speculative search):
    if ( _p.get_model_search_optimistic() )
    {
        NOMAD::Direction dir ( n , 0.0 , NOMAD::MODEL_SEARCH_DIR );
        dir.Point::operator = ( x - center );
        tk->set_direction  ( &dir );
    }
    
    tk->set ( n , _p.get_bb_nb_outputs() );
    tk->set_signature  ( &signature  );
    tk->Point::operator = ( x );
    
    // compute model f and h in order to accept or reject the trial point:
    NOMAD::Double h0 , f0; // model values of f and h at the center
    NOMAD::Double h1 , f1; // model values of f and h at the trial point
    
    const NOMAD::Double & h_min  = _p.get_h_min();
    NOMAD::hnorm_type     h_norm = _p.get_h_norm();
    
    model.scale ( x );
    
    model.eval_hf ( NOMAD::Point (n,0) , h_min , h_norm , h0 , f0 );
    model.eval_hf ( x                  , h_min , h_norm , h1 , f1 );
    
    if ( display_degree == NOMAD::FULL_DISPLAY )
    {
#ifdef DEBUG
        out << "model at center   : h=" << h0 << " f=" << f0 << std::endl;
#endif
        out << "model at candidate: h=" << h1 << " f=" << f1
        << std::endl << std::endl;
    }
    
    bool accept_point = true;
    
    if ( !f1.is_defined() || !h1.is_defined() )
        accept_point = false;
    else
    {
        if ( !f0.is_defined() || !h0.is_defined() )
            accept_point = true;
        else
            accept_point = (f1 <= f0) || (h1 <= h0);
    }
    
    // we check that the candidate does not correspond to another candidate:
    if ( accept_point )
    {
        const std::set<NOMAD::Priority_Eval_Point> & eval_lop
        = ev_control.get_eval_lop();
        
        std::set<NOMAD::Priority_Eval_Point>::const_iterator it , end = eval_lop.end();
        
        for ( it = eval_lop.begin() ; it != end ; ++it )
            if ( it->get_point()->NOMAD::Point::operator == ( *tk ) )
            {
                accept_point = false;
                break;
            }
    }
    
    // add the new point to the list of search trial points:
    if ( accept_point ) {
        ev_control.add_eval_point ( tk                      ,
                                   display_degree          ,
                                   _p.get_snap_to_bounds() ,
                                   NOMAD::Double()         ,
                                   NOMAD::Double()         ,
                                   f1                      ,
                                   h1                        );
#ifdef MODEL_STATS
        if ( tk ) {
            tk->set_mod_use ( 1                ); // 1 for model search
            tk->set_cond    ( model.get_cond() );
            tk->set_Yw      ( model.get_Yw  () );
            tk->set_nY      ( model.get_nY  () );
            tk->set_mh      ( h1               );
            tk->set_mf      ( f1               );
        }
#endif
        
    }
    else {
        if ( display_degree == NOMAD::FULL_DISPLAY )
            out << "candidate rejected" << std::endl;
        _one_search_stats.add_MS_rejected();
        delete tk;
    }
}

/*---------------------------------------------------------------*/
/*                    optimize a model (private)                 */
/*---------------------------------------------------------------*/
/*                                                               */
/*  . scaling puts points in [-1;1] but the model is optimized   */
/*    in [-1000;1000]                                            */
/*                                                               */
/*---------------------------------------------------------------*/
bool NOMAD::Quad_Model_Search::optimize_model
( const NOMAD::Quad_Model  & model          ,
 const NOMAD::Eval_Point ** xk             ,
 int                        i_inc          ,
 NOMAD::dd_type             display_degree ,
 const NOMAD::Display     & out            ,
 NOMAD::Point             & xf             ,
 NOMAD::Point             & xi             ,
 bool                     & stop           ,
 NOMAD::stop_type         & stop_reason      )
{
    xf.clear();
    xi.clear();
    
    int         n     = model.get_n();
    bool        error = false;
    std::string error_str;
    int         i;
    
    // initial displays:
    if ( display_degree == NOMAD::FULL_DISPLAY ) {
        std::ostringstream oss;
        oss << "model optimization";
        if ( xk[0] && xk[1] )
            oss << " (" << i_inc+1 << "/2)";
        out << std::endl << NOMAD::open_block ( oss.str() );
    }
    
    // parameters creation:
    NOMAD::Parameters model_param ( out );
    
    // number of variables:
    model_param.set_DIMENSION ( n );
    
    // blackbox outputs:
    model_param.set_BB_OUTPUT_TYPE ( _p.get_bb_output_type() );
    
    // barrier parameters:
    model_param.set_H_MIN  ( _p.get_h_min () );
    model_param.set_H_NORM ( _p.get_h_norm() );
    
    // starting points:
    {
        // 1/2: point (0,...,0):
        model_param.set_X0 ( NOMAD::Point ( n , 0.0 ) );
        
        // 2/2: model center if different than (0,..,0) and if in [-1;1]:
        NOMAD::Point x1 = model.get_center();
        
        if ( x1.size() == n && x1.is_complete() )
        {
            
            model.scale ( x1 );
            
            bool diff   = false;
            bool bnd_ok = true;
            
            for ( i = 0 ; i < n ; ++i )
            {
                
                if ( x1[i] != 0 )
                    diff = true;
                
                if ( x1[i].abs() > 1.0 )
                {
                    bnd_ok = false;
                    break;
                }
                
                x1[i] *= 1000.0;
            }
            
            if ( diff && bnd_ok )
                model_param.set_X0 ( x1 );
        }
    }
    
    // fixed variables:
    for ( i = 0 ; i < n ; ++i )
        if ( model.variable_is_fixed(i) || _p.variable_is_fixed(i) )
            model_param.set_FIXED_VARIABLE(i);
    
    // no model search and no model ordering:
    model_param.set_MODEL_SEARCH    ( false );
    model_param.set_MODEL_EVAL_SORT ( false );
    model_param.set_DIRECTION_TYPE(NOMAD::ORTHO_2N);
    
    
    // display:
    model_param.set_DISPLAY_DEGREE ( NOMAD::NO_DISPLAY );
    
    // mesh: use isotropic mesh
    model_param.set_ANISOTROPIC_MESH ( false );
    model_param.set_MESH_UPDATE_BASIS ( 4.0 );
    model_param.set_MESH_COARSENING_EXPONENT ( 1 );
    model_param.set_MESH_REFINING_EXPONENT ( -1 );
    model_param.set_INITIAL_MESH_INDEX ( 0 );
    model_param.set_INITIAL_MESH_SIZE ( NOMAD::Point ( n , 100.0 ) );
    
    // maximum number of evaluations:
    if (_p.get_nb_obj()==2)
        model_param.set_MULTI_OVERALL_BB_EVAL ( 50000 );
    else
        model_param.set_MAX_BB_EVAL ( 50000 );
    
    
    model_param.set_SNAP_TO_BOUNDS ( true );
    
    // disable user calls:
    model_param.set_USER_CALLS_ENABLED ( false );
    
    // set flags:
    bool flag_check_bimads , flag_reset_mesh , flag_reset_barriers , flag_p1_active;
    NOMAD::Mads::get_flags ( flag_check_bimads   ,
                            flag_reset_mesh     ,
                            flag_reset_barriers ,
                            flag_p1_active        );
    
    NOMAD::Mads::set_flag_check_bimads  (false  );		
    NOMAD::Mads::set_flag_reset_mesh     ( true  );
    NOMAD::Mads::set_flag_reset_barriers ( true  );
    NOMAD::Mads::set_flag_p1_active      ( false );
    
    // bounds:
    {
        NOMAD::Point lb ( n , -1.0 );
        NOMAD::Point ub ( n ,  1.0 );
        
        const NOMAD::Point & LB = _p.get_lb();
        const NOMAD::Point & UB = _p.get_ub();
        
        if ( LB.is_defined() || UB.is_defined() ) 
        {
            
            model.unscale ( lb );
            model.unscale ( ub );
            
            for ( i = 0 ; i < n ; ++i ) 
            {
                
                if ( LB[i].is_defined() && LB[i] > lb[i] )
                    lb[i] = LB[i];
                
                if ( UB[i].is_defined() && UB[i] < ub[i] )
                    ub[i] = UB[i];
            }
            
            model.scale ( lb );
            model.scale ( ub );
            
            for ( i = 0 ; i < n ; ++i ) 
            {
                if ( ub[i] < lb[i] || lb[i] > 0.0 || ub[i] < 0.0 ) 
                {
                    error     = true;
                    error_str = "optimization error: problem with bounds";
                    break; 
                }
                lb[i] *= 1000.0;
                ub[i] *= 1000.0;
            }
        }
        else 
        {
            lb *= 1000.0;
            ub *= 1000.0;
        }
        
        model_param.set_LOWER_BOUND ( lb );
        model_param.set_UPPER_BOUND ( ub );
    }
    
    if ( !error ) 
    {
        
        try
        {
            
            // parameters validation:
            model_param.check();
            
            // model evaluator creation:
            NOMAD::Evaluator *ev;
            if (model_param.get_nb_obj()==2)
                ev =new NOMAD::Multi_Obj_Quad_Model_Evaluator( model_param , model );
            else
                ev=new NOMAD::Single_Obj_Quad_Model_Evaluator(model_param, model);
            
            // algorithm creation and execution:
            NOMAD::Mads    mads ( model_param , ev );
            
            // Handle the case where nb_bb_obj>=2 but no-bimads (case of PhaseOneSearch) ---> need Phase_One_Evaluator for compute_f
            NOMAD::Phase_One_Evaluator * p1ev=NULL;
            if ( model_param.get_nb_obj() >= 2 && ! flag_check_bimads )
            {
                p1ev   = new NOMAD::Phase_One_Evaluator ( model_param , *ev );
                mads.get_evaluator_control().set_evaluator ( p1ev );
            }
            
            NOMAD::stop_type st = mads.run();
            
            delete ev;
            if (p1ev)
                delete p1ev;
            
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
            
            // display solution:
#ifdef DEBUG
            NOMAD::Display out_tmp = out;
            out_tmp.set_degrees ( NOMAD::NORMAL_DISPLAY );
            mads.display ( out_tmp );
#endif
            
            // update the stats on the number of model evaluations:
            _one_search_stats.update_MS_model_opt ( mads.get_stats().get_bb_eval() );
            
            // get the solution(s):
            const NOMAD::Eval_Point * best_feas   = mads.get_best_feasible  ();
            const NOMAD::Eval_Point * best_infeas = mads.get_best_infeasible();
            
            if ( best_feas ) {
                xf  = *best_feas;
                xf *= 0.001;
                
                if ( display_degree == NOMAD::FULL_DISPLAY ) {
                    out << "best feasible point after unscaling  : ( ";
                    xf.NOMAD::Point::display ( out );
                    out << " )" << std::endl;
                }
                
                model.unscale ( xf );
            }
            else if ( display_degree == NOMAD::FULL_DISPLAY )
                out << "no feasible solution" << std::endl;
            
            if ( best_infeas ) {
                xi  = *best_infeas;
                xi *= 0.001;
                
                if ( display_degree == NOMAD::FULL_DISPLAY ) {
                    out << "best infeasible point before unscaling: ( ";
                    xi.NOMAD::Point::display ( out );
                    out << " )" << std::endl;
                }
                
                model.unscale ( xi );
            }
            else if ( display_degree == NOMAD::FULL_DISPLAY )
                out << "no infeasible solution" << std::endl;
            
            if ( !xf.is_defined() && !xi.is_defined() ) {
                error     = true;
                error_str = "optimization error: no solution";
            }
        }
        catch ( std::exception & e ) {
            error     = true;
            error_str = std::string ( "optimization error: " ) + e.what();
        }
    }
    
    // error before run:
    else if ( display_degree == NOMAD::FULL_DISPLAY )
        out << "error before run" << std::endl;
    
    // reset flags:
    NOMAD::Mads::set_flag_check_bimads   ( flag_check_bimads   );
    NOMAD::Mads::set_flag_reset_mesh     ( flag_reset_mesh     );
    NOMAD::Mads::set_flag_reset_barriers ( flag_reset_barriers );
    NOMAD::Mads::set_flag_p1_active      ( flag_p1_active      );
    
    
    // close display block:
    if ( display_degree == NOMAD::FULL_DISPLAY ) {
        if ( error )
            out.close_block ( error_str );
        else
            out.close_block();
        out << std::endl;
    }
    
    return !error;
}
