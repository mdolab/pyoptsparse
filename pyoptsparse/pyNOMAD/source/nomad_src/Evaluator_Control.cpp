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
 \file   Evaluator_Control.cpp
 \brief  Control of the blackbox evaluations (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-15
 \see    Evaluator_Control.hpp
 */
#include "Evaluator_Control.hpp"
#include "Multi_Obj_Quad_Model_Evaluator.hpp"
#include "Single_Obj_Quad_Model_Evaluator.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
bool NOMAD::Evaluator_Control::_force_quit = false;
bool NOMAD::Evaluator_Control::_force_evaluation_failure = false;


/*---------------------------------------------------------*/
/*                       constructor                       */
/*---------------------------------------------------------*/
NOMAD::Evaluator_Control::Evaluator_Control
( const NOMAD::Parameters & p          ,
 NOMAD::Stats            & stats      ,
 NOMAD::Evaluator        * ev         ,   // can be NULL
 NOMAD::Cache            * cache      ,   // can be NULL
 NOMAD::Cache            * sgte_cache   ) // can be NULL
: _p                ( p          ) ,
_ev               ( ev         ) ,
_cache            ( cache      ) ,
_sgte_cache       ( sgte_cache ) ,
_model_eval_sort  ( true       ) ,
_del_ev           ( false      ) ,
_del_cache        ( false      ) ,
_del_sgte_cache   ( false      ) ,
#ifdef USE_MPI
_eval_in_progress ( NULL       ) ,
_nb_in_progress   ( 0          ) ,
_elop_tag         ( 0          ) ,
_slaves_elop_tags ( NULL       ) ,
_slave            ( NULL       ) ,
#endif
#ifdef USE_TGP
_last_TGP_model   ( NULL       ) ,
#endif
_stats            ( stats      ) ,
_last_stats_tag   ( -1         ) ,
_last_stats_bbe   ( -1         ) ,
_last_history_bbe ( -1         )
{
    NOMAD::Evaluator_Control::_force_quit = false;
    
    // Evaluator init:
    if ( !_ev ) {
        _ev = ( _p.get_index_obj().size() > 1 ) ? new NOMAD::Multi_Obj_Evaluator ( p ):
        new NOMAD::Evaluator           ( p );
        _del_ev = true;
    }
    
    if ( NOMAD::Slave::is_master() ) {
        
#ifdef USE_MPI
        
        int np = NOMAD::Slave::get_nb_processes();
        
        _eval_in_progress = new NOMAD::Eval_Point * [np];
        _slaves_elop_tags = new int                 [np];
        for ( int i = 0 ; i < np ; ++i ) {
            _eval_in_progress[i] = NULL;
            _slaves_elop_tags[i] = -1;
        }
        
        _slave = new NOMAD::Slave ( _p , _ev );
        
#endif
        
        const NOMAD::Display & out = _p.out();
        
        // caches creation:
        if ( !_cache ) {
            _cache     = new NOMAD::Cache ( out , NOMAD::TRUTH );
            _del_cache = true;
        }
        if ( !_sgte_cache ) {
            _sgte_cache     = new NOMAD::Cache ( out , NOMAD::SGTE );
            _del_sgte_cache = true;
        }
        
        // caches init (we only load cache file points with m blackbox outputs):
        std::string    file_name;
        int            m              = p.get_bb_nb_outputs();
        NOMAD::dd_type display_degree = out.get_gen_dd();
        
        if ( !_p.get_cache_file().empty() ) {
            file_name = _p.get_problem_dir() + _p.get_cache_file();
            if ( !_cache->load ( file_name , &m , display_degree == NOMAD::FULL_DISPLAY )
                && display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY)
                out << std::endl
                << "Warning (" << "Evaluator_Control.cpp" << ", " << __LINE__
                << "): could not load (or create) the cache file " << file_name
                << std::endl << std::endl;
        }
        
        if ( !_p.get_sgte_cache_file().empty() ) {
            file_name = _p.get_problem_dir() + _p.get_sgte_cache_file();
            if ( !_sgte_cache->load ( file_name , &m , display_degree==NOMAD::FULL_DISPLAY ) &&
                display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY )
                out << std::endl << "Warning (" << "Evaluator_Control.cpp" << ", " << __LINE__
                << "): could not load (or create) the surrogate cache file "
                << file_name << std::endl << std::endl;
        }
        
#ifdef MODEL_STATS
        if ( _p.has_model_search() ||
            ( _model_eval_sort &&
             _p.get_model_eval_sort() != NOMAD::NO_MODEL ) ) {
                out << std::endl
                << "MODEL_STATS is active. Displayed model stats are:"
                << std::endl
                << "mode ell nY wY cond";
                if ( _p.has_constraints() )
                    out << " h mh eh";
                out << " f mf ef" << std::endl
                << NOMAD::open_block()
                << "mode: model search (1) or model ordering (2)" << std::endl
                << "ell : mesh_index"                             << std::endl
                << "nY  : cardinality of the interpolation set Y" << std::endl
                << "wY  : width of Y"                             << std::endl
                << "cond: Y condition number"                     << std::endl;
                if ( _p.has_constraints() )
                    out << "h   : h value"           << std::endl
                    << "mh  : model value for h" << std::endl
                    << "eh  : relative error(%)" << std::endl;
                out << "f   : f value"             << std::endl
                << "mf  : model value for f"   << std::endl
                << "ef  : relative error(%)"   << std::endl
                << NOMAD::close_block()        << std::endl;
            }
#endif
    }
}

/*---------------------------------------------------------*/
/*                        destructor                       */
/*---------------------------------------------------------*/
NOMAD::Evaluator_Control::~Evaluator_Control ( void )
{
    if ( _del_ev )
        delete _ev;
    
    if ( _del_cache )
        delete _cache;
    
    if ( _del_sgte_cache )
        delete _sgte_cache;
    
    clear_eval_lop();
    
#ifdef USE_MPI
    
    if ( _eval_in_progress ) {
        int np = NOMAD::Slave::get_nb_processes();
        for ( int i = 0 ; i < np ; ++i )
            if ( _eval_in_progress[i] && !_eval_in_progress[i]->is_in_cache() )
                delete _eval_in_progress[i];
        delete [] _eval_in_progress;
    }
    if ( _slaves_elop_tags )
        delete [] _slaves_elop_tags;
    
    delete _slave;
    
#endif
}

/*---------------------------------------------------------*/
/*                           reset                         */
/*---------------------------------------------------------*/
void NOMAD::Evaluator_Control::reset ( void )
{
    _last_stats_tag = _last_stats_bbe = -1;
#ifdef USE_TGP
    _last_TGP_model = NULL;
#endif
}

/*---------------------------------------------------------*/
/*                     save the caches                     */
/*---------------------------------------------------------*/
bool NOMAD::Evaluator_Control::save_caches ( bool overwrite )
{
    const NOMAD::Display    & out = _p.out();
    NOMAD::dd_type display_degree = out.get_gen_dd();
    
    bool b1 = _cache->save      ( overwrite , display_degree == NOMAD::FULL_DISPLAY );
    bool b2 = _sgte_cache->save ( overwrite , display_degree == NOMAD::FULL_DISPLAY );
    
    if ( !b1 && display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY)
        out << std::endl << "Warning (" << "Evaluator_Control.cpp" << ", " << __LINE__
        << "): could not save the cache file "
        << _p.get_problem_dir() << _p.get_cache_file()
        << std::endl << std::endl;
    
    if ( !b2 && display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY)
        out << std::endl
        << "Warning (" << "Evaluator_Control.cpp" << ", " << __LINE__
        << "): could not save the surrogate cache file "
        << _p.get_problem_dir() << _p.get_sgte_cache_file()
        << std::endl << std::endl;
    return b1 && b2;
}

/*---------------------------------------------------------*/
/*    process an already evaluated Eval_Point (private)    */
/*---------------------------------------------------------*/
void NOMAD::Evaluator_Control::process_eval_point
( const NOMAD::Eval_Point & x            ,
 NOMAD::Barrier          & barrier      ,
 NOMAD::Pareto_Front     * pareto_front ) const
{
    // insertion of the Eval_Point in the barriers:
    barrier.insert(x);
    
    if ( x.get_eval_type() == NOMAD::TRUTH || _p.get_opt_only_sgte() )
    {
        
        // multi-objective:
        if ( pareto_front )
        {
            
            // insertion of the Eval_Point in the Pareto front:
            if ( x.is_feasible ( _p.get_h_min() ) &&
                pareto_front->insert ( x )       &&
                _p.get_user_calls_enabled()         )
                _ev->update_success ( _stats , x );
            
        }
        
        // single-objective: call virtual method Evaluator::update_success():
        else if ( _p.get_user_calls_enabled() &&
                 barrier.get_one_eval_succ() == NOMAD::FULL_SUCCESS )
            _ev->update_success ( _stats , x );
    }
}

/*---------------------------------------------------------*/
/*  update barrier b1 from points in barrier b2 and treat  */
/*  these points as evaluations (used in VNS search)       */
/*---------------------------------------------------------*/
NOMAD::success_type NOMAD::Evaluator_Control::process_barrier_points
( NOMAD::Barrier       & b1             ,
 const NOMAD::Barrier & b2             ,
 NOMAD::Pareto_Front  * pareto_front   ,
 NOMAD::dd_type         display_degree ,
 NOMAD::search_type     search           ) const
{
    b1.reset_success();
    
    NOMAD::Eval_Point                       *	modifiable_x;
    NOMAD::success_type							one_eval_succ;
    const NOMAD::Eval_Point                 *	last_success  = NULL;
    const std::list<const NOMAD::Eval_Point *>	& all_inserted  = b2.get_all_inserted();
    std::list<const NOMAD::Eval_Point *>::const_iterator it , end = all_inserted.end();
    for ( it = all_inserted.begin() ; it != end ; ++it )
    {
        
        // insertion in barrier:
        modifiable_x = &NOMAD::Cache::get_modifiable_point ( **it );
        
        modifiable_x->set_direction          ( NULL                              );
        modifiable_x->set_poll_center_type   ( NOMAD::UNDEFINED_POLL_CENTER_TYPE );
        modifiable_x->set_user_eval_priority ( NOMAD::Double()                   );
        modifiable_x->set_rand_eval_priority ( NOMAD::Double()                   );
        
        // process evaluation point:
        process_eval_point ( **it , b1 , pareto_front );
        
        one_eval_succ = b1.get_one_eval_succ();
        if ( one_eval_succ != NOMAD::UNSUCCESSFUL && one_eval_succ >= b1.get_success() )
            last_success = *it;
        
    }
    
    NOMAD::success_type success = b1.get_success();
    
    // display and save only the last success:
    if ( last_success && display_degree == NOMAD::FULL_DISPLAY)
        display_eval_result ( *last_success  ,
                             display_degree ,
                             search         ,
                             success        ,
                             success          );
    
    // barrier update:
    b1.update_and_reset_success();
    
    return success;
}

/*---------------------------------------------------------*/
/*      count the output stats (STAT_SUM and STAT_AVG)     */
/*      (private)                                          */
/*---------------------------------------------------------*/
void NOMAD::Evaluator_Control::count_output_stats ( const NOMAD::Eval_Point & x )
{
    const NOMAD::Point & bbo   = x.get_bb_outputs();
    int                  i_sum = _p.get_index_stat_sum();
    int                  i_avg = _p.get_index_stat_avg();
    
    // STAT_SUM:
    if ( i_sum >= 0 )
        _stats.update_stat_sum ( bbo[i_sum] );
    
    // STAT_AVG:
    if ( i_avg >= 0 )
        _stats.update_stat_avg ( bbo[i_avg] );
}

/*-------------------------------------------------------------------*/
/*                file displays for parameter STATS_FILE             */
/*-------------------------------------------------------------------*/
void NOMAD::Evaluator_Control::stats_file ( const std::string       & file_name ,
                                           const NOMAD::Eval_Point * x         ,
                                           bool                      feasible  ,
                                           const NOMAD::Point      * multi_obj   ) const
{
    std::string   fn = _p.get_problem_dir() + file_name;
    std::ofstream fout ( fn.c_str() , std::ios::app );
    
    if ( !fout.fail() )
    {
        fout.setf      ( std::ios::fixed             );
        fout.precision ( NOMAD::DISPLAY_PRECISION_BB );
        display_stats  ( false , fout , _p.get_stats_file() , x , feasible , multi_obj );
    }
    else
    {
        const NOMAD::Display & out = _p.out();
        if ( out.get_gen_dd() != NOMAD::NO_DISPLAY && out.get_gen_dd() != NOMAD::MINIMAL_DISPLAY)
            out << std::endl
            << "Warning (" << "Evaluator_Control.cpp" << ", " << __LINE__
            << "): could not save information in stats file \'"
            << file_name << "\'" << std::endl << std::endl;
    }
    fout.close();
}

/*-------------------------------------------------------------------*/
/*  display stats during Mads::run() for minimal and normal display  */
/*-------------------------------------------------------------------*/
void NOMAD::Evaluator_Control::display_stats
( bool                           header    ,
 const NOMAD::Display         & out       ,
 const std::list<std::string> & stats     ,
 const NOMAD::Eval_Point      * x         ,
 bool                           feasible  ,
 const NOMAD::Point           * multi_obj   ) const
{
    if ( stats.empty() )
    {
#ifndef R_VERSION
        out << std::endl;
#endif
        return;
    }
    
    if ( header )
    {
#ifndef R_VERSION
        out << std::endl;
#endif
    }
    
    NOMAD::Double            f;
    const NOMAD::Point     * sol		= NULL;
    const NOMAD::Point     * bbo		= NULL;
    const NOMAD::Signature * signature	= NULL;
    int                      bbe		= _stats.get_bb_eval();
    int						real_time	= _stats.get_real_time();
    int						blk_bbe		= _stats.get_block_eval();
    int                      i;
    
    // this integer is used for the default width display
    // of the various stats on the number of evaluations:
    int max_bbe = _p.get_max_bb_eval();
    if ( _p.get_max_sgte_eval() > max_bbe )
        max_bbe = _p.get_max_sgte_eval();
    if ( _p.get_max_sim_bb_eval() > max_bbe )
        max_bbe = _p.get_max_sim_bb_eval();
    if ( _p.get_max_eval() > max_bbe )
        max_bbe = _p.get_max_eval();
    
    if ( x )
    {
        signature       = x->get_signature();
        f               = (feasible) ? x->get_f() : NOMAD::INF;
        sol             = x;
        bbo             = &(x->get_bb_outputs());
        
        if (bbe < _last_stats_bbe && ! multi_obj)
            return;
        
        real_time = x->get_real_time();
        _last_stats_tag = x->get_tag();
        _last_stats_bbe = bbe;
    }
    
    
    
    std::string s1 , format;
    std::list<std::string>::const_iterator it , end = stats.end();
    for ( it = stats.begin() ; it != end ; ++it )
    {
        
        if ( it->empty() )
        {
#ifndef R_VERSION
            out << "\t";
#endif
        }
        else {
            
            if ( header )
            {
#ifndef R_VERSION
                s1 = *it;
                NOMAD::Display::extract_display_format ( s1 , format );
                out << s1;
#endif
            }
            
            else
            {
                
                // get the stats type:
                NOMAD::display_stats_type dst
                = NOMAD::Display::get_display_stats_type ( *it );
                
                // some stats types are disables in the multi-objective case:
                if ( multi_obj &&
                    ( dst == NOMAD::DS_SIM_BBE  ||
                     dst == NOMAD::DS_BBE      ||
                     dst == NOMAD::DS_SGTE     ||
                     dst == NOMAD::DS_EVAL     ||
                     dst == NOMAD::DS_TIME     ||
                     dst == NOMAD::DS_STAT_SUM ||
                     dst == NOMAD::DS_STAT_AVG    ) )
                    dst = NOMAD::DS_UNDEFINED;
                
                // display the stats:
                switch ( dst )
                {
                    case NOMAD::DS_UNDEFINED:
                        s1 = *it;
                        NOMAD::Display::extract_display_format ( s1 , format );
                        out << s1;
                        break;
                    case NOMAD::DS_OBJ:
                        if ( multi_obj )
                            display_stats_point ( out , stats , it , multi_obj );
                        else
                        {
#ifdef R_VERSION
                            {
                                std::ostringstream oss;
                                display_stats_real ( oss , f , format );
                                Rprintf ( "%s" , oss.str().c_str() );
                            }
#else
                            display_stats_real ( out , f , format );
#endif
                            format.clear();
                        }
                        break;
                    case NOMAD::DS_MESH_INDEX:
                    {
                        // display_stats_int ( out                           ,
                        //				   NOMAD::Mesh::get_mesh_index() ,
                        //				   10*L_LIMITS                   ,
                        //				   format                          );
                        
                        // format.clear();
                        
                        if ( signature )
                        {
                            NOMAD::Point mesh_indices=signature->get_mesh()->get_mesh_indices();
                            display_stats_point ( out , stats , it , &mesh_indices );
                        }
                        else
                            out << "-";
                        
                        
                        break;
                    }
                    case NOMAD::DS_DELTA_M:
                    case NOMAD::DS_MESH_SIZE:
                    {
                        if ( signature )
                        {
                            NOMAD::Point delta;
                            signature->get_mesh()->get_delta ( delta );
                            display_stats_point ( out , stats , it , &delta );
                        }
                        else
                            out << "-";
                    }
                        break;
                    case NOMAD::DS_DELTA_P:
                    case NOMAD::DS_POLL_SIZE:
                    {
                        if ( signature )
                        {
                            NOMAD::Point Delta;
                            signature->get_mesh()->get_Delta ( Delta );
                            display_stats_point ( out , stats , it , &Delta );
                            
                        }
                        else
                            out << "-";
                    }
                        break;
                    case NOMAD::DS_SIM_BBE:
                        display_stats_int ( out , _stats.get_sim_bb_eval() , max_bbe , format );
                        format.clear();
                        break;
                    case NOMAD::DS_BBE:
                        
#ifdef R_VERSION
                    {
                        std::ostringstream oss;
                        display_stats_int ( oss , bbe , max_bbe , format );
                        Rprintf ( "\t%s " , oss.str().c_str() );
                    }
#else
                    {
                        display_stats_int ( out , bbe , max_bbe , format );
                    }
#endif
                        format.clear();
                        break;
                    case NOMAD::DS_BLK_EVA:
                    {
                        display_stats_int ( out , blk_bbe , max_bbe , format );
                    }
                        format.clear();
                        break;
                        
                    case NOMAD::DS_SGTE:
                        //display_stats_int ( out , sgte_bbe , max_bbe , format );
                        display_stats_int ( out , _stats.get_sgte_eval() , max_bbe , format );
                        format.clear();
                        break;
                    case NOMAD::DS_EVAL:
                        display_stats_int ( out , _stats.get_eval() , max_bbe , format );
                        format.clear();
                        break;
                    case NOMAD::DS_TIME:
                        display_stats_int ( out , real_time , 3600 , format );
                        format.clear();
                        break;
                    case NOMAD::DS_STAT_SUM:
                        display_stats_real ( out , _stats.get_stat_sum() , format );
                        format.clear();
                        break;
                    case NOMAD::DS_STAT_AVG:
                        display_stats_real ( out , _stats.get_stat_avg() , format );
                        format.clear();
                        break;
                    case NOMAD::DS_BBO:
                        display_stats_point ( out , stats , it , bbo );
                        break;
                    case NOMAD::DS_SOL:
                        display_stats_point ( out , stats , it , sol , signature->get_input_type() );
                        break;
                    case NOMAD::DS_VAR:
                        ++it;
                        NOMAD::atoi ( *it , i );
                        if ( sol )
                            if (format.empty())
                                display_stats_type ( out , (*sol)[i] , (signature->get_input_type())[i] );
                            else
                                display_stats_real ( out , (*sol)[i] , format );
                            else
                                out << "-";
                        format.clear();
                        break;
                }
            }
        }
    }
    
    if ( !header )
#ifdef R_VERSION
        Rprintf("\n");
#else
    out << std::endl;
#endif
}

/*-----------------------------------------------------*/
/*  display a number with type                         */
/*-----------------------------------------------------*/
void NOMAD::Evaluator_Control::display_stats_type
( const NOMAD::Display        & out    ,
 const NOMAD::Double         & d      ,
 const NOMAD::bb_input_type  & bbType ) const
{
    
    // Default based on bbType
    std::string format2;
    switch (bbType)
    {
        case NOMAD::CONTINUOUS:
            format2 = "%0." + NOMAD::itos(DISPLAY_PRECISION_STD) + "g";
            break;
        case NOMAD::INTEGER || NOMAD::BINARY || NOMAD::CATEGORICAL:
            format2 = "%i";
            break;
        default:
            break;
    }
    d.display ( out , format2 );
    
}

/*-----------------------------------------------------*/
/*  display a real with DISPLAY_STATS (or STATS_FILE)  */
/*-----------------------------------------------------*/
void NOMAD::Evaluator_Control::display_stats_real
( const NOMAD::Display & out    ,
 const NOMAD::Double  & d      ,
 const std::string    & format ) const
{
    if ( format.empty() )
    {
        std::string format2 = "%0." + NOMAD::itos(DISPLAY_PRECISION_STD) + "g";
        d.display ( out , format2 );
    }
    else
        d.display ( out , format );
}


/*---------------------------------------------------------*/
/*  display an integer with DISPLAY_STATS (or STATS_FILE)  */
/*---------------------------------------------------------*/
void NOMAD::Evaluator_Control::display_stats_int
( const NOMAD::Display & out    ,
 int                    i      ,
 int                    max_i  ,
 const std::string    & format   ) const
{
    if ( format.empty() )
        out.display_int_w ( i , max_i );
    else {
        NOMAD::Double d = i;
        d.display ( out , format );
    }
}

/*---------------------------------------------------------*/
/*    display a point with DISPLAY_STATS (or STATS_FILE)   */
/*---------------------------------------------------------*/
void NOMAD::Evaluator_Control::display_stats_point
( const NOMAD::Display                      & out           ,
 const std::list<std::string>              & display_stats ,
 std::list<std::string>::const_iterator    & it            ,
 const NOMAD::Point                        * x             ,
 const std::vector<NOMAD::bb_input_type>   & bbType        ) const
{
    if ( x )
    {
        
        unsigned int n = x->size() , bbn = static_cast<int>(bbType.size());
        
        if ( bbn!=0 && n != bbn )
            throw NOMAD::Exception ( "Evaluator_Control.cpp" , __LINE__ ,
                                    "Evaluator_Control::display_stats_point(): bbType and x have different size" );
        
        
        // s1 is the string displayed befores and after
        // one coordinate (it may include format):
        std::string s1;
        if ( it != display_stats.begin() )
        {
            s1 = *(--it);
            ++it;
        }
        
        // extract the display format from s1:
        std::string format;
        if ( !s1.empty() )
            NOMAD::Display::extract_display_format ( s1 , format );
        
        // s2 is the string displayed between two coordinates:
        std::string s2;
        ++it;
        if ( it != display_stats.end() )
            s2 = *it;
        else if ( s2.empty() )
            --it;
        
        for ( unsigned int i = 0 ; i < n ; ++i )
        {
            if ( !s1.empty() && i > 0 )
                out << s1;
            
            if (bbn!=0 && format.empty())
                display_stats_type ( out , (*x)[i] , bbType[i]);
            else
                display_stats_real (out, (*x)[i] , format );
            
            if ( !s1.empty() )
                out << s1;
            if ( !s2.empty() && i < n-1  && s2.find("(VNS)")==std::string::npos && s2.find("(PhaseOne)")==std::string::npos && s2.find("(LH)")==std::string::npos && s2.find("(ExtendedPoll)")==std::string::npos )
                out << " " << s2;
            out << " ";
        }
        if ( !s2.empty() && (s2.find("(VNS)")!=std::string::npos || s2.find("(PhaseOne)")!=std::string::npos || s2.find("(LH)")!=std::string::npos || s2.find("(ExtendedPoll)")!=std::string::npos))
            out << s2;
    }
}

/*------------------------------------------*/
/*  save the solution file (SOLUTION_FILE)  */
/*------------------------------------------*/
void NOMAD::Evaluator_Control::write_solution_file ( const NOMAD::Eval_Point & x,
                                                    bool display_bimv) const
{
    const std::string & sol_file = _p.get_solution_file();
    if ( !sol_file.empty() && ( x.is_feasible ( _p.get_h_min() ) || display_bimv ) )
        write_sol_or_his_file ( _p.get_problem_dir() + sol_file , x , true , display_bimv );
}

/*----------------------------------------------*/
/*     save the solution file  (SOLUTION_FILE)  */
/*  or update the history file (HISTORY_FILE )  */
/*  (private)                                   */
/*----------------------------------------------*/
void NOMAD::Evaluator_Control::write_sol_or_his_file
( const std::string       & file_name ,
 const NOMAD::Eval_Point & x         ,
 bool                      is_sol    ,
 bool						display_bimv ) const
{
    // if is_sol == true: save the solution file
    //              else: update the history file
    bool          failed = false;
    std::ofstream fout;
    
    if ( is_sol )
        fout.open ( file_name.c_str() );
    else
        fout.open ( file_name.c_str() , std::ios::app );
    
    if ( !fout.fail() ) {
        
        fout.setf      ( std::ios::fixed             );
        fout.precision ( NOMAD::DISPLAY_PRECISION_BB );
        
        // solution display:
        if ( is_sol )
        {
            if ( _p.get_bb_input_include_seed() )
                fout << _p.get_seed() << std::endl;
            if ( _p.get_bb_input_include_tag() )
                fout << x.get_tag() << std::endl;
            x.Point::display ( fout , "\n" , -1 , -1 );
            if (display_bimv)
                fout << std::endl << "warning: best infeasible solution (min. violation)";
            fout << std::endl;
        }
        
        // history display:
        else {
            x.Point::display ( fout , " " , -1 , -1 );
            fout << " ";
            x.get_bb_outputs().Point::display ( fout , " " , -1 , -1 );
            fout << std::endl;
        }
        
        if ( fout.fail() )
            failed = true;
    }
    else
        failed = true;
    
    fout.close();
    
    if ( failed && _p.out().get_gen_dd() != NOMAD::NO_DISPLAY &&  _p.out().get_gen_dd() != NOMAD::MINIMAL_DISPLAY)
        _p.out() << std::endl
        << "Warning (" << "Evaluator_Control.cpp" << ", " << __LINE__
        << "): could not "
        << ( is_sol ? "save the current solution" :
            "update the history" )
        << " in \'"
        << file_name << "\'" << std::endl << std::endl;
}


/*---------------------------------------------------------*/
/*             display evaluation result (private)         */
/*---------------------------------------------------------*/
void NOMAD::Evaluator_Control::display_eval_result
( const NOMAD::Eval_Point & x                ,
 NOMAD::dd_type            display_degree   ,
 NOMAD::search_type        search           ,
 NOMAD::success_type       one_eval_success ,
 NOMAD::success_type       success            ) const
{
    const NOMAD::Display & out = _p.out();
    int cur_bbe;
    
    // surrogate evaluation:
    if ( x.get_eval_type() == NOMAD::SGTE )
    {
        if ( display_degree == NOMAD::FULL_DISPLAY )
        {
            out << std::endl << "point #" << x.get_tag() << " sgte eval: ";
            if ( x.is_eval_ok() )
            {
                out << "h=";
                if ( x.get_h().is_defined() )
                    out << x.get_h();
                else
                    out << "inf (extr. barrier)";
                out << " f=" << x.get_f();
            }
            else
                out << "failed";
            out << std::endl;
        }
        if ( !_p.get_opt_only_sgte() )
            return;
        
        cur_bbe = _stats.get_sgte_eval();
    }
    else
        cur_bbe = _stats.get_eval();
    
    const std::string & stats_file_name = _p.get_stats_file_name();
    bool                feas_x          = x.is_feasible ( _p.get_h_min() );
    
    // update the history file:
    // (contains surrogate evaluations if opt_only_sgte==true)
    const std::string & his_file = _p.get_history_file();
    if ( !his_file.empty() && cur_bbe > _last_history_bbe)
    {
        write_sol_or_his_file ( _p.get_problem_dir() + his_file , x , false );
        _last_history_bbe = cur_bbe;
    }
    
    // success displays:
    if ( one_eval_success != NOMAD::UNSUCCESSFUL &&
        one_eval_success >= success )
    {
        
        // save the current solution in file:
        write_solution_file ( x );
        
        bool ds_ok = ( cur_bbe > _last_stats_bbe)	&&
        ( _p.get_display_all_eval()		||
         ( one_eval_success == NOMAD::FULL_SUCCESS && feas_x ) );
        
        // normal display and minimal:
        if ( (display_degree == NOMAD::NORMAL_DISPLAY || display_degree == NOMAD::MINIMAL_DISPLAY ) && ds_ok )
            display_stats ( false , out , _p.get_display_stats() , &x , feas_x , NULL );
        // detailed display:
        else if ( display_degree == NOMAD::FULL_DISPLAY )
            out << std::endl << search << " " << one_eval_success
            << " point " << x;
        
        // stats file:
        if ( ds_ok && !stats_file_name.empty() )
            stats_file ( stats_file_name , &x , feas_x , NULL );
        
    }
    else
    {
        
        if ( display_degree == NOMAD::FULL_DISPLAY )
        {
            out << search << " " << one_eval_success
            << " point #" << x.get_tag();
            if ( x.is_eval_ok() )
                out << " [ h=" << x.get_h()
                << " f=" << x.get_f() << " ]" << std::endl;
            else if (x.check_rejected())
                out << ": evaluation rejected by user (this may alter convergence properties!)" << std::endl;
            else
                out << ": evaluation failed (you may need to check the source of the problem)." << std::endl;
        }
        
        if ( _p.get_display_all_eval() && cur_bbe > _last_stats_bbe  )
        {
            
            if ( display_degree == NOMAD::NORMAL_DISPLAY || display_degree == NOMAD::MINIMAL_DISPLAY )
                display_stats ( false , out , _p.get_display_stats() , &x , feas_x , NULL );
            
            if ( !stats_file_name.empty() )
                stats_file ( stats_file_name , &x , feas_x , NULL );
        }
    }
}


/*-------------------------------------------*/
/*        search a point in the cache        */
/*-------------------------------------------*/
/* . return true if the point is in cache    */
/* . private method                          */
/*-------------------------------------------*/
bool NOMAD::Evaluator_Control::cache_check
( const NOMAD::Eval_Point *& x              ,
 NOMAD::Barrier           & true_barrier   ,
 NOMAD::Barrier           & sgte_barrier   ,
 NOMAD::Pareto_Front      * pareto_front   ,
 bool                     & count_eval     ,
 const NOMAD::Double      & h_max          ,
 NOMAD::dd_type             display_degree   ) const
{
    NOMAD::eval_type          x_eval_type = x->get_eval_type();
    const NOMAD::Eval_Point * cache_x     = NULL;
    
    // first cache check:
    if ( x->is_in_cache() )
        cache_x = x;
    // second cache check:
    else
        cache_x = ( ( x->get_eval_type() == NOMAD::TRUTH ) ?
                   _cache : _sgte_cache )->find ( *x );
    
    // cache hit: transfer some data from x to cache_x:
    if ( cache_x )
    {
        
        if ( x_eval_type != cache_x->get_eval_type() )
            throw NOMAD::Exception ( "Evaluator_Control.cpp" , __LINE__ ,
                                    "Evaluator_Control::cache_check(): eval and cache pts have different eval_type" );
        
        if ( cache_x->is_eval_ok() )
        {
            
            NOMAD::Eval_Point * modifiable_cache_x
            = &NOMAD::Cache::get_modifiable_point ( *cache_x );
            
            // if wrong number of outputs, reset cache_x._bb_outputs:
            {
                int m = _p.get_bb_nb_outputs();
                if ( cache_x->get_bb_outputs().size() != m )
                    modifiable_cache_x->set_bb_output ( NOMAD::Point ( m ) );
            }
            
            modifiable_cache_x->set_signature          ( x->get_signature         () );
            modifiable_cache_x->set_direction          ( x->get_direction         () );
            // The point in cache is updated for the poll center to correspond to the new poll center of x (important for poll reduction)
            modifiable_cache_x->set_poll_center        ( x->get_poll_center       () );
            modifiable_cache_x->set_poll_center_type   ( x->get_poll_center_type  () );
            modifiable_cache_x->set_user_eval_priority ( x->get_user_eval_priority() );
            modifiable_cache_x->set_rand_eval_priority ( x->get_rand_eval_priority() );
            
#ifdef MODEL_STATS
            modifiable_cache_x->set_model_data ( *x );
#endif
            
            // set_f, set_h, and set_EB_ok:
            _ev->compute_f ( *modifiable_cache_x );
            _ev->compute_h ( *modifiable_cache_x );
        }
    }
    
    // point in cache but evaluation is to be made again:
    if ( cache_x && cache_x->is_eval_ok() &&
        ( !cache_x->get_f().is_defined() ||
         ( cache_x->is_EB_ok()                      &&
          !cache_x->get_bb_outputs().is_complete() &&
          cache_x->get_h().is_defined()            &&
          cache_x->get_h() < h_max                    ) ) )
    {
        x       = cache_x;
        cache_x = NULL;
    }
    
    // point in cache:
    if ( cache_x )
    {
        
        _stats.add_cache_hit();
        
        // displays:
        if ( display_degree == NOMAD::FULL_DISPLAY )
        {
            const NOMAD::Display & out = _p.out();
            if ( cache_x->get_eval_type() == NOMAD::SGTE )
                out << "surrogate ";
            out << "cache hit: #" << x->get_tag()
            << " --> #" << cache_x->get_tag() << std::endl;
        }
        
        // process the Eval_Point taken in cache:
        process_eval_point ( *cache_x ,
                            ( cache_x->get_eval_type() == NOMAD::TRUTH ) ?
                            true_barrier : sgte_barrier ,
                            pareto_front );
        
        // count the (simulated) bb eval ?
        int index_cnt_eval = _p.get_index_cnt_eval();
        if ( index_cnt_eval >= 0 && cache_x->get_bb_outputs()[index_cnt_eval] == 0.0 )
            count_eval = false;
        
        x = cache_x;
        
        return true;
    }
    
    return false;
}

/*----------------------------------------------------*/
/*                 eval a point (private)             */
/*----------------------------------------------------*/
void NOMAD::Evaluator_Control::eval_point ( NOMAD::Eval_Point       & x            ,
                                           NOMAD::Barrier          & true_barrier ,
                                           NOMAD::Barrier          & sgte_barrier ,
                                           NOMAD::Pareto_Front     * pareto_front ,
                                           bool                    & count_eval   ,
                                           bool                    & stop         ,
                                           NOMAD::stop_type        & stop_reason  ,
                                           const NOMAD::Double     & h_max          )
{
    int max_bb_eval   = _p.get_max_bb_eval();
    int max_sgte_eval = _p.get_max_sgte_eval();
    
    // blackbox or surrogate evaluations are allowed:
    if ( ( x.get_eval_type() == NOMAD::TRUTH && max_bb_eval   != 0 ) ||
        ( x.get_eval_type() == NOMAD::SGTE  && max_sgte_eval != 0 )    )
    {
        
        NOMAD::Eval_Point * eval_x = &NOMAD::Cache::get_modifiable_point ( x );
        
        // get the signature:
        NOMAD::Signature * signature = x.get_signature();
        if ( !signature )
            throw NOMAD::Exception ( "Evaluator_Control.cpp" , __LINE__ ,
                                    "Evaluator_Control::eval_point(): the point has no signature" );
        
        // evaluation of the point:
        // ------------------------
        bool eval_ok = true;
        
        NOMAD::Evaluator_Control::_force_evaluation_failure=false;
        {
            // 1. scaling:
            bool do_scaling = signature->get_scaling().is_defined();
            if ( do_scaling )
                eval_x->scale();
            
            // 2.1. evaluation:
            try
            {
                eval_ok = _ev->eval_x ( *eval_x , h_max , count_eval );
                
            }
            catch ( exception & e )
            {
                throw NOMAD::Exception ( "Evaluator_control.cpp" , __LINE__ , e.what() );
            }

            
            
            // 2.2. check the nan's:
            if ( eval_ok && eval_x->check_nan() )
                eval_ok = false;
            
            if ( _force_evaluation_failure )
                eval_ok = false;
            
            // 3. unscaling:
            if ( do_scaling )
                eval_x->unscale();
        }
        
        if ( eval_ok )
        {
            
            eval_x->set_eval_status ( NOMAD::EVAL_OK );
            
            // set_f, set_h and set_EB_ok:
            _ev->compute_f ( *eval_x );
            _ev->compute_h ( *eval_x );
            
        }
        else
        {
            eval_x->set_eval_status ( NOMAD::EVAL_FAIL );
            _stats.add_failed_eval();
        }
        
        // insertion in cache even if is_eval_ok == false:
        if ( !x.is_in_cache() )
        {
            
            int size_before , size_after;
            
            if ( x.get_eval_type() == NOMAD::SGTE ) {
                size_before = _sgte_cache->size();
                _sgte_cache->insert(x);
                size_after  = _sgte_cache->size();
            }
            else {
                size_before = _cache->size();
                _cache->insert(x);
                size_after  = _cache->size();
            }
            
            if ( size_after == size_before )
                x.set_in_cache ( false );
        }
        
    }
    
}


/*----------------------------------------------------*/
/*                 eval points in a list (private)    */
/*----------------------------------------------------*/
void NOMAD::Evaluator_Control::eval_points ( std::list<NOMAD::Eval_Point *>	& list_eval		,
                                            NOMAD::Barrier					& true_barrier	,
                                            NOMAD::Barrier					& sgte_barrier	,
                                            NOMAD::Pareto_Front				* pareto_front	,
                                            std::list<bool>					& count_list_eval,
                                            bool							& stop			,
                                            NOMAD::stop_type				& stop_reason	,
                                            const NOMAD::Double				& h_max          )
{
    int max_bb_eval   = _p.get_max_bb_eval();
    int max_sgte_eval = _p.get_max_sgte_eval();
    
    std::list<NOMAD::Eval_Point*>::iterator it_begin=list_eval.begin();
    
    // blackbox or surrogate evaluations are allowed:
    if ( ( (*it_begin)->get_eval_type() == NOMAD::TRUTH && max_bb_eval   != 0 ) ||
        ( (*it_begin)->get_eval_type() == NOMAD::SGTE  && max_sgte_eval != 0 )    )
    {
        
        // 1. Pre-evaluation tests and scaling
        for ( std::list<NOMAD::Eval_Point*>::iterator it=it_begin;it!=list_eval.end();++it)
        {
            // get the signature:
            NOMAD::Signature * signature = (*it)->get_signature();
            if ( !signature )
                throw NOMAD::Exception ( "Evaluator_Control.cpp" , __LINE__ ,
                                        "Evaluator_Control::eval_points(): the point has no signature" );
            
            // Scaling before evaluation of the points:
            bool do_scaling = signature->get_scaling().is_defined();
            if ( do_scaling )
                (*it)->scale();
            
        }
        
        // 2. list evaluation:
        bool eval_list_ok = true;
        NOMAD::Evaluator_Control::_force_evaluation_failure=false;
        
        try
        {
            _ev->eval_x ( list_eval , h_max,count_list_eval );
        }
        catch ( exception & e )
        {
            throw NOMAD::Exception ( "Evaluator_control.cpp" , __LINE__ , e.what() );
        }
        
        if ( _force_evaluation_failure )
            eval_list_ok = false;
        
        
        // One block of evaluations is counted
        if ( eval_list_ok )
            _stats.add_one_block_eval();
        
        
        
        // 3. Post list evaluation checks and operation
        std::list<bool>::iterator it_count=count_list_eval.begin();
        for ( std::list<NOMAD::Eval_Point*>::iterator it=it_begin;it!=list_eval.end();++it,++it_count)
        {
            bool eval_ok=true;
            bool eval_rejected=false;
            
            // 3.1. check the nan's and list evaluation failure:
            if ( !eval_list_ok || (*it)->check_nan() )
                eval_ok = false;
            
            if ((*it)->check_rejected())
            {
                eval_rejected=true;
                eval_ok=false;
            }
            
            // 3.2 unscaling:
            if ( (*it)->get_signature()->get_scaling().is_defined() )
                (*it)->unscale();
            
            
            if ( eval_ok )
            {
                (*it)->set_eval_status ( NOMAD::EVAL_OK );
                
                // set_f, set_h and set_EB_ok:
                _ev->compute_f ( *(*it) );
                _ev->compute_h ( *(*it));
                
            }
            else if (!eval_rejected)
            {
                (*it)->set_eval_status ( NOMAD::EVAL_FAIL );
                _stats.add_failed_eval();
            } // Do nothing if eval has been rejected
            
            // insertion in cache even if is_eval_ok == false. Exception: a point that has been rejected by user is not put in the cache.
            if ( !(*it)->is_in_cache() && !eval_rejected )
            {
                
                int size_before , size_after;
                
                if ( (*it)->get_eval_type() == NOMAD::SGTE )
                {
                    size_before = _sgte_cache->size();
                    _sgte_cache->insert(*(*it));
                    size_after  = _sgte_cache->size();
                }
                else
                {
                    size_before = _cache->size();
                    _cache->insert(*(*it));
                    size_after  = _cache->size();
                }
                
                if ( size_after == size_before )
                    (*it)->set_in_cache ( false );
            }
            
            
            // count the output stats (STAT_SUM and STAT_AVG):
            if ( (_p.check_stat_sum() || _p.check_stat_avg()) && !eval_rejected)
                count_output_stats(*(*it));
            
        }
    }
}



/*-------------------------------------------*/
/*      check stopping criteria (private)    */
/*-------------------------------------------*/
void NOMAD::Evaluator_Control::check_stopping_criteria
( NOMAD::search_type        search      ,
 bool                      count_eval  ,
 const NOMAD::Eval_Point & x           ,
 bool                    & stop        ,
 NOMAD::stop_type        & stop_reason   ) const
{
    // check the time:
    if ( !stop                 &&
        _p.get_max_time() > 0 &&
        _stats.get_real_time() >= _p.get_max_time() )
    {
        stop        = true;
        stop_reason = NOMAD::MAX_TIME_REACHED;
    }
    
    // count an evaluation or a simulated blackbox evaluation:
    if ( x.get_eval_type() == NOMAD::TRUTH )
    {
        _stats.add_eval();
        if ( count_eval && !x.get_current_run() )
            _stats.add_sim_bb_eval();
    }
    
    
    // check STAT_SUM_TARGET:
    if ( !stop	&&
        (_p.check_stat_sum() || _p.check_stat_avg()))
    {
        
        NOMAD::Double sum_target = _p.get_stat_sum_target();
        if ( sum_target.is_defined() )
        {
            NOMAD::Double sum = _stats.get_stat_sum();
            if ( sum.is_defined() && sum >= sum_target )
            {
                stop        = true;
                stop_reason = NOMAD::STAT_SUM_TARGET_REACHED;
            }
        }
    }
    
    // check the number of blackbox evaluations:
    if ( !stop )
    {
        int max_bb_eval   = _p.get_max_bb_eval();
        int max_sgte_eval = _p.get_max_sgte_eval();
        if ( max_bb_eval > 0 && _stats.get_bb_eval() >= max_bb_eval )
        {
            stop        = true;
            stop_reason = NOMAD::MAX_BB_EVAL_REACHED;
        }
        if ( max_sgte_eval > 0 && _stats.get_sgte_eval() >= max_sgte_eval )
        {
            stop        = true;
            stop_reason = NOMAD::MAX_SGTE_EVAL_REACHED;
        }
    }
    
    // check the stopping condition MAX_EVAL:
    if ( !stop                 &&
        _p.get_max_eval() > 0 &&
        _stats.get_eval() >= _p.get_max_eval() )
    {
        stop        = true;
        stop_reason = NOMAD::MAX_EVAL_REACHED;
    }
    
    // check the stopping condition MAX_SIM_BB_EVAL:
    if ( !stop                         &&
        _p.get_max_sim_bb_eval() >  0 &&
        _stats.get_sim_bb_eval() >= _p.get_max_sim_bb_eval() )
    {
        stop        = true;
        stop_reason = NOMAD::MAX_SIM_BB_EVAL_REACHED;
    }
    
    // check the stopping conditions F_TARGET and FEAS_REACHED
    // (for phase one: the evaluations must stop if all EB
    //  constraints are satisfied, but some PB constraints can
    //  be violated)
    if ( !stop          &&
        x.is_eval_ok() &&
        ( _p.get_opt_only_sgte() ||
         x.get_eval_type() == NOMAD::TRUTH ) )
    {
        
        bool feasible = x.is_feasible ( _p.get_h_min() );
        
        // check FEAS_REACHED:
        if ( feasible && _p.get_stop_if_feasible() )
        {
            stop        = true;
            stop_reason = NOMAD::FEAS_REACHED;
        }
        
        // check F_TARGET:
        {
            const NOMAD::Point           & f_target       = _p.get_f_target();
            const std::list<int>         & index_obj      = _p.get_index_obj();
            std::list<int>::const_iterator index_obj_end  = index_obj.end();
            bool                           check_f_target = f_target.is_defined();
            int                            nb_to_check    = (check_f_target) ?
            f_target.nb_defined() : 0;
            
            if ( check_f_target && ( feasible || search == NOMAD::LH_SEARCH_P1 ) )
            {
                const NOMAD::Point & bbo = x.get_bb_outputs();
                bool                 chk = true;
                int                  k   = 0;
                int                  cnt = 0;
                for ( std::list<int>::const_iterator it = index_obj.begin();
                     it != index_obj_end ; ++it , ++k )
                {
                    if ( bbo[*it].is_defined() && f_target[k].is_defined() )
                    {
                        if ( f_target[k] < bbo[*it] )
                        {
                            chk = false;
                            break;
                        }
                        cnt++;
                    }
                }
                
                if ( chk && cnt == nb_to_check )
                {
                    stop        = true;
                    stop_reason = NOMAD::F_TARGET_REACHED;
                }
            }
        }
    }
}

/*-------------------------------------------------------*/
/*  receive an evaluation result from a slave (private)  */
/*-------------------------------------------------------*/
#ifdef USE_MPI
void NOMAD::Evaluator_Control::receive_eval_result
( NOMAD::search_type    search       ,
 NOMAD::Eval_Point   * x            ,
 NOMAD::Barrier      & true_barrier ,
 NOMAD::Barrier      & sgte_barrier ,
 NOMAD::Pareto_Front * pareto_front ,
 int                   slave_rank   ,
 bool                & stop         ,
 NOMAD::stop_type    & stop_reason    )
{
    bool eval_ok , count_eval;
    
    
    // receive the evaluation result:
    _slave->receive_eval_result ( slave_rank , x , eval_ok , count_eval );
    
    // process the evaluation:
    if ( eval_ok ) {
        
        // set_f, set_h and set_EB_ok:
        _ev->compute_f ( *x );
        _ev->compute_h ( *x );
        
        // process the evaluated point:
        process_eval_point ( *x                                    ,
                            (x->get_eval_type()==NOMAD::TRUTH) ?
                            true_barrier : sgte_barrier           ,
                            pareto_front                            );
    }
    else
        _stats.add_failed_eval();
    
    // insertion in cache even if !eval_ok:
    if ( !x->is_in_cache() )
        ( ( x->get_eval_type() == NOMAD::SGTE ) ?
         _sgte_cache : _cache)->insert ( *x );
    
    // count the bb evaluation:
    if ( count_eval ) {
        if ( x->get_eval_type() == NOMAD::SGTE )
            _stats.add_sgte_eval();
        else
            _stats.add_bb_eval();
    }
    
    // count the output stats (STAT_SUM and STAT_AVG):
    if ( _p.check_stat_sum() || _p.check_stat_avg() ) {
        
        count_output_stats ( *x );
        
        // check STAT_SUM_TARGET:
        NOMAD::Double sum_target = _p.get_stat_sum_target();
        if ( sum_target.is_defined() ) {
            NOMAD::Double sum = _stats.get_stat_sum();
            if ( !stop && sum.is_defined() && sum >= sum_target ) {
                stop        = true;
                stop_reason = NOMAD::STAT_SUM_TARGET_REACHED;
            }
        }
    }
    
    // check stopping criteria:
    if ( !stop ) {
        
        int max_bb_eval   = _p.get_max_bb_eval();
        int max_sgte_eval = _p.get_max_sgte_eval();
        
        if ( max_bb_eval > 0 && _stats.get_bb_eval() >= max_bb_eval ) {
            stop        = true;
            stop_reason = NOMAD::MAX_BB_EVAL_REACHED;
        }
        if ( max_sgte_eval > 0 && _stats.get_sgte_eval() >= max_sgte_eval ) {
            stop        = true;
            stop_reason = NOMAD::MAX_SGTE_EVAL_REACHED;
        }
    }
    
    check_stopping_criteria ( search , count_eval , *x , stop , stop_reason );
}
#endif

/*----------------------------------------------------*/
/*          wait for evaluations in progress          */
/*----------------------------------------------------*/
#ifdef USE_MPI
void NOMAD::Evaluator_Control::wait_for_evaluations
( NOMAD::search_type                     search         ,
 NOMAD::Barrier                       & true_barrier   ,
 NOMAD::Barrier                       & sgte_barrier   ,
 NOMAD::Pareto_Front                  * pareto_front   ,
 bool                                 & stop           ,
 NOMAD::stop_type                     & stop_reason    ,
 NOMAD::success_type                  & success        ,
 std::list<const NOMAD::Eval_Point *> & evaluated_pts    )
{
    if ( _nb_in_progress == 0 )
        return;
    
    // display degree:
    const NOMAD::Display    & out = _p.out();
    NOMAD::dd_type display_degree = out.get_display_degree ( search );
    
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out << std::endl
        << NOMAD::open_block ( "wait for evaluations" );
    
    NOMAD::Barrier     & barrier = ( _p.get_opt_only_sgte() ) ?
    sgte_barrier : true_barrier;
    char                 signal;
    int                  source;
    NOMAD::Eval_Point  * eval_x;
    NOMAD::success_type  one_eval_success;
    
    while ( _nb_in_progress > 0 )
    {
        
        source = NOMAD::Slave::receive_signal ( signal );
        eval_x = _eval_in_progress[source];
        
        if ( eval_x )
        {
            
            if ( display_degree == NOMAD::FULL_DISPLAY )
                out << std::endl << "receive eval point #" << eval_x->get_tag()
                << " from slave " << source << std::endl << std::endl;
            
            receive_eval_result ( search       ,
                                 eval_x       ,
                                 true_barrier ,
                                 sgte_barrier ,
                                 pareto_front ,
                                 source       ,
                                 stop         ,
                                 stop_reason    );
            
            // list of processed points:
            if ( eval_x->is_in_cache() )
                evaluated_pts.push_back ( eval_x );
            
            // success:
            one_eval_success = barrier.get_one_eval_succ();
            success          = barrier.get_success();
            
            // asynchronous success count:
            if ( success == NOMAD::FULL_SUCCESS &&
                _elop_tag != _slaves_elop_tags[source] )
                _stats.add_asynchronous_success();
            
            // displays:
            display_eval_result ( *eval_x          ,
                                 display_degree   ,
                                 search           ,
                                 one_eval_success ,
                                 success            );
            
            if ( !_eval_in_progress[source]->is_in_cache() )
                delete _eval_in_progress[source];
            _eval_in_progress[source] = NULL;
            _slaves_elop_tags[source] = -1;
            --_nb_in_progress;
            
            // force quit (by pressing ctrl-c):
            if ( !stop && ( NOMAD::Evaluator_Control::_force_quit || NOMAD::Evaluator::get_force_quit() ) )
            {
                stop        = true;
                stop_reason = NOMAD::CTRL_C;
                break;
            }
            
            if ( stop && ( stop_reason==NOMAD::ERROR ||
                          stop_reason==NOMAD::UNKNOWN_STOP_REASON ) )
                break;
        }
        else
            NOMAD::Slave::send_signal ( NOMAD::WAIT_SIGNAL , source );
    }
    
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out.close_block();
}
#endif

/*----------------------------------------------------------------*/
/*  check if the evaluation at this point is already in progress  */
/*  (private)                                                     */
/*----------------------------------------------------------------*/
#ifdef USE_MPI
bool NOMAD::Evaluator_Control::already_in_progress
( const NOMAD::Eval_Point & x ) const
{
    if ( _eval_in_progress ) {
        
        int x_tag = x.get_tag();
        int np    = NOMAD::Slave::get_nb_processes();
        
        for ( int i = 0 ; i < np ; ++i )
            if ( _eval_in_progress[i] &&
                ( _eval_in_progress[i]->get_tag() == x_tag ||
                 _eval_in_progress[i]->Point::operator == ( x ) ) )
                return true;
    }
    return false;
}
#endif

/*----------------------------------------------------------------*/
/*     eval_list_of_points, private version (parallel version)    */
/*----------------------------------------------------------------*/
#ifdef USE_MPI
void NOMAD::Evaluator_Control::private_eval_list_of_points
( NOMAD::search_type              search         ,   // IN     : search type
 NOMAD::Barrier                & true_barrier   ,   // IN/OUT : the barrier
 NOMAD::Barrier                & sgte_barrier   ,   // IN/OUT : the surrogate barrier
 NOMAD::Pareto_Front           * pareto_front   ,   // IN/OUT : the Pareto front
 //          (can be NULL)
 bool                          & stop           ,   // IN/OUT : stopping criterion
 NOMAD::stop_type              & stop_reason    ,   // OUT    : stopping reason
 const NOMAD::Eval_Point      *& new_feas_inc   ,   // OUT    : new feasible incumbent
 const NOMAD::Eval_Point      *& new_infeas_inc ,   // OUT    : new infeas. incumbent
 NOMAD::success_type           & success        ,   // OUT    : type of success
 std::list<const NOMAD::Eval_Point *>
 & evaluated_pts    ) // OUT    : list of processed pts
{
    if ( stop || _eval_lop.empty() )
    {
        stop_reason = NOMAD::UNKNOWN_STOP_REASON;
        ++_elop_tag;
        return;
    }
    
    evaluated_pts.clear();
    
    
    // initial display:
    const NOMAD::Display    & out = _p.out();
    NOMAD::dd_type display_degree = out.get_display_degree ( search );
    
    if ( display_degree == NOMAD::FULL_DISPLAY )
    {
        std::ostringstream msg;
        msg << "list of points evaluation (" << search << ")";
        out << std::endl << NOMAD::open_block ( msg.str() );
    }
    
    // call the Evaluator (virtual) preprocessing of a list of points:
    _ev->list_of_points_preprocessing ( _eval_lop );
    
    const NOMAD::Eval_Point * old_feasible_incumbent   = NULL;
    const NOMAD::Eval_Point * old_infeasible_incumbent = NULL;
    
    // active barrier:
    NOMAD::Barrier & barrier = ( _p.get_opt_only_sgte() ) ?
    sgte_barrier : true_barrier;
    
    old_feasible_incumbent   = barrier.get_best_feasible();
    old_infeasible_incumbent = barrier.get_best_infeasible();
    
    NOMAD::Double f0;
    if ( _p.get_opportunistic_min_f_imprvmt().is_defined() &&
        old_feasible_incumbent )
        f0 = old_feasible_incumbent->get_f();
    
    new_feas_inc   = NULL;
    new_infeas_inc = NULL;
    stop           = false;
    success        = NOMAD::UNSUCCESSFUL;
    stop_reason    = NOMAD::NO_STOP;
    
    const NOMAD::Eval_Point  * x;
    NOMAD::check_failed_type   check_failed_reason;
    bool                       count_eval;
    std::vector<const NOMAD::Eval_Point *>
    to_be_evaluated;
    NOMAD::success_type        one_eval_success;
    bool                       one_for_luck = false;
    bool                       opp_stop     = false;
    int                        init_nb_eval = _stats.get_eval();
    int                        nb_success   = 0;
    int                        k            = 0;
    int                        nb_points    = static_cast<int> ( _eval_lop.size() );
    int                        max_bb_eval  = _p.get_max_bb_eval();
    
    // loop #1: search in cache:
    // -------------------------
    std::set<NOMAD::Priority_Eval_Point>::iterator
    it  = _eval_lop.begin() ,
    end = _eval_lop.end();
    while ( !stop && !opp_stop && it != end )
    {
        
        x = it->get_point();
        
        x->set_current_run ( true );
        
        // displays:
        if ( display_degree == NOMAD::FULL_DISPLAY )
        {
            
            // open the evaluation block:
            {
                std::ostringstream oss;
                if ( x->get_eval_type() == NOMAD::SGTE )
                    oss << "surrogate ";
                oss << "evaluation " << k+1 << "/" << nb_points;
                out << std::endl << NOMAD::open_block ( oss.str() );
            }
            
            out << std::endl << "point #" << x->get_tag() << "   ( ";
            x->Point::display ( out , " " , 2 , NOMAD::Point::get_display_limit() );
            out << " )" << std::endl;
            if ( x->get_direction() )
            {
                out << "direction    : " << *x->get_direction()  << std::endl;
                NOMAD::Point delta;
                x->get_signature()->get_mesh()->get_delta(delta);
                out << "direction d : ( " << *x->get_direction()/delta  << " )"  << std::endl;
            }
            if ( x->get_signature() )
                out << "mesh indices:  ( " << x->get_signature()->get_mesh()->get_mesh_indices() << " )" << std::endl;
            out << std::endl;
            
            
            
        }
        
        // check if the evaluation at this point is already in progress:
        if ( !already_in_progress ( *x ) )
        {
            
            // current point check (# of bb outputs, bounds, integer values, fixed-vars):
            if ( x->check ( _p.get_bb_nb_outputs() , check_failed_reason ) )
            {
                
                count_eval = true;
                
                // point in cache:
                if ( cache_check ( x                   ,
                                  true_barrier        ,
                                  sgte_barrier        ,
                                  pareto_front        ,
                                  count_eval          ,
                                  barrier.get_h_max() ,
                                  display_degree        ) )
                    
                {
                    
                    // list of processed points:
                    evaluated_pts.push_back ( x );
                    
                    // check stopping criteria:
                    check_stopping_criteria ( search , count_eval , *x , stop , stop_reason );
                    
                    // success:
                    one_eval_success = barrier.get_one_eval_succ();
                    success          = barrier.get_success();
                    
                    // displays:
                    display_eval_result ( *x               ,
                                         display_degree   ,
                                         search           ,
                                         one_eval_success ,
                                         success            );
                    
                    // stop the evaluations (opportunistic strategy) ?
                    if ( stop_evaluations ( *x               ,
                                           search           ,
                                           k                ,
                                           nb_points        ,
                                           stop             ,
                                           display_degree   ,
                                           one_eval_success ,
                                           success          ,
                                           init_nb_eval     ,
                                           f0               ,
                                           barrier          ,
                                           nb_success       ,
                                           one_for_luck       ) )
                    {
                        _stats.add_interrupted_eval();
                        opp_stop = true; // will break loop #1
                    }
                    
                    // close the evaluation block:
                    if ( display_degree == NOMAD::FULL_DISPLAY )
                        out.close_block();
                }
                
                // point not in cache (the point is saved for loop #2):
                else
                {
                    
                    // blackbox or surrogate evaluations are allowed:
                    if ( ( x->get_eval_type() == NOMAD::TRUTH && max_bb_eval != 0 ) ||
                        ( x->get_eval_type() == NOMAD::SGTE  && _p.get_max_sgte_eval() != 0 ) )
                        to_be_evaluated.push_back ( x );
                    
                    // close the evaluation block:
                    if ( display_degree == NOMAD::FULL_DISPLAY )
                        out.close_block();
                }
            }
            
            // points[k]->check() failed (close the evaluation block):
            else if ( display_degree == NOMAD::FULL_DISPLAY )
            {
                std::ostringstream oss;
                oss << "check failed (" << check_failed_reason << ")";
                out.close_block ( oss.str() );
            }
        }
        
        // evaluation already in progress (close the evaluation block):
        else if ( display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss;
            oss << "evaluation of point #" << x->get_tag()
            << " already in progress";
            out.close_block ( oss.str() );
        }
        
        ++it;
        ++k;
        
        // force quit (by pressing ctrl-c):
        if ( !stop && (NOMAD::Evaluator_Control::_force_quit || NOMAD::Evaluator::get_force_quit()) )
        {
            stop        = true;
            stop_reason = NOMAD::CTRL_C;
        }
        
    }  // end of loop #1
    // --------------
    
    // loop #2: evaluations:
    // ---------------------
    int                 nb_to_evaluate = static_cast<int> ( to_be_evaluated.size() );
    int                 nb_evaluated   = 0;
    int                 cur            = 0;
    int                 source;
    char                signal;
    NOMAD::Eval_Point * eval_x;
    
    while ( !stop && !opp_stop && nb_evaluated < nb_to_evaluate )
    {
        
        source = NOMAD::Slave::receive_signal ( signal );
        
        // 2.1: send the RESULT signal, receive and process the evaluation result:
        // -----------------------------------------------------------------------
        eval_x = _eval_in_progress[source];
        if ( eval_x )
        {
            
            if ( display_degree == NOMAD::FULL_DISPLAY )
                out << std::endl << "receive eval point #" << eval_x->get_tag()
                << " from slave " << source << std::endl << std::endl;
            
            receive_eval_result ( search       ,
                                 eval_x       ,
                                 true_barrier ,
                                 sgte_barrier ,
                                 pareto_front ,
                                 source       ,
                                 stop         ,
                                 stop_reason    );
            
            
            // list of processed points:
            if ( eval_x->is_in_cache() )
                evaluated_pts.push_back ( eval_x );
            
            // success:
            one_eval_success = barrier.get_one_eval_succ();
            success          = barrier.get_success();
            
            // asynchronous success count:
            if ( success == NOMAD::FULL_SUCCESS &&
                _elop_tag != _slaves_elop_tags[source] )
                _stats.add_asynchronous_success();
            
            // displays:
            display_eval_result ( *eval_x          ,
                                 display_degree   ,
                                 search           ,
                                 one_eval_success ,
                                 success            );
            
            // stop the evaluations (opportunistic strategy) ?
            if ( stop_evaluations ( *eval_x          ,
                                   search           ,
                                   nb_evaluated     ,
                                   nb_to_evaluate   ,
                                   stop             ,
                                   display_degree   ,
                                   one_eval_success ,
                                   success          ,
                                   init_nb_eval     ,
                                   f0               ,
                                   barrier          ,
                                   nb_success       ,
                                   one_for_luck       ) )
            {
                _stats.add_interrupted_eval();
                opp_stop = true; // will break loop #2
            }
            
            _eval_in_progress[source] = NULL;
            _slaves_elop_tags[source] = -1;
            --_nb_in_progress;
            ++nb_evaluated;
        }
        
        // 2.2: send the EVAL signal and launch a new evaluation:
        // ------------------------------------------------------
        else
        {
            
            // do not launch a new evaluation if...
            
            // there is no more points to be evaluated:
            if ( cur == nb_to_evaluate )
                NOMAD::Slave::send_signal ( NOMAD::WAIT_SIGNAL , source );
            
            // or if bbe+_nb_in_progress >= max_bb_eval:
            else if ( to_be_evaluated[cur]->get_eval_type() == NOMAD::TRUTH &&
                     max_bb_eval > 0 &&
                     _stats.get_bb_eval() + _nb_in_progress >= max_bb_eval    )
            {
                stop        = true;
                stop_reason = NOMAD::MAX_BB_EVAL_REACHED;
                NOMAD::Slave::send_signal ( NOMAD::WAIT_SIGNAL , source );
            }
            
            else
            {
                
                // get the signature:
                NOMAD::Signature * signature = to_be_evaluated[cur]->get_signature();
                
                // there is no signature (error):
                if ( !signature )
                {
                    stop        = true;
                    stop_reason = NOMAD::ERROR;
                    if ( display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY)
                        out << std::endl
                        << "Error in Evaluator_Control::private_eval_list_of_points():"
                        << " the point #" << to_be_evaluated[cur]->get_tag()
                        << " has no signature" << std::endl << std::endl;
                    NOMAD::Slave::send_signal ( NOMAD::WAIT_SIGNAL , source );
                }
                
                else
                {
                    
                    NOMAD::Slave::send_signal ( NOMAD::EVAL_SIGNAL , source );
                    
                    eval_x = &NOMAD::Cache::get_modifiable_point ( *to_be_evaluated[cur++] );
                    
                    if ( display_degree == NOMAD::FULL_DISPLAY )
                        out << std::endl
                        << "send eval point #" << eval_x->get_tag()
                        << " to slave " << source << std::endl;
                    
                    // 1. scaling:
                    bool do_scaling = signature->get_scaling().is_defined();
                    if ( do_scaling )
                        eval_x->scale();
                    
                    // 2. send the point:
                    _slave->send_eval_point ( eval_x , source , barrier.get_h_max() );
                    
                    // 3. unscaling:
                    if ( do_scaling )
                        eval_x->unscale();
                    
                    eval_x->set_eval_status ( NOMAD::EVAL_IN_PROGRESS );
                    
                    _eval_in_progress[source] = eval_x;
                    _slaves_elop_tags[source] = _elop_tag;
                    ++_nb_in_progress;
                }
            }
        }
        
        // force quit (by pressing ctrl-c):
        if ( !stop && ( NOMAD::Evaluator_Control::_force_quit || NOMAD::Evaluator::get_force_quit() ))
        {
            stop        = true;
            stop_reason = NOMAD::CTRL_C;
        }
        
    }  // end of loop #2
    // --------------
    
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out << std::endl
        << "number of evaluations in progress: " << _nb_in_progress
        << std::endl << std::endl;
    
    // the algorithm is not asynchronous: we have
    // to wait for all the evaluations in progress:
    if ( !_p.get_asynchronous() )
        
        wait_for_evaluations ( search        ,
                              true_barrier  ,
                              sgte_barrier  ,
                              pareto_front  ,
                              stop          ,
                              stop_reason   ,
                              success       ,
                              evaluated_pts   );
    
    // barriers update:
    if ( !stop )
    {
        true_barrier.update_and_reset_success();
        sgte_barrier.update_and_reset_success();
    }
    
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out << NOMAD::close_block ( "end of evaluations" ) << std::endl;
    
    // incumbents update:
    const NOMAD::Eval_Point * bf = barrier.get_best_feasible  ();
    const NOMAD::Eval_Point * bi = barrier.get_best_infeasible();
    if ( bf && bf != old_feasible_incumbent )
        new_feas_inc = bf;
    if ( bi && bi != old_infeasible_incumbent )
        new_infeas_inc = bi;
    
    // the list of eval. points is deleted (only points in the cache are kept):
    clear_eval_lop();
    
    // update the unique eval_lop() tag:
    ++_elop_tag;
    
} // end of eval_lop() parallel version

// C. Tribes may 28, 2014 --- method for points block evaluation of a given max size
/*----------------------------------------------------------------*/
/*       eval_list_of_points, private version (scalar version)    */
/*----------------------------------------------------------------*/
#else
void NOMAD::Evaluator_Control::private_eval_list_of_points
( NOMAD::search_type              search         ,   // IN     : search type
 NOMAD::Barrier                & true_barrier   ,   // IN/OUT : the barrier
 NOMAD::Barrier                & sgte_barrier   ,   // IN/OUT : the surrogate barrier
 NOMAD::Pareto_Front           * pareto_front   ,   // IN/OUT : the Pareto front
 //          (can be NULL)
 bool                          & stop           ,   // IN/OUT : stopping criterion
 NOMAD::stop_type              & stop_reason    ,   // OUT    : stopping reason
 const NOMAD::Eval_Point      *& new_feas_inc   ,   // OUT    : new feasible incumbent
 const NOMAD::Eval_Point      *& new_infeas_inc ,   // OUT    : new infeas. incumbent
 NOMAD::success_type           & success        ,   // OUT    : type of success
 std::list<const NOMAD::Eval_Point *>
 & evaluated_pts    ) // OUT    : list of processed pts
{
    if ( stop || _eval_lop.empty() )
    {
        stop_reason = NOMAD::UNKNOWN_STOP_REASON;
        return;
    }
    
    evaluated_pts.clear();
    
    // initial display:
    const NOMAD::Display    & out = _p.out();
    NOMAD::dd_type display_degree = out.get_display_degree ( search );
    
    
    if ( display_degree == NOMAD::FULL_DISPLAY )
    {
        std::ostringstream oss;
        oss << "list of points evaluation (" << search << ")";
        out << std::endl << NOMAD::open_block ( oss.str() );
    }
    
    // call the Evaluator (virtual) preprocessing of a list of points:
    _ev->list_of_points_preprocessing ( _eval_lop );
    
    const NOMAD::Eval_Point * old_feasible_incumbent   = NULL;
    const NOMAD::Eval_Point * old_infeasible_incumbent = NULL;
    
    // active barrier:
    NOMAD::Barrier & barrier = ( _p.get_opt_only_sgte() ) ?  sgte_barrier : true_barrier;
    
    old_feasible_incumbent   = barrier.get_best_feasible();
    old_infeasible_incumbent = barrier.get_best_infeasible();
    
    NOMAD::Double f0;
    if ( _p.get_opportunistic_min_f_imprvmt().is_defined() &&
        old_feasible_incumbent )
        f0 = old_feasible_incumbent->get_f();
    
    new_feas_inc   = NULL;
    new_infeas_inc = NULL;
    stop           = false;
    success        = NOMAD::UNSUCCESSFUL;
    stop_reason    = NOMAD::NO_STOP;
    
    const NOMAD::Eval_Point * x;
    NOMAD::check_failed_type  check_failed_reason;
    bool                      one_for_luck = false;
    bool                      stop_evals   = false;
    int                       init_nb_eval = _stats.get_eval();
    int                       nb_success   = 0;
    int                       k            = 0;
    int                       k_block      = 0;
    int                       nb_points    = get_nb_eval_points();
    int						  block_size   = _p.get_bb_max_block_size();
    int						  block_nb		= 1;
    
    // main loop (on the list of points):
    // ----------------------------------
    std::set<NOMAD::Priority_Eval_Point>::iterator it  = _eval_lop.begin() , end = _eval_lop.end();
    std::list<NOMAD::Eval_Point *> list_x,list_eval;
    std::list<bool> count_list_eval;
    
    while ( !stop_evals && !stop && it != end )
    {
        
        
        if ( block_size > 1 && display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss;
            oss << "Block of evaluations (" << block_nb <<")";
            out << std::endl << NOMAD::open_block ( oss.str() );
        }
        
        // Creation of a block of evaluations from the list
        //----------------------
        k_block=k;
        bool opportunistic_success_from_cache_point=false;
        while (list_eval.size()!=static_cast<size_t>(block_size) && it != end && ! stop_evals)
        {
            
            x = it->get_point();
            x->set_current_run ( true );
            
            // displays:
            if ( display_degree == NOMAD::FULL_DISPLAY )
            {
                {
                    // open the evaluation block:
                    std::ostringstream oss;
                    oss << "submitted ";
                    if ( x->get_eval_type() == NOMAD::SGTE )
                        oss << "surrogate ";
                    oss << "evaluation " << k+1 << "/" << nb_points;
                    out << std::endl << NOMAD::open_block ( oss.str() );
                }
                
                out << std::endl << "point #" << x->get_tag() << "   ( ";
                x->Point::display ( out , " " , 2 , NOMAD::Point::get_display_limit() );
                out << " )" << std::endl;
                if ( x->get_direction() )
                {
                    out << "direction   : " << *x->get_direction()  << std::endl;
                    NOMAD::Point delta;
                    x->get_signature()->get_mesh()->get_delta(delta);
                    out << "direction d : ( " << *x->get_direction()/delta  << " )"  << std::endl;
                }
                if ( x->get_signature() )
                    out << "mesh indices: ( " << x->get_signature()->get_mesh()->get_mesh_indices() << " )" << std::endl;
                out << std::endl;
                
            }
            
            // current point check (# of bb outputs, bounds, integer values, fixed-vars):
            if ( x->check ( _p.get_bb_nb_outputs() , check_failed_reason ) )
            {
                bool count_eval = true;
                bool has_been_in_cache=cache_check ( x                  ,
                                                    true_barrier        ,
                                                    sgte_barrier        ,
                                                    pareto_front        ,
                                                    count_eval			,
                                                    barrier.get_h_max() ,
                                                    display_degree        );
                
                
                // put the point in a block list for evaluation:
                if ( !has_been_in_cache )
                    list_eval.push_back(&NOMAD::Cache::get_modifiable_point ( *x ));
                else
                {
                    // check stopping criteria for points in cache
                    check_stopping_criteria ( search , count_eval , *x , stop , stop_reason );
                    
                    // process the evaluated point:
                    process_eval_point ( *x                                       ,
                                        ( x->get_eval_type() == NOMAD::TRUTH ) ? true_barrier : sgte_barrier ,
                                        pareto_front                                    );
                    
                    
                    // success:
                    NOMAD::success_type one_eval_success = barrier.get_one_eval_succ();
                    success                              = barrier.get_success();
                    
                    
                    opportunistic_success_from_cache_point = stop_evaluations ( *x               ,
                                                                               search           ,
                                                                               k                ,
                                                                               nb_points        ,
                                                                               stop             ,
                                                                               display_degree   ,
                                                                               one_eval_success ,
                                                                               success          ,
                                                                               init_nb_eval     ,
                                                                               f0               ,
                                                                               barrier          ,
                                                                               nb_success       ,
                                                                               one_for_luck       );
                    
                    
                }
                
                if (!stop)
                    list_x.push_back(&NOMAD::Cache::get_modifiable_point ( *x ));
                
                
                if ( opportunistic_success_from_cache_point )
                {
                    if ( display_degree == NOMAD::FULL_DISPLAY )
                        out << NOMAD::close_block();
                    
                    if ( block_size > 1 && display_degree == NOMAD::FULL_DISPLAY )
                        out << NOMAD::close_block ();
                    
                    stop_evals = true;
                    break;
                }
                
                
            }
            // points[k]->check() failed:
            else if ( display_degree == NOMAD::FULL_DISPLAY )
                out << "check failed (" << check_failed_reason << ")" << std::endl;
            
            if ( display_degree == NOMAD::FULL_DISPLAY )
                out << NOMAD::close_block();
            
            
            ++it;
            ++k;
        }
        if (list_eval.size()!=0)
        {
            
            count_list_eval.assign(list_eval.size(), false);
            
            if (_p.eval_points_as_block())
            {
                eval_points ( list_eval			,
                             true_barrier		,
                             sgte_barrier		,
                             pareto_front		,
                             count_list_eval	,
                             stop				,
                             stop_reason		,
                             barrier.get_h_max() );
                
                // check stopping criteria for points NOT in cache
                std::list<NOMAD::Eval_Point *>::iterator it_eval;
            }
            else
            {
                // bool count_eval=false;
                x=*(list_eval.begin());
                eval_point ( NOMAD::Cache::get_modifiable_point ( *x )	,
                            true_barrier								,
                            sgte_barrier								,
                            pareto_front								,
                            count_list_eval.front()					,
                            stop										,
                            stop_reason								,
                            barrier.get_h_max()						);
                
            }
            
        }
        
        // Stop evals and exit the loop
        if ( stop_evals )
            break;
        
        // Check all the points in the evaluation block
        std::list<NOMAD::Eval_Point *>::iterator it_x,it_eval;
        k=k_block;
        it_eval=list_eval.begin();
        for(it_x=list_x.begin();it_x!=list_x.end();++it_x)
        {
            
            x=(*it_x);
            
            // process the evaluated point:
            if ( x->is_eval_ok() && x->is_in_cache() )
                process_eval_point ( *x                                       ,
                                    ( x->get_eval_type() == NOMAD::TRUTH ) ?
                                    true_barrier : sgte_barrier                   ,
                                    pareto_front                                    );
            
            
            
            
            // success:
            NOMAD::success_type one_eval_success = barrier.get_one_eval_succ();
            success                              = barrier.get_success();
            
            // list of processed points:
            if ( x->is_in_cache() )
                evaluated_pts.push_back ( x );
            else
            {
                // this situation may occur on very thin meshes:
                // the point has not been found in the cache
                // and it failed to be inserted.
                one_eval_success = NOMAD::UNSUCCESSFUL;
            }
            
            // displays:
            if ( block_size > 0 && display_degree == NOMAD::FULL_DISPLAY )
            {
                // open the evaluation block:
                std::ostringstream oss;
                if ( x->get_eval_type() == NOMAD::SGTE )
                    oss << "surrogate ";
                oss << "evaluation " << k+1 << "/" << nb_points;
                out << std::endl << NOMAD::open_block ( oss.str() );
                out << std::endl << "point #" << x->get_tag() << std::endl;
            }
            
            
            std::list<bool>::iterator it_count=count_list_eval.begin();
            for(it_eval=list_eval.begin();it_eval!=list_eval.end();++it_eval,++it_count)
            {
                if ((*it_eval)==x)
                {
                    
                    // count the bb evaluation:
                    if ( *it_count )
                    {
                        if ( (*it_eval)->get_eval_type() == NOMAD::SGTE )
                            _stats.add_sgte_eval();
                        else
                        {
                            // current mads bbe evaluation
                            _stats.add_bb_eval();
                        }
                        
                        // count the output stats (STAT_SUM and STAT_AVG):
                        if ( _p.check_stat_sum() || _p.check_stat_avg() )
                            count_output_stats(*(*it_eval));
                    }
                    
                    check_stopping_criteria ( search , *it_count ,*(*it_eval) , stop , stop_reason );
                    
                    if ( *it_count )
                        display_eval_result ( *x, display_degree, search, one_eval_success, success );
                    
                    break;
                }
            }
            
            
            
            // close the evaluation block:
            if ( display_degree == NOMAD::FULL_DISPLAY )
                out << NOMAD::close_block ();
            
            // stop the evaluations (opportunistic strategy) ?
            if ( !stop_evals && stop_evaluations ( *x               ,
                                                  search           ,
                                                  k                ,
                                                  nb_points        ,
                                                  stop             ,
                                                  display_degree   ,
                                                  one_eval_success ,
                                                  success          ,
                                                  init_nb_eval     ,
                                                  f0               ,
                                                  barrier          ,
                                                  nb_success       ,
                                                  one_for_luck       ) )
            {
                _stats.add_interrupted_eval();
                stop_evals = true;
            }
            
            
            ++k;
            
        }
        
        if ( block_size > 1 && display_degree == NOMAD::FULL_DISPLAY )
            out << NOMAD::close_block ();
        
        
        // force quit (by pressing ctrl-c):
        if ( !stop && ( NOMAD::Evaluator_Control::_force_quit || NOMAD::Evaluator::get_force_quit()) )
        {
            stop        = true;
            stop_reason = NOMAD::CTRL_C;
        }
        
        list_x.clear();
        list_eval.clear();
        
        ++block_nb;
        
    }// end of test for list evaluation
    
    // barriers update:
    if ( !stop )
    {
        true_barrier.update_and_reset_success();
        sgte_barrier.update_and_reset_success();
    }
    
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out << std::endl << NOMAD::close_block ( "end of evaluations" )
        << std::endl;
    
    // incumbents update:
    const NOMAD::Eval_Point * bf = barrier.get_best_feasible  ();
    const NOMAD::Eval_Point * bi = barrier.get_best_infeasible();
    if ( bf && bf != old_feasible_incumbent )
        new_feas_inc = bf;
    if ( bi && bi != old_infeasible_incumbent )
        new_infeas_inc = bi;
    
    // the list of eval. points is deleted (only points in the cache are kept):
    clear_eval_lop();
    
} // end of eval_lop() scalar version


#endif

/*-------------------------------------------------*/
/*       reduce the list of evaluation points      */
/*-------------------------------------------------*/
void NOMAD::Evaluator_Control::reduce_eval_lop ( int n )
{
    int nb_eval_pts = get_nb_eval_points();
    
    if ( n < 0 || n >= nb_eval_pts )
        return;
    
    const NOMAD::Eval_Point * x;
    std::set<NOMAD::Priority_Eval_Point>::iterator it = _eval_lop.end();
    for( int i=0;i<nb_eval_pts-n;i++)
    {
        --it;
        x = it->get_point();
        if ( x && !x->is_in_cache() && x->get_eval_status() != NOMAD::EVAL_IN_PROGRESS )
            delete x;
    }
    _eval_lop.erase( it,_eval_lop.end());
}

/*-------------------------------------------------*/
/*            TGP model ordering (private)         */
/*-------------------------------------------------*/
void NOMAD::Evaluator_Control::TGP_model_ordering ( NOMAD::dd_type   display_degree ,
                                                   bool           & modified_list    )
{
    modified_list = false;
    
    if ( _p.get_opt_only_sgte() )
        return;
    
#ifdef USE_TGP
    
    // display:
    const NOMAD::Display & out = _p.out();
    
    // model stats:
    NOMAD::Model_Stats model_stats;
    NOMAD::Clock       clock;
    
#ifdef TGP_DEBUG
    out << std::endl << NOMAD::open_block ( "TGP model ordering") << std::endl;
#endif
    
    const std::vector<NOMAD::bb_output_type> & bbot = _p.get_bb_output_type();
    int i , j , n_XX = 0 , m = bbot.size();
    
    // construct prediction set (XX):
    // ------------------------------
    std::vector<NOMAD::Eval_Point *> XX;
    NOMAD::Point                     lb_XX , ub_XX;
    
    // save _eval_lop in XX and other_pts:
    const NOMAD::Eval_Point            * x;
    std::list<const NOMAD::Eval_Point *> other_pts;
    const NOMAD::Signature             * signature = NULL;
    int                                  n         = -1;
    
    std::set<NOMAD::Priority_Eval_Point>::const_iterator it , end = _eval_lop.end();
    for ( it = _eval_lop.begin() ; it != end ; ++it ) {
        x = it->get_point();
        if ( n < 0 ) {
            signature = x->get_signature();
            if ( !signature ) {
#ifdef TGP_DEBUG
                out << NOMAD::close_block ( "failure (no signature)" ) << std::endl;
#endif
                return;
            }
            n = signature->get_n();
            
            lb_XX = ub_XX = NOMAD::Point(n);
        }
        
        if ( x->size           () == n            &&
            x->get_m          () == m            &&
            x->get_eval_type  () == NOMAD::TRUTH &&
            !x->get_bb_outputs().is_defined()       ) {
            
            XX.push_back ( &NOMAD::Cache::get_modifiable_point ( *x ) );
            
            for ( i = 0 ; i < n ; ++i ) {
                if ( !lb_XX[i].is_defined() || (*x)[i] < lb_XX[i] )
                    lb_XX[i] = (*x)[i];
                if ( !ub_XX[i].is_defined() || (*x)[i] > ub_XX[i] )
                    ub_XX[i] = (*x)[i];
            }
        }
        else
            other_pts.push_back ( x );
    }
    
    n_XX = XX.size();
    
    if ( n_XX <= 1 ) {
#ifdef TGP_DEBUG
        out << NOMAD::close_block ( "failure (size(XX) <= 1)" ) << std::endl;
#endif
        return;
    }
    
    // the TGP model:
    NOMAD::TGP_Model * model;
    
    // Reuse the last TGP model from the TGP model search:
    if ( _last_TGP_model && _p.get_model_tgp_reuse_model() ) {
        
        model = _last_TGP_model;
        
        // individual predictions for XX points:
        for ( i = 0 ; i < n_XX ; ++i )
            if ( !model->predict ( *XX[i] , false ) ) // pred_outside_bnds = false
                for ( j = 0 ; j < m ; ++j )
                    XX[i]->set_bb_output ( j , NOMAD::Double() );
    }
    
    // creation of a new TGP model:
    else {
        
        model = new NOMAD::TGP_Model ( n , bbot , out , _p.get_model_tgp_mode() );
        
        NOMAD::Point center(n);
        for ( i = 0 ; i < n ; ++i )
            center[i] = ( lb_XX[i] + ub_XX[i] ) / 2.0;
        
        // construct interpolation set (X):
        // --------------------------------
        if ( !model->set_X ( *_cache       ,
                            &center       ,
                            _p.get_seed() ,
                            true            ) ) { // remove_fv = true
            
            if ( model->get_p() <= model->get_n() )
                model_stats.add_not_enough_pts();
#ifdef TGP_DEBUG
            out << NOMAD::close_block ( "failure: " + model->get_error_str() )
            << std::endl;
#endif
            
            delete model;
            
            return;
        }
        
        int p = model->get_p();
        
        // display sets X and XX:
        // ----------------------
#ifdef TGP_DEBUG
        {
            // max number of points displayed:
            const int set_display_limit = 15; // set to -1 for no limit
            
            // X:
            model->display_X ( out ,  set_display_limit );
            
            // XX:
            out << NOMAD::open_block ( "prediction points (XX)");
            for ( i = 0 ; i < n_XX ; ++i ) {
                out << "#";
                out.display_int_w ( i , n_XX );
                out << " x=(";
                XX[i]->NOMAD::Point::display ( out , " " , 15 , -1 );
                out << " )" << std::endl;
            }
            std::ostringstream oss;
            oss << "(size=" << n_XX << ")";
            out << NOMAD::close_block ( oss.str() ) << std::endl;
        }
#endif
        
        // TGP model construction:
        // -----------------------
#ifdef TGP_DEBUG
        out << "TGP model construction ...";
        out.flush();
#endif
        
        if ( !model->compute ( XX    ,
                              false ,       // compute_Ds2x      = false
                              false ,       // compute_improv    = false
                              false   ) ) { // pred_outside_bnds = false
            
            model_stats.add_construction_error();
            
#ifdef TGP_DEBUG
            out << "... error: " << model->get_error_str() << std::endl
            << NOMAD::close_block() << std::endl;
#endif
            
            // reset XX outputs:
            for ( i = 0 ; i < n_XX ; ++i )
                for ( j = 0 ; j < m ; ++j )
                    XX[i]->set_bb_output ( j , NOMAD::Double() );
            
            delete model;
            
            // check if ctrl-c has been pressed:
            if ( NOMAD::TGP_Output_Model::get_force_quit() )
                NOMAD::Evaluator_Control::_force_quit = true;
            
            return;
        }
#ifdef TGP_DEBUG
        out << "... OK" << std::endl << std::endl;
#endif
        
        // update model stats:
        model_stats.add_construction_time ( clock.get_CPU_time() );
        model_stats.update_nY             ( p                    );
        model_stats.update_ES_stats       ( n_XX , n_XX          );
        model_stats.add_nb_truth();
        model_stats.add_nb_TGP();
    }
    
    // open display block for model predictions:
#ifdef TGP_DEBUG
    out << NOMAD::open_block ( "TGP predictions (XX+ZZ)");
#endif
    
    // clear then fill _eval_lop again:
    // --------------------------------
    NOMAD::Double         f_model , h_model;
    const NOMAD::Double & h_min          = _p.get_h_min();
    NOMAD::hnorm_type     h_norm         = _p.get_h_norm();
    bool                  snap_to_bounds = _p.get_snap_to_bounds();
    
    modified_list = true;
    _eval_lop.clear();
    
    for ( i = 0 ; i < n_XX ; ++i ) {
        
        // compute model h and f values:
        model->eval_hf ( XX[i]->get_bb_outputs() ,
                        h_min                   ,
                        h_norm                  ,
                        h_model                 ,
                        f_model                   );
        
        // display model predictions:
#ifdef TGP_DEBUG
        out << "#";
        out.display_int_w ( i , n_XX );
        out << " x=(";
        XX[i]->NOMAD::Point::display ( out , " " , 15 , -1 );
        out << " ) m(x)=[";
        XX[i]->get_bb_outputs().display ( out , " " , 15 , -1 );
        out << " ]";
        
        if ( h_model.is_defined() && f_model.is_defined() )
            out << " hm=" << std::setw(15) << h_model
            << " fm=" << std::setw(15) << f_model;
        else
            out << " no model value";
        out << std::endl;
#endif
        
        // add the evaluation point:
        add_eval_point ( XX[i]           ,
                        display_degree  ,
                        snap_to_bounds  ,
                        NOMAD::Double() ,
                        NOMAD::Double() ,
                        f_model         ,
                        h_model           );
        
#ifdef MODEL_STATS
        if ( XX[i] && f_model.is_defined() && h_model.is_defined() ) {
            XX[i]->set_mod_use  ( 2                ); // 2 for model ordering
            XX[i]->set_Yw       ( model->get_Yw () );
            XX[i]->set_nY       ( p                );
            XX[i]->set_mh       ( h_model          );
            XX[i]->set_mf       ( f_model          );
        }
#endif
    }
    
#ifdef TGP_DEBUG
    {
        // close display block for model predictions:
        std::ostringstream oss;
        oss << "(size=" << n_XX << ")";
        out << NOMAD::close_block ( oss.str() ) << std::endl;
        
        // compute and display prediction errors:
        out << NOMAD::open_block ( "prediction relative errors on X(%)" );
        model->display_X_errors ( out );
        out << NOMAD::close_block() << std::endl;
    }
#endif
    
    // other points that have been previously discarded and have no model values:
    NOMAD::Eval_Point * y;
    std::list<const NOMAD::Eval_Point *>::const_iterator it2 , end2 = other_pts.end();
    for ( it2 = other_pts.begin() ; it2 != end2 ; ++it2 ) {
        y = &NOMAD::Cache::get_modifiable_point (**it2);
        add_eval_point ( y               ,
                        display_degree  ,
                        snap_to_bounds  ,
                        NOMAD::Double() ,
                        NOMAD::Double() ,
                        NOMAD::Double() ,
                        NOMAD::Double()   );
    }
    
    _stats.update_model_stats    ( model_stats );
    _model_ordering_stats.update ( model_stats );
    
    if ( model != _last_TGP_model )
        delete model;
    
#ifdef TGP_DEBUG
    out << NOMAD::close_block() << std::endl;
#else
    if ( display_degree == NOMAD::FULL_DISPLAY ) {
        out << std::endl << "model ordering";
        if ( !modified_list )
            out << " (no modification)";
        out << std::endl;
    }
#endif
#endif
}

/*-------------------------------------------------*/
/*         model_np1_quad_epsilon (private)      */
/*-------------------------------------------------*/
void NOMAD::Evaluator_Control::quad_model_ordering ( NOMAD::dd_type display_degree ,
                                                    bool         & modified_list    )
{
    const NOMAD::Display & out = _p.out();
    
#ifdef DEBUG
    out << std::endl << NOMAD::open_block ( "model_np1_quad_epsilon") << std::endl;
#endif
    
    // save _eval_lop in pts and other_pts:
    // ------------------------------------
    NOMAD::Point                         min , max , center , interpolation_radius;
    const NOMAD::Eval_Point *            y;
    std::list<const NOMAD::Eval_Point *> pts , other_pts;
    const NOMAD::Signature  *            signature     = NULL;
    const NOMAD::Double &                radius_factor = _p.get_model_quad_radius_factor();
    NOMAD::eval_type                     ev_type       = NOMAD::TRUTH;
    int                                  i , n = -1;
    
    std::set<NOMAD::Priority_Eval_Point>::const_iterator it , end = _eval_lop.end();
    for ( it = _eval_lop.begin() ; it != end ; ++it )
    {
        y = it->get_point();
        if ( n < 0 )
        {
            signature = y->get_signature();
            if ( !signature )
            {
#ifdef DEBUG
                out << NOMAD::close_block ( "failure (no signature)" ) << std::endl;
#endif
                modified_list = false;
                return;
            }
            n       = signature->get_n();
            ev_type = y->get_eval_type();
            min.resize                  ( n );
            max.resize                  ( n );
            center.resize               ( n );
            interpolation_radius.resize ( n );
        }
        
        if ( y->size() == n && y->get_eval_type() == ev_type ) {
            pts.push_back(y);
            for ( i = 0 ; i < n ; ++i ) {
                if ( !min[i].is_defined() || (*y)[i] < min[i] )
                    min[i] = (*y)[i];
                if ( !max[i].is_defined() || (*y)[i] > max[i] )
                    max[i] = (*y)[i];
            }
        }
        else
            other_pts.push_back ( y );
    }
    
    for ( i = 0 ; i < n ; ++i )
    {
        center              [i] = ( min[i] + max[i] ) / 2.0;
        interpolation_radius[i] = ( max[i] - min[i] ) * radius_factor / 2.0;
    }
    
#ifdef DEBUG
    out << NOMAD::open_block ( "points used to define interpolation radius")
    << "type of eval.   : " << ev_type    << std::endl
    << "number of points: " << pts.size() << std::endl
    << "min. coordinates: ( ";
    min.display ( out , " " , 2 );
    out << " )" << std::endl
    << "max. coordinates: ( ";
    max.display ( out , " " , 2 );
    out << " )" << std::endl
    << "center          : ( ";
    center.display ( out , " " , 2 );
    out << " )" << std::endl
    << "interp. radius  : ( ";
    interpolation_radius.display ( out , " " , 2 );
    out << " )" << std::endl
    << NOMAD::close_block() << std::endl;
#endif
    
    // create model:
    // -------------
    NOMAD::Clock       clock;
    NOMAD::Model_Stats model_stats;
    NOMAD::Quad_Model  model ( out                                              ,
                              _p.get_bb_output_type()                          ,
                              (ev_type==NOMAD::TRUTH) ? *_cache : *_sgte_cache ,
                              *signature                                         );
    
    int  max_Y_size = _p.get_model_quad_max_Y_size();
    int  min_Y_size = _p.get_model_quad_min_Y_size();
    bool use_WP     = _p.get_model_quad_use_WP    ();
    
    // construct interpolation set Y:
    model.construct_Y ( center               ,
                       interpolation_radius ,
                       max_Y_size             );
    
    int nY = model.get_nY();
    
#ifdef DEBUG
    out << "number of points in Y: " << nY
    << " (p=" << nY-1;
    if ( nY < 2 ) out << ", not enough";
    out << ")" << std::endl;
#endif
    
    // not enough points:
    if ( nY < 2 )
    {
        modified_list = false;
        model_stats.add_not_enough_pts();
    }
    else
    {
        
#ifdef DEBUG
        out << std::endl;
        model.display_Y ( out , "unscaled interpolation set Y" );
#endif
        
        // define scaling:
        model.define_scaling ( radius_factor );
        
#ifdef DEBUG
        out << std::endl;
        model.display_Y ( out , "scaled interpolation set Ys" );
#endif
        
        // model error flag:
        if ( model.get_error_flag() )
        {
            model_stats.add_construction_error();
            modified_list = false;
        }
        else
        {
            // not enough points:
            if ( nY < 2 || ( min_Y_size < 0 && nY <= model.get_nfree() ) )
            {
                model_stats.add_not_enough_pts();
                modified_list = false;
            }
            // enough points and no error:
            else
            {
                
                bool cautious         = _p.get_model_eval_sort_cautious();
                int  nb_inside_radius = 0;
                
                // check that there is at least two trial points inside the trust radius
                // (cautious strategy):
                nb_inside_radius = 0;
                std::list<const NOMAD::Eval_Point *>::const_iterator it2 , end2 = pts.end();
                if ( cautious )
                {
                    for ( it2 = pts.begin() ; it2 != end2 ; ++it2 )
                    {
                        NOMAD::Point scaled_pt ( **it2 );
                        model.scale ( scaled_pt );
                        if ( model.is_within_trust_radius ( scaled_pt ) )
                        {
                            if ( ++nb_inside_radius == 2 )
                                break;
                        }
                    }
                }
                
                // not enough points inside trust radius:
                if ( cautious && nb_inside_radius < 2 )
                    modified_list = false;
                
                // at least two trial points are inside trust radius:
                else
                {
                    // construct model:
                    // ----------------
                    model.construct ( use_WP , NOMAD::SVD_EPS , NOMAD::SVD_MAX_MPN , max_Y_size );
                    model_stats.add_construction_time ( clock.get_CPU_time() );
                    model_stats.update_nY ( model.get_nY() );
                    
                    // display model characteristics:
#ifdef DEBUG
                    out << std::endl;
                    model.display_model_coeffs ( out );
                    out << std::endl;
                    model.display_Y_error ( out );
#endif
                    
                    // count model:
                    if ( ev_type == NOMAD::TRUTH )
                        model_stats.add_nb_truth();
                    else
                        model_stats.add_nb_sgte();
                    
                    switch ( model.get_interpolation_type() )
                    {
                        case NOMAD::MFN:
                            model_stats.add_nb_MFN();
                            break;
                        case NOMAD::WP_REGRESSION:
                            model_stats.add_nb_WP_regression();
                            break;
                        case NOMAD::REGRESSION:
                            model_stats.add_nb_regression();
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
                        modified_list = false;
                        if ( model.get_error_flag() )
                            model_stats.add_construction_error();
                        else
                            model_stats.add_bad_cond();
                    }
                    else
                    {
                        // clear then fill _eval_lop again:
                        // --------------------------------
                        NOMAD::Double         f_model , h_model;
                        NOMAD::Eval_Point   * x;
                        bool                  snap_to_bounds = _p.get_snap_to_bounds();
                        
                        modified_list = true;
                        _eval_lop.clear();
                        
                        nb_inside_radius = 0;
                        
#ifdef DEBUG
                        out << std::endl << NOMAD::open_block ( "original trial points" );
#endif
                        
                        NOMAD::Quad_Model_Evaluator *quad_model_ev=new NOMAD::Quad_Model_Evaluator(_p , model);
                        
                        for ( it2 = pts.begin() ; it2 != end2 ; ++it2 )
                        {
                            NOMAD::Point scaled_pt ( **it2 );
                            model.scale ( scaled_pt );
                            
                            
                            
                            f_model.clear();
                            h_model.clear();
                            
                            
                            if ( !cautious || model.is_within_trust_radius ( scaled_pt ) )
                            {
                                
                                int m  = static_cast<int>(_p.get_bb_output_type().size());
                                NOMAD::Eval_Point x_eval(scaled_pt,m);
                                for (int i = 0 ; i < x_eval.size() ; ++i )
                                    x_eval[i] = scaled_pt[i].value() * 1000.0;
                                
                                bool count_eval;
                                
                                bool success=quad_model_ev->eval_x(x_eval,0.0,count_eval);
                                if (success)
                                {
                                    _ev->compute_f(x_eval);
                                    _ev->compute_h(x_eval);
                                    
                                    f_model=x_eval.get_f();
                                    h_model=x_eval.get_h();
                                }
                                
                                ++nb_inside_radius;
                            }
                            
                            x = &NOMAD::Cache::get_modifiable_point (**it2);
                            
#ifdef DEBUG
                            x->display_tag ( out );
                            out << ": ( ";
                            x->NOMAD::Point::display ( out , " " , 2 );
                            out << " ) scaled: (";
                            scaled_pt.NOMAD::Point::display ( out , " " , 2 );
                            out << ") ";
                            if ( h_model.is_defined() && f_model.is_defined() )
                                out << "hm=" << h_model << " fm=" << f_model;
                            else
                                out << "no model value";
                            out << std::endl;
#endif
                            
                            // add the evaluation point:
                            add_eval_point ( x               ,
                                            display_degree  ,
                                            snap_to_bounds  ,
                                            NOMAD::Double() ,
                                            NOMAD::Double() ,
                                            f_model         ,
                                            h_model           );
                            
#ifdef MODEL_STATS
                            if ( x && f_model.is_defined() && h_model.is_defined() )
                            {
                                x->set_mod_use  ( 2                ); // 2 for model ordering
                                x->set_cond     ( model.get_cond() );
                                x->set_Yw       ( model.get_Yw  () );
                                x->set_nY       ( model.get_nY  () );
                                x->set_mh       ( h_model          );
                                x->set_mf       ( f_model          );
                            }
#endif
                        }
                        
                        delete quad_model_ev;
                        
                        // other points that have been previously discarded
                        // and have no model values:
                        end2 = other_pts.end();
                        for ( it2 = other_pts.begin() ; it2 != end2 ; ++it2 )
                        {
                            
                            x = &NOMAD::Cache::get_modifiable_point (**it2);
#ifdef DEBUG
                            x->display_tag ( out );
                            out << ": ( ";
                            x->NOMAD::Point::display ( out , " " , 2 );
                            out << " ) no model value" << std::endl;
#endif	    
                            add_eval_point ( x               ,
                                            display_degree  ,
                                            snap_to_bounds  ,
                                            NOMAD::Double() ,
                                            NOMAD::Double() ,
                                            NOMAD::Double() ,
                                            NOMAD::Double()   );
                        }
#ifdef DEBUG
                        out.close_block();
#endif
                    }
                }
                model_stats.update_ES_stats ( nb_inside_radius , static_cast<int>(pts.size()) );
            }
        }
    }
    
    _stats.update_model_stats    ( model_stats );
    _model_ordering_stats.update ( model_stats );
    
#ifdef DEBUG
    out << NOMAD::close_block() << std::endl;
#else
    if ( display_degree == NOMAD::FULL_DISPLAY ) 
    {
        out << std::endl << "model ordering";
        if ( !modified_list )
            out << " (no modification)";
        out << std::endl;
    }
#endif
}

/*----------------------------------------------------------------------------------*/
/*  evaluation of a list of points (public version that calls the private version)  */
/*----------------------------------------------------------------------------------*/
void NOMAD::Evaluator_Control::eval_list_of_points
( NOMAD::search_type              search             , // IN    : search type
 NOMAD::Barrier                & true_barrier       , // IN/OUT: truth barrier
 NOMAD::Barrier                & sgte_barrier       , // IN/OUT: surrogate barrier
 NOMAD::Pareto_Front           * pareto_front       , // IN/OUT: Pareto front
 //         (can be NULL)
 bool                          & stop               , // IN/OUT: stopping criterion
 NOMAD::stop_type              & stop_reason        , // OUT   : stopping reason
 const NOMAD::Eval_Point      *& new_feas_inc       , // OUT   : new feas. incumbent
 const NOMAD::Eval_Point      *& new_infeas_inc     , // OUT   : new infeas. incumb.
 NOMAD::success_type           & success            , // OUT   : type of success
 std::list<const NOMAD::Eval_Point *>
 * evaluated_pts   )    // OUT   : list of processed
//         pts (can be NULL)
{
    
    bool del_evaluated_pts = false;
    if ( !evaluated_pts )
    {
        evaluated_pts     = new std::list<const NOMAD::Eval_Point *>;
        del_evaluated_pts = true;
    }
    
    bool sgte_eval_sort = _p.get_sgte_eval_sort() && _eval_lop.size() > 1;
    bool opt_only_sgte  = _p.get_opt_only_sgte ();
    bool snap_to_bounds = _p.get_snap_to_bounds();
    bool modified_list  = false;
    
    const NOMAD::Display    & out = _p.out();
    NOMAD::dd_type display_degree = out.get_display_degree ( search );
    
    // reset the success type:
    true_barrier.reset_success();
    sgte_barrier.reset_success();
    
    // define all points as surrogates:
    if ( opt_only_sgte || sgte_eval_sort ) 
    {
        for ( std::set<NOMAD::Priority_Eval_Point>::iterator it = _eval_lop.begin() ;
             it != _eval_lop.end() ; ++it )
            NOMAD::Cache::get_modifiable_point(*it->get_point()).set_eval_type(NOMAD::SGTE);
    }
    
    // use the surrogates to sort the eval. points:
    if ( !opt_only_sgte && sgte_eval_sort ) 
    {
        
        // evaluate the surrogate:
        private_eval_list_of_points ( search         ,
                                     true_barrier   ,
                                     sgte_barrier   ,
                                     NULL           , // Pareto front = NULL
                                     stop           ,
                                     stop_reason    ,
                                     new_feas_inc   ,
                                     new_infeas_inc ,
                                     success        ,
                                     *evaluated_pts   );
        if ( stop )
        {
            if ( del_evaluated_pts )
                delete evaluated_pts;
            return;
        }
        
        NOMAD::Eval_Point * x;
        
        // construct a new list of trial points that will be
        // ordered using surrogate values:   
        std::list<const NOMAD::Eval_Point *>::const_iterator
        end = evaluated_pts->end() , it2;
        for ( it2 = evaluated_pts->begin() ; it2 != end ; ++it2 ) {
            
            // Eval_Point construction:
            x = new NOMAD::Eval_Point;
            x->set ( (*it2)->size() , _p.get_bb_nb_outputs() );
            x->set_signature  ( (*it2)->get_signature () );
            x->set_direction  ( (*it2)->get_direction () );
            x->Point::operator = ( **it2 );
            
            modified_list = true;
            
            // add the new point to the ordered list of trial points:
            add_eval_point ( x               ,
                            display_degree  ,
                            snap_to_bounds  ,
                            (*it2)->get_f() ,
                            (*it2)->get_h() ,
                            NOMAD::Double() ,
                            NOMAD::Double()   );
        }
    }
    
    if ( stop ) {
        if ( del_evaluated_pts )
            delete evaluated_pts;
        return;
    }
    
    // model ordering:
    // ---------------
    if ( !modified_list && _model_eval_sort && _eval_lop.size() > 1 )
    {
        switch ( _p.get_model_eval_sort() ) {
            case NOMAD::TGP_MODEL:
                TGP_model_ordering ( display_degree , modified_list );
                if ( NOMAD::Evaluator_Control::_force_quit || NOMAD::Evaluator::get_force_quit() )
                {
                    stop        = true;
                    stop_reason = NOMAD::CTRL_C;
                }
                break;
            case NOMAD::QUADRATIC_MODEL:
                quad_model_ordering ( display_degree , modified_list );
                break;
            case NOMAD::NO_MODEL:;
        }
    }
    
    // this test is true if ctrl-c has been pressed:
    if ( stop ) {
        if ( del_evaluated_pts )
            delete evaluated_pts;
        return;
    }
    
    // display the re-ordered list of trial points:
    if ( modified_list && display_degree == NOMAD::FULL_DISPLAY ) {
        
        const NOMAD::Eval_Point * y;
        
        std::ostringstream oss;
        oss << "re-ordered list of " << _eval_lop.size()
        << " " << search << " trial points";
        
        out << NOMAD::open_block ( oss.str() ) << std::endl;
        
        std::set<NOMAD::Priority_Eval_Point>::const_iterator
        end = _eval_lop.end() , it;
        for ( it = _eval_lop.begin() ; it != end ; ++it ) {
            y =  it->get_point();
            y->display_tag ( out );
            out << ": ( ";
            y->Point::display ( out , " " , 2 , NOMAD::Point::get_display_limit() );
            out << " )";
            if ( y->get_direction() )
                out << " (dir " << y->get_direction()->get_index() << ")";
            out << std::endl;
        }
        out.close_block();
    }
    
    // evaluate the list of points on the 'true' function:
    private_eval_list_of_points ( search         ,
                                 true_barrier   ,
                                 sgte_barrier   ,
                                 pareto_front   ,
                                 stop           ,
                                 stop_reason    ,
                                 new_feas_inc   ,
                                 new_infeas_inc ,
                                 success        ,
                                 *evaluated_pts   ); 
    
#ifdef MODEL_STATS
    display_model_stats ( *evaluated_pts );
#endif
    
    if ( del_evaluated_pts )
        delete evaluated_pts;
}

/*------------------------------------------------------------------------------------*/
/*  ordering of a list of points based on surrogate (1st) or model (2nd) evaluations  */
/*------------------------------------------------------------------------------------*/
void NOMAD::Evaluator_Control::ordering_lop ( NOMAD::search_type             search             , // IN    : search type
                                             bool                          & stop               , // IN/OUT: stopping criterion
                                             NOMAD::stop_type              & stop_reason        , // OUT   : stopping reason
                                             NOMAD::Barrier                & true_barrier       , // IN/OUT: truth barrier
                                             NOMAD::Barrier                & sgte_barrier        // IN/OUT: surrogate barrier
)

{
    std::list<const NOMAD::Eval_Point *> * evaluated_pts     = new std::list<const NOMAD::Eval_Point *>;
    
    bool sgte_eval_sort = _p.get_sgte_eval_sort() && _eval_lop.size() > 1;
    bool opt_only_sgte  = _p.get_opt_only_sgte ();
    bool snap_to_bounds = _p.get_snap_to_bounds();
    bool modified_list  = false;
    
    const NOMAD::Display    & out = _p.out();
    NOMAD::dd_type display_degree = out.get_display_degree ( search );
    
    NOMAD::success_type success ;
    const NOMAD::Eval_Point *new_feas_inc ; 
    const NOMAD::Eval_Point *new_infeas_inc;
    
    
    // reset the success type:
    true_barrier.reset_success();
    sgte_barrier.reset_success();
    
    
    // use the surrogates to sort the eval. points:
    if ( !opt_only_sgte && sgte_eval_sort ) 
    {
        
        for ( std::set<NOMAD::Priority_Eval_Point>::iterator it = _eval_lop.begin() ; it != _eval_lop.end() ; ++it )
            NOMAD::Cache::get_modifiable_point(*it->get_point()).set_eval_type(NOMAD::SGTE);
        
        
        // evaluate the surrogate:
        private_eval_list_of_points ( search         ,
                                     true_barrier   ,
                                     sgte_barrier   ,
                                     NULL           , // Pareto front = NULL
                                     stop           ,
                                     stop_reason    ,
                                     new_feas_inc   ,
                                     new_infeas_inc ,
                                     success        ,
                                     *evaluated_pts   );
        if ( stop )	
        {
            delete evaluated_pts;
            return;
        }
        
        NOMAD::Eval_Point * x;
        
        // construct a new list of trial points that will be
        // ordered using surrogate values:   
        std::list<const NOMAD::Eval_Point *>::const_iterator
        end = evaluated_pts->end() , it2;
        for ( it2 = evaluated_pts->begin() ; it2 != end ; ++it2 ) 
        {
            
            // Eval_Point construction:
            x = new NOMAD::Eval_Point;
            x->set ( (*it2)->size() , _p.get_bb_nb_outputs() );
            x->set_signature  ( (*it2)->get_signature () );
            x->set_direction  ( (*it2)->get_direction () );
            x->set_poll_center( (*it2)->get_poll_center () );  // Poll center is needed for further testing (not needed when evaluating points)
            x->set_poll_center_type ( (*it2)->get_poll_center_type ()   );
            x->Point::operator = ( **it2 );
            
            modified_list = true;
            
            // add the new point to the ordered list of trial points:
            add_eval_point ( x               ,
                            display_degree  ,
                            snap_to_bounds  ,
                            (*it2)->get_f() ,
                            (*it2)->get_h() ,
                            NOMAD::Double() ,
                            NOMAD::Double()   );
        }
    }
    
    // model ordering:
    // ---------------
    if ( !modified_list && _model_eval_sort && _eval_lop.size() > 1 ) {
        switch ( _p.get_model_eval_sort() ) {
            case NOMAD::TGP_MODEL:
                TGP_model_ordering ( display_degree , modified_list );
                break;
            case NOMAD::QUADRATIC_MODEL:
                quad_model_ordering ( display_degree , modified_list );
                break;
            case NOMAD::NO_MODEL:;
        }
    }
    
    if ( NOMAD::Evaluator_Control::_force_quit || NOMAD::Evaluator::get_force_quit() )
    {
        stop        = true;
        stop_reason = NOMAD::CTRL_C;
    }
    
    delete evaluated_pts;
}



/*--------------------------------------------------------------*/
/*  return if a series of evaluations is opportunistic or not,  */
/*  depending on the search type (private)                      */
/*--------------------------------------------------------------*/
bool NOMAD::Evaluator_Control::is_opportunistic ( NOMAD::search_type t ) const
{
    switch ( t ) {
        case NOMAD::X0_EVAL:
            return false;
        case NOMAD::LH_SEARCH:
            return _p.get_opportunistic_LH();
        case NOMAD::CACHE_SEARCH:
            return _p.get_opportunistic_cache_search();
        default:
            return _p.get_opportunistic_eval();
    }
    return false;
}

/*----------------------------------------------------------------*/
/*                     stop the evaluations ?                     */
/*----------------------------------------------------------------*/
/* . check the opportunistic strategy stopping criterion          */
/* . private method                                               */
/*----------------------------------------------------------------*/
bool NOMAD::Evaluator_Control::stop_evaluations
( const NOMAD::Eval_Point & x                ,
 NOMAD::search_type        search           ,
 int                       k                ,
 int                       nb_points        ,
 bool                      stop             ,
 NOMAD::dd_type            display_degree   ,
 NOMAD::success_type       one_eval_success ,
 NOMAD::success_type       success          ,
 int                       init_nb_eval     ,
 const NOMAD::Double     & f0               ,
 const NOMAD::Barrier    & barrier          ,
 int                     & nb_success       ,
 bool                    & one_for_luck       ) const
{
    // opportunistic evaluation ?
    bool opportunistic = is_opportunistic ( search );
    
    if ( k < nb_points - 1 ) {
        
        if ( stop )
            return true;
        
        if ( opportunistic &&
            ( x.get_eval_type() == NOMAD::TRUTH || _p.get_opt_only_sgte() ) )
        {
            
            if ( one_for_luck && one_eval_success != NOMAD::FULL_SUCCESS ) 
            {
                if ( display_degree == NOMAD::FULL_DISPLAY )
                    _p.out() << std::endl
                    << "opportunistic termination of evaluations (lucky eval)"
                    << std::endl;
                return true;
            }
            
            if ( success == NOMAD::FULL_SUCCESS &&
                check_opportunistic_criterion ( display_degree   ,
                                               one_eval_success ,
                                               init_nb_eval     ,
                                               f0               ,
                                               barrier          ,
                                               nb_success       ,
                                               one_for_luck       ) )
                return true;
        }
    }
    return false;
}

/*-----------------------------------------------------------------*/
/*  check the opportunistic strategy stopping criterion (private)  */
/*            return true to stop the evaluations                  */
/*            return false to continue the evaluations             */
/*-----------------------------------------------------------------*/
bool NOMAD::Evaluator_Control::check_opportunistic_criterion
( NOMAD::dd_type         display_degree   ,
 NOMAD::success_type    one_eval_success ,
 int                    init_nb_eval     ,
 const NOMAD::Double  & f0               ,
 const NOMAD::Barrier & barrier          ,
 int                  & nb_success       ,
 bool                 & one_for_luck       ) const
{
    
    int                    min_nb_success = _p.get_opportunistic_min_nb_success();
    int                    min_eval       = _p.get_opportunistic_min_eval();
    NOMAD::Double          min_f_imprvmt  = _p.get_opportunistic_min_f_imprvmt();
    bool                   lucky_eval     = _p.get_opportunistic_lucky_eval();
    const NOMAD::Display & out            = _p.out();
    
    // min_nb_success:
    if ( min_nb_success > 0 )
    {
        
        if ( one_eval_success == NOMAD::FULL_SUCCESS )
            ++nb_success;
        
        if ( nb_success < min_nb_success ) 
        {
            
            if ( display_degree == NOMAD::FULL_DISPLAY )
                out << std::endl
                << "opport. strategy (nb_success=" << nb_success
                << " < min_nb_success=" << min_nb_success
                << "): continue evaluations"
                << std::endl;
            
            return false;
        }
    }
    
    // min_eval:
    if ( min_eval > 0 )
    {
        
        int eval = _stats.get_eval() - init_nb_eval;
        
        if ( eval < min_eval )
        {
            
            if ( display_degree == NOMAD::FULL_DISPLAY )
                out << std::endl
                << "opport. strategy (eval=" << eval
                << " < min_eval=" << min_eval
                << "): continue evaluations" << std::endl;
            return false;
        }
    }
    
    // min_f_imprvmt:
    if ( min_f_imprvmt.is_defined() ) 
    {
        
        const NOMAD::Eval_Point * bf = barrier.get_best_feasible();
        
        if ( f0.is_defined() && bf ) 
        {
            
            NOMAD::Double f = bf->get_f();
            
            if ( f.is_defined() ) 
            {
                
                NOMAD::Double f_imprvmt = f0.rel_err(f) * 100.0;
                
                if ( f_imprvmt < min_f_imprvmt )
                {
                    
                    if ( display_degree == NOMAD::FULL_DISPLAY )
                        out << std::endl
                        << "opport. strategy (f_improvement="
                        << f_imprvmt << " < min_f_imprvmt=" << min_f_imprvmt
                        << "): continue evaluations" << std::endl;
                    
                    return false;
                }
            }
        }
    }
    
    // lucky_eval:
    if ( lucky_eval && one_eval_success == NOMAD::FULL_SUCCESS )
    {
        one_for_luck = true;
        
        if ( display_degree == NOMAD::FULL_DISPLAY )
            out << std::endl
            << "opport. strategy: one more evaluation for luck"
            << std::endl;
        
        return false;
    }
    
    if ( display_degree == NOMAD::FULL_DISPLAY )
    {
        out << std::endl << "opport. strategy: stop evaluations" ;
        if (_p.get_bb_max_block_size() > 1)
            out << " at the end of the block evaluation";
        out << std::endl;
    }
    
    return true;
}

/*---------------------------------------------------------------*/
/*        display the list of evaluation points (_eval_lop)      */
/*---------------------------------------------------------------*/
void NOMAD::Evaluator_Control::display_eval_lop ( NOMAD::search_type t ) const
{
    const NOMAD::Display & out = _p.out();
    int cnt = 0 , nb = static_cast<int>(_eval_lop.size());
    
    if ( nb == 0 ) {
        out << std::endl << "no evaluation point" << std::endl;
        return;
    }
    
    // open indented block:
    std::ostringstream oss;
    if ( t != NOMAD::UNDEFINED_SEARCH )
        oss << t << " ";
    oss << "evaluation point";
    if ( nb > 1 )
        oss << "s";
    out << std::endl << NOMAD::open_block ( oss.str() ) << std::endl;
    
    // display the points:
    std::set<NOMAD::Priority_Eval_Point>::const_iterator it , end = _eval_lop.end();
    for ( it = _eval_lop.begin() ; it != end ; ++it ) {
        out << "point ";
        out.display_int_w ( ++cnt , nb );
        out << "/" << nb << ": ( ";
        it->get_point()->Point::display ( out                               ,
                                         " "                               ,
                                         2                                 ,
                                         NOMAD::Point::get_display_limit()   );
        out << " )" << std::endl;
    }
    
    // close indented block:
    out.close_block();
}

/*--------------------------------------------------------------*/
/*    add an Eval_Point to the list of points to be evaluated   */
/*--------------------------------------------------------------*/
/*  . x has to be a dynamic object                              */
/*  . it may be deleted into the method and be NULL after that  */
/*  . the point is also snapped to bounds                       */
/*  . periodic variables are checked                            */
/*--------------------------------------------------------------*/
void NOMAD::Evaluator_Control::add_eval_point( NOMAD::Eval_Point  *& x              ,
                                              NOMAD::dd_type        display_degree ,
                                              bool                  snap_to_bounds ,
                                              const NOMAD::Double & f_sgte         ,
                                              const NOMAD::Double & h_sgte         ,
                                              const NOMAD::Double & f_model        ,
                                              const NOMAD::Double & h_model         )
{
    if ( !x )
        return;
    
    const NOMAD::Display & out = _p.out();
    
    // treat the periodic variables:
    NOMAD::Direction * new_dir = NULL;
    
    if ( _p.has_periodic_variables() &&
        x->treat_periodic_variables ( new_dir ) ) 
    {
        
        if ( new_dir && new_dir->norm() == 0.0 ) 
        {
            
            if ( display_degree == NOMAD::FULL_DISPLAY )
                out << "point #" << x->get_tag()
                << " is flushed (||dir||==0)"
                << std::endl;
            
            delete x;
            x = NULL;
            
            delete new_dir;
            
            return;
        }
    }
    delete new_dir;
    
    if ( snap_to_bounds && x->snap_to_bounds() )
    {
        
        if ( display_degree == NOMAD::FULL_DISPLAY ) 
        {
            out << std::endl << "point #" << x->get_tag() << " ";
            if ( x->get_direction() && x->get_direction()->get_index() >= 0 )
                out << "(dir " << x->get_direction()->get_index() << ") ";
            out << "has been snapped to bounds" << std::endl;
        }
        
        if ( x->get_direction() && x->get_direction()->norm() == 0.0 )
        {
            
            if ( display_degree == NOMAD::FULL_DISPLAY )
                out << "point #" << x->get_tag()
                << " is flushed (||dir||==0)"
                << std::endl;
            delete x;
            x = NULL;
            
            return;
        }
    }
    
    // creation of the Priority_Eval_Point:
    NOMAD::Priority_Eval_Point pep ( x , _p.get_h_min() );
    
    // ordering elements of Priority_Eval_Point's:
    // -------------------------------------------
    
    // 1. surrogate values for f and h:
    pep.set_f_sgte ( f_sgte );
    pep.set_h_sgte ( h_sgte );
    
    // 2. model values for f and h:
    pep.set_f_model ( f_model );
    pep.set_h_model ( h_model );
    
    if ( x->get_direction() )
    {
        
        // get the signature:
        NOMAD::Signature * signature = x->get_signature();
        if ( !signature )
            throw NOMAD::Exception ( "Evaluator_Control.cpp" , __LINE__ ,
                                    "Evaluator_Control::add_eval_point(): the point has no signature" );
        
        // angle with last successful directions (feasible)
        const NOMAD::Direction & feas_success_dir = signature->get_feas_success_dir();
        if ( feas_success_dir.is_defined() &&
            x->get_poll_center_type() == NOMAD::FEASIBLE  )  
            pep.set_angle_success_dir ( feas_success_dir.get_angle ( *x->get_direction() ) );
        
        // angle with last infeasible success direction:
        const NOMAD::Direction & infeas_success_dir = signature->get_infeas_success_dir();
        if ( infeas_success_dir.is_defined() &&
            x->get_poll_center_type() == NOMAD::INFEASIBLE  )
            pep.set_angle_success_dir ( infeas_success_dir.get_angle ( *x->get_direction() ) );
        
    }
    
    
    
    // insertion of the point in _eval_lop:
    // ------------------------------------
    size_t size_before = _eval_lop.size();
    
    _eval_lop.insert ( pep );
    
    if ( _eval_lop.size() == size_before )
    {
        delete x;
        x = NULL;
    }
}

#ifdef MODEL_STATS
/*------------------------------------------------------------------*/
/*  display stats on an evaluation for which a model has been used  */
/*------------------------------------------------------------------*/
/*  The displayed stats are:                                        */
/*                                                                  */
/*     use (1:model_search, 2:model_ordering)                       */
/*     mesh_index                                                   */
/*     cardinality of Y                                             */
/*     width of Y                                                   */
/*     Y condition number                                           */
/*     h value, model for h, relative error (if constraints)        */
/*     f value, model for f, relative error                         */
/*                                                                  */
/*------------------------------------------------------------------*/
void NOMAD::Evaluator_Control::display_model_stats
( const std::list<const NOMAD::Eval_Point *> & evaluated_pts ) const
{
    const NOMAD::Display & out = _p.out();
    
    NOMAD::Double h , mh , eh , f , mf , ef;
    
    std::list<const NOMAD::Eval_Point *>::const_iterator it , end = evaluated_pts.end();
    for ( it = evaluated_pts.begin() ; it != end ; ++it ) {
        if ( *it && (*it)->get_mod_use() >= 0 ) {
            
            if ( _p.has_constraints() ) {
                h  = (*it)->get_h ();
                mh = (*it)->get_mh();
            }
            else
                h = mh = 0.0;
            
            f  = (*it)->get_f ();
            mf = (*it)->get_mf();
            
            if ( h.is_defined() && mh.is_defined() && f.is_defined() && mf.is_defined() ) {
                
                ef = f.rel_err ( mf ) * 100.0;
                
                out << (*it)->get_mod_use()
                << " " << std::setw(3) << NOMAD::Mesh::get_mesh_index()
                << " " << std::setw(4) << (*it)->get_nY()
                << " ";
                
                (*it)->get_Yw  ().display ( out , "%12.3g" ); out << " ";
                (*it)->get_cond().display ( out , "%12.3g" ); out << " ";
                if ( _p.has_constraints() ) {
                    eh = h.rel_err ( mh ) * 100.0;
                    h.display  ( out , "%14.3g" ); out << " ";
                    mh.display ( out , "%14.3g" ); out << " ";
                    eh.display ( out , "%14.3g" ); out << " ";
                }
                f.display  ( out , "%14.3g" ); out << " ";
                mf.display ( out , "%14.3g" ); out << " ";
                ef.display ( out , "%14.3g" );
                
                out << std::endl;
            }
            
            (*it)->clear_model_data();
        }
    }
}
#endif
