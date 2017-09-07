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
 \file   Display.cpp
 \brief  Custom class for display (implementation)
 \author Sebastien Le Digabel
 \date   2010-03-30
 \see    Display.hpp
 */
#include "Display.hpp"

/*---------------------------------------------------------*/
/*                    open an indented block               */
/*---------------------------------------------------------*/
void NOMAD::Display::open_block ( const std::string & msg ) const
{
    // open braces:
    if ( _newline )
        _out << _indent_str;
    if ( !msg.empty() )
        _out << msg << " ";
    _out << _open_brace << std::endl;
    
    _newline = true;
    
    // add a tabulation to the indentation string:
    _indent_str += '\t';
}

/*---------------------------------------------------------*/
/*                   close an indented block               */
/*---------------------------------------------------------*/
void NOMAD::Display::close_block ( const std::string & msg ) const
{
    _newline = true;
    
    if ( _indent_str.empty() )
        return;
    
    // remove a tabulation to the indentation string:
    _indent_str.erase ( 0 , 1 );
    
    // close braces:
    _out << _indent_str << _closed_brace << " " << msg << std::endl;
}

/*---------------------------------------------------------*/
/*                   set the display degrees               */
/*---------------------------------------------------------*/
void NOMAD::Display::set_degrees ( NOMAD::dd_type gen_dd    ,
                                  NOMAD::dd_type search_dd ,
                                  NOMAD::dd_type poll_dd   ,
                                  NOMAD::dd_type iter_dd     )
{
    // max = max { gen_dd , search_dd , poll_dd , iter_dd } :
    NOMAD::dd_type max = gen_dd;
    if ( search_dd > gen_dd )
        max = search_dd;
    if ( poll_dd > max )
        max = poll_dd;
    if ( iter_dd > max )
        max = iter_dd;
    
    // max=0: all to 0:
    if ( max == NOMAD::NO_DISPLAY )
        _gen_dd = _search_dd = _poll_dd = _iter_dd = NOMAD::NO_DISPLAY;
    
    // max=1: all to 1:
    else if ( max == NOMAD::MINIMAL_DISPLAY )
        _gen_dd = _search_dd = _poll_dd = _iter_dd = NOMAD::MINIMAL_DISPLAY;
    
    // max=2: all to 2:
    else if ( max == NOMAD::NORMAL_DISPLAY )
        _gen_dd = _search_dd = _poll_dd = _iter_dd = NOMAD::NORMAL_DISPLAY;
    
    // max=3: 0-->0, 1-->0, 2->0 and 3-->3:
    else
    {
        _gen_dd    = (gen_dd   ==NOMAD::FULL_DISPLAY)? NOMAD::FULL_DISPLAY:NOMAD::NO_DISPLAY;
        _search_dd = (search_dd==NOMAD::FULL_DISPLAY)? NOMAD::FULL_DISPLAY:NOMAD::NO_DISPLAY;
        _poll_dd   = (poll_dd  ==NOMAD::FULL_DISPLAY)? NOMAD::FULL_DISPLAY:NOMAD::NO_DISPLAY;
        _iter_dd   = (iter_dd  ==NOMAD::FULL_DISPLAY)? NOMAD::FULL_DISPLAY:NOMAD::NO_DISPLAY;
    }
}

/*---------------------------------------------------*/
/*  get the display degrees as a string (of size 4)  */
/*---------------------------------------------------*/
void NOMAD::Display::get_display_degree ( std::string & dd ) const
{
    dd.resize ( 4 );
    dd[0] = NOMAD::Display::dd_to_char ( _gen_dd    );
    dd[1] = NOMAD::Display::dd_to_char ( _search_dd );
    dd[2] = NOMAD::Display::dd_to_char ( _poll_dd   );
    dd[3] = NOMAD::Display::dd_to_char ( _iter_dd   );
}

/*---------------------------------------------*/
/*  get the display degree from a search type  */
/*---------------------------------------------*/
NOMAD::dd_type NOMAD::Display::get_display_degree ( NOMAD::search_type search ) const
{
    if ( search == NOMAD::X0_EVAL )
        return _gen_dd;
    if ( search == NOMAD::POLL || search == NOMAD::EXTENDED_POLL )
        return _poll_dd;
    return _search_dd;
}

/*------------------------------------------*/
/*  convert a dd_type into a char (static)  */
/*------------------------------------------*/
char NOMAD::Display::dd_to_char ( NOMAD::dd_type dd )
{
    if ( dd == NOMAD::NO_DISPLAY )
        return '0';
    if ( dd == NOMAD::MINIMAL_DISPLAY)
        return '1';
    if ( dd == NOMAD::NORMAL_DISPLAY )
        return '2';
    return '3';
}

/*------------------------------------------*/
/*  convert a dd_type into an int (static)  */
/*------------------------------------------*/
int NOMAD::Display::dd_to_int ( NOMAD::dd_type dd )
{
    if ( dd == NOMAD::NO_DISPLAY )
        return 0;
    if ( dd == NOMAD::MINIMAL_DISPLAY )
        return 1;
    if ( dd == NOMAD::NORMAL_DISPLAY )
        return 2;
    return 3;
}

/*------------------------------------------*/
/*  convert an int into a dd_type (static)  */
/*------------------------------------------*/
NOMAD::dd_type NOMAD::Display::int_to_dd ( int dd )
{
    if ( dd <= 0 )
        return NOMAD::NO_DISPLAY;
    if ( dd == 1 )
        return NOMAD::MINIMAL_DISPLAY;
    if ( dd == 2 )
        return NOMAD::NORMAL_DISPLAY;
    
    return NOMAD::FULL_DISPLAY;
}

/*---------------------------------------------------------*/
/*          display an integer with a specific width       */
/*---------------------------------------------------------*/
void NOMAD::Display::display_int_w ( int i , int max_i ) const
{
    (*this) << std::setw ( (max_i > 0) ?
                          (1+int(log(static_cast<double>(max_i))/NOMAD::LOG10)) : 1 )
    << i;
}

/*--------------------------------------------------------*/
/*           extract display format from a string         */
/*--------------------------------------------------------*/
/*  example: s="$%5.2f" becomes s="$" and format="%5.2f"  */
/*--------------------------------------------------------*/
void NOMAD::Display::extract_display_format ( std::string & s , std::string & format )
{
    format.clear();
    if ( s.empty() )
        return;
    size_t k = s.find("%");
    size_t n = s.size();
    if ( k < n )
    {
        if ( k > 0 && s[k-1]=='\\' )
        {
            std::string s1 = s.substr ( 0 , k-1 );
            std::string s2 = s.substr ( k , n-k );
            s = s1 + s2;
        }
        else
        {
            format = s.substr ( k , n-k );
            s      = s.substr ( 0 , k   );
        }
    }
}

/*---------------------------------------------------------*/
/*         to display a duration with a smart format       */
/*              (t is in seconds and integer)              */
/*---------------------------------------------------------*/
void NOMAD::Display::display_time ( int t ) const
{
    if ( t > 0 )
    {
        int h = t / 3600;
        t = t % 3600;
        int m = t / 60;
        t = t % 60;
        if ( h > 0 )
            (*this) << h << "h ";
        if ( m > 0 || h > 0)
            (*this) << m << "m ";
    }
    else
        t = 0;
    (*this) << t << "s";
}

/*-----------------------------------------------------------------*/
/*                   to display a size in memory                   */
/*-----------------------------------------------------------------*/
void NOMAD::Display::display_size_of ( float size ) const
{
    if ( size < 1024 )
        (*this) << static_cast<int>(size) << " B";
    else if ( size < 1048576 )
        (*this) << static_cast<int> ( 10 * size / 1024.0 ) / 10.0 << " KB";
    else if ( size < 1073741824 )
        (*this) << static_cast<int> ( 10 * size / 1048576.0 ) / 10.0 << " MB";
    else
        (*this) << static_cast<int> ( 10 * size / 1073741824.0 ) / 10.0 << " GB";
}

/*-------------------------------------------------------------------*/
/*          get the display_stats_type from a string (static)        */
/*-------------------------------------------------------------------*/
NOMAD::display_stats_type NOMAD::Display::get_display_stats_type
( const std::string & s )
{
    int         idst;
    std::string ss = s , keyword;
    NOMAD::toupper ( ss );
    NOMAD::display_stats_type dst = NOMAD::DS_OBJ;
    while ( dst < NOMAD::DS_UNDEFINED )
    {
        keyword = get_display_stats_keyword ( dst );
        if ( keyword == ss )
            return dst;
        idst = dst;
        ++idst;
        dst = static_cast<display_stats_type> ( idst );
    }
    return NOMAD::DS_UNDEFINED;
}

/*-------------------------------------------------------------------*/
/*     to get the keyword associated with a display_stats_type       */
/*  (used in Parameters::read() or Parameters::set_DISPLAY_STATS())  */
/*  (static)                                                         */
/*-------------------------------------------------------------------*/
std::string NOMAD::Display::get_display_stats_keyword ( NOMAD::display_stats_type dst )
{
    std::string s;
    switch ( dst )
    {
        case NOMAD::DS_OBJ:
            s = "OBJ";
            break;
        case NOMAD::DS_MESH_INDEX:
            s = "MESH_INDEX";
            break;
        case NOMAD::DS_DELTA_M:
        case NOMAD::DS_MESH_SIZE:
            s = "MESH_SIZE";
            break;
        case NOMAD::DS_DELTA_P:
        case NOMAD::DS_POLL_SIZE:
            s = "POLL_SIZE";
            break;
        case NOMAD::DS_EVAL:
            s = "EVAL";
            break;
        case NOMAD::DS_SIM_BBE:
            s = "SIM_BBE";
            break;
        case NOMAD::DS_BBE:
            s = "BBE";
            break;
        case NOMAD::DS_BLK_EVA:
            s = "BLK_EVA";
            break;
        case NOMAD::DS_SGTE:
            s = "SGTE";
            break;
        case NOMAD::DS_BBO:
            s = "BBO";
            break;
        case NOMAD::DS_SOL:
            s = "SOL";
            break;
        case NOMAD::DS_VAR:
            s = "VAR";
            break;
        case NOMAD::DS_TIME:
            s = "TIME";
            break;
        case NOMAD::DS_STAT_SUM:
            s = "STAT_SUM";
            break;
        case NOMAD::DS_STAT_AVG:
            s = "STAT_AVG";
            break;
        case NOMAD::DS_UNDEFINED:
            s = "undefined";
            break;
    }
    return s;
}

/*-----------------------------------------------------------------*/
/*                     to display a display degree                 */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream & out , NOMAD::dd_type dd )
{
    switch ( dd )
    {
        case NOMAD::NO_DISPLAY:
            out << "no display (0)";
            break;
        case NOMAD::MINIMAL_DISPLAY:
            out << "minimal display (1)";
            break;
        case NOMAD::NORMAL_DISPLAY:
            out << "normal (2)";
            break;
        case NOMAD::FULL_DISPLAY:
        default:
            out << "full (3)";
    }
    return out;
}

/*-----------------------------------------------------*/
/*              to display a success_type              */
/*-----------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream & out , NOMAD::success_type st )
{
    switch ( st )
    {
        case NOMAD::FULL_SUCCESS:
            out << "dominating";
            break;
        case NOMAD::PARTIAL_SUCCESS:
            out << "improving";
            break;
        case NOMAD::UNSUCCESSFUL:
            out << "unsuccessful";
    }
    return out;
}

/*-----------------------------------------------------------------*/
/*                 to display a vector of bb-input types           */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream                            & out   ,
                                   const std::vector<NOMAD::bb_input_type> & bbits   )
{
    if ( bbits.empty() )
        return out;
    size_t n = bbits.size()-1;
    for ( size_t k = 0 ; k < n ; ++k )
        out << std::setw(8) << bbits[k] << " ";
    out << std::setw(8) << bbits[n];
    return out;
}

/*-----------------------------------------------------------------*/
/*                     to display a bb-input type                  */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream & out , NOMAD::bb_input_type bi )
{
    switch ( bi ) {
        case NOMAD::CONTINUOUS:
            out << "cont(R)";
            break;
        case NOMAD::CATEGORICAL:
            out << "cat(C)";
            break;
        case NOMAD::INTEGER:
            out << "int(I)";
            break;
        case NOMAD::BINARY:
            out << "bin(B)";
    }
    return out;
}

/*-----------------------------------------------------------------*/
/*                     to display a search type                    */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream & out , NOMAD::search_type st )
{
    switch ( st )
    {
        case NOMAD::POLL:
            out << "poll";
            break;
        case NOMAD::EXTENDED_POLL:
            out << "extended poll";
            break;
        case NOMAD::SEARCH:
            out << "search";
            break;
        case NOMAD::CACHE_SEARCH:
            out << "cache search";
            break;
        case NOMAD::USER_SEARCH:
            out << "user search";
            break;
        case NOMAD::SPEC_SEARCH:
            out << "speculative search";
            break;
        case NOMAD::LH_SEARCH:
            out << "LH search";
            break;
        case NOMAD::LH_SEARCH_P1:
            out << "LH search - Phase one";
            break;
        case NOMAD::P1_SEARCH:
            out << "Phase one search";
            break;
        case NOMAD::MODEL_SEARCH:
            out << "model search";
            break;
        case NOMAD::VNS_SEARCH:
            out << "VNS search";
            break;
        case NOMAD::X0_EVAL:
            out << "x0 evaluation";
            break;
        case NOMAD::ASYNCHRONOUS:
            out << "asynchronous final evaluations";
            break;
        case NOMAD::UNDEFINED_SEARCH:
            out << "undefined";
    }
    return out;
}

/*-----------------------------------------------------------------*/
/*                     to display a model type                    */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream & out , NOMAD::model_type mt )
{
    switch ( mt )
    {
        case NOMAD::QUADRATIC_MODEL:
            out << "quadratic";
            break;
        case NOMAD::TGP_MODEL:
            out << "TGP";
            break;
        case NOMAD::NO_MODEL:
            out << "no models";
    }
    return out;
}

/*-----------------------------------------------------------------*/
/*                      to display a TGP mode                      */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream & out , NOMAD::TGP_mode_type m )
{
    switch ( m )
    {
        case NOMAD::TGP_FAST:
            out << "fast";
            break;
        case NOMAD::TGP_PRECISE:
            out << "precise";
            break;
        case NOMAD::TGP_USER:
            out << "user";
            break;
    }
    return out;
}

/*-----------------------------------------------------------------*/
/*                  to display an evaluation type                  */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream & out , NOMAD::eval_type et )
{
    switch ( et )
    {
        case NOMAD::TRUTH:
            out << "truth";
            break;
        case NOMAD::SGTE:
            out << "surrogate";
    }
    return out;
}

/*-----------------------------------------------------------------*/
/*               to display an evaluation status type              */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream & out , NOMAD::eval_status_type est )
{
    switch ( est )
    {
        case NOMAD::EVAL_FAIL:
            out << "fail";
            break;
        case NOMAD::EVAL_OK:
            out << "ok";
            break;
        case NOMAD::EVAL_IN_PROGRESS:
            out << "in progress";
            break;
        case NOMAD::UNDEFINED_STATUS:
            out << "undefined";
            break;
        case NOMAD::EVAL_USER_REJECT:
            out << "rejected";
            break;
            
    }
    return out;
}

/*-----------------------------------------------------------------*/
/*                      to display a stop reason                   */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream & out , NOMAD::stop_type st )
{
    switch ( st )
    {
        case NOMAD::DELTA_M_MIN_REACHED:
            out << "min mesh size";
            break;
        case NOMAD::DELTA_P_MIN_REACHED:
            out << "min poll size";
            break;
        case NOMAD::L_MAX_REACHED:
            out << "max mesh index";
            break;
        case NOMAD::L_MIN_REACHED:
            out << "min mesh index";
            break;
        case NOMAD::L_LIMITS_REACHED:
            out << "mesh index limits";
            break;
        case NOMAD::XL_LIMITS_REACHED:
            out << "mesh index limits";
            break;
        case NOMAD::MAX_TIME_REACHED:
            out << "max time";
            break;
        case NOMAD::MAX_BB_EVAL_REACHED:
            out << "max number of blackbox evaluations";
            break;
        case NOMAD::MAX_SGTE_EVAL_REACHED:
            out << "max number of sgte evaluations";
            break;
        case NOMAD::MAX_EVAL_REACHED:
            out << "max number of evaluations";
            break;
        case NOMAD::MAX_SIM_BB_EVAL_REACHED:
            out << "max number of sim. bb evaluations";
            break;
        case NOMAD::MAX_ITER_REACHED:
            out << "max number of iterations";
            break;
        case NOMAD::MAX_CONS_FAILED_ITER:
            out << "max number of consecutive failed iterations";
            break;
        case NOMAD::FEAS_REACHED:
            out << "feasibility achieved";
            break;
        case NOMAD::F_TARGET_REACHED:
            out << "objective target reached";
            break;
        case NOMAD::L_CURVE_TARGET_REACHED:
            out << "L-curve target reached";
            break;
        case NOMAD::P1_FAIL:
            out << "phase one failed";
            break;
        case NOMAD::STAT_SUM_TARGET_REACHED:
            out << "stat sum target reached";
            break;
        case NOMAD::X0_FAIL:
            out << "problem with starting point evaluation";
            break;
        case NOMAD::MESH_PREC_REACHED:
            out << "mesh size reached NOMAD precision";
            break;
        case NOMAD::MULTI_MAX_BB_REACHED:
            out << "max number of bb evaluations";
            break;
        case NOMAD::MULTI_NB_MADS_RUNS_REACHED:
            out << "max number of MADS runs";
            break;
        case NOMAD::MULTI_STAGNATION:
            out << "stagnation of the multi-obj. algo.";
            break;
        case NOMAD::MULTI_NO_PARETO_PTS:
            out << "initial runs cannot find Pareto points";
            break;
        case NOMAD::MAX_CACHE_MEMORY_REACHED:
            out << "max cache memory reached";
            break;
        case NOMAD::CTRL_C:
            out << "terminated by ctrl-c";
            break;
        case NOMAD::USER_STOPPED:
            out << "terminated by the user inside Evaluator::update_iteration()";
            break;
        case NOMAD::UNKNOWN_STOP_REASON:
        case NOMAD::NO_STOP:
            out << "unknown";
            break;
        case NOMAD::ERROR:
            out << "error";
            break;
    }
    return out;
}

/*-----------------------------------------------------------------*/
/*                    to display a bb_output_type                  */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream & out , NOMAD::bb_output_type bbt )
{
    switch ( bbt )
    {
        case NOMAD::OBJ:
            out << "OBJ";
            break;
        case NOMAD::EB:
            out << "EB";
            break;
        case NOMAD::PB:
            out << "PB";
            break;
        case NOMAD::PEB_P:
            out << "PEB(P)";
            break;
        case NOMAD::PEB_E:
            out << "PEB(E)";
            break;
        case NOMAD::FILTER:
            out << "F";
            break;
        case NOMAD::CNT_EVAL:
            out << "CNT_EVAL";
            break;
        case NOMAD::STAT_AVG:
            out << "STAT_AVG";
            break;
        case NOMAD::STAT_SUM:
            out << "STAT_SUM";
            break;
        case NOMAD::UNDEFINED_BBO:
            out << "-";
    }
    return out;
}

/*-----------------------------------------------------------------*/
/*                   to display a interpolation_type               */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream              & out ,
                                   NOMAD::interpolation_type   it    )
{
    switch ( it )
    {
        case NOMAD::MFN:
            out << "Minimum Frobenius Norm interpolation";
            break;
        case NOMAD::REGRESSION:
            out << "regression";
            break;
        case NOMAD::WP_REGRESSION:
            out << "well-poised regression";
            break;
        case NOMAD::UNDEFINED_INTERPOLATION_TYPE:
            out << "undefined";
    }
    return out;
}

/*-----------------------------------------------------------------*/
/*                       to display a hnorm_type                   */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream & out , NOMAD::hnorm_type hnorm )
{
    switch ( hnorm )
    {
        case NOMAD::L1:
            out << "L1";
            break;
        case NOMAD::L2:
            out << "L2";
            break;
        case NOMAD::LINF:
            out << "Linf";
    }
    return out;
}

/*-----------------------------------------------------------------*/
/*              to display a multi-obj. formulation type           */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream & out , NOMAD::multi_formulation_type f )
{
    switch ( f )
    {
        case NOMAD::NORMALIZED:
            out << "normalized";
            break;
        case NOMAD::PRODUCT:
            out << "product";
            break;
        case NOMAD::DIST_L1:
            out << "distance L1";
            break;
        case NOMAD::DIST_L2:
            out << "distance L2";
            break;
        case NOMAD::DIST_LINF:
            out << "distance Linf";
            break;
        case NOMAD::UNDEFINED_FORMULATION:
            out << "undefined";
            break;
    }
    return out;
}

/*-----------------------------------------------------------------*/
/*                      to display a direction_type                */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream          & out ,
                                   NOMAD::direction_type   dt    )
{
    switch ( dt )
    {
        case NOMAD::ORTHO_1:
            out << "Ortho-MADS 1";
            break;
        case NOMAD::ORTHO_2:
            out << "Ortho-MADS 2";
            break;
        case NOMAD::ORTHO_NP1_QUAD:
            out << "Ortho-MADS n+1 QUAD";
            break;
        case NOMAD::ORTHO_NP1_NEG:
            out << "Ortho-MADS n+1 NEG";
            break;
        case NOMAD::DYN_ADDED:
            out << "Dynamic n+1th direction";
            break;
        case NOMAD::ORTHO_2N:
            out << "Ortho-MADS 2n";
            break;
        case NOMAD::LT_1:
            out << "LT-MADS 1";
            break;
        case NOMAD::LT_2:
            out << "LT-MADS 2";
            break;
        case NOMAD::LT_2N:
            out << "LT-MADS 2n";
            break;
        case NOMAD::LT_NP1:
            out << "LT-MADS n+1";
            break;
        case NOMAD::GPS_BINARY:
            out << "GPS n, binary";
            break;
        case NOMAD::GPS_2N_STATIC:
            out << "GPS 2n, static";
            break;
        case NOMAD::GPS_2N_RAND:
            out << "GPS 2n, random";
            break;
        case NOMAD::GPS_NP1_STATIC_UNIFORM:
            out << "GPS n+1, static, uniform angles";
            break;
        case NOMAD::GPS_NP1_STATIC:
            out << "GPS n+1, static";
            break;
        case NOMAD::GPS_NP1_RAND_UNIFORM:
            out << "GPS n+1, random, uniform angles";
            break;
        case NOMAD::GPS_NP1_RAND:
            out << "GPS n+1, random";
            break;
        case NOMAD::NO_DIRECTION:
            out << "none";
            break;
        case NOMAD::MODEL_SEARCH_DIR:
            out << "model search direction";
            break;
        case NOMAD::UNDEFINED_DIRECTION:
            out << "undefined";
            break;
        case NOMAD::PROSPECT_DIR:
            out << "Prospect direction";
            break;
    }
    return out;
}

/*-----------------------------------------------------------------*/
/*                  to display a display_stats_type                */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream & out , NOMAD::display_stats_type dst )
{
    out << NOMAD::Display::get_display_stats_keyword ( dst );
    return out;
}

/*-----------------------------------------------------------------*/
/*                    to display a check_failed_type               */
/*-----------------------------------------------------------------*/
std::ostream & NOMAD::operator << ( std::ostream & out , NOMAD::check_failed_type cf )
{
    switch ( cf )
    {
        case NOMAD::CHECK_OK:
            out << "ok";
            break;
        case NOMAD::LB_FAIL:
            out << "lower bound";
            break;
        case NOMAD::UB_FAIL:
            out << "upper bound";
            break;
        case NOMAD::FIX_VAR_FAIL:
            out << "fixed variable";
            break;
        case NOMAD::BIN_FAIL:
            out << "binary variable";
            break;
        case NOMAD::CAT_FAIL:
            out << "categorical variable";
            break;
        case NOMAD::INT_FAIL:
            out << "integer variable";
            break;
    }
    return out;
}
