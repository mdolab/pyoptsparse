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
 \file   Eval_Point.cpp
 \brief  Evaluation point (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-14
 \see    Eval_Point.hpp
 */
#include "Cache_File_Point.hpp"
#include "Eval_Point.hpp"
#include "Slave.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
int NOMAD::Eval_Point::_current_tag = 0;
int NOMAD::Eval_Point::_current_bbe = 0;
int NOMAD::Eval_Point::_current_sgte_bbe = 0;
/*---------------------------------------------------------------------*/
/*                            constructor 1                            */
/*---------------------------------------------------------------------*/
NOMAD::Eval_Point::Eval_Point ( void )
: _tag              ( NOMAD::Eval_Point::_current_tag++ ) ,
_signature        ( NULL                              ) ,
_in_cache         ( false                             ) ,
_current_run      ( false                             ) ,
_eval_type        ( NOMAD::TRUTH                      ) ,
_direction        ( NULL                              ) ,
_poll_center_type ( NOMAD::UNDEFINED_POLL_CENTER_TYPE ) ,
_eval_status      ( NOMAD::UNDEFINED_STATUS           ) ,
_EB_ok            ( true                              )
{
#ifdef MODEL_STATS
    _mod_use = -1;
    _nY      = -1;
#endif
}

/*---------------------------------------------------------------------*/
/*                            constructor 2                            */
/*---------------------------------------------------------------------*/
NOMAD::Eval_Point::Eval_Point ( int n , int m )
: NOMAD::Point      ( n                                 ) ,
_tag              ( NOMAD::Eval_Point::_current_tag++ ) ,
_signature        ( NULL                              ) ,
_in_cache         ( false                             ) ,
_current_run      ( false                             ) ,
_eval_type        ( NOMAD::TRUTH                      ) ,
_direction        ( NULL                              ) ,
_poll_center_type ( NOMAD::UNDEFINED_POLL_CENTER_TYPE ) ,
_eval_status      ( NOMAD::UNDEFINED_STATUS           ) ,
_EB_ok            ( true                              ) ,
_bb_outputs       ( m                                 )
{
#ifdef MODEL_STATS
    _mod_use = -1;
    _nY      = -1;
#endif
}

/*---------------------------------------------------------------------*/
/*            constructor 3 (used for model evalaution)                */
/*---------------------------------------------------------------------*/
NOMAD::Eval_Point::Eval_Point ( const NOMAD::Point & x , int m )
: NOMAD::Point      ( x                                ) ,
_direction        ( NULL                              ) ,
_bb_outputs       ( m                                 )
{
#ifdef MODEL_STATS
    _mod_use = -1;
    _nY      = -1;
#endif
}




/*---------------------------------------------------------------------*/
/*                  constructor 4 ( used in Cache::load() )            */
/*---------------------------------------------------------------------*/
NOMAD::Eval_Point::Eval_Point ( const NOMAD::Cache_File_Point & x , NOMAD::eval_type et )
: NOMAD::Point      ( x.get_n()                         ) ,
_tag              ( NOMAD::Eval_Point::_current_tag++ ) ,
_signature        ( NULL                              ) ,
_in_cache         ( false                             ) ,
_current_run      ( false                             ) ,
_eval_type        ( et                                ) ,
_direction        ( NULL                              ) ,
_poll_center_type ( NOMAD::UNDEFINED_POLL_CENTER_TYPE ) ,
_EB_ok            ( true                              ) ,
_bb_outputs       ( x.get_bb_outputs()                )
{
    int n = size();
    for ( int i = 0 ; i < n ; ++i )
        (*this)[i] = x.get_coord(i);
    
    switch ( x.get_eval_status() ) {
        case 0:
            _eval_status = NOMAD::EVAL_FAIL;
            break;
        case 1:
            _eval_status = NOMAD::EVAL_OK;
            break;
        case 2:
            _eval_status = NOMAD::EVAL_IN_PROGRESS;
            break;
        case 3:
            _eval_status = NOMAD::UNDEFINED_STATUS;
            break;
    }
    
#ifdef MODEL_STATS
    _mod_use = -1;
    _nY      = -1;
#endif
}


/*---------------------------------------------------------------------*/
/*                           copy constructor                          */
/*---------------------------------------------------------------------*/
NOMAD::Eval_Point::Eval_Point ( const Eval_Point & x )
: NOMAD::Point        ( x.get_n()                         ) ,
_tag                ( NOMAD::Eval_Point::_current_tag++ ) ,
_signature          ( x._signature                      ) ,
_f                  ( x._f                              ) ,
_h                  ( x._h                              ) ,
_in_cache           ( x._in_cache                       ) ,
_current_run        ( x._current_run                    ) ,
_eval_type          ( x._eval_type                      ) ,
_direction          ( NULL                              ) ,
_poll_center_type   ( x._poll_center_type               ) ,
_eval_status        ( x._eval_status                    ) ,
_EB_ok              ( x._EB_ok                          ) ,
_bb_outputs         ( x.get_bb_outputs()                ) ,
_user_eval_priority ( x._user_eval_priority             ) ,
_rand_eval_priority ( x._rand_eval_priority             )
{
    // point coordinates:
    int n = size();
    for ( int i = 0 ; i < n ; ++i )
        (*this)[i] = x[i];
    
    // _direction:
    if ( x._direction )
        _direction = new Direction ( *x._direction );
    
#ifdef MODEL_STATS
    set_model_data ( x );
#endif
}

/*---------------------------------------------------------------------*/
/*                               destructor                            */
/*---------------------------------------------------------------------*/
NOMAD::Eval_Point::~Eval_Point ( void )
{
    if ( _direction )
        delete _direction;
}

/*-------------------------------------------------------*/
/*  SET methods used to complete a default construction  */
/*-------------------------------------------------------*/
void NOMAD::Eval_Point::set ( int n , int m )
{
    reset ( n );
    _bb_outputs.reset ( m );
}

void NOMAD::Eval_Point::set ( const NOMAD::Point & x , int m )
{
    NOMAD::Point::operator = ( x );
    _bb_outputs.reset ( m );
}

/*-------------------------------------------------------*/
/*           manually set the tag of a point             */
/*   (used in parallel version so that all points have   */
/*    still unique tags                                  */
/*-------------------------------------------------------*/
void NOMAD::Eval_Point::set_tag ( int tag )
{
    _tag = tag;
    NOMAD::Eval_Point::_current_tag = tag+1;
}

/*-------------------------------------------------------*/
/*           increment counted black box evaluations     */
/* All points have  unique tags                          */
/*-------------------------------------------------------*/
void NOMAD::Eval_Point::increment_bbe ( void )
{
    NOMAD::Eval_Point::_current_bbe ++;
    _bbe=_current_bbe;
}


/*-------------------------------------------------------*/
/*           increment counted black box evaluations     */
/* All points have  unique tags                          */
/*-------------------------------------------------------*/
void NOMAD::Eval_Point::increment_sgte_bbe ( void )
{
    NOMAD::Eval_Point::_current_sgte_bbe ++;
    _sgte_bbe=_current_sgte_bbe;
}


void NOMAD::Eval_Point::set_direction ( const NOMAD::Direction * dir )
{
    delete _direction;
    _direction = ( dir ) ? new NOMAD::Direction ( *dir ) : NULL;
}

void NOMAD::Eval_Point::set_poll_center ( const NOMAD::Eval_Point * pc )
{
    _poll_center=pc;
}

/*-------------------------------------------------------*/
/*                    set the signature                  */
/*-------------------------------------------------------*/
void NOMAD::Eval_Point::set_signature ( NOMAD::Signature * s )
{
    if ( !s ) {
        _signature = NULL;
        return;
    }
    
    if ( !s->is_compatible(*this) )
        throw NOMAD::Exception ( "Eval_Point.cpp" , __LINE__ ,
                                "x.Eval_Point::set_signature(s): x and s are incompatible" );
    
    _signature = s;
}

/*------------------------------------------*/
/*             get the signature            */
/*------------------------------------------*/
NOMAD::Signature * NOMAD::Eval_Point::get_signature ( void ) const
{
#ifdef USE_MPI
    if ( !NOMAD::Slave::is_master() )
        throw NOMAD::Exception ( "Eval_Point.cpp" , __LINE__ ,
                                "Eval_Point::get_signature(): cannot be invoked by slave processes" );
#endif
    return _signature;
}

/*-------------------------------------------------------*/
/*                          sizeof                       */
/*-------------------------------------------------------*/
int NOMAD::Eval_Point::size_of ( void ) const
{
    return NOMAD::Point::size_of () +
    _bb_outputs.size_of        () +
    _f.size_of                 () +
    _h.size_of                 () +
    _user_eval_priority.size_of() +
    _rand_eval_priority.size_of() +
    sizeof (_tag                ) +
    sizeof (NOMAD::Signature *  ) +
    sizeof (_current_run        ) +
    sizeof (_in_cache           ) +
    sizeof (_eval_type          ) +
    sizeof (_eval_status        ) +
    sizeof (_EB_ok              ) +
    sizeof (_direction          ) +
    ((_direction     ) ? _direction->size_of() : 0);
}

/*-------------------------------------------------------*/
/*                        scaling                        */
/*-------------------------------------------------------*/
void NOMAD::Eval_Point::scale ( void )
{
    if ( !_signature )
        throw NOMAD::Exception ( "Eval_Point.cpp" , __LINE__ ,
                                "x.Eval_Point::scale(): x has no signature" );
    _signature->scale ( *this );
}

/*-------------------------------------------------------*/
/*                       unscaling                       */
/*-------------------------------------------------------*/
void NOMAD::Eval_Point::unscale ( void )
{
    if ( !_signature )
        throw NOMAD::Exception ( "Eval_Point.cpp" , __LINE__ ,
                                "x.Eval_Point::unscale(): x has no signature" );
    _signature->unscale ( *this );
}

/*-------------------------------------------------------*/
/*                     snap to bounds                    */
/*       returns true if the point has been modified     */
/*-------------------------------------------------------*/
bool NOMAD::Eval_Point::snap_to_bounds ( void )
{
    if ( !_signature )
        throw NOMAD::Exception ( "Eval_Point.cpp" , __LINE__ ,
                                "x.Eval_Point::snap_to_bounds(): x has no signature" );
    return _signature->snap_to_bounds ( *this , _direction );
}

/*-------------------------------------------------------*/
/*             treat the periodic variables              */
/*       returns true if the point has been modified     */
/*-------------------------------------------------------*/
bool NOMAD::Eval_Point::treat_periodic_variables ( NOMAD::Direction *& new_dir )
{
    if (!_signature)
        throw NOMAD::Exception ( "Eval_Point.cpp" , __LINE__ ,
                                "x.Eval_Point::treat_periodic_variables(): x has no signature" );
    
    return _signature->treat_periodic_variables ( *this , _direction , new_dir );
}

/*--------------------------------------------------*/
/*                 Eval_Point::check                */
/*--------------------------------------------------*/
bool NOMAD::Eval_Point::check ( int m , NOMAD::check_failed_type & cf ) const
{
    if ( size() <= 0 || !_signature || m != _bb_outputs.size() ) {
        std::string err = "Eval_Point::check() could not procede";
        if ( !_signature )
            err += " (no signature)";
        else if ( m != _bb_outputs.size() )
            err += " (wrong number of blackbox outputs)";
        else
            err += " (point size <= 0 !)";
        throw NOMAD::Exception ( "Eval_Point.cpp" , __LINE__ , err );
    }
    
    cf = NOMAD::CHECK_OK;
    
    const std::vector<NOMAD::bb_input_type>
    & input_types = _signature->get_input_types();
    const NOMAD::Point & lb          = _signature->get_lb();
    const NOMAD::Point & ub          = _signature->get_ub();
    const NOMAD::Point & fv          = _signature->get_fixed_variables();
    int                  n           = size();
    NOMAD::bb_input_type iti;
    
    for ( int i = 0 ; i < n ; ++i )
    {
        
        const NOMAD::Double xi = (*this)[i];
        
        // undefined coordinates ?
        if ( !xi.is_defined() )
            throw NOMAD::Exception ( "Eval_Point.cpp" , __LINE__ ,
                                    "Eval_Point::check() could not procede (undefined coordinates)" );
        
        // check the bounds:
        const NOMAD::Double & lbi = lb[i];
        if ( lbi.is_defined() && xi < lbi )
        {
            cf = NOMAD::LB_FAIL;
            return false;
        }
        
        const NOMAD::Double & ubi = ub[i];
        if ( ubi.is_defined() && xi > ubi )
        {
            cf = NOMAD::UB_FAIL;
            return false;
        }
        
        // check the integer/categorical/binary variables:
        iti = input_types[i];
        if ( iti == NOMAD::BINARY && !xi.is_binary() ) {
            cf = NOMAD::BIN_FAIL;
            return false;
        }
        if ( ( iti == NOMAD::INTEGER || iti == NOMAD::CATEGORICAL )
            && !xi.is_integer() ) {
            cf = ( iti == NOMAD::INTEGER ) ? NOMAD::INT_FAIL : NOMAD::CAT_FAIL;
            return false;
        }
        
        // check the fixed-variables:
        const NOMAD::Double & fvi = fv[i];
        if ( fvi.is_defined() && fvi != xi ) {
            cf = NOMAD::FIX_VAR_FAIL;
            return false;
        }
    }
    return true;
}

/*--------------------------------------------------*/
/*                  display the tag                 */
/*--------------------------------------------------*/
void NOMAD::Eval_Point::display_tag ( const NOMAD::Display & out ) const
{
    out << "#";
    out.display_int_w ( _tag , NOMAD::Eval_Point::_current_tag );
}

/*--------------------------------------------------*/
/*                      display                     */
/*--------------------------------------------------*/
void NOMAD::Eval_Point::display_eval( const NOMAD::Display & out , bool in_block ) const
{
    if ( in_block ) {
        
        std::ostringstream oss;
        oss << "#" << _tag;
        out << NOMAD::open_block ( oss.str() )
        << "x    = ( ";
        NOMAD::Point::display ( out , " " , 2  , NOMAD::Point::get_display_limit() );
        out << " )" << std::endl
        << "F(x) = [ ";
        _bb_outputs.display ( out , " " , 2 , NOMAD::Point::get_display_limit() );
        out << " ]" << std::endl;
        if ( _h.is_defined() )
            out << "h    = " << _h << std::endl;
        if ( _f.is_defined() )
            out << "f    = " << _f << std::endl;
        out.close_block();
    }
    else {
        display_tag ( out );
        out << " x=( ";
        NOMAD::Point::display ( out , " " , 2  , NOMAD::Point::get_display_limit() );
        out << " ) F(x)=[ ";
        _bb_outputs.display ( out , " " , 2 , NOMAD::Point::get_display_limit() );
        out << " ]";
        if ( _h.is_defined() )
            out << " h=" << _h;
        if ( _f.is_defined() )
            out << " f=" << _f;
    }
}

/*--------------------------------------------------------------*/
/*  comparison operator '<': used to find and store the points  */
/*                           in the filter                      */
/*--------------------------------------------------------------*/
bool NOMAD::Eval_Point::operator < ( const NOMAD::Eval_Point & x ) const
{
    if ( this == &x || !is_eval_ok() || !_EB_ok )
        return false;
    
    double h  = _h.value();
    double f  = _f.value();
    double hx = x._h.value();
    double fx = x._f.value();
    
    if ( h < hx )
        return ( f <= fx );
    
    if ( h == hx )
        return ( f < fx );
    
    return false;
}

/*--------------------------------------------------------------*/
/*        check if there are nan's in the blackbox outputs      */
/*--------------------------------------------------------------*/
bool NOMAD::Eval_Point::check_nan ( void ) const
{
    int m = _bb_outputs.size();
    for ( int i = 0 ; i < m ; ++i ) {
        if ( _bb_outputs[i].is_defined() ) {
#ifdef WINDOWS
            if ( isnan ( _bb_outputs[i].value() ) )
                return true;
#else
            if ( std::isnan ( _bb_outputs[i].value() ) )
                return true;
#endif
        }
    }
    return false;
}

#ifdef MODEL_STATS

/*--------------------------------------------------------------*/
/*             set model data when debugging models             */
/*--------------------------------------------------------------*/
void NOMAD::Eval_Point::set_model_data ( const NOMAD::Eval_Point & x ) const
{
    _mod_use = x._mod_use;
    _nY      = x._nY;
    _cond    = x._cond;
    _Yw      = x._Yw;
    _mh      = x._mh;
    _mf      = x._mf;
}

/*--------------------------------------------------------------*/
/*                       clear model data                       */
/*--------------------------------------------------------------*/
void NOMAD::Eval_Point::clear_model_data ( void ) const
{
    _mod_use = -1;
    _nY      = -1;
    _cond.clear();
    _Yw.clear();
    _mh.clear();
    _mf.clear();
}

#endif
