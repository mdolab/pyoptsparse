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
 \file   Quad_Model.cpp
 \brief  Quadratic regression or MFN interpolation model (implementation)
 \author Sebastien Le Digabel
 \date   2010-08-31
 \see    Quad_Model.hpp
 */
#include "Quad_Model.hpp"

/*-----------------------------------------------------------*/
/*                         constructor                       */
/*-----------------------------------------------------------*/
NOMAD::Quad_Model::Quad_Model
( const NOMAD::Display                     & out       ,
 const std::vector<NOMAD::bb_output_type> & bbot      ,
 const NOMAD::Cache                       & cache     ,
 const NOMAD::Signature                   & signature   )
: _out                  ( out                                 ) ,
_bbot                 ( bbot                                ) ,
_interpolation_type   ( NOMAD::UNDEFINED_INTERPOLATION_TYPE ) ,
_n                    ( signature.get_n()                   ) ,
_nfree                ( _n                                  ) ,
_fixed_vars           ( new bool [_n]                       ) ,
_index                ( NULL                                ) ,
_alpha                ( NULL                                ) ,
_cache                ( cache                               ) ,
_signature            ( signature                           ) ,
_error_flag           ( true                                )
{
    for ( int i = 0 ; i < _n ; ++i )
        _fixed_vars[i] = false;
    init_alpha();
}

/*-----------------------------------------------------------*/
/*                           destructor                      */
/*-----------------------------------------------------------*/
NOMAD::Quad_Model::~Quad_Model ( void )
{
    int m = static_cast<int> ( _bbot.size() );
    for ( int i = 0 ; i < m ; ++i )
        delete _alpha[i];
    delete [] _alpha;
    delete [] _fixed_vars;
    delete [] _index;
    
    // clear the interpolation set Y:
    for ( size_t k = 0 ; k < _Y.size() ; ++k )
        delete _Y[k];
}

/*-----------------------------------------------------------*/
/*      initialize alpha, the model parameters (private)     */
/*-----------------------------------------------------------*/
void NOMAD::Quad_Model::init_alpha ( void )
{
    _n_alpha = ( _nfree + 1 ) * ( _nfree + 2 ) / 2;
    
    int i , m = static_cast<int> ( _bbot.size() );
    
    // initialize _alpha:
    // ------------------
    if ( _alpha )
    {
        for ( i = 0 ; i < m ; ++i )
            delete _alpha[i];
        delete [] _alpha;
    }
    
    _alpha = new NOMAD::Point * [m];
    
    for ( i = 0 ; i < m ; ++i )
        _alpha[i] = ( _bbot[i] == NOMAD::OBJ || NOMAD::bbot_is_constraint(_bbot[i]) ) ?
        new NOMAD::Point ( _n_alpha ) : NULL;
    
    // initialize _index:
    // ------------------
    
    // example: with 3 variables (X,Y,Z) with Y fixed.
    // --------
    //   the problem is reduced to the two variables x=X and y=Z,
    //   and _index is corresponds to:
    //
    // 0     1      0     1     : index[0] = 0
    // 1     X      1     x     : index[1] = 1
    // 2     Y      2     y     : index[2] = 3
    // 3     Z      3  .5 x^2   : index[3] = 4
    // 4  .5 X^2    4  .5 y^2   : index[4] = 6
    // 5  .5 Y^2    5     xy    : index[5] = 8
    // 6  .5 Z^2
    // 7     XY
    // 8     XZ
    // 9     YZ
    //
    // If there are no fixed variables, index is of size (n+1)(n+2)/2
    //   with index[i] = i.
    
    delete [] _index;
    
    _index = new int [_n_alpha];
    
    int nm1 = _n - 1;
    int c1  = 2*_n + 1;
    int c2  = 1;
    int k1 , k2;
    
    _index[0] = 0;
    for ( i = 0 ; i < _n ; ++i )
    {
        if ( !_fixed_vars[i] )
        {
            _index[c2       ] = i+1;
            _index[c2+_nfree] = i+1+_n;
            ++c2;
        }
    }
    
    c2 += _nfree;
    
    for ( k1 = 0 ; k1 < nm1 ; ++k1 )
        for ( k2 = k1+1 ; k2 < _n ; ++k2 )
        {
            if ( !_fixed_vars[k1] && !_fixed_vars[k2] )
                _index[c2++] = c1;
            ++c1;
        }
}

/*---------------------------------------------------------*/
/*  check evaluation point outputs before the integration  */
/*  into an interpolation set (private)                    */
/*---------------------------------------------------------*/
bool NOMAD::Quad_Model::check_outputs ( const NOMAD::Point & bbo , int m ) const {
    
    if ( bbo.size() != m )
        return false;
    
    for ( int i = 0 ; i < m ; ++i )
        if ( !bbo[i].is_defined() || bbo[i].value() > NOMAD::MODEL_MAX_OUTPUT )
            return false;
    
    return true;
}

/*-----------------------------------------------------------*/
/*               construct the interpolation set Y           */
/*-----------------------------------------------------------*/
void NOMAD::Quad_Model::construct_Y ( const NOMAD::Point & center               ,
                                     const NOMAD::Point & interpolation_radius ,
                                     int                  max_Y_size             )
{
    _error_flag = true;
    
    if ( center.size()               != _n   ||
        interpolation_radius.size() != _n   ||
        !center.is_complete()               ||
        !interpolation_radius.is_complete()    )
        return;
    
    _error_flag = false;
    _center     = center;
    
    int m = static_cast<int> ( _bbot.size() );
    
    // browse the cache:
    const NOMAD::Eval_Point * cur = _cache.begin();
    while ( cur )
    {
        
        if ( cur->get_eval_status() == NOMAD::EVAL_OK &&
            cur->get_n          () == _n             &&
            _signature.is_compatible (*cur)             )
        {
            
            const NOMAD::Point & bbo = cur->get_bb_outputs();
            
            if ( check_outputs ( bbo , m ) )
            {
                
                // the center point has been found
                // (it is put in first position):
                if ( _center == *cur )
                {
                    _Y.push_back ( new NOMAD::Eval_Point ( *cur ) );
                    int nYm1 = get_nY() - 1;
                    if ( nYm1 > 0 )
                    {
                        NOMAD::Eval_Point * tmp = _Y[0];
                        _Y[0   ] = _Y[nYm1];
                        _Y[nYm1] = tmp;
                    }
                }
                
                // other points must within the interpolation radius:
                else if ( is_within_radius ( *cur , interpolation_radius ) )
                {
                    _Y.push_back ( new NOMAD::Eval_Point ( *cur ) );
                }
            }
        }
        cur = _cache.next();
    }
    
    // respect the limit on the number of points:
    if ( get_nY() > max_Y_size )
        reduce_Y ( center , max_Y_size );
}

/*-----------------------------------------------------------------*/
/*             reduce the number of interpolation points           */
/*-----------------------------------------------------------------*/
void NOMAD::Quad_Model::reduce_Y ( const NOMAD::Point & center     ,
                                  int                  max_Y_size   )
{
    int nY = get_nY();
    
    if ( nY <= max_Y_size )
        return;
    
    std::multiset<NOMAD::Model_Sorted_Point> Ys;
    for ( int k = 0 ; k < nY ; ++k )
        Ys.insert ( NOMAD::Model_Sorted_Point ( _Y[k] , center ) );
    
    _Y.clear();
    
    std::multiset<NOMAD::Model_Sorted_Point>::const_iterator it , end = Ys.end();
    for ( it = Ys.begin() ; it != end ; ++it ) {
        if ( get_nY() < max_Y_size )
            _Y.push_back ( static_cast<NOMAD::Eval_Point *> ( it->get_point() ) );
        else
            delete it->get_point();
    }
}

/*-----------------------------------------------------------*/
/*  check if an unscaled point is in B(center,radius) for a  */
/*  given radius (private)                                   */
/*-----------------------------------------------------------*/
bool NOMAD::Quad_Model::is_within_radius ( const NOMAD::Point & x      ,
                                          const NOMAD::Point & radius   ) const
{
    if ( x.size() != _n || radius.size() != _n )
        return false;
    
    for ( int i = 0 ; i < _n ; ++i )
    {
        if ( !x[i].is_defined()                     ||
            !radius[i].is_defined()                ||
            radius[i] < ( x[i] - _center[i]).abs()    )
            return false;
    }
    return true;
}

/*------------------------------------------------------*/
/*  check if a scaled point is inside the trust radius  */
/*------------------------------------------------------*/
bool NOMAD::Quad_Model::is_within_trust_radius ( const NOMAD::Point & x ) const
{
    // check that all scaled coordinates are in [-1;1] and
    // that fixed variables are equal to zero:
    for ( int i = 0 ; i < _n ; ++i )
        if ( !_ref    [i].is_defined()         ||
            !_scaling[i].is_defined()         ||
            !       x[i].is_defined()         ||
            x[i].abs() > 1.0                  ||
            ( _fixed_vars[i] && x[i] != 0.0 )    )
            return false;
    return true;
}

/*--------------------------------------------------------------*/
/*  . define scaling to put all coordinates centered in [-r;r]  */
/*  . looks also for fixed variables                            */
/*--------------------------------------------------------------*/
void NOMAD::Quad_Model::define_scaling ( const NOMAD::Double & r )
{
    if ( _error_flag || _Y.empty() ) {
        _error_flag = true;
        return;
    }
    
    int           i , k;
    int           nY = get_nY();
    NOMAD::Point  min(_n) , max(_n);
    NOMAD::Double tmp;
    
    // The parameters defining the scaling with rotation (see define_scaling_by_direction) are cleared.
    // Only the parameters for the basic scaling are set
    _dirP.clear();
    _epsilon.clear();
    _delta_m.clear();
    
    
    _scaling.clear();
    _ref.clear    ();
    _ref.reset     ( _n );
    _scaling.reset ( _n );
    
    // compute the reference (center of Y):
    for ( k = 0 ; k < nY ; ++k )
    {
        
        if ( !_Y[k] || _n != _Y[k]->size() )
        {
            _error_flag = true;
            return;
        }
        
        for ( i = 0 ; i < _n ; ++i )
        {
            tmp = (*_Y[k])[i];
            if ( !min[i].is_defined() || tmp < min[i] )
                min[i] = tmp;
            if ( !max[i].is_defined() || tmp > max[i] )
                max[i] = tmp;
        }
    }
    
    for ( i = 0 ; i < _n ; ++i )
        _ref[i] = ( min[i] + max[i] ) / 2.0;
    
#ifdef MODEL_STATS
    _Yw = NOMAD::Double();
    for ( i = 0 ; i < _n ; ++i )
    {
        tmp = max[i]-min[i];
        if ( !_Yw.is_defined() || tmp > _Yw )
            _Yw = tmp;
    }
#endif
    
#ifdef DEBUG
    _out << std::endl
    << "define_scaling(): reference = ( " << _ref << " )" << std::endl;
#endif
    
    // compute the scaling (and detect fixed variables):
    for ( k = 0 ; k < nY ; ++k ) {
        
        for ( i = 0 ; i < _n ; ++i ) {
            tmp = ( (*_Y[k])[i] - _ref[i] ).abs();
            if ( !_scaling[i].is_defined() || _scaling[i] < tmp )
                _scaling[i] = tmp;
        }
    }
    
    _nfree = _n;
    
    for ( i = 0 ; i < _n ; ++i )
    {
        if ( _scaling[i] == 0.0 )
        {
            _scaling   [i] = 0.0;
            _fixed_vars[i] = true;
            --_nfree;
            if ( _nfree == 0 )
            {
                _scaling.clear();
                _ref.clear();
                _error_flag = true;
                return;
            }
        }
        else
            _scaling[i] *= 1.0/r; // all coordinates in [-r;r]
    }
    
    if ( _nfree < _n )
        init_alpha();
    
    for ( k = 0 ; k < nY ; ++k )
    {
        if ( !scale ( *_Y[k] ) )
        {
            _scaling.clear();
            _error_flag = true;
            return;
        }
    }
    
#ifdef DEBUG
    _out << "define_scaling(): scaling   = ( " << _scaling << " )" << std::endl;
#endif
    
    _error_flag = false;
}

/*-------------------------------------------------------------------*/
/*  . Scaling with rotation based on a set of directions.             */
/*     See paper:                                                    */
/*       Reducing the number of function evaluations in              */
/*       Mesh Adaptive Direct Search algorithms, Audet, Ianni,       */
/*       LeDigabel, Tribes, 2014                                     */
/*  . looks also for fixed variables                                 */
/*-------------------------------------------------------------------*/
void NOMAD::Quad_Model::define_scaling_by_directions ( const std::list<NOMAD::Direction> & dirs, const NOMAD::Point & delta_m, const NOMAD::Double & epsilon  )
{
    if ( _error_flag || _Y.empty() ) {
        _error_flag = true;
        return;
    }
    
    int           i , k;
    int           nY = get_nY();
    NOMAD::Point  min(_n) , max(_n);
    NOMAD::Double tmp;
    
    
    // The parameters defining the basic scaling (see define_scaling) are cleared.
    // Only the parameters for the direction scaling are set
    _scaling.clear();
    _ref.clear    ();
    
    // For direction scaling
    if (static_cast<int> (dirs.size())!=_n || static_cast<int>(delta_m.size()) != _n || epsilon<=0.0 || epsilon>=1)
    {
        _error_flag = true;
        return;
    }
    _delta_m=delta_m;
    // Get D' from dirs (scaling with delta_m
    std::list<NOMAD::Direction>::const_iterator itDir;
    for (itDir=dirs.begin(); itDir != dirs.end(); itDir++)
    {
        NOMAD::Direction dir_i(_n,0.0,itDir->get_type());
        dir_i.set_index(itDir->get_index());
        for ( int i = 0 ; i < _n ; ++i )
        {
            if (_delta_m[i]<=0.0)
            {
                _error_flag=true;
                return;
            }
            dir_i[i]=(*itDir)[i]/_delta_m[i];
        }
        _dirP.push_back(dir_i);
    }
    
    _epsilon=epsilon;
    
    
    // compute the min and the max:
    for ( k = 0 ; k < nY ; ++k )
    {
        
        if ( !_Y[k] || _n != _Y[k]->size() )
        {
            _error_flag = true;
            return;
        }
        
        for ( i = 0 ; i < _n ; ++i )
        {
            tmp = (*_Y[k])[i];
            if ( !min[i].is_defined() || tmp < min[i] )
                min[i] = tmp;
            if ( !max[i].is_defined() || tmp > max[i] )
                max[i] = tmp;
        }
    }
    
#ifdef MODEL_STATS
    _Yw = NOMAD::Double();
    for ( i = 0 ; i < _n ; ++i )
    {
        tmp = max[i]-min[i];
        if ( !_Yw.is_defined() || tmp > _Yw )
            _Yw = tmp;
    }
#endif
    
    // Detect fixed variables:
    _nfree = _n;
    for ( i = 0 ; i < _n ; ++i )
    {
        bool fixed_var_i=true;
        for ( k = 0 ; k < nY ; ++k )
        {
            if ( ( (*_Y[k])[i] - _center[i] ).abs() > 0.0 )
            {
                fixed_var_i=false;
                break;
            }
        }
        _fixed_vars[i]=fixed_var_i;
        if (fixed_var_i)
            --_nfree;
        
        if ( _nfree == 0 )
        {
            _scaling.clear();
            _ref.clear();
            _dirP.clear();
            _error_flag = true;
            return;
        }
    }
    if ( _nfree < _n )
        init_alpha();
    
    // Perform scaling of Y
    for ( k = 0 ; k < nY ; ++k )
    {
        if ( !scale ( *_Y[k] ) )
        {
            _scaling.clear();
            _dirP.clear();
            _error_flag = true;
            return;
        }
    }
    
#ifdef DEBUG
    _out << "define_scaling_by_direction(): " << std::endl;
    for ( itDir = _dirP.begin() ; itDir != _dirP.end() ; ++itDir )
    {
        _out << "dirPrime ";
        _out.display_int_w ( (*itDir).get_index() , static_cast<int>(_dirP.size()) );
        _out << " : " << *itDir << std::endl;
    }
#endif
    
    _error_flag = false;
}


/*--------------------------------------------------------------*/
/*                          scale a point                       */
/*--------------------------------------------------------------*/
bool NOMAD::Quad_Model::scale ( NOMAD::Point & x ) const
{
    if ( _error_flag || _n != x.size() )
        return false;
    
    if (_dirP.size()==0)
    {
        // Scale without rotation
        for ( int i = 0 ; i < _n ; ++i )
        {
            if ( !_ref [i].is_defined() ||
                !_scaling[i].is_defined()  ||
                !       x[i].is_defined()    )
                return false;
            x[i] -= _ref[i];
            if ( _scaling[i] != 0 )
                x[i] /= _scaling[i];
        }
    }
    else
    {
        if (! _epsilon.is_defined() || !_delta_m.is_complete())
            return false;
        // Scale with rotation based on direction and center (see paper Reducing the number of function evaluations in Mesh Adaptive Direct Search algorithms, Audet, Ianni, LeDigabel, Tribes, 2014
        // T(y)=(D')^-1*(center-x)/delta_m/(1-epsilon) - epsilon*1/(1-epsilon)
        // (D')^-1=(D')^T/normCol^2
        NOMAD::Point temp(_n,0.0);
        NOMAD::Double normCol2=0.0;
        std::list<NOMAD::Direction>::const_iterator itDir=_dirP.begin();
        for ( int i = 0 ; i < _n ; ++i )
        {
            normCol2+=pow((*itDir)[i].value(),2.0);
            
            if (_delta_m[i] !=0.0)
                temp[i]=(_center[i].value()-x[i].value())/_delta_m[i].value()/(1.0-_epsilon.value());
            else
                return false;
            x[i]=0.0;
        }
        int j=0;
        for (itDir=_dirP.begin(); itDir != _dirP.end(); itDir++,j++)
        {
            for ( int i = 0 ; i < _n ; ++i )
            {
                x[j]+=temp[i].value()*(*itDir)[i].value()/normCol2.value();
            }
            x[j]-=_epsilon.value()/(1.0-_epsilon.value());
        }
    }
    
    return true;
}

/*-----------------------------------------------------------*/
/*                       unscale a point                     */
/*-----------------------------------------------------------*/
bool NOMAD::Quad_Model::unscale ( NOMAD::Point & x ) const
{
    if ( _error_flag || _n != x.size() )
        return false;
    
    if (_dirP.size()==0)
    {
        // Scale without rotation
        for ( int i = 0 ; i < _n ; ++i )
        {
            if ( !_ref    [i].is_defined() ||
                !_scaling[i].is_defined() ||
                !       x[i].is_defined()    )
                return false;
            
            x[i] *= _scaling[i];
            x[i] += _ref    [i];
        }
    }
    else
    {
        
        if (! _epsilon.is_defined() || !_delta_m.is_complete())
            return false;
        
        // UnScale with rotation see paper Reducing the number of function evaluations in Mesh Adaptive Direct Search algorithms, Audet, Ianni, LeDigabel, Tribes, 2014
        //T^−1(x) = center+ _delta_m Dp ((ε−1)x−ε1)
        NOMAD::Point temp(_n,0.0);
        for ( int i = 0 ; i < _n ; ++i )
        {
            temp[i]=(x[i]*(_epsilon-1.0)-_epsilon)*_delta_m[i];
            x[i]=0.0;
        }
        std::list<NOMAD::Direction>::const_iterator itDir;
        int j=0;
        for (itDir=_dirP.begin(); itDir != _dirP.end(); itDir++,j++)
        {
            for (int i=0 ; i< _n ; i++)
            {
                x[i]+=temp[j]*(*itDir)[i];
            }
        }
        for ( int i = 0 ; i < _n ; ++i )
        {
            x[i]+=_center[i];
        }
    }
    
    return true;
}

/*-----------------------------------------------------------*/
/*                       unscale the slope at a point                     */
/*-----------------------------------------------------------*/
bool NOMAD::Quad_Model::unscale_grad ( NOMAD::Point & x ) const
{
    if ( _error_flag || _n != x.size() )
        return false;
    
    for ( int i = 0 ; i < _n ; ++i )
    {
        
        if (!_scaling[i].is_defined() || !x[i].is_defined()    )
            return false;
        
        x[i] *= _scaling[i];
    }
    
    return true;
}

/*------------------------------------------------------------------*/
/*  compute the element (i,j) of the interpolation matrix M(phi,Y)  */
/*  (private)                                                       */
/*------------------------------------------------------------------*/
double NOMAD::Quad_Model::compute_M ( int i , int j ) const {
    
    if ( _error_flag )
        return 0.0;
    
    if ( j == 0 )
        return 1.0;
    
    if ( j <= _nfree )
        return (*_Y[i])[_index[j]-1].value();
    
    if ( j <= 2 * _nfree )
        return 0.5 * pow ( (*_Y[i])[_index[j-_nfree]-1].value() , 2.0 );
    
    int nm1 = _nfree - 1;
    int jm2n , dec , r , i1 , i2;
    
    jm2n = j - 2 * _nfree;
    dec  = nm1;
    r    = jm2n;
    i1   = -1;
    
    while ( r > 0 ) {
        r -= dec;
        ++i1;
        --dec;
    }
    
    i2 = r + nm1;
    
    return (*_Y[i])[_index[i1+1]-1].value() * (*_Y[i])[_index[i2+1]-1].value();
}

/*-----------------------------------------------------------*/
/*               construct m models (one by output)          */
/*-----------------------------------------------------------*/
void NOMAD::Quad_Model::construct ( bool   use_WP     ,
                                   double eps        ,
                                   int    max_mpn    ,
                                   int    max_Y_size   )
{
    if ( _error_flag )
        return;
    
    int p1 = get_nY();
    
    
    // MFN interpolation:
    if ( p1 < _n_alpha ) {
        _interpolation_type = NOMAD::MFN;
        _error_flag = !construct_MFN_model ( eps , max_mpn , max_Y_size );
    }
    else {
        
        _error_flag = true;
        
        // well-poised regression:
        if ( use_WP && p1 > _n_alpha ) {
            _interpolation_type = NOMAD::WP_REGRESSION;
            _error_flag = !construct_WP_model ( max_Y_size );
        }
        
        // regression:
        if ( _error_flag ) {
            _interpolation_type = NOMAD::REGRESSION;
            _error_flag = !construct_regression_model ( eps , max_mpn , max_Y_size );
        }
    }
}

/*---------------------------------------------------------------*/
/*  find interpolation point with max Lagrange polynomial value  */
/*---------------------------------------------------------------*/
/*    . ji = argmax |li(x)| for x in Y                           */
/*    . used in construct_WP_model()                             */
/*    . private                                                  */
/*---------------------------------------------------------------*/
int NOMAD::Quad_Model::find_max_lix
( const NOMAD::Point                     & li      ,
 const std::vector<NOMAD::Eval_Point *> & Y       ,
 int                                      i1      ,
 int                                      i2      ,
 NOMAD::Double                          & max_lix   ) const
{
    max_lix = -1.0;
    int  ji = -1;
    NOMAD::Double tmp;
    for ( int j = i1 ; j <= i2 ; ++j ) {
        tmp = eval ( *Y[j] , li );
        if ( tmp.is_defined() ) {
            tmp = tmp.abs();
            if ( tmp > max_lix ) {
                max_lix = tmp;
                ji      = j;
            }
        }
    }
    if ( ji < 0 )
        max_lix.clear();
    return ji;
}

/*-----------------------------------------------------------*/
/*          construct well-poised (WP) model (private)       */
/*-----------------------------------------------------------*/
bool NOMAD::Quad_Model::construct_WP_model ( int max_Y_size )
{
    
#ifdef DEBUG
    _out << std::endl
    << NOMAD::open_block ( "NOMAD::Quad_Model::construct_WP_model()" );
#endif
    
    // check the set Y:
    if ( !check_Y() )
        return false;
    
    int i , j , k , p1 = get_nY();
    
    // the number of points (p+1) must be in [1+(n+1)(n+2)/2;MS_MAX_Y_SIZE]:
    if ( p1 <= _n_alpha || p1 > max_Y_size ) {
#ifdef DEBUG
        _out << std::endl
        << "NOMAD::Quad_Model::construct_WP_model(): "
        << "(p+1) not in [1+(n+1)(n+2)/2;" << max_Y_size << "]"
        << std::endl << NOMAD::close_block() << std::endl;
#endif
        return false;
    }
    
    // Lagrange polynomials:
    std::vector<NOMAD::Point *> l;
    for ( i = 0 ; i < _n_alpha ; ++i ) {
        l.push_back ( new NOMAD::Point ( _n_alpha ) );
        for ( j = 0 ; j < _n_alpha ; ++j )
            (*l[i])[j] = (i==j) ? 1.0 : 0.0;
    }
    
    // creation of sets Y1 and Y2; Y2 contains all available points
    // of _Y and Y1 will be the 'well-poised' set with n_alpha points:
    std::vector<NOMAD::Eval_Point *> Y1 , Y2 = _Y;
    int                              iy2 , ny2m1 = p1-1;
    NOMAD::Double                    max_lix , liyi , ljyi;
    
    // we init Y1 with the first point of Y:
    Y1.push_back ( Y2[0] );
    Y2[0] = Y2[ny2m1];
    Y2.resize ( ny2m1 );
    
    // use algo 6.2 p.95 of the DFO book in order to construct Lagrange polynomials:
    // -----------------------------------------------------------------------------
    for ( i = 0 ; i < _n_alpha ; ++i ) {
        
        // 1. point selection (select a point in Y2: Y2[iy2]):
        // -------------------
        if ( i > 0 ) {
            
            ny2m1 = static_cast<int>(Y2.size())-1;
            iy2   = find_max_lix ( *l[i] , Y2 , 0 , ny2m1 , max_lix );
            
            if ( iy2 < 0 ) {
#ifdef DEBUG
                _out << std::endl
                << "NOMAD::Quad_Model::construct_WP_model(): "
                << "cannot find candidate in Y"
                << std::endl << NOMAD::close_block() << std::endl;
#endif
                for ( i = 0 ; i < _n_alpha ; ++i )
                    delete l[i];
                return false;
            }
            
            // add Y2[iy2] in Y1:
            Y1.push_back ( Y2[iy2] );
            Y2[iy2] = Y2[ny2m1];
            Y2.resize (ny2m1);
        }
        
        // 2. normalization:
        // -----------------
        liyi = eval ( *Y1[i] , *l[i] );
        
        if ( liyi.abs().value() < 1e-15 ) {
#ifdef DEBUG
            _out << std::endl
            << "NOMAD::Quad_Model::construct_WP_model(): set Y is not poised"
            << std::endl << NOMAD::close_block() << std::endl;
#endif
            for ( i = 0 ; i < _n_alpha ; ++i )
                delete l[i];
            return false;
        }
        
        for ( k = 0 ; k < _n_alpha ; ++k ) {
            (*l[i])[k] /= liyi;
            if ( (*l[i])[k].abs().value() < 1e-15 )
                (*l[i])[k] = 0.0;
        }
        
        // 3. orthogonalization:
        // ---------------------
        for ( j = 0 ; j < _n_alpha ; ++j )
            if ( j != i ) {
                ljyi = eval ( *Y1[i] , *l[j] );
                for ( k = 0 ; k < _n_alpha ; ++k ) {
                    (*l[j])[k] = (*l[j])[k] - ljyi * (*l[i])[k];
                    if ( (*l[j])[k].abs().value() < 1e-15 )
                        (*l[j])[k] = 0.0;
                }
            }
    }
    
#ifdef DEBUG
    display_lagrange_polynomials ( l , Y1 );
#endif
    
    // compute alpha:
    // --------------
    int m = static_cast<int> ( _bbot.size() );
    for ( i = 0 ; i < m ; ++i )
        if ( _alpha[i] ) {
            for ( j = 0 ; j < _n_alpha ; ++j ) {
                (*_alpha[i])[j] = 0.0;
                for ( k = 0 ; k < _n_alpha ; ++k )
                    (*_alpha[i])[j] += Y1[k]->get_bb_outputs()[i] * (*l[k])[j];
            }
        }
    
    // poisedness improvement using algorithm 6.3 page 95:
    // ---------------------------------------------------
    
    // old alpha:
    NOMAD::Point ** old_alpha = new NOMAD::Point * [m] , ** tmp_alpha;
    for ( i = 0 ; i < m ; ++i )
        old_alpha[i] = ( _alpha[i] ) ?
        new NOMAD::Point ( _n_alpha ) : NULL;
    
    int           ik;
    NOMAD::Double ljyk , lkyk , lix , new_rel_err ,
    cur_rel_err = compute_max_rel_err();
    
    if ( cur_rel_err.is_defined() && cur_rel_err.value() > 1e-15 ) {
        
        for ( int niter = 0 ; niter < 10 ; ++niter ) {
            
            ny2m1 = static_cast<int>(Y2.size())-1;
            
            if ( ny2m1 < 0 )
                break;
            
            max_lix = -1.0;
            iy2     = -1;
            ik      = -1;
            
            for ( i = 0 ; i < _n_alpha ; ++i ) {
                
                j = find_max_lix ( *l[i] , Y2 , 0 , ny2m1 , lix );
                if ( j >= 0 && lix > max_lix ) {
                    max_lix = lix;
                    iy2     = j;
                    ik      = i;
                }
            }
            
            if ( ik < 0 )
                break;
            
            // set Y1[ik] = Y2[iy2]:
            Y1[ik ] = Y2[iy2];
            Y2[iy2] = Y2[ny2m1];
            Y2.resize ( ny2m1 );
            
            lkyk = eval ( *Y1[ik] , *l[ik] );
            
            if ( lkyk.abs() <= 1e-15 )
                break;
            
            // update Lagrange polynomials:
            // ----------------------------
            
            // normalization and orthogonalization:
            for ( i = 0 ; i < _n_alpha ; ++i )
                (*l[ik])[i] /= lkyk;
            
            for ( j = 0 ; j < _n_alpha ; ++j ) {
                if ( j != ik ) {
                    ljyk = eval ( *Y1[ik] , *l[j] );
                    for ( i = 0 ; i < _n_alpha ; ++i )
                        (*l[j])[i] = (*l[j])[i] - ljyk * (*l[ik])[i];
                }
            }
            
            // save old alpha and compute new one:
            for ( i = 0 ; i < m ; ++i )
                if ( _alpha[i] ) {
                    *(old_alpha[i]) = *(_alpha[i]);
                    for ( j = 0 ; j < _n_alpha ; ++j ) {
                        (*_alpha[i])[j] = 0.0;
                        for ( k = 0 ; k < _n_alpha ; ++k )
                            (*_alpha[i])[j] += Y1[k]->get_bb_outputs()[i] * (*l[k])[j];
                    }
                }
            
            // compute new error:
            new_rel_err = compute_max_rel_err();
            
            // if no better error, restore old alpha and exit loop:
            if ( !new_rel_err.is_defined() || new_rel_err >= cur_rel_err ) {
                tmp_alpha = _alpha;
                _alpha    = old_alpha;
                old_alpha = tmp_alpha;
                break;
            }
            
            cur_rel_err = new_rel_err;
        }
    }
    
    for ( i = 0 ; i < m ; ++i )
        delete old_alpha[i];
    delete [] old_alpha;
    for ( i = 0 ; i < _n_alpha ; ++i )
        delete l[i];
    
#ifdef DEBUG
    _out.close_block();
#endif
    
    return true;
}

/*-----------------------------------------------------------*/
/*             construct regression model (private)          */
/*-----------------------------------------------------------*/
bool NOMAD::Quad_Model::construct_regression_model ( double eps        ,
                                                    int    max_mpn    ,
                                                    int    max_Y_size   )
{
#ifdef DEBUG
    _out << std::endl
    << NOMAD::open_block ( "NOMAD::Quad_Model::construct_regression_model()" );
#endif
    
    _error_flag = false;
    
    // check the set Y:
    if ( !check_Y() )
        return false;
    
    int p1 = get_nY();
    
    // the number of points (p+1) must be in [(n+1)(n+2)/2;MS_MAX_Y_SIZE]:
    if ( p1 < _n_alpha || p1 > max_Y_size ) {
#ifdef DEBUG
        _out << std::endl
        << "NOMAD::Quad_Model::construct_regression_model(): "
        << "(p+1) not in [(n+1)(n+2)/2;"
        << max_Y_size << "]"
        << std::endl << NOMAD::close_block() << std::endl;
#endif
        return false;
    }
    
    // for this procedure, the number of points is limited to 500
    // (because of the SVD decomposition):
    if ( p1 > 500 ) {
        reduce_Y ( NOMAD::Point ( _n , 0.0 ) , 500 );
        p1 = 500;
    }
    
    // construct the matrix F=M'M (_n_alpha,_n_alpha):
    // -----------------------------------------------
    int       i , j , k;
    double ** F = new double *[_n_alpha];
    double ** M = new double *[p1];
    for ( i = 0 ; i < p1 ; ++i ) {
        M[i] = new double[_n_alpha];
        for ( j = 0 ; j < _n_alpha ; ++j )
            M[i][j] = compute_M ( i , j );
    }
    
    for ( i = 0 ; i < _n_alpha ; ++i ) {
        F[i] = new double[_n_alpha];
        for ( j = 0 ; j <= i ; ++j ) {
            F[i][j] = 0.0;
            for ( k = 0 ; k < p1 ; ++k )
                F[i][j] += M[k][i] * M[k][j];
            if ( i != j )
                F[j][i] = F[i][j];
        }
    }
    
#ifdef DEBUG
    _out << std::endl << "F=";
    for ( i = 0 ; i < _n_alpha ; ++i ) {
        _out << "\t";
        for ( j = 0 ; j < _n_alpha ; ++j )
            _out << std::setw(12) << F[i][j] << " ";
        _out << std::endl;
    }
#endif
    
    bool error = false;
    
    // SVD decomposition of the F matrix (F=U.W.V'):
    // ---------------------------------------------
    // (F will be transformed in U)
    
    double  * W = new double  [_n_alpha];
    double ** V = new double *[_n_alpha];
    for ( i = 0 ; i < _n_alpha ; ++i )
        V[i] = new double[_n_alpha];
    
    std::string error_msg;
    if ( NOMAD::SVD_decomposition ( error_msg , F , W , V , _n_alpha , _n_alpha , max_mpn ) ) {
        
        // compute condition number:
        compute_cond ( W , _n_alpha , eps );
        
#ifdef DEBUG
        _out << std::endl << "F=";
        for ( i = 0 ; i < _n_alpha ; ++i ) {
            _out << "\t";
            for ( j = 0 ; j < _n_alpha ; ++j )
                _out << std::setw(12) << F[i][j] << " ";
            _out << std::endl;
        }
        
        _out << std::endl << "W=\t";
        for ( i = 0 ; i < _n_alpha ; ++i )
            _out << std::setw(12) << W[i] << " ";
        _out << std::endl << std::endl << "cond=" << _cond << std::endl;
        
        _out << std::endl << "V=";
        for ( i = 0 ; i < _n_alpha ; ++i ) {
            _out << "\t";
            for ( j = 0 ; j < _n_alpha ; ++j )
                _out << std::setw(12) << V[i][j] << " ";
            _out << std::endl;
        }
#endif
        
    }
    else {
#ifdef DEBUG
        _out << std::endl << "NOMAD::Quad_Model::construct_regression_model(): "
        << "SVD decomposition (" << error_msg << ")"
        << std::endl << NOMAD::close_block() << std::endl;
#endif
        error = true;
        _cond.clear();
    }
    
    // resolution of system F.alpha = M'.f(Y):
    // ---------------------------------------
    if ( !error ) {
        int m = static_cast<int> ( _bbot.size() );
        for ( i = 0 ; i < m ; ++i )
            if ( _alpha[i] )
                solve_regression_system ( M , F , W  , V , i , *_alpha[i] , eps );
    }
    
    // free memory:
    for ( i = 0 ; i < _n_alpha ; ++i ) {
        delete [] F[i];
        delete [] V[i];
    }
    for ( i = 0 ; i < p1 ; ++i )
        delete [] M[i];
    delete [] M;
    delete [] F;
    delete [] V;
    delete [] W;
    
#ifdef DEBUG
    _out.close_block();
#endif
    
    return !error;
}

/*-------------------------------------------------------------*/
/*              compute condition number (private)             */
/*-------------------------------------------------------------*/
void NOMAD::Quad_Model::compute_cond ( const double * W , int n , double eps )
{
    double min = NOMAD::INF;
    double max = -min;
    for ( int i = 0 ; i < n ; ++i ) {
        if ( W[i] < min )
            min = W[i];
        if ( W[i] > max )
            max = W[i];
    }
    if ( min < eps )
        min = eps;
    _cond = max / min;
}

/*-------------------------------------------------------------*/
/*  resolution of system F.alpha = M'.f(Y) for the regression  */
/*  (private)                                                  */
/*-------------------------------------------------------------*/
void NOMAD::Quad_Model::solve_regression_system ( double      ** M         ,
                                                 double      ** F         ,
                                                 double       * W         ,
                                                 double      ** V         ,
                                                 int            bbo_index ,
                                                 NOMAD::Point & alpha     ,
                                                 double         eps	      ) const
{
    // resize the alpha vector:
    if ( alpha.size() != _n_alpha )
        alpha.reset ( _n_alpha , 0.0 );
    
    double * alpha_tmp  = new double [_n_alpha];
    int      i , k , p1 = get_nY();
    
    // solve the system:
    for ( i = 0 ; i < _n_alpha ; ++i ) {
        alpha_tmp[i] = 0.0;
        for ( k = 0 ; k < p1 ; ++k )
            alpha_tmp[i] += M[k][i] * ( _Y[k]->get_bb_outputs()[bbo_index].value() );
    }
    
    double * alpha_tmp2 = new double [_n_alpha];
    
    // some W values will be zero (or near zero);
    // each value that is smaller than eps is ignored
    
    for ( i = 0 ; i < _n_alpha ; ++i ) {
        alpha_tmp2[i] = 0.0;
        for ( k = 0 ; k < _n_alpha ; ++k )
            if ( W[i] > eps )
                alpha_tmp2[i] += F[k][i] * alpha_tmp[k] / W[i];
    }
    
    delete [] alpha_tmp;
    
    for ( i = 0 ; i < _n_alpha ; ++i ) {
        alpha[i] = 0.0;
        for ( k = 0 ; k < _n_alpha ; ++k )
            alpha[i] += V[i][k] * alpha_tmp2[k];
    }
    
    delete [] alpha_tmp2;
}

/*----------------------------------------------------------*/
/*  construct Minimum Frobenius Norm (MFN) model (private)  */
/*----------------------------------------------------------*/
bool NOMAD::Quad_Model::construct_MFN_model ( double eps        ,
                                             int    max_mpn    ,
                                             int    max_Y_size   )
{
#ifdef DEBUG
    _out << std::endl
    << NOMAD::open_block ( "NOMAD::Quad_Model::construct_MFN_model()" );
#endif
    
    // check the set Y:
    if ( !check_Y() )
        return false;
    
    int p1 = get_nY();
    
    // the number of points (p+1) must be in [n+1;(n+1)(n+2)/2-1]:
    if ( p1 <= _nfree || p1 >= _n_alpha ) {
#ifdef DEBUG
        _out << std::endl
        << "NOMAD::Quad_Model::construct_MFN_model(): "
        << "(p+1) not in [n+1;(n+1)(n+2)/2-1]"
        << std::endl << NOMAD::close_block() << std::endl;
#endif
        return false;
    }
    
    // for this procedure, the number of points is limited to 250
    // (because of the SVD decomposition):
    if ( p1 > 250 ) {
        reduce_Y ( NOMAD::Point ( _n , 0.0 ) , 250 );
        p1 = 250;
    }
    
    // construct the matrix F (4 parts):
    // ---------------------------------
    // [ 1 | 2 ]
    // [ --+-- ]
    // [ 3 | 4 ]
    
    int       i , j , k;
    int      np1  = _nfree + 1;
    int       nF  = np1 + p1;
    double ** F   = new double *[nF];
    double ** M   = new double *[p1];
    for ( i = 0 ; i < nF ; ++i )
        F[i] = new double[nF];
    
    // 1/4: MQ.MQ' (p+1,p+1):
    {
        for ( i = 0 ; i < p1 ; ++i )
        {
            M[i] = new double[_n_alpha];
            for ( j = 0 ; j < _n_alpha ; ++j )
                M[i][j] = compute_M ( i , j );
            for ( j = 0 ; j <= i ; ++j )
            {
                F[i][j] = 0.0;
                for ( k = np1 ; k < _n_alpha ; ++k )
                    F[i][j] += M[i][k] * M[j][k];
                if ( i != j )
                    F[j][i] = F[i][j];
            }
        }
    }
    
    // 2/4: ML (p+1,n+1):
    for ( i = 0 ; i < p1 ; ++i )
    {
        F[i][p1] = 1.0;
        for ( j = p1+1 ; j < nF ; ++j )
            F[i][j] = M[i][j-p1];
    }
    
    // 3/4: ML' (n+1,p+1):
    for ( j = 0 ; j < p1 ; ++j )
    {
        F[p1][j] = 1.0;
        for ( i = p1+1 ; i < nF ; ++i )
            F[i][j] = M[j][i-p1];
    }
    
    // 4/4: 0 (n+1,n+1):
    for ( i = p1 ; i < nF ; ++i )
        for ( j = p1 ; j < nF ; ++j )
            F[i][j] = 0.0;
    
    
#ifdef DEBUG
    _out << std::endl << "F=";
    for ( i = 0 ; i < nF ; ++i )
    {
        _out << "\t";
        for ( j = 0 ; j < nF ; ++j )
            _out << std::setw(12) << F[i][j] << " ";
        _out << std::endl;
    }
#endif
    
    for ( i = 0 ; i < p1 ; ++i )
        delete [] M[i];
    delete [] M;
    
    bool error = false;
    
    // SVD decomposition of the F matrix (F = U.W.V'):
    // -----------------------------------------------
    // (F will be transformed in U)
    
    double  * W = new double  [nF];
    double ** V = new double *[nF];
    for ( i = 0 ; i < nF ; ++i )
        V[i] = new double[nF];
    
    std::string error_msg;
    
    if ( NOMAD::SVD_decomposition ( error_msg , F , W , V , nF , nF , max_mpn ) ) {
        
        // compute condition number:
        compute_cond ( W , nF , eps );
        
#ifdef DEBUG
        _out << std::endl << "F=";
        for ( i = 0 ; i < nF ; ++i )
        {
            _out << "\t";
            for ( j = 0 ; j < nF ; ++j )
                _out << std::setw(12) << F[i][j] << " ";
            _out << std::endl;
        }
        
        _out << std::endl << "W=\t";
        for ( i = 0 ; i < nF ; ++i )
            _out << std::setw(12) << W[i] << " ";
        _out << std::endl << std::endl << "cond=" << _cond << std::endl;
        
        _out << std::endl << "V=";
        for ( i = 0 ; i < nF ; ++i )
        {
            _out << "\t";
            for ( j = 0 ; j < nF ; ++j )
                _out << std::setw(12) << V[i][j] << " ";
            _out << std::endl;
        }
#endif
        
    }
    else
    {
#ifdef DEBUG
        _out << std::endl << "NOMAD::Quad_Model::construct_MFN_model(): "
        << "SVD decomposition (" << error_msg << ")"
        << std::endl << std::endl;
#endif
        error = true;
        _cond.clear();
    }
    
    // resolution of system F.[mu alpha_L]'=[f(Y) 0]' :
    // ------------------------------------------------
    if ( !error )
    {
        int m = static_cast<int> ( _bbot.size() );
        for ( i = 0 ; i < m ; ++i )
            if ( _alpha[i] )
                solve_MFN_system ( F , W  , V , i , *_alpha[i] , eps );
    }
    
    // free memory:
    for ( i = 0 ; i < nF ; ++i )
    {
        delete [] F[i];
        delete [] V[i];
    }
    delete [] F;
    delete [] V;
    delete [] W;
    
#ifdef DEBUG
    _out.close_block();
#endif
    
    return !error;
}

/*--------------------------------------------------*/
/*  resolution of system F.[mu alpha_L]'=[f(Y) 0]'  */
/*  for MFN interpolation (private)                 */
/*--------------------------------------------------*/
void NOMAD::Quad_Model::solve_MFN_system ( double      ** F         ,
                                          double       * W         ,
                                          double      ** V         ,
                                          int            bbo_index ,
                                          NOMAD::Point & alpha     ,
                                          double         eps	      ) const
{
    // resize the alpha vector:
    if ( alpha.size() != _n_alpha )
        alpha.reset ( _n_alpha , 0.0 );
    
    int i , k , k1 , k2 ,
    np1 = _nfree + 1  ,
    nm1 = _nfree - 1  ,
    p1  = get_nY()    ,
    nF  = np1 + p1;
    
    // step 1/2: find alpha_L and mu:
    // ---------
    double * alpha_tmp = new double [np1];
    double * mu_tmp    = new double [ p1];
    double * mu        = new double [ p1];
    
    // if F is singular, some W values will be zero (or near zero);
    // each value that is smaller than eps is ignored:
    for ( i = 0 ; i < p1 ; ++i )
    {
        mu_tmp[i] = 0.0;
        if ( W[i] > eps )
            for ( k = 0 ; k < p1 ; ++k )
                mu_tmp[i] += F[k][i] *
                ( _Y[k]->get_bb_outputs()[bbo_index].value() ) / W[i];
    }
    
    for ( i = p1 ; i < nF ; ++i )
    {
        alpha_tmp[i-p1] = 0.0;
        if ( W[i] > eps )
            for ( k = 0 ; k < p1 ; ++k )
                alpha_tmp[i-p1] += F[k][i] *
                ( _Y[k]->get_bb_outputs()[bbo_index].value() ) / W[i];
    }
    
    for ( i = 0 ; i < p1 ; ++i )
    {
        mu[i] = 0.0;
        for ( k = 0 ; k < p1 ; ++k )
            mu[i] += V[i][k] * mu_tmp[k];
        for ( k = p1 ; k < nF ; ++k )
            mu[i] += V[i][k] * alpha_tmp[k-p1];
    }
    
    for ( i = p1 ; i < nF ; ++i )
    {
        alpha[i-p1] = 0.0;
        for ( k = 0 ; k < p1 ; ++k )
            alpha[i-p1] += V[i][k] * mu_tmp[k];
        for ( k = p1 ; k < nF ; ++k )
            alpha[i-p1] += V[i][k] * alpha_tmp[k-p1];
    }
    
    delete [] alpha_tmp;
    delete [] mu_tmp;
    
#ifdef DEBUG
    _out << std::endl << "output #" << bbo_index << ": mu=\t";
    for ( i = 0 ; i < p1 ; ++i )
        _out << std::setw(12) << mu[i] << " ";
    _out << std::endl;
    
    _out << std::endl << "output #" << bbo_index << ": alpha_intermediate=\t";
    for ( i = 0 ; i < alpha.size() ; ++i )
        _out << std::setw(12) << alpha[i] << " ";
    _out << std::endl;
    
#endif
    
    // step 2/2: find alpha_Q:
    // ---------
    for ( i = 0 ; i < _nfree ; ++i )
    {
        alpha[i+np1] = 0.0;
        for ( k = 0 ; k < p1 ; ++k )
            alpha[i+np1] += mu[k] * pow ( (*_Y[k])[_index[i+1]-1].value() , 2.0 ) / 2.0;
    }
    
    for ( k1 = 0 ; k1 < nm1 ; ++k1 )
        for ( k2 = k1+1 ; k2 < _nfree ; ++k2 )
        {
            alpha[i+np1] = 0.0;
            for ( k = 0 ; k < p1 ; ++k )
                alpha[i+np1] += mu[k] * (*_Y[k])[_index[k1+1]-1].value() * (*_Y[k])[_index[k2+1]-1].value();
            ++i;
        }
    
    delete [] mu;
}

/*-----------------------------------------------------------*/
/*            check the interpolation set Y (private)        */
/*-----------------------------------------------------------*/
bool NOMAD::Quad_Model::check_Y ( void ) const
{
    if ( _Y.empty() ) {
#ifdef DEBUG
        _out << std::endl << "NOMAD::Quad_Model::check_Y(): set Y is empty"
        << std::endl << std::endl;
#endif
        return false;
    }
    
    int nY = get_nY();
    int m  = static_cast<int> ( _bbot.size() );
    
    for ( int k = 0 ; k < nY ; ++k ) {
        
        if ( _Y[k] == NULL ) {
#ifdef DEBUG
            _out << std::endl
            << "NOMAD::Quad_Model::check_Y(): NULL pointer in the set Y"
            << std::endl << std::endl;
#endif
            return false;
        }
        
        if ( _Y[k]->get_eval_status() != NOMAD::EVAL_OK ) {
#ifdef DEBUG
            _out << std::endl
            << "NOMAD::Quad_Model::check_Y(): a point in Y failed to evaluate"
            << std::endl << std::endl;
#endif
            return false;
        }
        
        const NOMAD::Point & bbo = _Y[k]->get_bb_outputs();
        
        if ( !bbo.is_complete() ) {
#ifdef DEBUG
            _out << std::endl
            << "NOMAD::Quad_Model::check_Y(): some bb outputs in Y are not defined"
            << std::endl << std::endl;
#endif
            return false;
        }
        
        if ( bbo.size() != m ) {
#ifdef DEBUG
            _out << std::endl
            << "NOMAD::Quad_Model::check_Y(): "
            << "bb outputs in Y do not have the same dimension"
            << std::endl << std::endl;
#endif
            return false;
        }
        
        if ( _Y[k]->size() != _n ) {
#ifdef DEBUG
            _out << std::endl
            << "NOMAD::Quad_Model::check_Y(): "
            << "points in Y do not have the same dimension"
            << std::endl << std::endl;
#endif
            return false;
        }
    }
    return true;
}

/*----------------------------------------------------*/
/*     check if the model is ready for evaluations    */
/*----------------------------------------------------*/
bool NOMAD::Quad_Model::check ( void ) const
{
    if ( !_alpha )
        return false;
    
    int nalpha = (_nfree+1)*(_nfree+2)/2;
    int i , m  = static_cast<int> ( _bbot.size() );
    
    for ( int bbo_index = 0 ; bbo_index < m ; ++bbo_index )
    {
        
        if ( _alpha[bbo_index] )
        {
            
            if ( _alpha[bbo_index]->size() != nalpha )
                return false;
            
            for ( i = 0 ; i < nalpha ; ++i )
                if ( !(*_alpha[bbo_index])[i].is_defined() )
                    return false;
        }
    }
    
    return true;
}


/*--------------------------------------------------------------------------*/
/*                      evaluate a model at a given point                   */
/*--------------------------------------------------------------------------*/
/*  . the method assumes that x.size()==_n, alpha.is_complete(), and        */
/*    alpha.size()==(_nfree+1)*(_nfree+2)/2                                 */
/*  . a more efficient version is used in Quad_Model_Evaluator::eval_x() )  */
/*--------------------------------------------------------------------------*/
NOMAD::Double NOMAD::Quad_Model::eval ( const NOMAD::Point & x     ,
                                       const NOMAD::Point & alpha   ) const
{
    int i , j , k = 1 , nm1 = _n-1;
    NOMAD::Double z = alpha[0];
    
    for ( i = 0 ; i < _n ; ++i )
    {
        if ( !_fixed_vars[i] )
        {
            z += x[i] * ( alpha[k] + 0.5 * alpha[k+_nfree] * x[i] );
            ++k;
        }
    }
    
    k += _nfree;
    
    for ( i = 0 ; i < nm1 ; ++i )
        if ( !_fixed_vars[i] )
            for ( j = i+1 ; j < _n ; ++j )
                if ( !_fixed_vars[j] )
                    z += alpha[k++] * x[i] * x[j];
    
    return z;
}


/*----------------------------------------------------------------*/
/*             compute model h and f values at a point            */
/*----------------------------------------------------------------*/
void NOMAD::Quad_Model::eval_hf ( const NOMAD::Point  & x      ,
                                 const NOMAD::Double & h_min  ,
                                 NOMAD::hnorm_type     h_norm ,
                                 NOMAD::Double       & h      ,
                                 NOMAD::Double       & f        ) const
{
    f.clear();
    h = 0.0;
    
    int           m  = static_cast<int>(_bbot.size());
    NOMAD::Double bboi;
    
    for ( int i = 0 ; i < m ; ++i ) {
        
        if ( _alpha[i] ) {
            
            bboi = eval ( x , *_alpha[i] );
            
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
    }
    
    if ( h_norm == NOMAD::L2 )
        h = h.sqrt();
}




/*-----------------------------------------------------*/
/*  compute the maximal relative error of a model for  */
/*  the interpolation set (private)                    */
/*-----------------------------------------------------*/
NOMAD::Double NOMAD::Quad_Model::compute_max_rel_err ( void ) const
{
    NOMAD::Double truth_value , model_value , rel_err , max_rel_err;
    int           k , nY = get_nY() , m = static_cast<int> ( _bbot.size() );
    
    for ( int bbo_index = 0 ; bbo_index < m ; ++bbo_index ) {
        if ( _alpha[bbo_index] ) {
            for ( k = 0 ; k < nY ; ++k ) {
                if ( _Y[k] && _Y[k]->get_eval_status() == NOMAD::EVAL_OK ) {
                    truth_value = _Y[k]->get_bb_outputs()[bbo_index];
                    if ( truth_value.is_defined() ) {
                        model_value = eval ( *_Y[k] , *_alpha[bbo_index] );
                        if ( model_value.is_defined() ) {
                            if ( truth_value.abs() != 0.0 ) {
                                rel_err = (truth_value-model_value).abs() / truth_value.abs();
                                if ( !max_rel_err.is_defined() || rel_err > max_rel_err )
                                    max_rel_err = rel_err;
                            }
                        }
                    }
                }
            }
        }
    }
    return max_rel_err;
}

/*---------------------------------------------*/
/*  compute the cumulated error of a model     */
/*  for the points of the interpolation set Y  */
/*  and for one output in particular           */
/*  (private)                                  */
/*---------------------------------------------*/
void NOMAD::Quad_Model::compute_model_error ( int             bbo_index   ,
                                             NOMAD::Double & error       ,
                                             NOMAD::Double & min_rel_err ,
                                             NOMAD::Double & max_rel_err ,
                                             NOMAD::Double & avg_rel_err   ) const
{
    NOMAD::Double truth_value , model_value , rel_err;
    int  nY  = get_nY() , cnt = 0;
    bool chk = true;
    
    max_rel_err.clear();
    min_rel_err.clear();
    avg_rel_err = error = 0.0;
    
#ifdef DEBUG
    std::ostringstream msg;
    msg << "output #" << bbo_index;
    _out.open_block ( msg.str() );
#endif
    
    for ( int k = 0 ; k < nY ; ++k )
        if ( _Y[k] && _Y[k]->get_eval_status() == NOMAD::EVAL_OK )
        {
            truth_value = _Y[k]->get_bb_outputs()[bbo_index];
            
            if ( truth_value.is_defined() ) 
            {
                model_value = eval ( *_Y[k] , *_alpha[bbo_index] );
                if ( model_value.is_defined() )
                {
                    rel_err.clear();
                    if ( truth_value.abs() != 0.0 )
                        rel_err = (truth_value-model_value).abs() / truth_value.abs();
                    else 
                    {
                        if (truth_value.abs()==model_value.abs())
                            rel_err=0.0;
                        else 
                            rel_err=NOMAD::INF;
                    }
                    if ( !max_rel_err.is_defined() || rel_err > max_rel_err )
                        max_rel_err = rel_err;
                    if ( !min_rel_err.is_defined() || rel_err < min_rel_err )
                        min_rel_err = rel_err;
                    avg_rel_err += rel_err;
                    ++cnt;
                    
#ifdef DEBUG
                    _out << "Y[" << k << "]= ( ";
                    _Y[k]->NOMAD::Point::display ( _out );
                    _out << " )" << " f=" << truth_value
                    << " m=" << model_value << " error^2="
                    << ( model_value - truth_value ).pow2()
                    << " rel_err=" << rel_err
                    << std::endl;
#endif
                    error += ( model_value - truth_value ).pow2();
                }
                else
                {
                    chk = false;
                    break;
                }
            }
            else
            {
                chk = false;
                break;
            }
        }
    
#ifdef DEBUG
    _out.close_block();
#endif
    
    if ( chk)
    {  // Case where chk is true (at least one model_value and the corresponding thruth value were defined => cnt != 0)
        error       = error.sqrt();
        avg_rel_err = avg_rel_err / cnt;
    }
    else
    {
        error.clear();
        min_rel_err.clear();
        max_rel_err.clear();
        avg_rel_err.clear();
    }
}

/*-----------------------------------------------------------*/
/*                display the model coefficients             */
/*-----------------------------------------------------------*/
void NOMAD::Quad_Model::display_model_coeffs ( const NOMAD::Display & out ) const
{
    if ( _error_flag ) {
        out << "model coefficients: could not be constructed" << std::endl;
        return;
    }
    
    int m = static_cast<int> ( _bbot.size() );
    
    out << NOMAD::open_block ( "model coefficients" );
    for ( int i = 0 ; i < m ; ++i ) {
        out << "output #";
        out.display_int_w ( i , m );
        out << ": ";
        if ( _alpha[i] ) {
            out<< "[ ";
            _alpha[i]->display ( out , " " , 6 );
            out << " ]";
        }
        else
            out << "NULL";
        out << std::endl;
    }
    out.close_block();
}

/*-----------------------------------------------------------*/
/*                 display the interpolation set Y           */
/*-----------------------------------------------------------*/
void NOMAD::Quad_Model::display_Y ( const NOMAD::Display & out   ,
                                   const std::string    & title   ) const
{
    out << NOMAD::open_block ( title );
    int nY = get_nY();
    for ( int k = 0 ; k < nY ; ++k ) {
        out << "#";
        out.display_int_w ( k , nY );
        out << ": ";
        if ( _Y[k] ) {
            out << "( ";
            _Y[k]->NOMAD::Point::display    ( out , " " , 12 );
            out << " ) bbo=[ ";
            _Y[k]->get_bb_outputs().display ( out , " " , 12 );
            out << " ]";
        }
        else
            out << "NULL";
        out << std::endl;
    }
    out.close_block();
}

/*-------------------------------------------------------*/
/*  display cumulated error on the interpolation points  */
/*-------------------------------------------------------*/
void NOMAD::Quad_Model::display_Y_error ( const NOMAD::Display & out ) const
{
    if ( _error_flag ) {
        out << "model error on the interpolation set: cannot be computed"
        << std::endl;
        return;
    }
    
    int i ;
    int index = -1;
    int m = static_cast<int> ( _bbot.size() );
    
    for ( i = 0 ; i < m ; ++i )
        if ( _alpha[i] ) {
            if ( index >= 0 ) {
                index = -1;
                break;
            }
            else
                index = i;
        }
    
    NOMAD::Double error , min_rel_err , max_rel_err , avg_rel_err;
    
    // only one output:
    if ( index >= 0 ) {
        compute_model_error ( index , error , min_rel_err , max_rel_err , avg_rel_err );
        out << "model errors on the interpolation set: error="
        << error << " min_rel_err=" << min_rel_err
        << " max_rel_err=" << max_rel_err << " avg_rel_err=" << avg_rel_err
        << std::endl;
    }
    
    // several outputs:
    else {
        
        out.open_block ( "model error on the interpolation set" );
        
        NOMAD::Double error_i , min_rel_err_i , max_rel_err_i , avg_rel_err_i;
        
        error = avg_rel_err = 0.0;
        min_rel_err.clear();
        max_rel_err.clear();
        
        int cnt = 0;
        
        for ( i = 0 ; i < m ; ++i )
            if ( _alpha[i] ) {
                
                ++cnt;
                
                compute_model_error ( i             ,
                                     error_i       ,
                                     min_rel_err_i ,
                                     max_rel_err_i ,
                                     avg_rel_err_i   );
                
                if (error_i.is_defined())
                    error       += error_i;
                if (avg_rel_err_i.is_defined())
                    avg_rel_err += avg_rel_err_i;
                if ( !min_rel_err.is_defined() || min_rel_err_i < min_rel_err )
                    min_rel_err = min_rel_err_i;
                if ( !max_rel_err.is_defined() || max_rel_err_i > max_rel_err )
                    max_rel_err = max_rel_err_i;
                
                out << "output #" << i << ": error=" << error_i
                << " min_rel_err=" << min_rel_err_i << " max_rel_err="
                << max_rel_err_i << " avg_rel_err=" << avg_rel_err_i << std::endl;
            }
        
        out << std::endl << "global: error=" << error
        << " min_rel_err=" << min_rel_err
        << " max_rel_err=" << max_rel_err
        << " avg_rel_err=" << avg_rel_err / cnt
        << std::endl << NOMAD::close_block();
    }
}

/*-----------------------------------------------*/
/*     display Lagrange polynomials (private)    */
/*-----------------------------------------------*/
void NOMAD::Quad_Model::display_lagrange_polynomials
( const std::vector<NOMAD::Point      *> & l ,
 const std::vector<NOMAD::Eval_Point *> & Y   ) const
{
    int i , j , nY = static_cast<int> ( Y.size() );
    
    // display Lagrange polynomials:
    _out << std::endl << NOMAD::open_block ( "Lagrange polynomials" );
    for ( i = 0 ; i < _n_alpha ; ++i ) {
        _out << "l[";
        _out.display_int_w ( i , _n_alpha );
        _out << "] = [ ";
        l[i]->NOMAD::Point::display ( _out , " " , 14 , -1 );
        _out << "]" << std::endl;
    }
    _out.close_block();
    
    // display current set Y:
    _out << std::endl << NOMAD::open_block ( "current set Y" );
    for ( i = 0 ; i < nY ; ++i ) {
        _out << "Y[";
        _out.display_int_w ( i , nY );
        _out << "] = ";
        if ( Y[i] ) {
            _out << "( ";
            Y[i]->NOMAD::Point::display ( _out , " " , 6 , -1 );
            _out << " )";
        }
        else
            _out << "NULL";
        _out << std::endl;
    }
    
    _out.close_block();
    
    // display l(Y): should be the identity matrix:
    NOMAD::Double tmp , err = 0.0;
    _out << std::endl << NOMAD::open_block ( "l(Y)" );
    for ( i = 0 ; i < _n_alpha ; ++i ) {
        _out << "l[";
        _out.display_int_w ( i , _n_alpha );
        _out << "]: ";
        for ( j = 0 ; j < _n_alpha ; ++j ) {
            tmp.clear();
            if ( j < nY && Y[j] ) {
                tmp  = eval ( *Y[j] , *l[i] );
                err += (i==j) ? (tmp-1.0).abs() : tmp.abs();
            }
            tmp.display ( _out , "%6.3f" );
            _out << " ";
        }
        _out << std::endl;
    }
    _out << std::endl << "error (with identity) = "
    << err << std::endl << NOMAD::close_block() << std::endl;
}
