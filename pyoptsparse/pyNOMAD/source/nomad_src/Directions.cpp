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
 \file   Directions.cpp
 \brief  Set of polling directions (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-13
 \see    Directions.hpp
 */
#include "Directions.hpp"
#include <math.h>
#include "XMesh.hpp"
#include "SMesh.hpp"


/*---------------------------------------------------------*/
/*                       constructor                       */
/*---------------------------------------------------------*/
NOMAD::Directions::Directions
( int                                     nc                 ,
 const std::set<NOMAD::direction_type> & direction_types    ,
 const std::set<NOMAD::direction_type> & sec_poll_dir_types ,
 const NOMAD::Display                  & out                  )
: _nc                 ( nc                 ) ,
_direction_types    ( direction_types    ) ,
_sec_poll_dir_types ( sec_poll_dir_types ) ,
_is_binary          ( false              ) ,
_is_categorical     ( false              ) ,
_lt_initialized     ( false              ) ,
_out                ( out                )
{
    // check the directions:
    if ( _direction_types.find ( NOMAD::NO_DIRECTION ) != _direction_types.end()    )
        _direction_types.clear();
    
    if ( _sec_poll_dir_types.find ( NOMAD::NO_DIRECTION ) != _sec_poll_dir_types.end() )
        _sec_poll_dir_types.clear();
    
    // is_orthomads: true if at least one direction is of type Ortho-MADS:
    _is_orthomads = NOMAD::dirs_have_orthomads ( _direction_types );
    if ( !_is_orthomads )
        _is_orthomads = NOMAD::dirs_have_orthomads ( _sec_poll_dir_types );
    
    
}

/*---------------------------------------------------------*/
/*                         destructor                      */
/*---------------------------------------------------------*/
NOMAD::Directions::~Directions ( void )
{
    if ( _lt_initialized )
    {
        int n = 2 * NOMAD::L_LIMITS;
        for ( int i = 0 ; i <= n ; ++i )
            delete _bl[i];
    }
    
}

/*---------------------------------------------------------*/
/*              LT-MADS initializations (private)          */
/*---------------------------------------------------------*/
void NOMAD::Directions::lt_mads_init ( void )
{
    int n = 2 * NOMAD::L_LIMITS;
    for ( int i = 0 ; i <= n ; ++i )
    {
        _bl   [i] = NULL;
        _hat_i[i] = -1;
    }
    _lt_initialized = true;
}

/*---------------------------------------------------------*/
/*                         set_binary                      */
/*---------------------------------------------------------*/
void NOMAD::Directions::set_binary ( void )
{
    _is_binary      = true;
    _is_categorical = false;
    _is_orthomads   = false;
    _direction_types.clear();
    _direction_types.insert ( NOMAD::GPS_BINARY );
    if ( !_sec_poll_dir_types.empty() ) {
        _sec_poll_dir_types.clear();
        _sec_poll_dir_types.insert ( NOMAD::GPS_BINARY );
    }
}

/*---------------------------------------------------------*/
/*                       set_categorical                   */
/*---------------------------------------------------------*/
void NOMAD::Directions::set_categorical ( void )
{
    _is_categorical = true;
    _is_binary      = false;
    _is_orthomads   = false;
    _direction_types.clear();
    _sec_poll_dir_types.clear();
}

/*----------------------------------------------------------------------*/
/*  compute binary directions when all groups of variables are binary   */
/*  (private)                                                           */
/*----------------------------------------------------------------------*/
void NOMAD::Directions::compute_binary_directions
( std::list<NOMAD::Direction> & d ) const
{
    // _GPS_BINARY_ n directions:
    NOMAD::Direction * pd;
    for ( int i = 0 ; i < _nc ; ++i ) {
        d.push_back ( NOMAD::Direction ( _nc , 0.0 , NOMAD::GPS_BINARY ) );
        pd = &(*(--d.end()));
        (*pd)[i] = 1.0;
    }
}


/*----------------------------------------------------------------------------*/
/*            get the directions on a unit n-sphere for a given mesh          */
/*----------------------------------------------------------------------------*/
void NOMAD::Directions::compute ( std::list<NOMAD::Direction> & dirs		,
                                 NOMAD::poll_type              poll			,
                                 const NOMAD::OrthogonalMesh  & mesh		)
{
    
    dirs.clear();
    
    // categorical variables: do nothing:
    if ( _is_categorical )
        return;
    
    // binary variables: we use special directions:
    if ( _is_binary )
    {
        compute_binary_directions ( dirs );
        return;
    }
    
    NOMAD::Double                               pow_gps , cst;
    const NOMAD::Direction                    * bl;
    NOMAD::Direction                          * pd;
    int                                         i , j , k , hat_i;
    std::list<NOMAD::Direction>::iterator it2 , end2;
    
    const NOMAD::Point mesh_indices=mesh.get_mesh_indices();
    int mesh_index=int(mesh_indices[0].value()); // single mesh_index used only for isotropic mesh
    
    
    /*-----------------------------------*/
    /*  loop on the types of directions  */
    /*-----------------------------------*/
    const std::set<NOMAD::direction_type> & dir_types =
    (poll == NOMAD::PRIMARY) ? _direction_types : _sec_poll_dir_types;
    
    std::set<NOMAD::direction_type>::const_iterator it , end = dir_types.end() ;
    for ( it = dir_types.begin() ; it != end ; ++it )
    {
        
        if ( *it == NOMAD::UNDEFINED_DIRECTION ||
            *it == NOMAD::NO_DIRECTION        ||
            *it == NOMAD::MODEL_SEARCH_DIR       )
            continue;
        
        /*--------------------------------------------------*/
        /*  Ortho-MADS                                      */
        /* ---> dirs are on a unit n-sphere by construction */
        /*--------------------------------------------------*/
        if ( NOMAD::dir_is_orthomads (*it) )
        {
            bool success_dir;
            
            
            NOMAD::Direction dir ( _nc , 0.0 , *it );
            NOMAD::Double    alpha_t_l;
            
            success_dir=compute_dir_on_unit_sphere ( dir );
            
            if ( success_dir )
            {
#ifdef DEBUG
                if ( dynamic_cast<const NOMAD::XMesh*> ( &mesh ) )
                {
                    _out << std::endl
                    << NOMAD::open_block ( "compute Ortho-MADS directions with XMesh" )
                    << "type              = " << *it << std::endl;
                    NOMAD::Point mesh_indices=mesh.get_mesh_indices();
                    _out<< "Mesh indices  = ( ";
                    mesh_indices.NOMAD::Point::display( _out, " " , -1 , -1 );
                    _out << " )" << std::endl;
                    _out<< "Unit sphere Direction  = ( ";
                    dir.NOMAD::Point::display ( _out , " " , -1 , -1 );
                    _out << " )" << std::endl << NOMAD::close_block();
                }
                else if ( dynamic_cast<const NOMAD::SMesh*> ( &mesh ) )
                {
                    _out << std::endl
                    << NOMAD::open_block ( "compute Ortho-MADS directions with SMesh" )
                    << "type              = " << *it
                    << " on isotropic mesh " << std::endl
                    << "mesh index   (lk) = " << mesh_index   << std::endl
                    << "alpha     (tk,lk) = " << alpha_t_l    << std::endl
                    << "Unit sphere Direction  = ( ";
                    dir.NOMAD::Point::display ( _out , " " , -1 , -1 );
                    _out << " )" << std::endl << NOMAD::close_block();
                    
                }
                
#endif
                // Ortho-MADS 2n and n+1:
                // ----------------------
                if ( *it == NOMAD::ORTHO_2N || *it == NOMAD::ORTHO_NP1_QUAD || *it == NOMAD::ORTHO_NP1_NEG )
                {
                    
                    // creation of the 2n directions:
                    // For ORTHO_NP1 the reduction from 2N to N+1 is done in MADS::POLL
                    NOMAD::Direction ** H = new NOMAD::Direction * [2*_nc];
                    
                    // Ordering D_k alternates Hk and -Hk instead of [H_k -H_k]
                    for ( i = 0 ; i < _nc ; ++i )
                    {
                        dirs.push_back ( NOMAD::Direction ( _nc , 0.0 , *it ) );
                        H[i    ] = &(*(--dirs.end()));
                        dirs.push_back ( NOMAD::Direction ( _nc , 0.0 , *it ) );
                        H[i+_nc] = &(*(--dirs.end()));
                    }
                    
                    // Householder transformations on the 2n directions on a unit n-sphere:
                    householder ( dir , true , H );
                    
                    
                    delete [] H;
                    
                    
                }
                
                // Ortho-MADS 1 or Ortho-MADS 2:
                // -----------------------------
                else
                {
                    dirs.push_back ( dir );
                    if ( *it == NOMAD::ORTHO_2 )
                        dirs.push_back ( -dir );
                }
            }
        }
        
        /*--------------------------------*/
        /*  LT-MADS                       */
        /* ---> dirs NOT on unit n-sphere */
        /*--------------------------------*/
        else if ( NOMAD::dir_is_ltmads (*it) )
        {
            
            if ( !_lt_initialized)
                lt_mads_init();
            
            bl = get_bl ( mesh , *it , hat_i );
            
            // LT-MADS 1 or LT-MADS 2: -b(l) and/or b(l):
            // ------------------------------------------
            if ( *it == NOMAD::LT_1 || *it == NOMAD::LT_2 )
            {
                dirs.push_back ( - *bl );
                if ( *it == NOMAD::LT_2 )
                    dirs.push_back ( *bl );
            }
            
            // LT-MADS 2n or LT-MADS n+1:
            // --------------------------
            else {
                
                // row permutation vector:
                int * row_permutation_vector  = new int [_nc];
                row_permutation_vector[hat_i] = hat_i;
                
                NOMAD::Random_Pickup rp ( _nc );
                
                for ( i = 0 ; i < _nc ; ++i )
                    if ( i != hat_i )
                    {
                        j = rp.pickup();
                        if ( j == hat_i )
                            j = rp.pickup();
                        row_permutation_vector[i] = j;
                    }
                
                rp.reset();
                
                for ( j = 0 ; j < _nc ; ++j )
                {
                    i = rp.pickup();
                    if ( i != hat_i )
                    {
                        
                        dirs.push_back ( NOMAD::Direction ( _nc , 0.0 , *it ) );
                        pd = &(*(--dirs.end()));
                        
                        create_lt_direction ( mesh , *it , i , hat_i , pd );
                        permute_coords ( *pd , row_permutation_vector );
                    }
                    else
                        dirs.push_back ( *bl );
                    
                    // completion to maximal basis:
                    if ( *it == NOMAD::LT_2N )
                        dirs.push_back ( NOMAD::Direction ( - *(--dirs.end()) , NOMAD::LT_2N ) );
                    
                }
                
                delete [] row_permutation_vector;
                
                // completion to minimal basis:
                if ( *it == NOMAD::LT_NP1 )
                {
                    dirs.push_back ( NOMAD::Direction ( _nc , 0.0 , NOMAD::LT_NP1 ) );
                    pd = &(*(--dirs.end()));
                    
                    end2 = --dirs.end();
                    for ( it2 = dirs.begin() ; it2 != end2 ; ++it2 )
                    {
                        for ( i = 0 ; i < _nc ; ++i )
                            (*pd)[i] -= (*it2)[i];
                    }
                }
            }
            
            
        }
        
        /*--------------------------------*/
        /*  GPS                           */
        /* ---> dirs NOT on unit n-sphere */
        /*--------------------------------*/
        else
        {
            
            // GPS for binary variables: should'nt be here:
            if ( *it == NOMAD::GPS_BINARY )
            {
#ifdef DEBUG
                _out << std::endl << "Warning (" << "Directions.cpp" << ", " << __LINE__
                << "): tried to generate binary directions at the wrong place)"
                << std::endl << std::endl;
#endif
                dirs.clear();
                return;
            }
            
            // this value represents the non-zero values in GPS directions
            // (it is tau^{|ell_k|/2}, and it is such that delta_k * pow_gps = Delta_k):
            if ( !pow_gps.is_defined() )
                pow_gps = pow ( mesh.get_update_basis() , abs(mesh_index) / 2.0 );
            
            // GPS 2n, static:
            // ---------------
            
            if ( *it == NOMAD::GPS_2N_STATIC )
            {
                for ( i = 0 ; i < _nc ; ++i )
                {
                    dirs.push_back ( NOMAD::Direction ( _nc , 0.0 , NOMAD::GPS_2N_STATIC ) );
                    pd = &(*(--dirs.end()));
                    (*pd)[i] = pow_gps;
                    
                    dirs.push_back ( NOMAD::Direction ( _nc , 0.0 , NOMAD::GPS_2N_STATIC ) );
                    pd = &(*(--dirs.end()));
                    (*pd)[i] = -pow_gps;
                }
            }
            
            // GPS 2n, random:
            // ---------------
            else if ( *it == NOMAD::GPS_2N_RAND )
            {
                
                int v , cnt;
                
                std::list  <int>::const_iterator end3;
                std::list  <int>::iterator       it3;
                std::list  <int> rem_cols;
                std::vector<int> rem_col_by_row ( _nc );
                
                // creation of the 2n directions:
                std::vector<NOMAD::Direction *> pdirs ( 2 * _nc );
                
                for ( i = 0 ; i < _nc ; ++i )
                {
                    
                    dirs.push_back ( NOMAD::Direction ( _nc , 0.0 , NOMAD::GPS_2N_RAND ) );
                    pdirs[i] = &(*(--dirs.end()));
                    (*pdirs[i])[i] = pow_gps;
                    
                    dirs.push_back ( NOMAD::Direction ( _nc , 0.0 , NOMAD::GPS_2N_RAND ) );
                    pdirs[i+_nc] = &(*(--dirs.end()));
                    (*pdirs[i+_nc])[i] = -pow_gps;
                    
                    rem_cols.push_back(i);
                    rem_col_by_row[i] = i;
                }
                
                // random perturbations:
                for ( i = 1 ; i < _nc ; ++i )
                {
                    if ( rem_col_by_row[i] > 0 )
                    {
                        v = NOMAD::RNG::rand()%3 - 1; // v in { -1 , 0 , 1 }
                        if ( v )
                        {
                            
                            // we decide a (i,j) couple:
                            k    = NOMAD::RNG::rand()%(rem_col_by_row[i]);
                            cnt  = 0;
                            end3 = rem_cols.end();
                            it3  = rem_cols.begin();
                            while ( cnt != k )
                            {
                                ++it3;
                                ++cnt;
                            }
                            j = *it3;
                            
                            // the perturbation:
                            (*pdirs[i])[j] = (*pdirs[j+_nc])[i] = -v * pow_gps;
                            (*pdirs[j])[i] = (*pdirs[i+_nc])[j] =  v * pow_gps;
                            
                            // updates:
                            rem_cols.erase(it3);
                            it3 = rem_cols.begin();
                            while ( *it3 != i )
                                ++it3;
                            rem_cols.erase(it3);
                            for ( k = i+1 ; k < _nc ; ++k )
                                rem_col_by_row[k] -= j<k ? 2 : 1;
                        }
                    }
                }
            }
            
            // GPS n+1, static:
            // ----------------
            else if ( *it == NOMAD::GPS_NP1_STATIC )
            {
                
                // (n+1)^st direction:
                dirs.push_back ( NOMAD::Direction ( _nc , -pow_gps , NOMAD::GPS_NP1_STATIC ) );
                
                // first n directions:
                for ( i = 0 ; i < _nc ; ++i )
                {
                    dirs.push_back ( NOMAD::Direction ( _nc , 0.0 , NOMAD::GPS_NP1_STATIC ) );
                    pd = &(*(--dirs.end()));
                    (*pd)[i] = pow_gps;
                }
            }
            
            // GPS n+1, random:
            // ----------------
            else if ( *it == NOMAD::GPS_NP1_RAND )
            {
                
                // (n+1)^st direction:
                dirs.push_back ( NOMAD::Direction ( _nc , 0.0 , NOMAD::GPS_NP1_RAND ) );
                NOMAD::Direction * pd1 = &(*(--dirs.end()));
                
                NOMAD::Double d;
                
                // first n directions:
                for ( i = 0 ; i < _nc ; ++i )
                {
                    dirs.push_back ( NOMAD::Direction ( _nc , 0.0 , NOMAD::GPS_NP1_RAND ) );
                    pd = &(*(--dirs.end()));
                    
                    d = NOMAD::RNG::rand()%2 ? -pow_gps : pow_gps;
                    (*pd )[i] =  d;
                    (*pd1)[i] = -d;
                }
            }
            
            // GPS n+1, static, uniform angles:
            // --------------------------------
            else if ( *it == NOMAD::GPS_NP1_STATIC_UNIFORM )
            {
                
                cst = pow_gps * sqrt(static_cast<double>(_nc)*(_nc+1))/_nc;
                
                // n first directions:
                for ( j = _nc-1 ; j >= 0 ; --j )
                {
                    dirs.push_back ( NOMAD::Direction ( _nc , 0.0 , NOMAD::GPS_NP1_STATIC_UNIFORM ) );
                    pd = &(*(--dirs.end()));
                    
                    k = _nc-j-1;
                    
                    for ( i = 0 ; i < k ; ++i )
                        (*pd)[i] = -cst / sqrt(static_cast<double>(_nc-i)*(_nc-i+1));
                    
                    (*pd)[k] = (cst * (j+1)) / sqrt(static_cast<double>(j+1)*(j+2));
                    
                }
                
                // (n+1)^st direction:
                dirs.push_back ( NOMAD::Direction ( _nc , 0.0 , NOMAD::GPS_NP1_STATIC_UNIFORM ) );
                pd = &(*(--dirs.end()));
                for ( i = 0 ; i < _nc ; ++i )
                    (*pd)[i] = -cst / sqrt(static_cast<double>(_nc-i)*(_nc-i+1));
                
            }
            
            // GPS n+1, random, uniform angles:
            // --------------------------------
            // (we apply the procedure defined in
            // "Pattern Search Methods for user-provided points:
            // application to molecular geometry problems",
            // by Alberto, Nogueira, Rocha and Vicente,
            // SIOPT 14-4, 1216-1236, 2004, doi:10.1137/S1052623400377955)
            else if ( *it == NOMAD::GPS_NP1_RAND_UNIFORM )
            {
                
                cst = pow_gps * sqrt(static_cast<double>(_nc)*(_nc+1))/_nc;
                
                // n first directions:
                for ( j = _nc-1 ; j >= 0 ; --j )
                {
                    dirs.push_back ( NOMAD::Direction ( _nc , 0.0 , NOMAD::GPS_NP1_STATIC_UNIFORM ) );
                    pd = &(*(--dirs.end()));
                    
                    k = _nc-j-1;
                    
                    for ( i = 0 ; i < k ; ++i )
                        (*pd)[i] = -cst / sqrt(static_cast<double>(_nc-i)*(_nc-i+1));
                    
                    (*pd)[k] = (cst * (j+1)) / sqrt(static_cast<double>(j+1)*(j+2));
                    
                }
                
                // (n+1)^st direction:
                dirs.push_back ( NOMAD::Direction ( _nc , 0.0 , NOMAD::GPS_NP1_STATIC_UNIFORM ) );
                pd = &(*(--dirs.end()));
                for ( i = 0 ; i < _nc ; ++i )
                    (*pd)[i] = -cst / sqrt(static_cast<double>(_nc-i)*(_nc-i+1));
                
                // random perturbations:
                NOMAD::Random_Pickup                  rp   ( _nc );
                std::vector<bool>                     done ( _nc );
                bool                                  chg_sgn;
                std::list<NOMAD::Direction>::iterator it;
                NOMAD::Double                         tmp;
                
                end2 = dirs.end();
                for ( i = 0 ; i < _nc ; ++i )
                    done[i] = false;
                
                for ( i = 0 ; i < _nc ; ++i )
                {
                    k = rp.pickup();
                    if ( i != k && !done[i] )
                    {
                        chg_sgn = ( NOMAD::RNG::rand()%2 != 0 );
                        
                        for ( it = dirs.begin() ; it != end2 ; ++it )
                        {
                            tmp      = (*it)[i];
                            (*it)[i] = ( chg_sgn ? -1.0 : 1.0 ) * (*it)[k];
                            (*it)[k] = tmp;
                        }
                        
                        done[i] = done[k] = true;
                    }
                    else
                        if ( NOMAD::RNG::rand()%2 )
                            for ( it = dirs.begin() ; it != end2 ; ++it )
                                (*it)[i] = -(*it)[i];
                }
            }
        }
    } // end of loop on the types of directions
    
    // The direction are unscaled on a unit n-sphere when necessary (that is directions that are not orthomads)
    for ( it2 = dirs.begin() ; it2 != dirs.end() ; ++it2 )
        if ( ! NOMAD::dir_is_orthomads ( it2->get_type() ) )
        {
            NOMAD::Double norm_dir=it2->norm();
            for (int i=0 ; i < _nc ; ++i )
                (*it2)[i] /= norm_dir;
        }
    
}


/*-----------------------------------------------------------------------------*/
/*       compute a random direction	on a unit N-Sphere          			   */
/*  see http://en.wikipedia.org/wiki/N-sphere                                  */
/*-----------------------------------------------------------------------------*/
bool NOMAD::Directions::compute_dir_on_unit_sphere ( NOMAD::Direction	& random_dir ) const
{
    
    int           i;
    NOMAD::Double norm,d;
    
    for ( i = 0 ; i < _nc ; ++i )
        random_dir[i]=NOMAD::RNG::normal_rand(0,1);
    
    norm=random_dir.norm();
    
    if ( norm==0 )
        return false;
    
    
    for ( i = 0 ; i < _nc ; ++i )
        random_dir[i]/=norm;
    
    return true;
    
}


/*-----------------------------------------------------------------*/
/*  compute the squared norm of normalized(2u_t-e) for Ortho-MADS  */
/*  (private)                                                      */
/*-----------------------------------------------------------------*/
NOMAD::Double NOMAD::Directions::eval_ortho_norm ( const NOMAD::Double & x      ,
                                                  const NOMAD::Double & norm   ,
                                                  const NOMAD::Point  & b      ,
                                                  NOMAD::Point        & new_b    ) const
{
    NOMAD::Double fx = 0.0;
    
    for ( int i = 0 ; i < _nc ; ++i ) {
        new_b[i] = ( x * b[i] / norm ).round();
        fx += new_b[i]*new_b[i];
    }
    
    return fx;
}

/*--------------------------------------------------------*/
/*  get the expression of an integer t in inverse base p  */
/*  (private, static)                                     */
/*--------------------------------------------------------*/
NOMAD::Double NOMAD::Directions::get_phi ( int t , int p )
{
    int         div;
    int        size = int ( ceil ( log(static_cast<double>(t+1)) /
                                  log(static_cast<double>(p)) ) );
    int          ll = t;
    NOMAD::Double d = 0.0;
    
    for ( int i = 0 ; i < size ; ++i ) {
        div = NOMAD::Double ( pow ( p , size-i-1.0 ) ).round();
        d  += ( ll / div ) * pow ( static_cast<double>(p) , i-size );
        ll  = ll % div;
    }
    return d;
}

/*----------------------------------------------------------------*/
/*  . Householder transformation to generate _nc directions from  */
/*    a given direction                                        */
/*  . compute also H[i+nc] = -H[i] (completion to 2n directions)  */
/*  . private method                                              */
/*----------------------------------------------------------------*/
void NOMAD::Directions::householder ( const NOMAD::Direction  & dir     ,
                                     bool                      complete_to_2n ,
                                     NOMAD::Direction       ** H                ) const
{
    int i , j;
    
    NOMAD::Double norm2 = dir.squared_norm() , v , h2i;
    
    for ( i = 0 ; i < _nc ; ++i ) {
        
        h2i = 2 * dir[i];
        
        for ( j = 0 ; j < _nc ; ++j ) {
            
            // H[i]:
            (*H[i])[j] = v = (i==j) ? norm2 - h2i * dir[j] : - h2i * dir[j];
            
            // -H[i]:
            if ( complete_to_2n )
                (*H[i+_nc])[j] = -v;
        }
    }
}

/*---------------------------------------------------------*/
/*          get the LT-MADS b(l) direction (private)       */
/*---------------------------------------------------------*/
const NOMAD::Direction * NOMAD::Directions::get_bl ( const NOMAD::OrthogonalMesh & mesh,
                                                    NOMAD::direction_type   dtype      ,
                                                    int                   & hat_i       )
{
    
    NOMAD::Point  mesh_indices = mesh.get_mesh_indices();
    int mesh_index =static_cast<int>(mesh_indices[0].value());
    NOMAD::Direction * bl = _bl    [ mesh_index + NOMAD::L_LIMITS ];
    hat_i                 = _hat_i [ mesh_index + NOMAD::L_LIMITS ];
    
    if ( !bl )
    {
        hat_i = -1;
        create_lt_direction ( mesh , dtype , -1 , hat_i , bl );
    }
    
    return bl;
}

/*---------------------------------------------------------*/
/*          create a new LT-MADS direction (private)       */
/*---------------------------------------------------------*/
/*           (if hat_i == -1, a new b(l) direction         */
/*            is created and hat_i is set)                 */
/*---------------------------------------------------------*/
void NOMAD::Directions::create_lt_direction ( const NOMAD::OrthogonalMesh  & mesh ,
                                             NOMAD::direction_type dtype      ,
                                             int                   diag_i     ,
                                             int                 & hat_i      ,
                                             NOMAD::Direction   *& dir          )
{
    
    int mesh_index=static_cast<int>(mesh.get_mesh_indices()[0].value());
    int i_pow_tau =
    static_cast<int>
    ( ceil ( pow ( mesh.get_update_basis() , abs(mesh_index) / 2.0 ) ) );
    
    int j = diag_i+1;
    
    // new b(l) direction:
    if ( hat_i < 0 )
    {
        _hat_i [ mesh_index + NOMAD::L_LIMITS ] = diag_i = hat_i = NOMAD::RNG::rand()%_nc;
        _bl    [ mesh_index + NOMAD::L_LIMITS ] = dir
        = new NOMAD::Direction ( _nc, 0.0, dtype );
        
        j = 0;
    }
    
    (*dir)[diag_i] = (NOMAD::RNG::rand()%2) ? -i_pow_tau : i_pow_tau;
    
    for ( int k = j ; k < _nc ; ++k )
        if ( k != hat_i )
        {
            (*dir)[k] = NOMAD::RNG::rand()%i_pow_tau;
            if ( NOMAD::RNG::rand()%2 && (*dir)[k] > 0.0 )
                (*dir)[k] = -(*dir)[k];
        }
    
#ifdef DEBUG
    if ( j==0 )
        _out << "new LT-MADS b(l) direction: b(" << mesh_index << ") = "
        << *dir << std::endl << std::endl;
#endif
}

/*---------------------------------------------------------*/
/*      permute the coordinates of a direction (private)   */
/*---------------------------------------------------------*/
void NOMAD::Directions::permute_coords ( NOMAD::Direction & dir                ,
                                        const int        * permutation_vector   ) const
{
    NOMAD::Point tmp = dir;
    for ( int i = 0 ; i < _nc ; ++i )
        dir [ permutation_vector[i] ] = tmp[i];
}


/*---------------------------------------------------------*/
/*                           display                       */
/*---------------------------------------------------------*/
void NOMAD::Directions::display ( const NOMAD::Display & out ) const
{
    out << "n             : " << _nc << std::endl
    << "types         : { ";
    std::set<NOMAD::direction_type>::const_iterator it , end = _direction_types.end();
    for ( it = _direction_types.begin() ; it != end ; ++it )
        out << "[" << *it << "] ";
    out << "}" << std::endl
    << "sec poll types: { ";
    end = _sec_poll_dir_types.end();
    for ( it = _sec_poll_dir_types.begin() ; it != end ; ++it )
        out << "[" << *it << "] ";
    out << "}" << std::endl;
    if ( _is_orthomads )
    {
        out << "seed          : " << NOMAD::RNG::get_seed() << endl;
    }
}

/*---------------------------------------------------------*/
/*                     comparison operator                 */
/*---------------------------------------------------------*/
bool NOMAD::Directions::operator < ( const NOMAD::Directions & d ) const
{
    // number of variables:
    if ( _nc < d._nc )
        return true;
    if ( d._nc < _nc )
        return false;
    
    // boolean variables:
    if ( _is_binary != d._is_binary )
        return _is_binary;
    if ( _is_categorical != d._is_categorical )
        return _is_categorical;
    if ( _is_orthomads != d._is_orthomads )
        return _is_orthomads;
    
    // direction types:
    size_t nd = _direction_types.size();
    if ( nd < d._direction_types.size() )
        return true;
    if ( d._direction_types.size() < nd )
        return false;
    
    size_t ns = _sec_poll_dir_types.size();
    if ( ns < d._sec_poll_dir_types.size() )
        return true;
    if ( d._sec_poll_dir_types.size() < ns )
        return false;
    
    std::set<NOMAD::direction_type>::const_iterator
    it1 = _direction_types.begin()   ,
    it2 = d._direction_types.begin() ,
    end = _direction_types.end();
    
    while ( it1 != end ) {
        if ( *it1 < *it2 )
            return true;
        if ( *it2 < *it1 )
            return false;
        ++it1;
        ++it2;
    }
    
    it1 = _sec_poll_dir_types.begin();
    it2 = d._sec_poll_dir_types.begin();
    end = _sec_poll_dir_types.end();
    
    while ( it1 != end ) {
        if ( *it1 < *it2 )
            return true;
        if ( *it2 < *it1 )
            return false;
        ++it1;
        ++it2;
    }
    
    return false;
}

