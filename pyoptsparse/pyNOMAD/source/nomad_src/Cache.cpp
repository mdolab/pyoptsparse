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
 \file   Cache.cpp
 \brief  Cache memorizing all evaluations (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-12
 \see    Cache.hpp
 */
#include "Cache.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
std::set<std::string> NOMAD::Cache::_locked_files;

/*------------------------------------------------------*/
/*                init of _sizeof (private)             */
/*------------------------------------------------------*/
int NOMAD::Cache::sizeof_init ( void ) const
{
    return
    sizeof ( _cache1     ) +
    sizeof ( _cache2     ) +
    sizeof ( _cache3     ) +
    sizeof ( _locked_file) +
    sizeof ( _it         ) +
    sizeof ( _out        ) +
    sizeof ( _sizeof     );
}

/*---------------------------------------------------------------------*/
/*  insertion of a point into the list of extern points (_extern_pts)  */
/*---------------------------------------------------------------------*/
void NOMAD::Cache::insert_extern_point ( const NOMAD::Eval_Point & x ) const
{
    if ( !x.get_current_run() )
        _extern_pts.push_front ( &x );
}

/*----------------------------------------------------------*/
/*  get the first extern point and remove it from the list  */
/*  (returns NULL if there is no extern point)              */
/*----------------------------------------------------------*/
const NOMAD::Eval_Point * NOMAD::Cache::get_and_remove_extern_point ( void ) const
{
    if ( _extern_pts.empty() )
        return NULL;
    const NOMAD::Eval_Point * extern_point = *_extern_pts.begin();
    _extern_pts.pop_front();
    return extern_point;
}

/*------------------------------------------*/
/*  . search x in the cache                 */
/*  . return NULL if x is not in the cache  */
/*------------------------------------------*/
const NOMAD::Eval_Point * NOMAD::Cache::find ( const NOMAD::Eval_Point & x ) const
{
    // check the eval types:
    if ( x.get_eval_type() != _eval_type )
        throw NOMAD::Cache::Cache_Error ( "Cache.cpp" , __LINE__ ,
                                         "NOMAD::Cache:find(x): x.eval_type != cache.eval_type" );
    
    // find:
    std::set<NOMAD::Cache_Point>::const_iterator it;
    NOMAD::cache_index_type                      cache_index;
    const NOMAD::Eval_Point                    * cache_x
    = NOMAD::Cache::find ( x , it , cache_index );
    
    return cache_x;
}

/*----------------------------------------------------*/
/*  . search x in the cache                           */
/*  . return NULL if x is not in the cache            */
/*  . this version returns an iterator and indicates  */
/*    in which set the point has been found           */
/*  . private method                                  */
/*----------------------------------------------------*/
const NOMAD::Eval_Point * NOMAD::Cache::find
( const NOMAD::Eval_Point                      &  x          ,
 std::set<NOMAD::Cache_Point>::const_iterator & it          ,
 NOMAD::cache_index_type                      & cache_index   ) const
{
    // search in _cache2 (points to write in a cache file):
    NOMAD::Cache_Point cp ( &x );
    it = _cache2.find ( cp );
    if ( it != _cache2.end() ) {
        cache_index = NOMAD::CACHE_2;
        return it->get_point();
    }
    
    // search in _cache3 (points saved in a file):
    it = _cache3.find ( cp );
    if ( it != _cache3.end() ) {
        cache_index = NOMAD::CACHE_3;
        return it->get_point();
    }
    
    // search in _cache1 (points read in an initial cache file):
    it = _cache1.find ( cp );
    if ( it != _cache1.end() ) {
        cache_index = NOMAD::CACHE_1;
        return it->get_point();
    }
    
    cache_index = NOMAD::UNDEFINED_CACHE;
    return NULL;
}

/*-----------------------------------------------------*/
/*                     erase a point                   */
/*-----------------------------------------------------*/
/*  . return true if the point has been found and      */
/*      removed from the cache                         */
/*  . the point is not deleted from memory if its      */
/*    address does not match the address of the point  */
/*    in cache                                         */
/*-----------------------------------------------------*/
bool NOMAD::Cache::erase ( const NOMAD::Eval_Point & x )
{
    // check the eval types:
    if ( x.get_eval_type() != _eval_type )
        throw NOMAD::Cache::Cache_Error ( "Cache.cpp" , __LINE__ ,
                                         "NOMAD::Cache:erase(x): x.eval_type != cache.eval_type" );
    
    std::set<NOMAD::Cache_Point>::iterator it;
    NOMAD::cache_index_type       cache_index;
    
    // search in cache:
    const NOMAD::Eval_Point * cache_x = find ( x , it , cache_index );
    
    // the point has been found:
    if ( cache_x ) {
        
        // remove the point from the list of extern points:
        if ( cache_x->get_current_run() || x.get_current_run() ) {
            std::list<const NOMAD::Eval_Point*>::const_iterator end2 = _extern_pts.end();
            std::list<const NOMAD::Eval_Point*>::iterator       it2  = _extern_pts.begin();
            while ( it2 != end2 ) {
                if ( *it2 == cache_x || *it2 == &x ) {
                    _extern_pts.erase ( it2 );
                    break;
                }
                ++it2;
            }
        }
        
        // erase the point in cache if its address is different from &x:
        if ( cache_x != &x )
            delete cache_x;
        
        // remove the point from the cache:
        _sizeof -= x.size_of();
        
        switch ( cache_index ) {
            case NOMAD::CACHE_1:
                _cache1.erase ( it );
                break;
            case NOMAD::CACHE_2:
                _cache2.erase ( it );
                break;
            case NOMAD::CACHE_3:
                _cache3.erase ( it );
                break;
            case NOMAD::UNDEFINED_CACHE:
                break;
        }
        return true;
    }
    return false;
}

/*-----------------------------------------------------*/
/*  erase all points in cache and unlock _locked_file  */
/*-----------------------------------------------------*/
void NOMAD::Cache::clear ( void )
{
    const NOMAD::Eval_Point * x = begin();
    while ( x ) {
        delete x;
        x = next();
    }
    _cache1.clear();
    _cache2.clear();
    _cache3.clear();
    unlock();
    
    _extern_pts.clear();
    
    _sizeof = static_cast<float> ( sizeof_init() );
}

/*------------------------------------------------*/
/* . insertion of a point in the cache (_cache2)  */
/* . supposes that x is NOT already in cache      */
/*------------------------------------------------*/
void NOMAD::Cache::insert ( const NOMAD::Eval_Point & x )
{
    // check the eval types:
    if ( x.get_eval_type() != _eval_type )
        throw NOMAD::Cache::Cache_Error ( "Cache.cpp" , __LINE__ ,
                                         "NOMAD::Cache:insert(x): x.eval_type != cache.eval_type" );
    
    // insertion in _extern_pts:
    insert_extern_point ( x );
    
    // insertion in _cache2:
    NOMAD::Cache_Point cp  ( &x );
    _cache2.insert  (  cp  );
    x.set_in_cache ( true );
    _sizeof += x.size_of();
}

/*-------------------------------------------------------------------------*/
/*  . insert all points of 'c' in the current cache                        */
/*  . this empties 'c' points (to avoid Eval_Point copies)                 */
/*  . c._locked_file and this->_locked_file are different by construction  */
/*-------------------------------------------------------------------------*/
void NOMAD::Cache::insert ( Cache & c )
{
    if ( &c == this )
        return;
    
    // check the eval types:
    if ( c._eval_type != _eval_type )
        throw NOMAD::Cache::Cache_Error ( "Cache.cpp" , __LINE__ ,
                                         "NOMAD::Cache:insert(c): c._eval_type != this->_eval_type" );
    
    // insertion:
    NOMAD::Point              bbo_cache , bbo_cur;
    const NOMAD::Eval_Point * cache_x;
    const NOMAD::Eval_Point * cur = c.begin();
    
    while ( cur ) {
        
        cache_x = find ( *cur );
        
        // the current point is already in cache:
        if ( cache_x ) {
            update ( get_modifiable_point ( *cache_x ) , *cur );
            delete cur;
        }
        
        // point not in cache:
        else
            insert ( *cur );
        
        cur = c.next();
    }
    
    c._sizeof = static_cast<float> ( sizeof_init() );
    
    c._cache1.clear();
    c._cache2.clear();
    c._cache3.clear();
    c._extern_pts.clear();
}

/*------------------------------------------------------------------*/
/*  . begin() and next() methods, to browse the cache               */
/*                                                                  */
/*  . example of use:                                               */
/*                                                                  */
/*        const NOMAD::Eval_Point * cur = cache.begin();            */
/*        while ( cur ) {                                           */
/*          ...                                                     */
/*          cur = cache.next();                                     */
/*        }                                                         */
/*                                                                  */
/*  . we browse in this order: _cache2, _cache3, and _cache_1       */
/*------------------------------------------------------------------*/

// begin():
// --------
const NOMAD::Eval_Point * NOMAD::Cache::begin ( void ) const
{
    if ( !_cache2.empty() ) {
        _it = _cache2.begin();
        return _it->get_point();
    }
    if ( !_cache3.empty() ) {
        _it = _cache3.begin();
        return _it->get_point();
    }
    if ( !_cache1.empty() ) {
        _it = _cache1.begin();
        return _it->get_point();
    }
    return NULL;
}

// next() (supposes that begin() has been called)
// -------
const NOMAD::Eval_Point * NOMAD::Cache::next ( void ) const
{
    ++_it;
    
    if ( !_cache2.empty() && _it == _cache2.end() ) {
        if ( !_cache3.empty() ) {
            _it = _cache3.begin();
            return _it->get_point();
        }
        if ( !_cache1.empty() ) {
            _it = _cache1.begin();
            return _it->get_point();
        }
        return NULL;
    }
    
    
    if ( !_cache3.empty() && _it == _cache3.end() ) {
        if ( !_cache1.empty() ) {
            _it = _cache1.begin();
            return _it->get_point();
        }
        return NULL;
    }
    
    if ( !_cache1.empty() && _it == _cache1.end() )
        return NULL;
    
    return _it->get_point();
}

/*---------------------------------------------------------------------*/
/*                      check if a file is locked (private)            */
/*---------------------------------------------------------------------*/
bool NOMAD::Cache::is_locked ( const std::string & file_name )
{
    if ( file_name == _locked_file )
        return true;
    return ( Cache::_locked_files.find ( file_name ) != Cache::_locked_files.end() );
}

/*---------------------------------------------------------------------*/
/*                          lock a file (private)                      */
/*---------------------------------------------------------------------*/
bool NOMAD::Cache::lock ( const std::string & file_name )
{
    if ( is_locked ( file_name ) )
        return false;
    
    Cache::_locked_files.insert ( file_name );
    _locked_file = file_name;
    
    return true;
}

/*---------------------------------------------------------------------*/
/*                      unlock the locked file (private)               */
/*---------------------------------------------------------------------*/
void NOMAD::Cache::unlock ( void )
{
    if ( _locked_file.empty() )
        return;
    
    std::set<std::string>::iterator it = Cache::_locked_files.find ( _locked_file );
    if ( it != Cache::_locked_files.end() )
        _locked_files.erase(it);
    
    _locked_file.clear();
}

/*----------------------------------------------------*/
/*  . reads points from a cache file                  */
/*  . fills the set _cache1                           */
/*  . private method                                  */
/*----------------------------------------------------*/
bool NOMAD::Cache::read_points_from_cache_file ( std::ifstream & fin             ,
                                                const int     * p_nb_bb_outputs ,
                                                bool            display           )
{
    try {
        
        NOMAD::Clock c;
        
        // the stream is placed at the first point (after the CACHE_FILE_ID tag):
        fin.seekg ( sizeof ( NOMAD::CACHE_FILE_ID ) , std::ios::beg );
        
        NOMAD::Cache_File_Point   cfp;
        NOMAD::Eval_Point       * cur;
        const NOMAD::Eval_Point * cache_x;
        
        // main loop:
        while ( !fin.eof() ) {
            
            // reading of the Cache_File_Point:
            if ( !cfp.read ( fin ) ) {
                if ( fin.eof() )
                    break;
                return false;
            }
            
            // we ignore this cache file point if it has a different
            // number of blackbox outputs than *p_nb_bb_outputs:
            if ( p_nb_bb_outputs && cfp.get_m() != *p_nb_bb_outputs )
                continue;
            
            // creation of the Eval_Point:
            cur = new NOMAD::Eval_Point ( cfp , _eval_type );
            
            // we look if the current point is already in cache:
            cache_x = find ( *cur );
            
            // the current point is already in cache:
            if ( cache_x ) {
                update ( get_modifiable_point ( *cache_x ) , *cur );
                delete cur;
            }
            
            // point not in cache: insertion:
            else {
                
                // insertion in _extern_pts:
                insert_extern_point ( *cur );
                
                // insertion in _cache1:
                NOMAD::Cache_Point cp ( cur );
                _cache1.insert    ( cp   );
                cur->set_in_cache ( true );
                _sizeof += cur->size_of();
            }
            
        } // end of main loop
        
        // display stats on the cache load:
        if ( display ) {
            _out << "number of points: " << static_cast<int>(_cache1.size()) << std::endl
            << "size            : ";
            _out.display_size_of ( _sizeof );
            _out << std::endl
            << "load time       : " << c.get_real_time() << 's' << std::endl;
        }
    }
    catch ( ... ) {
        return false;
    }
    return true;
}

/*---------------------------------------------------------------------*/
/* . load a cache file (fill the set _cache1)                          */
/* . lock the file 'file_name'                                         */
/* . return false if 'file_name' exists and is not a valid cache file  */
/* . return false if 'file_name' is already locked                     */
/* . return false if the object already locked a file                  */
/* . return true if this object is already locked with this file       */
/*   (and just do nothing)                                             */
/* . return true and create a valid cache file if 'file_name'          */
/*   does not exist                                                    */
/* . return false if the previous step did not work                    */
/*                                                                     */
/* . parameter p_nb_bb_outputs (default=NULL): to indicate a number of */
/*   blackbox outputs; points in file with a different value are       */
/*   ignored                                                           */
/*---------------------------------------------------------------------*/
bool NOMAD::Cache::load ( const std::string & file_name       ,
                         const int         * p_nb_bb_outputs ,
                         bool                display           )
{
    if ( !file_name.empty() && file_name == _locked_file )
        return true;
    
    if ( file_name.empty() || !_locked_file.empty() || is_locked(file_name) )
        return false;
    
    // the file exists:
    if ( NOMAD::check_read_file ( file_name ) ) {
        
        int           id;
        std::ifstream fin ( file_name.c_str() , std::ios::binary );
        
        fin.read ( (char *) &id , sizeof(int) );
        
        // it is a valid cache file:
        if ( !fin.fail() && id == NOMAD::CACHE_FILE_ID ) {
            
            // display:
            if ( display )
                _out << std::endl
                << NOMAD::open_block ( "loading of \'" + file_name + "\'" );
            
            // read the points:
            if ( !read_points_from_cache_file ( fin , p_nb_bb_outputs , display ) ) {
                fin.close();
                return false;  // it is not a valid cache file
            }
            
            // lock the file:
            lock ( file_name );
            
            fin.close();
            
            if ( display )
                _out.close_block();
            
            return true;
        }
        
        // it is not a valid cache file:
        else {
            fin.close();
            return false;
        }
    }
    
    // the file does not exist:
    else {
        
        // display:
        if ( display )
            _out << std::endl << "creating cache file \'" << file_name << "\'" << std::endl;
        
        // create the file as a valid cache file:
        std::ofstream fout ( file_name.c_str() , std::ios::binary );
        
        if ( fout.fail() ) {
            fout.close();
            return false;
        }
        
        fout.write ( (char *) &NOMAD::CACHE_FILE_ID , sizeof ( NOMAD::CACHE_FILE_ID ) );
        fout.close();
        
        // lock:
        lock ( file_name );
    }
    
    return true;
}

/*------------------------------------------------------------------*/
/*  if overwrite == false:                                          */
/*         . write points at the end of the locked cache file       */
/*         . write only points from _cache2                         */
/*         . transfer points from _cache2 to _cache3                */
/*  else:                                                           */
/*         . write all points in the locked cache file              */
/*------------------------------------------------------------------*/
bool NOMAD::Cache::save ( bool overwrite , bool display )
{
    if ( _locked_file.empty() )
        return true;
    
    // display:
    if ( display )
        _out << std::endl << "saving cache file \'" << _locked_file << "\'" << std::endl;
    
    std::ofstream fout;
    
    if ( overwrite ) {
        
        // open:
        fout.open ( _locked_file.c_str() , std::ios::binary );
        if ( fout.fail() ) {
            fout.close();
            return false;
        }
        
        // cache file tag:
        fout.write ( (char *) &NOMAD::CACHE_FILE_ID , sizeof ( NOMAD::CACHE_FILE_ID ) );
        
        // save all cache points:
        const NOMAD::Eval_Point * cur = begin();
        while ( cur ) {
            NOMAD::Cache_File_Point cfp ( *cur );
            if ( !cfp.write ( fout ) ) {
                fout.close();
                return false;
            }
            cur = next();
        }
    }
    
    else {
        
        // open and go at the end of the file:
        fout.open ( _locked_file.c_str() , std::ios::binary | std::ios::app );
        if ( fout.fail() ) {
            fout.close();
            return false;
        }
        
        std::set<NOMAD::Cache_Point>::iterator it = _cache2.begin();
        while ( it != _cache2.end() ) {
            
            // write it->get_point() in 'fout':
            NOMAD::Cache_File_Point cfp ( *it->get_point() );
            if ( !cfp.write ( fout ) ) {
                fout.close();
                return false;
            }
            
            // transfer the point from _cache2 to _cache3:
            NOMAD::Cache_Point cp = *it; // ( it->get_point() );
            _cache3.insert ( cp   );
            _cache2.erase  ( it++ );
        }
    }
    
    // close the file:
    fout.close();
    
    return true;
}

/*-----------------------------------------------------------------*/
/*  . update a point already in cache from another point with the  */
/*    same coordinates                                             */
/*  . if both points have a different number of blackbox outputs,  */
/*    they are not comparable and we set cache_x = x               */
/*  . private method                                               */
/*-----------------------------------------------------------------*/
void NOMAD::Cache::update ( NOMAD::Eval_Point       & cache_x ,
                           const NOMAD::Eval_Point & x         ) const
{
    const NOMAD::Point & bbo_x = x.get_bb_outputs();
    
    if ( &cache_x == &x         ||
        !x.is_eval_ok()        ||
        !cache_x.is_in_cache() ||
        bbo_x.empty()          ||
        cache_x != x              )
        return;
    
    // check the eval types:
    if ( x.get_eval_type      () != _eval_type ||
        cache_x.get_eval_type() != _eval_type    )
        throw NOMAD::Cache::Cache_Error ( "Cache.cpp" , __LINE__ ,
                                         "NOMAD::Cache:update(): problem with the eval. types" );
    
    const NOMAD::Point & bbo_cache_x = cache_x.get_bb_outputs();
    int                  m           = bbo_cache_x.size();
    
    _sizeof -= cache_x.size_of();
    
    // if the current point could not evaluate and x did,
    // or if the number of bb outputs is different, we set cache_x = x:
    if ( !cache_x.is_eval_ok() || m != bbo_x.size() )
    {
        cache_x.set_eval_status ( NOMAD::EVAL_OK     );
        cache_x.set_bb_output   ( bbo_x              );
        cache_x.set_signature   ( x.get_signature () );
        cache_x.set_direction   ( x.get_direction () );
        
        _sizeof += cache_x.size_of();
        return;
    }
    
    // we complete _bb_outputs:
    int c1 = 0;
    int c2 = 0;
    
    for ( int i = 0 ; i < m ; ++i ) {
        
        if ( bbo_cache_x[i].is_defined() )
            ++c1;
        
        if ( bbo_x[i].is_defined() )
            ++c2;
        
        if ( !bbo_cache_x[i].is_defined() && bbo_x[i].is_defined() )
            cache_x.set_bb_output ( i , bbo_x[i] );
    }
    
    // the two points are 'eval_ok' and have comparable outputs:
    // we select the best as the one with the more defined bb_outputs:
    if ( c2 > c1 ) {
        cache_x.set_signature  ( x.get_signature () );
        cache_x.set_direction  ( x.get_direction () );
        
    }
    
    _sizeof += cache_x.size_of();
}

/*---------------------------------------------------------*/
/*              display the list of extern points          */
/*---------------------------------------------------------*/
void NOMAD::Cache::display_extern_pts ( const NOMAD::Display & out ) const
{
    int  nb = static_cast<int>(_extern_pts.size());
    int cnt = 0;
    std::list<const NOMAD::Eval_Point *>::const_iterator it , end = _extern_pts.end();
    for ( it = _extern_pts.begin() ; it != end ; ++it ) {
        out << "point ";
        out.display_int_w ( ++cnt , nb );
        out << "/" << nb << ": ";
        (*it)->display_eval ( out , false );
        out << std::endl;
    }
}

/*---------------------------------------------------------*/
/*                         display                         */
/*---------------------------------------------------------*/
void NOMAD::Cache::display ( const NOMAD::Display & out ) const
{
    out << "number of cache points: " << size() << std::endl
    << "size in memory        : ";
    out.display_size_of ( _sizeof );
    out << std::endl << "cache file            : ";
    if ( _locked_file.empty() )
        out << "-" << std::endl;
    else
        out << _locked_file << std::endl;
    
#ifdef DEBUG
    int  nb = size();
    int cnt = 0;
    out << NOMAD::open_block ( "cache points" ) << std::endl;
    const NOMAD::Eval_Point * cur = begin();
    while ( cur ) {
        out << "point ";
        out.display_int_w ( ++cnt , nb );
        out << "/" << nb << ": ";
        cur->display_eval ( out , false );
        out << std::endl;
        cur = next();
    }
    out.close_block();
    
    // this display is repeated here as there can
    // be a lot of points displayed just before:
    out << "number of cache points: " << size() << std::endl
    << "size in memory        : ";
    out.display_size_of ( _sizeof );
    out << std::endl << "cache file            : ";
    if ( _locked_file.empty() )
        out << "-";
    else
        out << _locked_file;
    out << std::endl;
#endif
}
