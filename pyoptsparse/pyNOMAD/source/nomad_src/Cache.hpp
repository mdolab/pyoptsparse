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
  \file   Cache.hpp
  \brief  Cache memorizing all evaluations (headers)
  \author Sebastien Le Digabel
  \date   2010-04-12
  \see    Cache.cpp
*/

//  To browse the cache:
//  ---------------------
//    const NOMAD::Eval_Point * cur = cache.begin();
//    while ( cur ) {
//      ...
//      cur = cache.next();
//    }

#ifndef __CACHE__
#define __CACHE__

#include "Cache_Point.hpp"
#include "Clock.hpp"

namespace NOMAD {

  /// Cache memorizing all evaluations.
  /**
     - The evaluation are stored as NOMAD::Eval_Point objects in
       three \c std::set containers
       ( \c _cache1 , \c _cache2 and \c _cache3 ).
     - Browse the points inside the cache with the following instructions:
       \code
       const NOMAD::Eval_Point * cur = cache.begin();
       while ( cur ) {
         ...
         cur = cache.next();
       }
       \endcode
  */
  class Cache : private NOMAD::Uncopyable {

  protected:

    /// Display:
    const NOMAD::Display & _out;

  private:
  
    /// Static list of locked files.
    static std::set<std::string> _locked_files;

    /// Cache file that is currently in use for load/save (this file is locked).
    std::string _locked_file;
    
    /// Type of cache (truth or surrogate).
    NOMAD::eval_type _eval_type;

    /// Points read in the (optional) initial cache file.
    std::set<NOMAD::Cache_Point> _cache1;

    /// Points of the cache to write in the cache file.
    std::set<NOMAD::Cache_Point> _cache2;

    /// Points of the cache already saved in the file.
    std::set<NOMAD::Cache_Point> _cache3;

    /// Points not from the current run (from a cache file or from a user).
    mutable std::list<const NOMAD::Eval_Point*> _extern_pts;
    
    /// Size in memory.
    mutable float _sizeof;

    /// Iterator to browse the cache with \c begin() and \c next().
    mutable std::set<NOMAD::Cache_Point>::const_iterator _it;
    
    /*---------------------------------------------------------------------------*/

    /// Read points in a cache file.
    /**
       - Reads points from a cache file.
       - Fills the set \c _cache1.
       \param fin             The already opened input file -- \b IN.
       \param p_nb_bb_outputs A pointer on an integer with the number
                              of blackbox outputs; can be \c NULL  -- \b IN.
       \param display         A boolean equal to \c true if displays
                              are authorized -- \b IN.
       \return A boolean equal to \c true if the file could be read.
    */
    bool read_points_from_cache_file ( std::ifstream & fin             ,
				       const int     * p_nb_bb_outputs ,
				       bool            display           );
    
    /// Check if a file is locked.
    /**
       \param file_name Name of the file -- \b IN.
       \return A boolean equal to \c true if the file is locked.
    */
    bool is_locked ( const std::string & file_name );

    /// Lock a file.
    /**
       \param file_name Name of the file -- \b IN.
       \return A boolean equal to \c false if the file was already locked.
    */
    bool lock ( const std::string & file_name );

    /// Unlock the locked file.
    void unlock ( void );

    /// Update a point already in cache.
    /**
       - The point already in cache is updated from another
         point with the same coordinates.
       - If both points have a different number of blackbox outputs,
         they are not comparable and we set \c cache_x \c = \c x .
       \param cache_x The point already in cache -- \b IN/OUT.
       \param x       The other point            -- \b IN.
    */
    void update ( NOMAD::Eval_Point & cache_x , const NOMAD::Eval_Point & x ) const;
    
    /// Find a point in the cache.
    /**
       \param x           The point -- \b IN.
       \param it          An iterator indicating the position of the point
                          -- \b OUT.
       \param cache_index Index of the std::set in which the point has been found
                          -- \b OUT.
       \return A pointer to the point in cache.
       \return \c NULL if \c x is not in the cache.
    */
    const NOMAD::Eval_Point * find
    ( const NOMAD::Eval_Point                      & x           ,
      std::set<NOMAD::Cache_Point>::const_iterator & it          ,
      NOMAD::cache_index_type                      & cache_index   ) const;
    
    /// Initialization of \c _sizeof.
    /**
       \return The size of an empty cache.
    */
    int sizeof_init ( void ) const;  

    /*---------------------------------------------------------------------------*/
    
  public:

    /// Exception class for a cache error.
    /**
       Occurs when points have not the same \c eval_type as \c this->_eval_type.
    */
    class Cache_Error : public NOMAD::Exception
    {
    public:
      /// Constructor.
      Cache_Error ( const std::string & file ,
		    int                 line ,
		    const std::string & msg    )
	: NOMAD::Exception ( file , line , msg ) {}
    };
    
    /*---------------------------------------------------------------------------*/

    /// Constructor.
    /**
       \param out  The NOMAD::Display object -- \b IN.
       \param type Type of the cache
                   -- \b IN -- \b optional (default = NOMAD::TRUTH).
    */
    explicit Cache ( const NOMAD::Display & out ,
		     NOMAD::eval_type       type = NOMAD::TRUTH )
      : _out       ( out                               ) ,
	_eval_type ( type                              ) ,
	_sizeof    ( static_cast<float>(sizeof_init()) )   {}
    
    /// Destructor.
    virtual ~Cache ( void ) { clear(); }
    
    /// Access to the size of cache in memory.
    /**
       \return The size of cache in memory, in bytes.
    */
    float size_of ( void ) const { return _sizeof; }

    /// Access to the number of points.
    /**
       \return The number of points in the cache.
    */
    int size ( void ) const
    {
      return static_cast<int> ( _cache1.size() + _cache2.size() + _cache3.size() );
    }
    
    /// check if the cache is empty.
    /**
       \return A boolean equal to \c true if the cache is empty.
    */
    bool empty ( void ) const
    {
      return _cache1.empty() && _cache2.empty() && _cache3.empty();
    }
    
    /// Access to the evaluation type (truth or surrogate).
    /**
       \return The evaluation type.
    */
    NOMAD::eval_type get_eval_type ( void ) const { return _eval_type; }

    /// Const cast for an evaluation point.
    /**
       Transforms a \c const NOMAD::Eval_Point \c & into a NOMAD::Eval_Point \c &.
       \param x The const point -- \b IN.
       \return  The non-const point.
    */
    static NOMAD::Eval_Point & get_modifiable_point ( const NOMAD::Eval_Point & x )
    {
      return const_cast<NOMAD::Eval_Point&> ( x );
    }
  
    /// Insertion of a point into the list of extern points.
    /**
       \param x The extern point -- \b IN.
    */
    void insert_extern_point  ( const NOMAD::Eval_Point & x ) const;

    /// Access to the number of extern points.
    /**
       \return The number of extern points.
    */
    virtual int get_nb_extern_points ( void ) const { return static_cast<int>(_extern_pts.size()); }

    /// Access to an extern point.
    /**
       Get the first extern point and remove it from the list.
       \return A pointer to the extern point;
               \c NULL if there is no extern point.
    */
    virtual const NOMAD::Eval_Point * get_and_remove_extern_point ( void ) const;

    /// Find a point in the cache.
    /**
       \param x The point -- \b IN.
       \return A pointer to the point in cache.
       \return \c NULL if the point is not in cache.
    */
    virtual const NOMAD::Eval_Point * find ( const NOMAD::Eval_Point & x ) const;
    
    /// Access to the first point in cache.
    /**
       We browse in this order: \c _cache2, \c _cache3, and \c _cache_1 .
       \return A pointer to the first point in cache;
               \c NULL if the cache is empty.
    */
    const NOMAD::Eval_Point * begin ( void ) const;

    /// Access to the next point when browsing the cache.
    /**
       Supposes that \c begin() has already been called.
       \return A pointer to the next point; \c NULL if there
               is no more point.
    */
    const NOMAD::Eval_Point * next ( void ) const;

    /// Erase a point.
    /**
       The point is not deleted from memory if its
       address does not match the address of the point
       in cache.
       \param x The point -- \b IN.
       \return A boolean equal to \c true
               if the point has been found and removed from the cache.
    */
    bool erase ( const NOMAD::Eval_Point & x );

    /// Insertion of a point in the cache ( \c _cache2 ).
    /**
       Supposes that \c x is not already in cache.
       \param x The point -- \b IN.
    */
    virtual void insert ( const NOMAD::Eval_Point & x );

    /// Insert all points of another cache into the current cache.
    /**
       This empties the points in \c c in order to avoid NOMAD::Eval_Point copies.
       \c c._locked_file and \c this->_locked_file are different by construction.
       \param c The other cache -- \b IN/OUT.
    */
    void insert ( Cache & c );

    /// Load a cache file.
    /**
       - fills the set \c _cache1.
       - locks the file \c file_name.

       \param file_name The file -- \b IN.

       \param p_nb_bb_outputs A pointer to the number of blackbox outputs.
                              It is ignored if equal to \c NULL.
			      Points from the file with a different number
			      of outputs are ignored
			      -- \b IN -- \b optional (default = \c NULL).

       \param display A boolean equal to \c true if displays
                      are authorized
		      -- \b IN -- \b optional (default = \c true).

       \return \c false if the file exists and is not a valid cache file.
       \return \c false if the file is already locked.
       \return \c false if the object already locked another file.
       \return \c true if the object is already locked with this file (and do nothing).
       \return \c true and create a valid cache file if the file does not exist.
       \return \c false if the previous anything else did not work.
    */
    bool load ( const std::string & file_name              ,
		const         int * p_nb_bb_outputs = NULL ,
		bool                display         = true   );

    /// Save a cache file.
    /**
       \param overwrite A boolean.
			If \c false, points are written
			  at the end of the locked file and only points
			  from \c _cache2 are considered.
			  These points are then transfered to \c _cache3.
			If \c true, all the points are considered.
		        -- \b IN -- \b optional (default = \c false).
       \param display   A boolean equal to \c true if displays
                        are authorized
		        -- \b IN -- \b optional (default = \c true).
       \return A boolean equal to \c true if the save could complete.
    */
    bool save ( bool overwrite = false ,
		bool display   = true    );
    
    /// Erase all points in cache and unlock the file.
    void clear ( void );

    /// Display the list of extern points.
    /**
       \param out The NOMAD::Display object -- \b IN.
    */
    void display_extern_pts ( const NOMAD::Display & out ) const;

    /// Display the list of extern points.
    /**
       Uses the \c this->_out member as NOMAD::Display object.
    */
    void display_extern_pts ( void ) const { display_extern_pts ( _out ); }

    /// Display.
    /**
       \param out The NOMAD::Display object -- \b IN.
    */
    void display ( const NOMAD::Display & out ) const;

    /// Display.
    /**
       Uses the \c this->_out member as NOMAD::Display object.
    */
    void display ( void ) const { display ( _out ); }

  };
  
  /// Display a NOMAD::Cache object.
  /**
     \param out The NOMAD::Display object -- \b IN.
     \param c   The NOMAD::Cache object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
  */
  inline const NOMAD::Display & operator << ( const NOMAD::Display & out ,
					      const NOMAD::Cache   & c     )
  {
    c.display ( out );
    return out;
  }
}

#endif
