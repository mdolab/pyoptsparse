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
  \file   Cache_File_Point.hpp
  \brief  Class for points in binary files (headers)
  \author Sebastien Le Digabel
  \date   2010-04-06
  \see    Cache_File_Point.cpp
*/

#ifndef __CACHE_FILE_POINT__
#define __CACHE_FILE_POINT__

#include "Display.hpp"
#include "Uncopyable.hpp"

namespace NOMAD {
  
  // forward declarations:
  class Point;
  class Eval_Point;

  /// Class to represent NOMAD::Eval_Point objects in binary files.
  /**
   - All point coordinates are saved.
   - Only defined blackbox outputs are saved.
   - To get all bb_outputs:
   \code
   NOMAD::Point bbo ( _m );
   for ( int i = 0 ; i < _m_def ; ++i )
     bbo [ _bbo_index[i] ] = _bbo_def[i];
   \endcode
  */
  class Cache_File_Point : private NOMAD::Uncopyable {

  private:
    
#ifdef MEMORY_DEBUG
    /// Number of NOMAD::Cache_File_Point objects in memory.
    static int _cardinality;

    /// Max number of NOMAD::Cache_File_Point objects in memory.
    static int _max_cardinality;
#endif

    int           _n;      ///< Dimension of the point.
    int           _m;      ///< Number of blackbox outputs (both defined and undefined).
    int           _m_def;  ///< Number of defined blackbox outputs.

    /// Evaluation status.
    /**
	- 0: fail.
	- 1: ok.
	- 2: in progress.
	- 3: undefined.
    */
    unsigned char _eval_status;

    double      * _coords;     ///< The \c _n coordinates.
    double      * _bbo_def;    ///< The \c _m_def defined blackbox output values.
    int         * _bbo_index;  ///< The index for the blackbox output values.

    /// Reset.
    void reset ( void );

  public:

#ifdef MEMORY_DEBUG
    /// Access to the number of NOMAD::Cache_File_Point objects in memory.
    /**
       \return The number of NOMAD::Cache_File_Point objects in memory.
    */
    static int get_cardinality ( void ) { return Cache_File_Point::_cardinality; }

    /// Access to the max number of NOMAD::Cache_File_Point objects in memory.
    /**
       \return The max number of NOMAD::Cache_File_Point objects in memory.
    */
    static int get_max_cardinality ( void )
    {
      return Cache_File_Point::_max_cardinality;
    }
#endif

    /// Constructor #1.
    explicit Cache_File_Point ( void );
    
    /// Constructor #2.
    /**
       From a NOMAD::Eval_Point object.
       \param x The evaluation point.
    */
    explicit Cache_File_Point ( const NOMAD::Eval_Point & x );
    
    /// Destructor.
    virtual ~Cache_File_Point ( void );

    /// Access to the dimension of the point.
    /**
       \return The dimension of the point.
    */    
    int get_n ( void  ) const { return _n; }

    /// Access to the number of blackbox outputs.
    /**
       \return The number of blackbox outputs.
    */  
    int get_m ( void  ) const { return _m; }

    /// Access to the evaluation status.
    /**
       Evaluation status as a character:
       - 0: fail.
       - 1: ok.
       - 2: in progress.
       - 3: undefined.
       \return A character as the evaluation status.
    */  
    unsigned char get_eval_status ( void ) const { return _eval_status; }

    /// Access to the coordinates.
    /**
       \param i The index (0 for the first element) -- \b IN.
       \return The \c (i+1)th coordinate.
    */
    double get_coord ( int i ) const;

    /// Access to the blackbox outputs.
    /**
       \return A NOMAD::Point as the \c _m outputs.
    */
    const NOMAD::Point get_bb_outputs ( void ) const;
    
    /// Write in a binary file.
    /**
       \param fout The output file -- \b IN/OUT.
       \return A boolean equal to \c true if the file could be written.
    */
    bool write ( std::ofstream & fout ) const;
    
    /// Read in a binary file.
    /**
       \param fin The input file -- \b IN/OUT.
       \return A boolean equal to \c true if the file could be read.
    */
    bool read ( std::ifstream & fin );
    
    /// Display.
    /**
       \param out The NOMAD::Display object -- \b IN.
    */
    void display ( const NOMAD::Display & out ) const;
  };

  /// Display a NOMAD::Cache_File_Point object.
  /**
     \param out The NOMAD::Display object -- \b IN.
     \param x   The NOMAD::Cache_File_Point object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
  */
  inline const NOMAD::Display & operator << ( const NOMAD::Display          & out ,
					      const NOMAD::Cache_File_Point & x     )
  {
    x.display ( out );
    return out;
  }
}

#endif
