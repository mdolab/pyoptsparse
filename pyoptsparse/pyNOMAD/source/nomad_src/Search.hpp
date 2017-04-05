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
  \file   Search.hpp
  \brief  Generic class for search strategies (headers)
  \author Sebastien Le Digabel
  \date   2010-04-08
*/
#ifndef __SEARCH__
#define __SEARCH__

#include "Evaluator_Control.hpp"

namespace NOMAD {
  
  // Forward declarations.
  class Mads;

  /// Generic class for search strategies.
  /**
     This is an abstract class (it is not possible to create NOMAD::Search objects).
  */
  class Search {

  protected:

    NOMAD::Parameters  & _p;    ///< Parameters.
    NOMAD::search_type   _type; ///< Search type.

  public:

    /// Constructor.
    /**
      \param p Parameters  -- \b IN.
      \param t Search type -- \b IN.
    */
    Search ( NOMAD::Parameters  & p ,
	     NOMAD::search_type   t   )
      : _p    ( p ) ,
	_type ( t ) {}

    /// Destructor.
    virtual ~Search ( void ) {}

    /// The search.
    /**
      - Has to be implemented by every NOMAD::Search subclass.
      - Pure virtual method.
      \param mads           NOMAD::Mads object invoking this search -- \b IN/OUT.
      \param nb_search_pts  Number of generated search points       -- \b OUT.
      \param stop           Stop flag                               -- \b IN/OUT.
      \param stop_reason    Stop reason                             -- \b OUT.
      \param success        Type of success                         -- \b OUT.
      \param count_search   Count or not the search                 -- \b OUT.
      \param new_feas_inc   New feasible incumbent                  -- \b IN/OUT.
      \param new_infeas_inc New infeasible incumbent                -- \b IN/OUT.
    */
    virtual void search
    ( NOMAD::Mads              & mads           ,
      int                      & nb_search_pts  ,
      bool                     & stop           ,
      NOMAD::stop_type         & stop_reason    ,
      NOMAD::success_type      & success        ,
      bool                     & count_search   ,
      const NOMAD::Eval_Point *& new_feas_inc   ,
      const NOMAD::Eval_Point *& new_infeas_inc   ) = 0;

    /// Reset.
    virtual void reset ( void ) {}

    /// Display.
    /**
       \param out The NOMAD::Display object -- \b IN.
    */
    virtual void display ( const NOMAD::Display & out ) const {}

  };

  /// Display a NOMAD::Search object.
  /**
     \param out The NOMAD::Display object                -- \b IN.
     \param s   The NOMAD::Search object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
  */
  inline const NOMAD::Display & operator << ( const NOMAD::Display & out ,
 					      const NOMAD::Search  & s     ) {
    s.display ( out );
    return out;
  }

}
#endif
