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
  \file   VNS_Search.hpp
  \brief  VNS search (headers)
  \author Sebastien Le Digabel
  \date   2010-04-12
  \see    VNS_Search.cpp
*/
#ifndef __VNS_SEARCH__
#define __VNS_SEARCH__

#include "Mads.hpp"
#include "Search.hpp"

namespace NOMAD {

  /// Variable Neighborhood Search (VNS) search.
  class VNS_Search : public NOMAD::Search , private NOMAD::Uncopyable {

  private:

    int                           _k; ///< VNS neighborhood parameter.
    int                       _k_max; ///< Maximum value of \c _k.
    const NOMAD::Eval_Point * _old_x; ///< Previous reference point (updates \c _k).
    

  public:

    /// Constructor.
    /**
       \param p Parameters -- \b IN.
    */
    VNS_Search ( NOMAD::Parameters & p )
      : NOMAD::Search ( p , NOMAD::VNS_SEARCH     ) ,
	_k            ( 1                         ) ,
	_k_max        ( 1                         ) ,
	_old_x        ( NULL                      )   {}
    
    /// Destructor.
    virtual ~VNS_Search ( void ) {}

    /// Reset.
    virtual void reset ( void );

    /// The VNS search.
    /**
       Principle: x --[shaking(k)]--> x' --[descent]--> x'' .
       \param mads           NOMAD::Mads object invoking this search -- \b IN/OUT.
       \param nb_search_pts  Number of generated search points       -- \b OUT.
       \param stop           Stop flag                               -- \b IN/OUT.
       \param stop_reason    Stop reason                             -- \b OUT.
       \param success        Type of success                         -- \b OUT.
       \param count_search   Count or not the search                 -- \b OUT.
       \param new_feas_inc   New feasible incumbent                  -- \b IN/OUT.
       \param new_infeas_inc New infeasible incumbent                -- \b IN/OUT.
    */
    virtual void search ( NOMAD::Mads              & mads           ,
			  int                      & nb_search_pts  ,
			  bool                     & stop           ,
			  NOMAD::stop_type         & stop_reason    ,
			  NOMAD::success_type      & success        ,
			  bool                     & count_search   ,
			  const NOMAD::Eval_Point *& new_feas_inc   ,
			  const NOMAD::Eval_Point *& new_infeas_inc   );

  };
}

#endif
