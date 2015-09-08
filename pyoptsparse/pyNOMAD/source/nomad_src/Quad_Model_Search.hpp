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
  \file   Quad_Model_Search.hpp
  \brief  Quadratic Model search (headers)
  \author Sebastien Le Digabel
  \date   2010-08-30
  \see    Quad_Model_Search.cpp
*/
#ifndef __QUAD_MODEL_SEARCH__
#define __QUAD_MODEL_SEARCH__

#include "Mads.hpp"
#include "Single_Obj_Quad_Model_Evaluator.hpp"
#include "Multi_Obj_Quad_Model_Evaluator.hpp"

namespace NOMAD {

  /// Model search.
  class Quad_Model_Search : public NOMAD::Search , private NOMAD::Uncopyable {

  private:

    NOMAD::Model_Stats _one_search_stats;   ///< Stats for one search.
    NOMAD::Model_Stats _all_searches_stats; ///< Stats for all searches.

    /// Model optimization.
    /**
       \param model          The model                          -- \b IN.
       \param xk             The two model centers              -- \b IN.
       \param i_inc          Model center index (\c 0 or \c 1 ) -- \b IN.
       \param display_degree Display degree                     -- \b IN.
       \param out            The NOMAD::Display object          -- \b IN.
       \param xf             Feasible solution \c xf            -- \b OUT.
       \param xi             Infeasible solution \c xi          -- \b OUT.
       \param stop           Stop flag                          -- \b OUT.
       \param stop_reason    Stop reason                        -- \b OUT.
    */
    bool optimize_model ( const NOMAD::Quad_Model  & model          ,
			  const NOMAD::Eval_Point ** xk             ,
			  int                        i_inc          ,
			  NOMAD::dd_type             display_degree ,
			  const NOMAD::Display     & out            ,
			  NOMAD::Point             & xf             ,
			  NOMAD::Point             & xi             ,
			  bool                     & stop           ,
			  NOMAD::stop_type         & stop_reason      );

    /// Project to mesh and create a trial point.
    /**
       \param ev_control     The NOMAD::Evaluator_Control object -- \b IN.
       \param x              The point coordinates               -- \b IN.
       \param model          The model                           -- \b IN.
       \param signature      Signature                           -- \b IN.
       \param mesh_indices   Mesh indic                          -- \b IN.
       \param delta			 Mesh size parameter                 -- \b IN.
       \param display_degree Display degree                      -- \b IN.
       \param out            The NOMAD::Display object           -- \b IN.
    */
    void create_trial_point
    ( NOMAD::Evaluator_Control & ev_control     ,
      NOMAD::Point               x              ,
      const NOMAD::Quad_Model  & model          ,
      NOMAD::Signature         & signature      ,
	  const NOMAD::Point        & mesh_indices  ,
      const NOMAD::Point       & delta          ,
      NOMAD::dd_type             display_degree ,
      const NOMAD::Display     & out              );
      
    /*----------------------------------------------------------------------*/

  public:

    /// Constructor.
    /**
       \param p Parameters -- \b IN.
    */
    Quad_Model_Search ( NOMAD::Parameters & p )
      : NOMAD::Search ( p , NOMAD::MODEL_SEARCH ) {}
    
    /// Destructor.
    virtual ~Quad_Model_Search ( void ) {}

    /// Reset.
    virtual void reset ( void ) {}

    /// The quadratic model search.
    /**
       Based on quadratic regression/MFN interpolation models.
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
    //// Display stats.
    /**
       \param out The NOMAD::Display object -- \b IN.
    */
    virtual void display ( const NOMAD::Display & out ) const
    {
      out << _all_searches_stats;
    }
  };
}

#endif
