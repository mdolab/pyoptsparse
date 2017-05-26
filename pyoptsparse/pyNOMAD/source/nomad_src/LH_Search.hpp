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
  \file   LH_Search.hpp
  \brief  Latin-Hypercube search (headers)
  \author Sebastien Le Digabel
  \date   2010-04-09
  \see    LH_Search.cpp
*/
#ifndef __LH_SEARCH__
#define __LH_SEARCH__

#include "Search.hpp"
#include "Mads.hpp"
#include "RNG.hpp"  

namespace NOMAD {

  /// Class for the Latin-Hypercube search.
  class LH_Search : public NOMAD::Search , private NOMAD::Uncopyable {

  private:

    bool _initial_search;  /// Initial search flag (for display only).

    /// Decide \c p values for one variable.
    /**
       If no bounds, values are scaled with the largest Delta^m_k value obtained so far.
       \param p           Number of values to decide                 -- \b IN.
       \param delta_m     Delta^m_k (for the projection to the mesh) -- \b IN.
       \param delta_m_max Largest Delta^m_k value                    -- \b IN.
       \param bbit        Black-box input type                       -- \b IN.
       \param lb          Lower bound                                -- \b IN.
       \param ub          Upper bound                                -- \b IN.
       \param x           The \p values                              -- \b OUT.
    */
    void values_for_var_i ( int                          p           ,
			    const NOMAD::Double        & delta_m     ,
			    const NOMAD::Double        & delta_m_max ,
			    const NOMAD::bb_input_type & bbit        ,
			    const NOMAD::Double        & lb          ,
			    const NOMAD::Double        & ub          ,
			    NOMAD::Point               & x             ) const;
  public:
    
    /// Constructor.
    /**
       \param p              Parameters          -- \b IN.
       \param initial_search Initial search flag -- \b IN.
       \param phase_one      Phase one flag      -- \b IN.
    */
    LH_Search ( NOMAD::Parameters        & p              ,
		bool                       initial_search ,
		bool                       phase_one        )
      : NOMAD::Search   ( p , phase_one ? NOMAD::LH_SEARCH_P1 : NOMAD::LH_SEARCH ) ,
	_initial_search ( initial_search                                         )   {}

    /// Destructor.
    virtual ~LH_Search ( void ) {}

    /// The Latin-Hypercube search.
    /**
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

    /// Simpler method used to generate a list of LH points.
    /**
       \param n   Dimension           -- \b IN.
       \param m   Number of outputs   -- \b IN.
       \param p   Number of LH points -- \b IN.
       \param lb  Lower bounds        -- \b IN.
       \param ub  Upper bounds        -- \b IN.
       \param pts LH points           -- \b OUT.
       \return A boolean equal to \c true if no error occured.
    */
    static bool LH_points ( int                                n   ,
			    int                                m   ,
			    int                                p   ,
			    const NOMAD::Point               & lb  ,
			    const NOMAD::Point               & ub  ,
			    std::vector<NOMAD::Eval_Point *> & pts   );
  };
}

#endif
