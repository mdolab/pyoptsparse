/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search - version 3.7.0.beta        */
/*                                                                                     */
/*  Copyright (C) 2001-2014  Mark Abramson        - the Boeing Company, Seattle        */
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
  \file   TGP_Model_Search.hpp
  \brief  TGP Model search (headers)
  \author Sebastien Le Digabel
  \date   2011-02-17
  \see    TGP_Model_Search.cpp
*/
#ifdef USE_TGP

#ifndef __TGP_MODEL_SEARCH__
#define __TGP_MODEL_SEARCH__

#include "LH_Search.hpp"
#include "TGP_Model_Evaluator.hpp"

namespace NOMAD {

  /// Model search.
  class TGP_Model_Search : public NOMAD::Search , private NOMAD::Uncopyable {

  private:

    NOMAD::TGP_Model * _model;

    NOMAD::Model_Stats _one_search_stats;   ///< Stats for one search.
    NOMAD::Model_Stats _all_searches_stats; ///< Stats for all searches.

    /// Delete a list of points.
    /**
       \param pts The points -- \b IN/OUT.
    */
    static void clear_pts ( std::vector<NOMAD::Point *> & pts );

    /// Delete a list of evaluation points.
    /**
       \param pts The points -- \b IN/OUT.
    */
    static void clear_pts ( std::vector<NOMAD::Eval_Point *> & pts );

    /// Model construction.
    /**
       \param  cache          Cache of true evaluations                 -- \b IN.
       \param  incumbent      The incumbent                             -- \b IN.
       \param  delta_m        Mesh size parameter                       -- \b IN.
       \param  out            The NOMAD::Display object                 -- \b IN.
       \param  display_degree Display degree                            -- \b IN.
       \param  display_lim    Max number of pts when sets are displayed -- \b IN.
       \param  stats          Model search stats                        -- \b IN/OUT.
       \param  compute_Ds2x   Flag to enable/disable Ds2x computation   -- \b OUT.
       \param  XX             The set of prediction points              -- \b OUT.
       \param  stop           Stop flag                                 -- \b OUT.
       \param  stop_reason    Stop reason                               -- \b OUT.
       \param  error_std      Error string                              -- \b OUT.
       \return A boolean equal to \c true if the model has been constructed.
    */
    bool model_construction ( const NOMAD::Cache               & cache          ,
			      const NOMAD::Point               & incumbent      ,
			      const NOMAD::Point               & delta_m        ,
			      const NOMAD::Display             & out            ,
			      NOMAD::dd_type                     display_degree ,
			      int                                display_lim    ,
			      NOMAD::Stats                     & stats          ,
			      bool                             & compute_Ds2x   ,
			      std::vector<NOMAD::Eval_Point *> & XX             ,
			      bool                             & stop           ,
			      NOMAD::stop_type                 & stop_reason    ,
			      std::string                      & error_str        );

    /// Create a list of prediction points.
    /**
       \param  cache     Cache of true evaluations    -- \b IN.
       \param  n         Number of variables          -- \b IN.
       \param  m         Number of outputs            -- \b IN.
       \param  incumbent The incumbent                -- \b IN.
       \param  delta_m   Mesh size parameter          -- \b IN.
       \param  XX        The set of prediction points -- \b OUT.
    */
    void set_XX ( const NOMAD::Cache               & cache     ,
		  int                                n         ,
		  int                                m         ,
		  const NOMAD::Point               & incumbent ,
		  const NOMAD::Point               & delta_m   ,
		  std::vector<NOMAD::Eval_Point *> & XX          ) const;

    /// Create the complete list of trial points (oracle + Ds2x + improv).
    /**
       \param oracle_pts     Oracle points              -- \b IN.
       \param Ds2x_pts       Ds2x points                -- \b IN.
       \param improv_pts     Improv points              -- \b IN.
       \param incumbent      The incumbent              -- \b IN.
       \param max_pts        Max number of trial points -- \b IN.
       \param out            The NOMAD::Display object  -- \b IN.
       \param display_degree Display degree             -- \b IN.
       \param trial_pts      The list of trial points   -- \b OUT.
     */
    void create_trial_pts
    ( const std::vector<NOMAD::Point *> & oracle_pts     ,
      const std::vector<NOMAD::Point *> & Ds2x_pts       ,
      const std::vector<NOMAD::Point *> & improv_pts     ,
      const NOMAD::Point                & incumbent      ,
      int                                 max_pts        ,
      const NOMAD::Display              & out            ,
      NOMAD::dd_type                      display_degree ,
      std::vector<NOMAD::Point *>       & trial_pts        ) const;

    /// Create oracle points by optimizing the model.
    /**
       \param  cache          Cache of true evaluations                 -- \b IN.
       \param  incumbent      The incumbent                             -- \b IN.
       \param  delta_m        Mesh size parameter                       -- \b IN.
       \param  out            The NOMAD::Display object                 -- \b IN.
       \param  display_degree Display degree                            -- \b IN.
       \param  display_lim    Max number of pts when sets are displayed -- \b IN.
       \param  XX             The set of prediction points              -- \b IN.
       \param  oracle_pts     Oracle candidates points                  -- \b OUT.
       \param  stop           Stop flag                                 -- \b OUT.
       \param  stop_reason    Stop reason                               -- \b OUT.
       \return A boolean equal to \c true oracle points are proposed.
    */
    bool create_oracle_pts
    ( const NOMAD::Cache                     & cache          ,
      const NOMAD::Point                     & incumbent      ,
      const NOMAD::Point                     & delta_m        ,
      const NOMAD::Display                   & out            ,
      NOMAD::dd_type                           display_degree ,
      int                                      display_lim    ,
      const std::vector<NOMAD::Eval_Point *> & XX             ,
      std::vector<NOMAD::Point *>            & oracle_pts     ,
      bool                                   & stop           ,
      NOMAD::stop_type                       & stop_reason      );

    /// Model optimization.
    /**
       \param x0s            The three starting points -- \b IN.
       \param out            The NOMAD::Display object -- \b IN.
       \param display_degree Display degree            -- \b IN.
       \param xf             Feasible solution \c xf   -- \b OUT.
       \param xi             Infeasible solution \c xi -- \b OUT.
       \param stop           Stop flag                 -- \b OUT.
       \param stop_reason    Stop reason               -- \b OUT.
    */
    bool optimize_model ( const NOMAD::Eval_Point * x0s[3]         ,
			  const NOMAD::Display    & out            ,
			  NOMAD::dd_type            display_degree ,
			  NOMAD::Point           *& xf             ,
			  NOMAD::Point           *& xi             ,
			  bool                    & stop           ,
			  NOMAD::stop_type        & stop_reason      );

    /// Project and accept or reject an oracle trial point.
    /**
       \param  cache          Cache of true evaluations -- \b IN.
       \param  incumbent      The incumbent             -- \b IN.
       \param  delta_m        Mesh size parameter       -- \b IN.
       \param  out            The NOMAD::Display object -- \b IN.
       \param  display_degree Display degree            -- \b IN.
       \param  x              The oracle point         -- \b IN/OUT.
       \return A boolean equal to \c true if the point is accepted.
    */
    bool check_oracle_point
    ( const NOMAD::Cache   & cache          ,
      const NOMAD::Point   & incumbent      ,
      const NOMAD::Point   & delta_m        ,
      const NOMAD::Display & out            ,
      NOMAD::dd_type         display_degree ,
      NOMAD::Point         & x                );

    /// Insert a trial point in the evaluator control object.
    /**
       \param x              The point coordinates               -- \b IN.
       \param signature      Signature                           -- \b IN.
       \param incumbent      The incumbent                       -- \b IN.
       \param display_degree Display degree                      -- \b IN.
       \param ev_control     The NOMAD::Evaluator_Control object -- \b IN/OUT.
    */
    void register_point ( NOMAD::Point               x              ,
			  NOMAD::Signature         & signature      ,
			  const NOMAD::Point       & incumbent      ,
			  // C.Tribes august 26, 2014 --- not needed
              // int                        mesh_index     ,
			  NOMAD::dd_type             display_degree ,
			  NOMAD::Evaluator_Control & ev_control       ) const;

    /// Create the list of improv points.
    /**
       These points (from the set \c XX) maximize
         the expected improvement of the objective.
       Priority is given to predicted feasible points.
       \param  XX            The set of prediction points              -- \b IN.
       \param incumbent      The incumbent                             -- \b IN.
       \param max_pts        Max number of points                      -- \b IN.
       \param out            The NOMAD::Display object                 -- \b IN.
       \param display_degree Display degree                            -- \b IN.
       \param display_lim    Max number of pts when sets are displayed -- \b IN.
       \param Ds2x_pts       The list of improv points                 -- \b OUT.
    */
    void create_improv_pts
    ( const std::vector<NOMAD::Eval_Point *> & XX             ,
      const NOMAD::Point                     & incumbent      ,
      int                                      max_pts        ,
      const NOMAD::Display                   & out            ,
      NOMAD::dd_type                           display_degree ,
      int                                      display_lim    ,
      std::vector<NOMAD::Point *>            & improv_pts       ) const;

    /// Create the list of Ds2x points.
    /**
       These points (from the set \c XX) maximize the expected reduction in
         predictive variance for each output.
       \param  XX            The set of prediction points              -- \b IN.
       \param out            The NOMAD::Display object                 -- \b IN.
       \param display_degree Display degree                            -- \b IN.
       \param display_lim    Max number of pts when sets are displayed -- \b IN.
       \param Ds2x_pts       The list of Ds2x points                   -- \b OUT.
    */
    void create_Ds2x_pts
    ( const std::vector<NOMAD::Eval_Point *> & XX             ,
      const NOMAD::Display                   & out            ,
      NOMAD::dd_type                           display_degree ,
      int                                      display_lim    ,
      std::vector<NOMAD::Point *>            & Ds2x_pts         ) const;

    /// Prediction at one point.
    /**
       \param  x The point     -- \b IN.
       \param  h Value of \c h -- \b OUT.
       \param  f Value of \c f -- \b OUT.
       \return A boolean equal to \c true if the prediction was possible.
    */
    bool predict ( const NOMAD::Point & x ,
		   NOMAD::Double      & h ,
		   NOMAD::Double      & f   ) const;

    /// Display the prediction error for the evaluated points.
    /**
       \param evaluated_pts List of evaluated points  -- \b IN.
       \param out           The NOMAD::Display object -- \b IN.
    */
    void display_eval_pred_errors
    ( const std::list<const NOMAD::Eval_Point *> & evaluated_pts ,
      const NOMAD::Display                       & out             );

    /*----------------------------------------------------------------------*/

  public:

    /// Constructor.
    /**
       \param p Parameters -- \b IN.
    */
    TGP_Model_Search ( NOMAD::Parameters & p )
      : NOMAD::Search ( p , NOMAD::MODEL_SEARCH ) , _model ( NULL ) {}
    
    /// Destructor.
    virtual ~TGP_Model_Search ( void ) { reset(); }

    /// Reset.
    virtual void reset ( void );

    /// The TGP model search.
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

    /// Access to the model.
    /**
       \return The model.
    */
    NOMAD::TGP_Model * get_model ( void ) const { return _model; }

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
#endif
