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
  \file   TGP_Model.hpp
  \brief  TGP (Bayesian treed Gaussian process) model for all outputs (headers)
  \author Sebastien Le Digabel
  \date   2011-02-07
  \see    TGP_Model.cpp
*/
#ifdef USE_TGP

#ifndef __TGP_MODEL__
#define __TGP_MODEL__

#include "TGP_Output_Model.hpp"
#include "Model_Sorted_Point.hpp"
#include "Cache.hpp"

namespace NOMAD {

  /// TGP models for all outputs.
  class TGP_Model : private NOMAD::Uncopyable {

  private:

    const NOMAD::Display & _out;  ///< Display.

    std::vector<NOMAD::bb_output_type> _bbot; ///< Blackbox output types.

    std::string  _error_str; ///< Error string.

    int          _p;    ///< Number of interpolation points.
    
    int          _n0;   ///< Original dimension (all variables).
    int          _n;    ///< Number of free variables (\c n \c <= \c n0).

    int          _n_XX; ///< Number of points (rows) in \c XX.

    NOMAD::Point _lb;   ///< Lower bounds for \c X    (size=\c n0).
    NOMAD::Point _ub;   ///< Upper bounds for \c X    (size=\c n0).
    NOMAD::Point _fv;   ///< Fixed variables for \c X (size=\c n0).

    int   * _av_index; ///< All variables index.
    int   * _fv_index; ///< Free variables index.
    //
    // avi[j] in [0;n -1], j in [0;n0-1]; avi[j]=-1: var j fixed
    // fvi[j] in [0;n0-1], j in [0;n-1]
    //
    // example with n0=3 and n=2 (second variable is fixed):
    //   _av_index: [ 0 -1 1 ]
    //   _fv_index: [ 0 2 ]

    double ** _X;      ///< Interpolation matrix \c X (size = \c p \c x \c n).
    
    double ** _XX;     ///< Points at which evaluate the model (prediction points).
    
    double ** _Ds2x;   ///< Expected reduction in predictive variance (\c Ds2x).

    int     * _improv; ///< Expected improvement (ranks).

    bool _model_computed;            ///< Flag to check if the model has been computed.
    bool _nep_flag;                  ///< 'Not enought points' (nep) flag.

    NOMAD::TGP_mode_type       _tgp_mode;   ///< TGP mode (fast, precise, or user).

    NOMAD::TGP_Output_Model ** _tgp_models; ///< TGP models (one for each output).
    //
    // tgp_models[i] may be NULL if bbot[i]
    // is not an objective or a constraint

    double ** _tgp_rect;       ///< TGP rectangle (bounds).

    int       _usr_BTE[3];     ///< User \c BTE parameters (TGP default: 2000,7000,2).

    bool      _tgp_linburn;    ///< TGP \c linburn parameter.

    int       _usr_pmax;       ///< User max number of interpolation points in \c X.

    bool      _tgp_verb;       ///< Display (\c verb) parameter.

    /// Clear memory.
    void clear ( void );

    /// Filter and sort the intepolation set \c X.
    /**
       \param  X          Original interpolation set (typically the cache) -- \b IN.
       \param  center     Center of the model (may be NULL)                -- \b IN.
       \param  filtered_X Filtered and sorted interpolation set            -- \b OUT.
       \return The number of filtered interpolation points.
    */
    int filter_and_sort_X
    ( const std::vector<const NOMAD::Eval_Point *> & X          ,
      const NOMAD::Point                           * center     ,
      std::list<const NOMAD::Eval_Point *>         & filtered_X   ) const;

    /**
       Compute the ranges, check the fixed variables, set the
       different indexes, and return the number of free variables.
       \param  X          Filtered and sorted interpolation set       -- \b IN.
       \param  remove_fv  Set to \c true to eliminate fixed variables -- \b IN.
       \return The number of free variables (member \c _n).
    */
    int check_fixed_variables ( const std::list<const NOMAD::Eval_Point *> & X         ,
				bool                                         remove_fv   );

    /// Tests to check if an interpolation point is valid for interpolation.
    /**
       \param  x The interpolation point -- \b IN.
       \return A boolean equal to \c true if \c x is valid.
    */
    bool test_interpolation_point ( const NOMAD::Eval_Point * x ) const;

    /// Get the limits on the interpolation set size \c p.
    /**
       \param  n     Number of free variables                              -- \b IN.
       \param  pmin  Inf limit on \c p the number of interpolation points  -- \b OUT.
       \param  pmax  Sup limit on \c p                                     -- \b OUT.
       \return A boolean equal to \c true if the limits are valid.
    */
    bool get_p_limits ( int n , int & pmin , int & pmax );

    /// Get the TGP \c BTE parameters.
    /**
       \param  BTE The BTE parameters -- \b OUT.
       \return A boolean equal to \c true if the BTE parameters are valid.
    */
    bool get_BTE ( int BTE[3] );

    /// Compute the TGP \c dparam array.
    /**
       \param  n Dimension -- \b IN.
       \return A pointer to the \c dparam array.
    */
    static double * get_TGP_dparams ( int n );

    /// Check if the interpolation matrix is of full rank.
    /**
       \param  X The \c X matrix    -- \b IN.
       \param  p Number of rows     -- \b IN.
       \param  n Number of columns  -- \b IN.
       \return A boolean equal to \c true if \c X is of full rank.
    */
    static bool check_full_rank ( double ** X , int p , int n );

  public:

    /// Constructor 1/2.
    /**
       \param n0         Dimension including fixed variables -- \b IN.
       \param bbot       Blackbox output types               -- \b IN.
       \param out        Display object                      -- \b IN.
       \param mode       TGP mode (fast or precise only)     -- \b IN.
    */
    explicit TGP_Model ( int                                        n0   ,
			 const std::vector<NOMAD::bb_output_type> & bbot ,
			 const NOMAD::Display                     & out  ,
			 NOMAD::TGP_mode_type                       mode   );

    /// Constructor 2/2 (with no blackbox output types; a value of \c m=1 is taken).
    /**
       This is the user TGP mode.
       \param n0         Dimension including fixed variables        -- \b IN.
       \param out        Display object                             -- \b IN.
       \param BTE        TGP \c B, \c T, and \c E parameters        -- \b IN.
       \param linburn    TGP \c linburn parameter                   -- \b IN.
       \param pmax       Max number of interpolation points in \c X -- \b IN.
       \param verb       TGP \c verb parameter                      -- \b IN.
    */
    explicit TGP_Model ( int                    n0      ,
			 const NOMAD::Display & out     ,
			 int                    BTE[3]  ,
			 bool                   linburn ,
			 int                    pmax    ,
			 bool                   verb      );

    /// Destructor.
    virtual ~TGP_Model ( void );

    /// Set the interpolation set \c X from a cache.
    /**
       \param  cache      Cache of true evaluations                           -- \b IN.
       \param  center     Center of the model (the incumbent; may be \c NULL) -- \b IN.
       \param  seed       Random seed                                         -- \b IN.
       \param  remove_fv  Set to \c true to eliminate fixed variables         -- \b IN.
       \return A boolean equal to \c true if no error occured.
    */
    bool set_X ( const NOMAD::Cache & cache     ,
		 const NOMAD::Point * center    ,
		 int                  seed      ,
		 bool                 remove_fv   );

    /// Set the interpolation set \c X from a set of points.
    /**
       \param  X          The set of points                           -- \b IN.
       \param  center     Center of the model ( may be \c NULL)       -- \b IN.
       \param  seed       Random seed                                 -- \b IN.
       \param  remove_fv  Set to \c true to eliminate fixed variables -- \b IN.
       \return A boolean equal to \c true if no error occured.
    */
    bool set_X ( const std::vector<const NOMAD::Eval_Point *> & X         ,
		 const NOMAD::Point                           * center    ,
		 int                                            seed      ,
		 bool                                           remove_fv   );

    /// Compute the models (one for each output).
    /**
       \param  XX_pts            Prediction points                         -- \b IN/OUT.
       \param  compute_Ds2x      Flag to activate the Ds2x computation     -- \b IN.
       \param  compute_improv    Flag to activate the improv computation   -- \b IN.
       \param  pred_outside_bnds If \c false, no prediction outside bounds -- \b IN.
       \return A boolean equal to \c true if the model computation worked.
    */
    bool compute ( std::vector<NOMAD::Eval_Point *> & XX_pts            ,
		   bool                               compute_Ds2x      ,
		   bool                               compute_improv    ,
		   bool                               pred_outside_bnds   );

    /// Prediction at one point.
    /**
       \param  x                 The point (size \c _n or \c _n0)              -- \b IN.
       \param  pred_outside_bnds Set to \c false to not predict outside bounds -- \b IN.
       \return A boolean equal to \c true if the prediction worked.
    */
    bool predict ( NOMAD::Eval_Point & x , bool pred_outside_bnds );

    /// Compute model \c h and \c f values given one blackbox output.
    /**
       \param  bbo    Blackbox output           -- \b IN.
       \param  h_min  Value of \c h_min         -- \b IN..
       \param  h_norm Norm used to compute \c h -- \b IN..
       \param  h      Value of \c h             -- \b OUT.
       \param  f      Value of \c f             -- \b OUT.
    */
    void eval_hf ( const NOMAD::Point  & bbo    ,
		   const NOMAD::Double & h_min  ,
		   NOMAD::hnorm_type     h_norm ,
		   NOMAD::Double       & h      ,
		   NOMAD::Double       & f        ) const;

    /// Get the \c XX points with the largest expected improvement of the objective.
    /**
       \param pts_index The \c XX points indexes -- \b OUT.
    */
    void get_improv_points ( std::list<int> & pts_indexes ) const;

    /// Get the \c XX points with the largest expected reduction in predictive variance.
    /**
       \param pts_index The \c XX points indexes -- \b OUT.
    */
    void get_Ds2x_points ( std::set<int> & pts_indexes ) const;

    /// Get error string.
    /**
       \return The error string.
    */
    const std::string & get_error_str ( void ) const { return _error_str; }

    /// Get \c nep flag (nep=not enough points).
    /**
       \return The nep flag.
    */
    bool get_nep_flag ( void ) const { return _nep_flag; }

    /// Access to the number of interpolation points \c p.
    /**
       \return The number of interpolation points.
    */
    int get_p ( void ) const { return _p; }

    /// Access to the number of free variables \c n.
    /**
       \return The number of free variables.
    */
    int get_n ( void ) const { return _n; }

    /// Access to the total number of variables \c n0.
    /**
       \return The total number of variables.
    */
    int get_n0 ( void ) const { return _n0; }

    /// Access to the lower bounds.
    /**
       \return The lower bounds.
    */
    const NOMAD::Point & get_lb ( void ) const { return _lb; }

    /// Access to the upper bounds.
    /**
       \return The upper bounds.
    */
    const NOMAD::Point & get_ub ( void ) const { return _ub; }

#ifdef MODEL_STATS
    /// Access to the width of the interpolation set \c X (or \c Y).
    /**
       \return The width of the interpolation set.
    */
    NOMAD::Double get_Yw ( void ) const;
#endif

    /// Access to improv.
    /**
       \param  i Index in the \c XX matrix -- \b IN.
       \return   The improv value of the ith point in \c XX.
    */
    int get_improv ( int i ) const
    {
      return ( _improv == NULL || i < 0 || i >= _n_XX ) ? -1 : _improv[i];
    }

    /// Display the expected reduction in predictive variance (\c Ds2x).
    /**
       \param out Display                                                 -- \b IN.
       \param display_limit Max number of pts displayed (-1 for no limit) -- \b IN.
    */
    void display_Ds2x ( const NOMAD::Display & out , int display_limit = -1 ) const;

    /// Display the expected improvement ranks (\c improv).
    /**
       \param out Display                                                 -- \b IN.
       \param display_limit Max number of pts displayed (-1 for no limit) -- \b IN.
    */
    void display_improv ( const NOMAD::Display & out , int display_limit = -1 ) const;

    /// Display error stats for the interpolation set \c X.
    /**
       \param out Display -- \b IN.
    */
    void display_X_errors ( const NOMAD::Display & out );

    /// Display interpolation set \x X.
    /**
       \param out Display                                                 -- \b IN.
       \param display_limit Max number of pts displayed (-1 for no limit) -- \b IN.
    */
    void display_X ( const NOMAD::Display & out , int display_limit = -1 ) const;

    /// Default display.
    void display ( void ) { display ( _out ); }

    /// Display.
    /**
       \param out Display -- \b IN.
    */
    void display ( const NOMAD::Display & out ) const;   
  };

  /// Display a NOMAD::TGP_Model object.
  /**
     \param out The NOMAD::Display object                   -- \b IN.
     \param s   The NOMAD::TGP_Model object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
  */
  inline const NOMAD::Display & operator << ( const NOMAD::Display   & out ,
					      const NOMAD::TGP_Model & s     ) {
    s.display ( out );
    return out;
  }
}

#endif

#endif
