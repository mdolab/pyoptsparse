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
  \file   Quad_Model.hpp
  \brief  Quadratic regression or MFN interpolation model (headers)
  \author Sebastien Le Digabel
  \date   2010-08-31
  \see    Quad_Model.cpp
*/
#ifndef __QUAD_MODEL__
#define __QUAD_MODEL__

#include "Cache.hpp"
#include "Model_Sorted_Point.hpp"
#include "Evaluator.hpp"

namespace NOMAD {

  /// Class for quadratic regression or MFN interpolation model.
  class Quad_Model : private NOMAD::Uncopyable {

    /*-------------------------------------------------------------------------*/
  private:

    const NOMAD::Display                     & _out;  ///< Display.

    std::vector<NOMAD::Eval_Point *>           _Y;    ///< Interpolation points.

    const std::vector<NOMAD::bb_output_type> & _bbot; ///< Blackbox output types.

    NOMAD::interpolation_type _interpolation_type;    ///< Interpolation type.

    int                       _n;                     ///< Dimension.
    int                       _nfree;                 ///< Number of free variables.
    int                       _n_alpha;               ///< Number of model coefficients.
    bool                    * _fixed_vars;            ///< Fixed variables.
    int                     * _index;                 ///< Var. indexes with fixed var.
    NOMAD::Point           ** _alpha;                 ///< Model coefficients.
    NOMAD::Point              _center;                ///< Model center.
    NOMAD::Point              _ref;                   ///< Reference for scaling.
    NOMAD::Point              _scaling;               ///< Scaling.
    const NOMAD::Cache      & _cache;                 ///< Cache.
    const NOMAD::Signature  & _signature;             ///< Signature.
    bool                      _error_flag;            ///< Error flag.
	std::list<NOMAD::Direction>  _dirP;               ///< Directions used for scaling (may be empty)
	NOMAD::Point _delta_m;                           ///< Mesh size used for scaling 
	NOMAD::Double _epsilon;                           ///< Offset for direction scaling 
	  
	  
    NOMAD::Double             _cond;                  ///< Condition number.

    /// Initialize alpha (model parameters).
    void init_alpha ( void );

    /// Check if an unscaled point is in \c B(center,radius) for a given radius.
    /**
       \param x      The unscaled point -- \b IN.
       \param radius The radius         -- \b IN.
       \return \c true is \c x is in \c B(center,radius).
    */
    bool is_within_radius ( const NOMAD::Point & x      ,
			    const NOMAD::Point & radius   ) const;

    /// Check the interpolation set \c Y.
    /**
       \return \c true if the interpolation set is valid.
    */
    bool check_Y ( void ) const;

    /// Check outputs before the integration into \c Y.
    /**
       \param bbo The outputs       -- \b IN.
       \param m   Number of outputs -- \b IN.
       return \c true if the \c m outputs are valid.
    */
    bool check_outputs ( const NOMAD::Point & bbo , int m ) const;

    /// Reduce the number of interpolation points.
    /**
       The points are sorted accorded to their distance to the model center.
       \param center      Center of the model                -- \b IN.
       \param max_Y_size  Max number of interpolation points -- \b IN.
    */
    void reduce_Y  ( const NOMAD::Point & center , int max_Y_size );

    /// Compute condition number.
    /**
       \param W   Matrix W given as a vector -- \b IN.
       \param n   Size of \c W               -- \b IN
       \param eps Epsilon                    -- \b IN.
    */
    void compute_cond ( const double * W , int n , double eps );

    /// Compute the cumulated error of a model for one output.
    /**
       The errors are computed on the interpolation set \c Y.
       \param bbo_index   Blackbox output index  -- \b IN.
       \param error       Cumulated error        -- \b OUT. 
       \param min_rel_err Min relative error     -- \b OUT.
       \param max_rel_err Max relative error     -- \b OUT.
       \param avg_rel_err Average relative error -- \b OUT.
    */
    void compute_model_error ( int             bbo_index   ,
			       NOMAD::Double & error       ,
			       NOMAD::Double & min_rel_err ,
			       NOMAD::Double & max_rel_err ,
			       NOMAD::Double & avg_rel_err   ) const;

    /// Compute the maximal relative error of a model.
    /**
       The error is computed for the interpolation set \c Y.
       \return The maximal relative error.
    */
    NOMAD::Double compute_max_rel_err ( void ) const;

    /// Compute the element \c (i,j) of the interpolation matrix \c M(phi,Y).
    /**
       \param i Row index    -- \b IN.
       \param j Column index -- \b IN.
    */
    double compute_M ( int i , int j ) const;

    /// Construct Minimum Frobenius Norm (MFN) model.
    /**
       - This occurs when \c p+1 \c < \c (n+1)(n+2)/2.
       \param eps        Epsilon                               -- \b IN.
       \param max_mpn    Maximum \c m+n value for SVD matrices -- \b IN.
       \param max_Y_size Maximum number of elements in \c Y    -- \b IN.
       \return \c true if the construction succeeded
    */
    bool construct_MFN_model ( double eps , int max_mpn , int max_Y_size );

    /// Construct regression model.
    /**
       - This occurs when \c p+1 \c >= \c (n+1)(n+2)/2.
       \param eps        Epsilon                               -- \b IN.
       \param max_mpn    Maximum \c m+n value for SVD matrices -- \b IN.
       \param max_Y_size Maximum number of elements in \c Y    -- \b IN.
       \return \c true if the construction succeeded
    */
    bool construct_regression_model ( double eps        ,
				      int    max_mpn    ,
				      int    max_Y_size   );

    /// Construct well-poised (WP) model.
    /**
       \param max_Y_size Maximum number of elements in \c Y -- \b IN.
       \return \c true if the construction succeeded
    */
    bool construct_WP_model ( int max_Y_size );

    /// Find interpolation point with max Lagrange polynomial value.
    /**
       \param  li      Lagrange polynomial             -- \b IN.
       \param  Y       Interpolation points candidates -- \b IN.
       \param  i1      Initial index in \c Y           -- \b IN.
       \param  i2      Final index in \c Y             -- \b IN.
       \param  max_lix Absolute value of the max value -- \b OUT.
       \return Index of interpolation point.
    */
    int find_max_lix ( const NOMAD::Point                     & li      ,
		       const std::vector<NOMAD::Eval_Point *> & Y       ,
		       int                                      i1      ,
		       int                                      i2      ,
		       NOMAD::Double                          & max_lix   ) const;

    /// Resolution of system \c F.[mu alpha_L]'=[f(Y) 0]' for MFN interpolation.
    /**
       \param U         Matrix \c F=U from the SVD decomposition \c U.W.V' -- \b IN.
       \param W         Matrix \c W from the SVD decomposition \c U.W.V'   -- \b IN.
       \param V         Matrix \c V from the SVD decomposition \c U.W.V'   -- \b IN.
       \param bbo_index Blackbox output index                              -- \b IN.
       \param alpha     Model parameters                                   -- \b IN.
       \param eps       Epsilon                                            -- \b IN.
    */
    void solve_MFN_system ( double      ** U         ,
			    double       * W         , 
			    double      ** V         ,
			    int            bbo_index ,
			    NOMAD::Point & alpha     ,
			    double         eps	      ) const;

    /// Resolution of system \c F.alpha=M'.f(Y) for the regression.
    /**
       \param M         Matrix \c M                                        -- \b IN.
       \param U         Matrix \c F=U from the SVD decomposition \c U.W.V' -- \b IN.
       \param W         Matrix \c W from the SVD decomposition \c U.W.V'   -- \b IN.
       \param V         Matrix \c V from the SVD decomposition \c U.W.V'   -- \b IN.
       \param bbo_index Blackbox output index                              -- \b IN.
       \param alpha     Model parameters                                   -- \b IN.
       \param eps       Epsilon                                            -- \b IN.
    */
    void solve_regression_system ( double      ** M         ,
				   double      ** U         ,
				   double       * W         , 
				   double      ** V         ,
				   int            bbo_index ,
				   NOMAD::Point & alpha     ,
				   double         eps	      ) const;

    /// Display Lagrange polynomials.
    /**
       \param l Lagrange polynomials -- \b IN.
       \param Y Interpolation set    -- \b IN.
    */
    void display_lagrange_polynomials 
    ( const std::vector<NOMAD::Point      *> & l ,
      const std::vector<NOMAD::Eval_Point *> & Y   ) const;

#ifdef MODEL_STATS
    mutable NOMAD::Double _Yw; ///< Width of the interpolation set \c Y.

  public:

    /// Access to the width of the interpolation set \c X (or \c Y).
    /**
       \return The width of the interpolation set.
    */
    const NOMAD::Double & get_Yw ( void ) const { return _Yw; }
#endif

    /*-------------------------------------------------------------------------*/
  public:

    /// Constructor.
    /**
       \param out           The NOMAD::Display object   -- \b IN.
       \param bbot          Output types                -- \b IN.
       \param cache         Cache                       -- \b IN.
       \param signature     Signature                   -- \b IN.
    */
    Quad_Model ( const NOMAD::Display                     & out       ,
		 const std::vector<NOMAD::bb_output_type> & bbot      ,
		 const NOMAD::Cache                       & cache     ,
		 const NOMAD::Signature                   & signature   );

    /// Destructor.
    virtual ~Quad_Model ( void );

    /// Evaluate a model at a given point.
    /**
       \param x     The point        -- \b IN.
       \param alpha Model parameters -- \b IN.
       \return Model value.
    */
    NOMAD::Double eval ( const NOMAD::Point & x     ,
			 const NOMAD::Point & alpha   ) const;

	  /// Compute model \c h and \c f values at a point.
    /**
       \param x      The point                 -- \b IN.
       \param h_min  Value of \c h_min         -- \b IN..
       \param h_norm Norm used to compute \c h -- \b IN..
       \param h      Value of \c h             -- \b OUT.
       \param f      Value of \c f             -- \b OUT.
    */
    void eval_hf ( const NOMAD::Point  & x      ,
		   const NOMAD::Double & h_min  ,
		   NOMAD::hnorm_type     h_norm ,
		   NOMAD::Double       & h      ,
		   NOMAD::Double       & f        ) const;

	  
    /// Access to the interpolation type.
    /**
       \return The interpolation type.
    */
    const NOMAD::interpolation_type & get_interpolation_type ( void ) const
    {
      return _interpolation_type;
    }

    /// Access to the center of the model.
    /**
       \return The center.
    */
    const NOMAD::Point & get_center ( void ) const { return _center; }

    /// Access to the dimension.
    /**
       \return The dimension \c n.
    */
    int get_n ( void ) const { return _n; }

    /// Access to the number of free variables.
    /**
       \return The number of free variables \c n.
    */
    int get_nfree ( void ) const { return _nfree; }

    /// Access to the model parameters.
    /**
       \return The model parameters \c alpha.
    */
    NOMAD::Point ** get_alpha ( void ) const { return _alpha; }

    /// Check if the model is ready for evaluations.
    /**
       \return A boolean equal to \c true if the model is ready.
    */
    bool check ( void ) const;

    /// Access to the fixed variables.
    /**
       \param i Variable index -- \b IN.
       \return \c true if variable \c i is fixed.
    */
    bool variable_is_fixed ( int i ) const { return _fixed_vars[i]; }

    /// Access to the number of interpolation points.
    /**
       \return The number of interpolation points \c nY=p+1.
    */
    int get_nY ( void ) const { return static_cast<int> ( _Y.size() ); }

    /// Access to the condition number.
    /**
       \return The condition number.
    */
    const NOMAD::Double & get_cond ( void ) const { return _cond; }

    /// Access to the error flag.
    /**
       \return The error flag.
    */
    bool get_error_flag ( void ) const { return _error_flag; }

    /// Construct the interpolation set \c Y.
    /**
       \param center               Model center                       -- \b IN.
       \param interpolation_radius Interpolation radius               -- \b IN.
       \param max_Y_size           Maximum number of elements in \c Y -- \b IN.
    */
    void construct_Y ( const NOMAD::Point & center               ,
		       const NOMAD::Point & interpolation_radius ,
		       int                  max_Y_size             );

    /// Construct \c m models (one by output).
    /**
       \param use_WP     Use or not well-poisedness            -- \b IN.
       \param eps        Epsilon                               -- \b IN.
       \param max_mpn    Maximum \c m+n value for SVD matrices -- \b IN.
       \param max_Y_size Maximum number of elements in \c Y    -- \b IN.
    */
    void construct ( bool   use_WP     ,
		     double eps        ,
		     int    max_mpn    ,
		     int    max_Y_size   );

    /// Define scaling to put all coordinates centered in \c [-r;r].
    /**
       - Looks also for fixed variables.
       \param r The \c r parameter corresponds to \c MODEL_RADIUS_FACTOR -- \b IN.
    */
    void define_scaling ( const NOMAD::Double & r );

	  /// Define scaling based on directions. See paper: Reducing the number of function evaluations in Mesh Adaptive Direct Search algorithms, Audet, Ianni, LeDigabel, Tribes, 2014     
	  /**
       - Looks also for fixed variables.
       \param dirP    The \c dirP parameter corresponds to set of directions formin a hyper-cube centered on poll center -- \b IN.
	   \param delta_m The \c delta_m parameter is the dimension of the mesh -- \b IN.
	   \param epsilon The \c epsilon parameter is the hyper-cube offset from the poll center -- \b IN.
	   */	  
	void define_scaling_by_directions ( const std::list<NOMAD::Direction> & dirP, const NOMAD::Point & delta_m, const NOMAD::Double &epsilon  );

	  
    /// Scale a point.
    /**
       \param x The point to scale -- \b IN/OUT.
       \return \c true if the scaling worked.
    */
    bool scale ( NOMAD::Point & x ) const;

    /// Unscale a point.
    /**
       \param x The point to unscale -- \b IN/OUT.
       \return \c true if the unscaling worked.
    */
    bool unscale ( NOMAD::Point & x ) const;
	  
	  /// Unscale the gradient at a point.
	  /**
       \param x The grad to unscale -- \b IN/OUT.
       \return \c true if the unscaling worked.
	   */
	  bool unscale_grad ( NOMAD::Point & x ) const;	  
	  

    /// Check if a caled point is inside the trust radius.
    /**
       \param x The scaled point -- \b IN.
       \return  \c true is \c x is in the trust radius.
    */
    bool is_within_trust_radius ( const NOMAD::Point & x ) const;

    /// Display the model coefficients.
    /**
       \param out   The NOMAD::Display object -- \b IN.
    */
    void display_model_coeffs ( const NOMAD::Display & out ) const;

    /// Display the interpolation set \c Y.
    /**
       \param out   The NOMAD::Display object  -- \b IN.
       \param title Title of the display block -- \b IN
                    --\b optional (default="interpolation set Y").
    */
    void display_Y ( const NOMAD::Display & out ,
		     const std::string    & title = "interpolation set Y" ) const;

    /// Display cumulated error on the interpolation points.
    /**
       \param out   The NOMAD::Display object -- \b IN.
    */
    void display_Y_error ( const NOMAD::Display & out ) const;
  };
}

#endif
