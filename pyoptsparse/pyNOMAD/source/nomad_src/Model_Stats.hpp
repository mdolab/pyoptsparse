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
  \file   Model_Stats.hpp
  \brief  Model stats (headers)
  \author Sebastien Le Digabel
  \date   2010-09-24
  \see    Model_Stats.cpp
*/
#ifndef __MODEL_STATS__
#define __MODEL_STATS__

#include "Double.hpp"

namespace NOMAD {

  /// Model stats.
  class Model_Stats {

  private:
    
    int    _nb_truth;           ///< Number of truth models.
    int    _nb_sgte;            ///< Number of surrogate models.
    int    _nb_MFN;             ///< Number of quadr. MFN models.
    int    _nb_WP_regression;   ///< Number of quadr. well-poised regressions models.
    int    _nb_regression;      ///< Number of quadr. regression_max_nY models.
    int    _nb_TGP;             ///< Number of TGP models.
    int    _not_enough_pts;     ///< Number of too small \c Y sets.
    int    _nb_Y_sets;          ///< Number of \c Y sets.
    float  _sum_nY;             ///< Total number of \c Y points.
    int    _min_nY;             ///< Minimal \c Y size.
    int    _max_nY;             ///< Maximal \c Y size.
    int    _construction_error; ///< Number of construction errors.
    double _construction_time;  ///< Models construction CPU time.
    double _optimization_time;  ///< Models optimization CPU time.

    /// Number of times that \c cond exceeded \c SVD_MAX_COND.
    int    _bad_cond;         

    // model search (MS):
    int           _MS_nb_searches;    ///< Number of searches.
    int           _MS_pts;            ///< Number of search points.
    int           _MS_success;        ///< Number of search successes.
    int           _MS_bb_eval;        ///< Number of search blackbox evaluations.
    int           _MS_sgte_eval;      ///< Number of search surrogate evaluations.
    int           _MS_cache_hits;     ///< Number of search cache hits.
    int           _MS_rejected;       ///< Number of rejected trial points.
    int           _MS_max_search_pts; ///< Max number of trial points for one search.
    int           _MS_nb_opt;         ///< Number of optimizations.
    int           _MS_opt_error;      ///< Number of optimization errors.
    int           _MS_avg_model_eval; ///< Avg number of evaluations for a model opt.
    int           _MS_max_model_eval; ///< Max number of evaluations for a model opt.
    int           _MS_max_bbe;        ///< Number of times model eval. reached limit.

    // eval sort (ES):
    int           _ES_nb_inside_radius; ///< Total number of points inside radius.
    int           _ES_nb_pts;           ///< Total number of points.

public:

    /// Constructor.
    explicit Model_Stats ( void ) { reset(); }

    /// Copy constructor.
    /**
       \param s The copied object -- \b IN.
    */
    explicit Model_Stats ( const Model_Stats & s )
      : _nb_truth            ( s._nb_truth            ) ,
	_nb_sgte             ( s._nb_sgte             ) ,
	_nb_MFN              ( s._nb_MFN              ) ,
	_nb_WP_regression    ( s._nb_WP_regression    ) ,
	_nb_regression       ( s._nb_regression       ) ,
	_nb_TGP              ( s._nb_TGP              ) ,
	_not_enough_pts      ( s._not_enough_pts      ) ,
	_nb_Y_sets           ( s._nb_Y_sets           ) ,
	_sum_nY              ( s._sum_nY              ) ,
	_min_nY              ( s._min_nY              ) ,
	_max_nY              ( s._max_nY              ) ,
	_construction_error  ( s._construction_error  ) ,
	_construction_time   ( s._construction_time   ) ,
	_optimization_time   ( s._optimization_time   ) ,
	_bad_cond            ( s._bad_cond            ) ,
	_MS_nb_searches      ( s._MS_nb_searches      ) ,
	_MS_pts              ( s._MS_pts              ) ,
	_MS_success          ( s._MS_success          ) ,
	_MS_bb_eval          ( s._MS_bb_eval          ) ,
	_MS_sgte_eval        ( s._MS_sgte_eval        ) ,
	_MS_cache_hits       ( s._MS_cache_hits       ) ,
	_MS_rejected         ( s._MS_rejected         ) ,
	_MS_max_search_pts   ( s._MS_max_search_pts   ) ,
	_MS_nb_opt           ( s._MS_nb_opt           ) ,
	_MS_opt_error        ( s._MS_opt_error        ) ,
	_MS_avg_model_eval   ( s._MS_avg_model_eval   ) ,
	_MS_max_model_eval   ( s._MS_max_model_eval   ) ,
	_MS_max_bbe          ( s._MS_max_bbe          ) ,
	_ES_nb_inside_radius ( s._ES_nb_inside_radius ) ,
	_ES_nb_pts           ( s._ES_nb_pts           ) {}

    /// Affectation operator.
    /**
       \param s The right-hand side object -- \b IN.
    */
    Model_Stats & operator = ( const Model_Stats & s );
    
    /// Destructor.
    virtual ~Model_Stats ( void ) {}

    /// Reset the stats.
    void reset ( void );

    /// Update stats from another NOMAD::Model_Stats object.
    /**
       \param s The other NOMAD::Model_Stats object -- \b IN.
    */
    void update ( const Model_Stats & s );

    /// Update stats on the size of interpolation set \c Y.
    /**
       \param nY Size of an interpolation set \c Y -- \b IN.
    */
    void update_nY ( int nY );

    /// Update stat \c _MS_max_search_pts.
    /**
       \param sp Number of trial points for a model search -- \b IN.
    */
    void update_MS_max_search_pts ( int sp );

    /// Update stats \c _MS_nb_opt, \c _MS_max_model_eval, and \c _MS_avg_model_eval.
    /**
       \param eval Number of model evaluations during one model optimization -- \b IN.
    */
    void update_MS_model_opt ( int eval );

    /// Update model ordering stats.
    /**
       \param nb_inside_radius Number of points inside radius -- \b IN.
       \param nb_pts           Number of points               -- \b IN.
    */
    void update_ES_stats ( int nb_inside_radius , int nb_pts )
    {
      _ES_nb_inside_radius += nb_inside_radius;
      _ES_nb_pts           += nb_pts;
    }

    /// Add a real to stat \c _construction_time.
    /**
       \param t Time -- \b IN.
    */
    void add_construction_time ( double t ) { _construction_time += t; }

    /// Add a real to stat \c _optimization_time.
    /**
       \param t Time -- \b IN.
    */
    void add_optimization_time ( double t ) { _optimization_time += t; }

    /// Add \c 1 to stat \c _nb_truth.
    void add_nb_truth ( void ) { ++_nb_truth; }

    /// Add \c 1 to stat \c _nb_sgte.
    void add_nb_sgte ( void ) { ++_nb_sgte; }

    /// Add \c 1 to stat \c _nb_MFN.
    void add_nb_MFN ( void ) { ++_nb_MFN; }

    /// Add \c 1 to stat \c _nb_WP_regression.
    void add_nb_WP_regression ( void ) { ++_nb_WP_regression; }

    /// Add \c 1 to stat \c _nb_regression.
    void add_nb_regression ( void ) { ++_nb_regression; }

    /// Add \c 1 to stat \c _nb_TGP.
    void add_nb_TGP ( void ) { ++_nb_TGP; }

    /// Add \c 1 to stat \c _not_enough_pts.
    void add_not_enough_pts ( void ) { ++_not_enough_pts; }

    /// Add \c 1 to stat \c _construction_error.
    void add_construction_error ( void ) { ++_construction_error; }

    /// Add \c 1 to stat \c _bad_cond.
    void add_bad_cond ( void ) { ++_bad_cond; }

    /// Add \c 1 to stat \c _MS_nb_searches.
    void add_MS_nb_searches ( void ) { ++_MS_nb_searches; }

    /// Add \c 1 to stat \c _MS_success.
    void add_MS_success ( void ) { ++_MS_success; }

    /// Add \c 1 to stat \c _MS_opt_error.
    void add_MS_opt_error ( void ) { ++_MS_opt_error; }

    /// Add \c 1 to stat \c _MS_max_bbe.
    void add_MS_max_bbe ( void ) { ++_MS_max_bbe; }

    /// Add \c 1 to stat \c _MS_rejected.
    void add_MS_rejected ( void ) { ++_MS_rejected; }

    /// Add an integer to stat \c _MS_pts.
    /**
       \param i The integer -- \b IN.
    */
    void add_MS_pts ( int i ) { _MS_pts += i; }

    /// Add an integer to stat \c _MS_bb_eval.
    /**
       \param i The integer -- \b IN.
    */
    void add_MS_bb_eval ( int i ) { _MS_bb_eval += i; }

    /// Add an integer to stat \c _MS_sgte_eval.
    /**
       \param i The integer -- \b IN.
    */
    void add_MS_sgte_eval ( int i ) { _MS_sgte_eval += i; }

    /// Add an integer to stat \c _MS_cache_hits.
    /**
       \param i The integer -- \b IN.
    */
    void add_MS_cache_hits ( int i ) { _MS_cache_hits += i; }

    /// Access to the number of model searches.
    /*
      \return The number of model searches.
    **/
    int get_MS_nb_searches ( void ) const { return _MS_nb_searches; }

    /// Access to the number of model optimizations.
    /*
      \return The number of model optimizations.
    **/
    int get_MS_nb_opt ( void ) const { return _MS_nb_opt; }

    /// Access to the number of model search points.
    /*
      \return The number of model search points.
    **/
    int get_MS_pts ( void ) const { return _MS_pts; }

    /// Access to the number of model blackbox evaluations.
    /*
      \return The number of model blackbox evaluations.
    **/
    int get_MS_bb_eval ( void ) const { return _MS_bb_eval; }

    /// Access to the number of model search surrogate evaluations.
    /*
      \return The number of model search surrogate evaluations.
    **/
    int get_MS_sgte_eval ( void ) const { return _MS_sgte_eval; }

    /// Access to the number of model search successes.
    /*
      \return The number of model search successes.
    **/
    int get_MS_success ( void ) const { return _MS_success; }

    /// Access to the number of models.
    /*
      \return The number of models.
    **/
    int get_nb_models ( void ) const { return _nb_truth + _nb_sgte; }

    /// Access to stat \c _min_nY.
    /**
       \return The stat \c _min_nY.
    */
    int get_min_nY ( void ) const { return _min_nY; }

    /// Access to stat \c _max_nY.
    /**
       \return The stat \c _max_nY.
    */
    int get_max_nY ( void ) const { return _max_nY; }

    /// Access to the average size of interpolation sets.
    /**
       \return The average size of interpolation sets.
    */
    float get_avg_nY ( void ) const
    {
      return ( _nb_Y_sets == 0 ) ? 0 : _sum_nY / _nb_Y_sets;
    }

    /// Access to stat \c _nb_MFN.
    /**
       \return \c The stat _nb_MFN.
    */
    int get_nb_MFN ( void ) const { return _nb_MFN; }
    
    /// Access to stat \c _nb_WP_regression.
    /**
       \return \c The stat _nb_WP_regression.
    */
    int get_nb_WP_regression ( void ) const { return _nb_WP_regression; }

    /// Access to stat \c _nb_regression.
    /**
       \return \c The stat _nb_regression.
    */
    int get_nb_regression ( void ) const { return _nb_regression; }

    /// Access to stat \c _nb_TGP.
    /**
       \return \c The stat _nb_TGP.
    */
    int get_nb_TGP ( void ) const { return _nb_TGP; }

    /// Display.
    /**
       \param out The NOMAD::Display object -- \b IN.
    */
    void display ( const NOMAD::Display & out ) const;
  };

  /// Display a NOMAD::Model_Stats object.
  /**
     \param out The NOMAD::Display object                     -- \b IN.
     \param s   The NOMAD::Model_Stats object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
  */
  inline const NOMAD::Display & operator << ( const NOMAD::Display     & out ,
					      const NOMAD::Model_Stats & s     ) {
    s.display ( out );
    return out;
  }
}

#endif
