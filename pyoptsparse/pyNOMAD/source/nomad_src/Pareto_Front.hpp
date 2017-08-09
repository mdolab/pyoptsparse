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
  \file   Pareto_Front.hpp
  \brief  Pareto front (headers)
  \author Sebastien Le Digabel
  \date   2010-04-09
  \see    Pareto_Front.cpp
*/
#ifndef __PARETO_FRONT__
#define __PARETO_FRONT__

#include "Pareto_Point.hpp"

namespace NOMAD {

  /// Pareto front for two objective functions.
  /**
     Browse the front with the following instructions:
     \code
     const Eval_Point * cur = pareto_front.begin();
     while ( cur ) {
       ...
       cur = pareto_front.next();
     }
     \endcode
  */
  class Pareto_Front : private NOMAD::Uncopyable {

  private:

    /// The set of Pareto points.
    std::set<NOMAD::Pareto_Point> _pareto_pts;

    /// Iterator to browse the front with begin() and next().
    mutable std::set<NOMAD::Pareto_Point>::const_iterator _it;

  public:

    /// Constructor.
    Pareto_Front ( void ) {}

    /// Destructor.
    virtual ~Pareto_Front ( void ) {}

    /// Access to the first Pareto point.
    /**
       Used to initialize a loop on the Pareto points.
       \return A pointer to the first Pareto point and \c NULL if the front
               is empty.
    */
    const NOMAD::Eval_Point * begin ( void ) const;

    /// Access to the next Pareto point.
    /**
       Used to increment a loop on the Pareto points.
       \return A pointer to the next Pareto point and \c NULL if
               the current point is the last point in the front.
    */
    const NOMAD::Eval_Point * next  ( void ) const;

    /// Access to the number of Pareto points.
    /**
       \return The number of Pareto points.
    */
    int size ( void ) const { return static_cast<int>(_pareto_pts.size()); }

    /// Check if the front is empty.
    /**
       \return A boolean equal to \c true if the front is empty.
    */
    bool empty ( void ) const { return _pareto_pts.empty(); }

    /// Computation and access to the reference point.
    /**
       \param xj      A pointer to the reference point; Is equal to \c NULL if
                      no reference exists -- \b OUT.
       \param delta_j The \c delta stats measuring the front repartition -- \b OUT.
       \return A pointer to the reference point and
               \c NULL if there is no reference point.
    */
    NOMAD::Point * get_ref ( const NOMAD::Pareto_Point *& xj      ,
			     NOMAD::Double              & delta_j   ) const;

    /// Access to the Pareto point minimizing f2(x).
    /**
       \return A pointer to the Pareto point minimizing f2 and
               \c NULL if such a point does not exist.
    */
    const NOMAD::Eval_Point * get_best_f2 ( void ) const;

    /// Compute the stats \c delta and \c surf.
    /**
       - \c delta measures the front repartition (lower is best).
       - \c surf measures the front quality (lower is best).
       \param delta_j   The \c delta stat -- \b OUT.
       \param surf      The \c surf stat  -- \b OUT.
       \param f_bounds  NOMAD::Point with 4 values (f1_min, f1_max, f2_min, and f2_max)
                        defining bounds for f1 and f2 for the computation
			of the \c surf stat -- \b IN.
    */
    void get_delta_surf ( NOMAD::Double      & delta_j ,
			  NOMAD::Double      & surf    ,
			  const NOMAD::Point & f_bounds  ) const;       
    
    /// Insertion of a point in the Pareto front.
    /**
       \param x The point to be inserted -- \b IN.
       \return A boolean equal to \c true if the point is a new Pareto point.
    */
    bool insert ( const NOMAD::Eval_Point & x );

    /// Display the Pareto points.
    /**
       \param out The NOMAD::Display object -- \b IN.
    */
    void display ( const NOMAD::Display & out ) const;
  };

  /// Display a NOMAD::Pareto_Front object.
  /**
     \param out The NOMAD::Display object -- \b IN.
     \param pf  The NOMAD::Pareto_Front object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
  */
  inline const NOMAD::Display & operator << ( const NOMAD::Display      & out ,
					      const NOMAD::Pareto_Front & pf    ) {
    pf.display ( out );
    return out;
  }
}

#endif
