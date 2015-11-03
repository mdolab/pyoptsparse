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
  \file   Pareto_Point.hpp
  \brief  Pareto point (headers)
  \author Sebastien Le Digabel
  \date   2010-04-22
  \see    Pareto_Point.cpp
*/
#ifndef __PARETO_POINT__
#define __PARETO_POINT__

#include "Multi_Obj_Evaluator.hpp"

namespace NOMAD {

  /// Pareto point for two objective functions.
  class Pareto_Point : public NOMAD::Set_Element<NOMAD::Eval_Point> {
  
  private:

    int _w; ///< Weight.

    /// Affectation operator.
    /**
       \param p The right-hand side object -- \b IN.
    */
    Pareto_Point & operator = ( const Pareto_Point & p );

  public:

    /// Constructor.
    /**
       \param ep A pointer to an evaluation point -- \b IN.
    */
    Pareto_Point ( const NOMAD::Eval_Point * ep )
      : NOMAD::Set_Element<NOMAD::Eval_Point> ( ep ) ,
	_w                                    ( 0  )   {}

    /// Copy constructor.
    /**
       \param pp The copied object -- \b IN.
    */
    explicit Pareto_Point ( const Pareto_Point & pp )
      : NOMAD::Set_Element<NOMAD::Eval_Point> ( pp.get_element() ) ,
	_w                                    ( 0                )   {}

    /// Destructor.
    virtual ~Pareto_Point ( void ) {}

    /// Update the weight.
    /**
       A more evolved formula than \c ++w is used in order
       to avoid stagnation with large number of evaluations.
    */
    void update_w ( void ) { _w = 2 * _w + 2; }

    /// Access to the weight.
    /**
       \return The weight.
    */
    int get_w ( void ) const { return _w; }

    /// Access to the value of the first objective function.
    /**
       \return The value of the first objective function.
    */
    const NOMAD::Double & get_f1 ( void ) const
    {
      return get_element()->get_bb_outputs()[NOMAD::Multi_Obj_Evaluator::get_i1()];
    }

    /// Access to the value of the second objective function.
    /**
       \return The value of the second objective function.
    */
    const NOMAD::Double & get_f2 ( void ) const
    {
      return get_element()->get_bb_outputs()[NOMAD::Multi_Obj_Evaluator::get_i2()];
    }

    /// Comparison operator.
    /**
       - Supposes that \c y is a Pareto point.
       - Used for the insertion in a set (the Pareto front).
       - \c f1(*this) and \c f1(y) are compared.
       \param y The right-hand side of the comparison -- \b IN.
       \return A boolean equal to \c true if \c *this \c < \c y.
    */
    virtual bool operator < ( const NOMAD::Set_Element<NOMAD::Eval_Point> & y ) const;

    /// Dominance operator.
    /**
       Used for the comparison (dominance) of two Pareto points
       before they are inserted into the Pareto front.
       \param y The right-hand side of the comparison -- \b IN.
       \return A boolean equal to \c true if \c *this dominates \c y.
    */
    bool dominates ( const Pareto_Point & y ) const;
    
    /// Display the Pareto point.
    /**
       \param out The NOMAD::Display object -- \b IN.
    */
    void display ( const NOMAD::Display & out ) const;
  };

  /// Display a NOMAD::Pareto_Point object.
  /**
     \param out The NOMAD::Display object -- \b IN.
     \param pp  The NOMAD::Pareto_Point object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
  */
  inline const NOMAD::Display & operator << ( const NOMAD::Display      & out ,
					      const NOMAD::Pareto_Point & pp    )
  {
    pp.display ( out );
    return out;
  }
}

#endif
