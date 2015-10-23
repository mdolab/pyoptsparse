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
  \file   TGP_Model_Evaluator.hpp
  \brief  NOMAD::Evaluator subclass for TGP model optimization (headers)
  \author Sebastien Le Digabel
  \date   2011-02-17
  \see    TGP_Model_Evaluator.cpp
*/
#ifdef USE_TGP

#ifndef __TGP_MODEL_EVALUATOR__
#define __TGP_MODEL_EVALUATOR__

#include "Search.hpp"
#include "Evaluator.hpp"

namespace NOMAD {

  /// NOMAD::Evaluator subclass for quadratic model optimization.
  class TGP_Model_Evaluator : public NOMAD::Evaluator {

  private:

//     int       _n;           ///< Number of variables.
//     int       _nm1;         ///< Number of variables minus one.
//     int       _m;           ///< Number of blackbox outputs.
//     double  * _x;           ///< An evaluation point.
//     double ** _alpha;       ///< Model parameters.
//     bool      _model_ready; ///< \c true if model ready to evaluate.

    NOMAD::TGP_Model & _model; ///< The TGP model.

  public:

    /// Constructor.
    /**
       \param p     Parameters -- \b IN.
       \param model Model      -- \b IN.
    */
    TGP_Model_Evaluator ( const NOMAD::Parameters & p     ,
			  NOMAD::TGP_Model        & model   )
      : NOMAD::Evaluator ( p     ) ,        
	_model           ( model )   {}

    /// Destructor.
    virtual ~TGP_Model_Evaluator ( void ) {}

    /// Evaluate the blackboxes at a given trial point.
    /**
       \param x The trial point -- \b IN/OUT.
       \param h_max      Maximal feasibility value \c h_max -- \b IN.
       \param count_eval Flag indicating if the evaluation has to be counted
                         or not -- \b OUT.
       \return A boolean equal to \c false if the evaluation failed.
     */
    virtual bool eval_x ( NOMAD::Eval_Point   & x          ,
			  const NOMAD::Double & h_max      ,
			  bool                & count_eval   ) const;
  };
}

#endif

#endif
