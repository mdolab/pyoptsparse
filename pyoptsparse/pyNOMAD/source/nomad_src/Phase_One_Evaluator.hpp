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
 \file   Phase_One_Evaluator.hpp
 \brief  NOMAD::Evaluator subclass for the phase one (headers)
 \author Sebastien Le Digabel
 \date   2010-04-22
 \see    Phase_One_Evaluator.cpp
 */
#ifndef __PHASE_ONE_EVALUATOR__
#define __PHASE_ONE_EVALUATOR__

#include "Evaluator.hpp"

namespace NOMAD {
    
    /// NOMAD::Evaluator subclass for the phase one.
    class Phase_One_Evaluator : public NOMAD::Evaluator {
        
    private:
        
        NOMAD::Evaluator & _basic_ev; ///< The original evaluator.
        
    public:
        
        /// Constructor.
        /**
         \param p  Parameters         -- \b IN.
         \param ev Original evaluator -- \b IN.
         */
        Phase_One_Evaluator ( const NOMAD::Parameters & p , NOMAD::Evaluator & ev )
        : NOMAD::Evaluator ( p  ) ,
        _basic_ev        ( ev )   {}
        
        /// Destructor.
        virtual ~Phase_One_Evaluator ( void ) {}
        
        /// User updates after a success.
        /**
         This virtual method is called every time a new (full) success is made.
         \param s Stats -- \b IN.
         \param x Successful point -- \b IN.
         */
        virtual void update_success ( const NOMAD::Stats & s , const NOMAD::Eval_Point & x )
        {
            _basic_ev.update_success ( s , x );
        }
        
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
                             bool                & count_eval   ) const
        {
            return _basic_ev.eval_x ( x , h_max , count_eval );
        }
        
        
        /// Evaluate the blackboxes at a given list of trial points.
        /**
         \param list_x The list of trial points -- \b IN/OUT.
         \param h_max      Maximal feasibility value \c h_max -- \b IN.
         \param list_count_eval Flags indicating if the evaluations have to be counted
         or not -- \b OUT.
         \return A boolean equal to \c false if the evaluation failed.
         */
        virtual bool eval_x ( std::list<NOMAD::Eval_Point *>   & list_x          ,
                             const NOMAD::Double & h_max      ,
                             std::list<bool>                & list_count_eval   ) const
        {
            return _basic_ev.eval_x ( list_x , h_max , list_count_eval );
        }
        
        
        /// User preprocessing of points before evaluations.
        /**
         This method is called before the evaluation of a list of points.
         \param pts List of points to preprocess -- \b IN/OUT.
         */
        virtual void list_of_points_preprocessing
        ( std::set<NOMAD::Priority_Eval_Point> & pts ) const
        {
            _basic_ev.list_of_points_preprocessing ( pts );
        }
        
        
        /// Access to the model evaluator flag.
        /**
         \return The model evaluator flag.
         */
        virtual bool is_model_evaluator ( void ) const
        {
            return _basic_ev.is_model_evaluator();
        }
        
        /// User updates after an iteration.
        /**
         This virtual method is called every time a MADS iteration is terminated.
         \param success      Success of the iteration              -- \b IN.
         \param stats        Stats                                 -- \b IN.
         \param ev_control   The NOMAD::Evaluator_Control object   -- \b IN.
         \param true_barrier Barrier for true evaluations          -- \b IN.
         \param sgte_barrier Barrier for surrogate evaluations     -- \b IN.
         \param pareto_front Pareto front                          -- \b IN.
         \param stop         Allows the user to stop the algorithm -- \b OUT.
         */
        virtual void update_iteration ( NOMAD::success_type              success      ,
                                       const NOMAD::Stats             & stats        ,
                                       const NOMAD::Evaluator_Control & ev_control   ,
                                       const NOMAD::Barrier           & true_barrier ,
                                       const NOMAD::Barrier           & sgte_barrier ,
                                       const NOMAD::Pareto_Front      & pareto_front ,
                                       bool                           & stop           ) {}
        
        
        
        /// Objective computation.
        /**
         - Compute \c f(x) from the blackbox outputs of a point.
         - Special objective for MADS phase one.
         \param x The trial point -- \b IN/OUT.
         */
        virtual void compute_f ( NOMAD::Eval_Point & x ) const;
    };
}

#endif
