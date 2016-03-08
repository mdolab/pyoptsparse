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
  \file   Priority_Eval_Point.hpp
  \brief  Evaluation point with a priority (headers)
  \author Sebastien Le Digabel
  \date   2010-04-22
  \see    Priority_Eval_Point.cpp
*/
#ifndef __PRIORITY_EVAL_POINT__
#define __PRIORITY_EVAL_POINT__

#include "Set_Element.hpp"
#include "Eval_Point.hpp"

namespace NOMAD {
	
	/// Evaluation point with a priority.
	class Priority_Eval_Point : public NOMAD::Set_Element<NOMAD::Eval_Point> {
		
	private:
		
		NOMAD::Double _h_min;              ///< \c h_min value for comparison operator.
		NOMAD::Double _f_sgte;             ///< Objective surrogate value.
		NOMAD::Double _h_sgte;             ///< Feasibility surrogate value.
		NOMAD::Double _f_model;            ///< Objective model value.
		NOMAD::Double _h_model;            ///< Feasibility model value.
		NOMAD::Double _angle_success_dir;  ///< Angle with last successful direction.
		NOMAD::Double _angle_simplex_grad; ///< Angle with simplex gradient.
		static bool	  _lexicographic_order; ///< Use lexicographic order for comparison 
		
		/// Affectation operator.
		/**
		 \param x The right-hand side object -- \b IN.
		 */
		Priority_Eval_Point & operator = ( const Priority_Eval_Point & x );
		
		/// Compare the \c h values of two points.
		/**
		 The two points to compare are \c x1 and \c x2.
		 \param hx1 \c h(x1) -- \b IN.
		 \param hx2 \c h(x2) -- \b IN.
		 \return \c h(x1) \c < \c h(x2)
		 with the following format:
		 -  1: \c x1 best than \c x2.
		 - -1: \c x2 best than \c x1.
		 -  0: undetermined.
		 */
		int compare_h_values ( const NOMAD::Double & hx1 ,
							  const NOMAD::Double & hx2   ) const;
		
		/// Compare the \c h and \c f values of two points.
		/**
		 The two points to compare are \c x1 and \c x2.
		 \param hx1 \c h(x1) -- \b IN.
		 \param fx1 \c f(x1) -- \b IN.
		 \param hx2 \c h(x2) -- \b IN.
		 \param fx2 \c f(x2) -- \b IN.
		 \return \c (h(x1),f(x1)) \c < \c (h(x2),f(x2))
		 with the following format:
		 -  1: \c x1 best than \c x2.
		 - -1: \c x2 best than \c x1.
		 -  0: undetermined.
		 */
		int compare_hf_values ( const NOMAD::Double & hx1 ,
							   const NOMAD::Double & fx1 ,
							   const NOMAD::Double & hx2 ,
							   const NOMAD::Double & fx2   ) const;
	public:
		
		/// Constructor.
		/**
		 \param x A pointer to the evaluation point -- \b IN.
		 \param h_min \c h_min value                -- \b IN.
		 */
		explicit Priority_Eval_Point ( const NOMAD::Eval_Point * x     ,
									  const NOMAD::Double     & h_min   )
		: NOMAD::Set_Element<NOMAD::Eval_Point> ( x     ) ,
		_h_min                                ( h_min )  {}
		
		/// Copy constructor.
		/**
		 \param pep The copied object -- \b IN.
		 */
		explicit Priority_Eval_Point ( const Priority_Eval_Point & pep )
		: NOMAD::Set_Element<NOMAD::Eval_Point> ( pep.get_element()       ) ,
		_h_min                                ( pep._h_min              ) ,
		_f_sgte                               ( pep._f_sgte             ) ,
		_h_sgte                               ( pep._h_sgte             ) ,
		_f_model                              ( pep._f_model            ) ,
		_h_model                              ( pep._h_model            ) ,
		_angle_success_dir                    ( pep._angle_success_dir  ) ,
		_angle_simplex_grad                   ( pep._angle_simplex_grad )   {}
		
		/// Destructor.
		virtual ~Priority_Eval_Point ( void ) {}
		
		/// Access to specific elements of comparison.
		/**
		 - This method is defined virtual in NOMAD::Set_Element so that
         \c operator \c < \c (Set_Element x) can invoke
         it on \c x (which is in fact a \c Priority_Eval_Point).
		 - This avoids an expensive downcast in \c operator \c < .
		 \param f_sgte              Objective surrogate value            -- \b OUT.
		 \param h_sgte              Feasibility surrogate value          -- \b OUT.
		 \param f_model             Objective model value                -- \b OUT.
		 \param h_model             Feasibility model value              -- \b OUT.
		 \param angle_success_dir   Angle with last successful direction -- \b OUT.
		 \param angle_simplex_grad  Angle with simplex gradient          -- \b OUT.
		 */
		virtual void get_priority_criteria ( NOMAD::Double & f_sgte             ,
											NOMAD::Double & h_sgte             ,
											NOMAD::Double & f_model            ,
											NOMAD::Double & h_model            ,
											NOMAD::Double & angle_success_dir  ,
											NOMAD::Double & angle_simplex_grad   ) const
		{
			f_sgte             = _f_sgte;
			h_sgte             = _h_sgte;
			f_model            = _f_model;
			h_model            = _h_model;
			angle_success_dir  = _angle_success_dir;
			angle_simplex_grad = _angle_simplex_grad;
		}
		
		/// Comparison operator.
		/**
		 This virtual function directly call \c dominates().
		 \param x The right-hand side object -- \b IN.
		 \return A boolean equal to \c true if \c *this \c < \c x.
		 */
		virtual bool operator < ( const NOMAD::Set_Element<NOMAD::Eval_Point> & x ) const
		{ return dominates ( x ); }
		
		/// Comparison operator.
		/**
		 \param x The right-hand side object -- \b IN.
		 \return A boolean equal to \c true if \c *this \c < \c x.
		 */
		bool dominates ( const NOMAD::Set_Element<NOMAD::Eval_Point> & x ) const;
		
		/// Access to the evaluation point.
		/**
		 \return A pointer to the evaluation point.
		 */
		const NOMAD::Eval_Point * get_point ( void ) const { return get_element(); }
		
		/// Set the angle with last successful direction.
		/**
		 \param a The angle with last successful direction -- \b IN.
		 */
		void set_angle_success_dir ( const NOMAD::Double & a ) { _angle_success_dir = a; }
		
		/// Set the angle with simplex gradient .
		/**
		 \param a The angle with simplex gradient -- \b IN.
		 */
		void set_angle_simplex_grad ( const NOMAD::Double & a ) { _angle_simplex_grad = a; }
		
		
		/// Set the objective surrogate value.
		/**
		 \param f The objective surrogate value -- \b IN.
		 */
		void set_f_sgte ( const NOMAD::Double & f ) { _f_sgte = f; }
		
		/// Set the feasibility surrogate value.
		/**
		 \param h The feasibility surrogate value -- \b IN.
		 */
		void set_h_sgte ( const NOMAD::Double & h ) { _h_sgte = h; }
		
		/// Set the objective model value.
		/**
		 \param f The objective model value -- \b IN.
		 */
		void set_f_model ( const NOMAD::Double & f ) { _f_model = f; }
		
		/// Set the feasibility model value.
		/**
		 \param h The feasibility model value -- \b IN.
		 */
		void set_h_model ( const NOMAD::Double & h ) { _h_model = h; }
		
		/// Set the lexicographic order for sorting.
		/**
		 */
		static void set_lexicographic_order ( bool order ) { _lexicographic_order = order; }
		
		
	};
}

#endif
