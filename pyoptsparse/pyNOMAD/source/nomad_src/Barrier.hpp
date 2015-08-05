/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search - version 3.6.2        */
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
 \file   Barrier.hpp
 \brief  Barrier for constraints handling (headers)
 \author Sebastien Le Digabel
 \date   2010-04-12
 \see    Barrier.cpp
 */
#ifndef __BARRIER__
#define __BARRIER__

#include "Cache.hpp"
#include "Filter_Point.hpp"

namespace NOMAD {
	
	/// Barrier class for constraints handling.
	/**
     The barrier is basically a structure that stores evaluation points.
	 */
	class Barrier : private NOMAD::Uncopyable {
		
	private:
		
		const NOMAD::Parameters & _p;              ///< Parameters.
		NOMAD::eval_type          _eval_type;      ///< Truth or surrogate.
		NOMAD::Double             _h_max;          ///< Maximal value for \c h.
		const NOMAD::Eval_Point * _best_feasible;  ///< Best feasible solution.
		
		/// Progressive barrier reference point.
		const NOMAD::Eval_Point * _ref;
		
		std::set<NOMAD::Filter_Point> _filter;     ///< Filter.
		std::set<int>                 _prefilter;  ///< Pre-filter based on point tags.
		
		int                           _rho_leaps;  ///< Number of trigger (rho) leaps.
		
		const NOMAD::Eval_Point     * _poll_center;      ///< Primary poll center.
		const NOMAD::Eval_Point     * _sec_poll_center;  ///< Secondary poll center.
		
		/**
		 Number of constraints that have been changed to NOMAD::PEB_E
		 (Extreme status for the PEB strategy).
		 */
		int _peb_changes;
		
		/**
		 Number of times that the filter has been reseted following a PEB
		 change of status.
		 */
		int _peb_filter_reset;
		
		/**
		 List of all points that we tried to insert into the filter,
		 with PEB constraints.
		 */
		std::list<const NOMAD::Eval_Point *> _peb_lop;
		
		/// List of all points that we tried to insert into the filter.
		std::list<const NOMAD::Eval_Point *> _all_inserted;
		
		NOMAD::success_type _one_eval_succ; ///< Success for one evaluation.
		NOMAD::success_type _success;       ///< Success for a list of evaluations.
		
		/*----------------------------------------------------*/
		
		/// Insertion of a feasible point.
		/**
		 \param  x The feasible point to insert -- \b IN.
		 \return Success type of the insertion
		 (NOMAD::FULL_SUCCESS or NOMAD::UNSUCCESSFUL).
		 */
		NOMAD::success_type insert_feasible ( const NOMAD::Eval_Point & x );
		
		/// Insertion of an infeasible point.
		/**
		 \param  x The infeasible point to insert -- \b IN.
		 \return Success type of the insertion
		 (NOMAD::FULL_SUCCESS, NOMAD::UNSUCCESSFUL, or NOMAD::PARTIAL_SUCCESS).
		 */
		NOMAD::success_type insert_infeasible ( const NOMAD::Eval_Point & x );
		
		/// Change the value of \c h_max.
		/**
		 \param h_max The new value of \c h_max -- \b IN.
		 */
		void set_h_max ( const NOMAD::Double & h_max );
		
		/// Insertion into the filter.
		/**
		 \param x      The point to insert -- \b IN.
		 \param insert A boolean equal to \c true if the point has been
		 inserted into the barrier -- \b OUT.
		 */
		void filter_insertion ( const NOMAD::Eval_Point & x , bool & insert );
		
		/*----------------------------------------------------*/
		
	public:
		
		/// Exception class for barrier update error.
		class Update_Error : public NOMAD::Exception
		{
		public:
			/// Constructor.
			Update_Error ( const std::string & file ,
						  int                 line ,
						  const std::string & msg    )
			: NOMAD::Exception ( file , line , msg ) {}
		};
		
		/// Exception class for an insertion error.
		class Insert_Error : public NOMAD::Exception
		{
		public:
			/// Constructor.
			Insert_Error ( const std::string & file ,
						  int                 line ,
						  const std::string & msg    )
			: NOMAD::Exception ( file , line , msg ) {}
		};
		
		/*----------------------------------------------------*/
		
		/// Constructor.
		/**
		 \param p  Parameters         -- \b IN.
		 \param et Truth or surrogate -- \b IN.
		 */
		Barrier ( const NOMAD::Parameters & p , NOMAD::eval_type et )
		: _p               ( p                   ) ,
		_eval_type       ( et                  ) ,
		_h_max           ( p.get_h_max_0()     ) ,
		_best_feasible   ( NULL                ) ,
		_ref             ( NULL                ) ,
		_rho_leaps       ( 0                   ) ,
		_poll_center     ( NULL                ) , 
		_sec_poll_center ( NULL                ) ,
		_peb_changes     ( 0                   ) ,
		_peb_filter_reset( 0                   ) ,
		_one_eval_succ   ( NOMAD::UNSUCCESSFUL ) ,    
		_success         ( NOMAD::UNSUCCESSFUL )   {}
		
		/// Destructor.
		virtual ~Barrier ( void ) {}
		
		/// Reset the barrier.
		void reset ( void );
		
		/// Reset the success types.
		void reset_success ( void ) { _one_eval_succ = _success = NOMAD::UNSUCCESSFUL; }
		
		/// Poll center selection.
		/**
		 \param last_it_success Success of the last iteration -- \b IN.
		 */
		void select_poll_center ( NOMAD::success_type last_it_success );
		
		/// Barrier updates.
		void update_and_reset_success ( void );
		
		/// Access to the best infeasible point.
		/**
		 \return A pointer to the best infeasible point.
		 */
		const NOMAD::Eval_Point * get_best_infeasible ( void ) const;
		
		/// Access to the best feasible point.
		/**
		 \return A pointer to the best feasible point.
		 */
		const NOMAD::Eval_Point * get_best_feasible ( void ) const { return _best_feasible; }
		
		/// Access to the best infeasible point with min. violation.
		/**
		 \return A pointer to the best infeasible point with min. viol.
		 */
		const NOMAD::Eval_Point * get_best_infeasible_min_viol ( void ) const;
		
		/// Access to \c h_max.
		/**
		 \return The value of \c h_max.
		 */
		const NOMAD::Double & get_h_max ( void ) const { return _h_max; }
		
		/// Access to the poll center.
		/**
		 \return A pointer to the poll center.
		 */
		const NOMAD::Eval_Point * get_poll_center ( void ) const { return _poll_center; }
		
		/// Access to the success type of one evaluation.
		/**
		 \return The success type of one evaluation.
		 */
		NOMAD::success_type get_one_eval_succ ( void ) const { return _one_eval_succ;}
		
		/// Access to the success type of a list of evaluations.
		/**
		 \return The success type of a list of evaluations.
		 */
		NOMAD::success_type get_success ( void ) const { return _success; }
		
		/// Access to the secondary poll center.
		/**
		 \return A pointer to the secondary poll center.
		 */
		const NOMAD::Eval_Point * get_sec_poll_center ( void ) const
		{
			return _sec_poll_center;
		}
		
		/// Access the list of all inserted points.
		/**
		 \return A list with pointers to the inserted points.
		 */
		const std::list<const NOMAD::Eval_Point *> & get_all_inserted ( void ) const
		{
			return _all_inserted;
		}
		
		/// Insertion of a point into the barrier.
		/**
		 \param x The point -- \b IN.
		 */
		void insert ( const NOMAD::Eval_Point & x );
		
		/// Update the barrier from another barrier.
		/**
		 Used by the VNS search.
		 \param b The other barrier -- \b IN.
		 */
		void insert ( const Barrier & b );
		
		/// Check the PEB constraints.
		/**
		 This will change eventually the status of PEB constraints
		 from NOMAD::PEB_P to NOMAD::PEB_E.
		 \param x       Point at which the constraints are checked -- \b IN.
		 \param display A boolean equal to \c true if displays are authorized -- \b IN.
		 */
		void check_PEB_constraints ( const NOMAD::Eval_Point & x , bool display );
		
		/// Display.
		/**
		 \param out The NOMAD::Display object -- \b IN.
		 */
		void display ( const NOMAD::Display & out ) const;
		
		/// Display.
		/**
		 Uses the NOMAD::Display object of the NOMAD::Parameters class.
		 */
		void display ( void ) const { display ( _p.out() ); }
		
	};
	
	/// Display a NOMAD::Barrier object.
	/**
     \param out The NOMAD::Display object -- \b IN.
     \param b   The NOMAD::Barrier object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
	 */
	inline const NOMAD::Display & operator << ( const NOMAD::Display & out ,
											   const NOMAD::Barrier & b     )
	{
		b.display ( out );
		return out;
	}
}

#endif
