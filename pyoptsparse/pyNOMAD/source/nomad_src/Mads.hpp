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
 \file   Mads.hpp
 \brief  MADS algorithm (headers)
 \author Sebastien Le Digabel
 \date   2010-04-20
 \see    Mads.cpp
 */
#ifndef __MADS__
#define __MADS__

#include "Quad_Model_Search.hpp"
#include "Speculative_Search.hpp"
#include "LH_Search.hpp"
#include "VNS_Search.hpp"
#include "Cache_Search.hpp"
#include "Phase_One_Search.hpp"
#include "L_Curve.hpp"
#include "Extended_Poll.hpp"

#include "XMesh.hpp"
#include "SMesh.hpp"

#ifdef USE_TGP
#include "TGP_Model_Search.hpp"
#endif

namespace NOMAD {
    
    // forward declaration of Extended_Poll:
    class Extended_Poll;
    
    /// The MADS algorithm.
    class Mads : private NOMAD::Uncopyable {
        
    private:
        
        static bool _force_quit; ///< Forces NOMAD to terminate if Ctrl-C is pressed.
        
        NOMAD::Parameters      & _p;             ///< Parameters.
        NOMAD::Stats             _stats;         ///< Statistics.
        NOMAD::Evaluator_Control _ev_control;    ///< Evaluator control.
        NOMAD::Evaluator_Control _ev_control_for_sorting;    ///< Evaluator control.
        NOMAD::Barrier           _true_barrier;  ///< Barrier for true function evaluations.
        NOMAD::Barrier           _sgte_barrier;  ///< Barrier for surrogate evaluations.
        
        NOMAD::OrthogonalMesh * _mesh;   ///< Access to the OrthogonalMesh
        
        /// Pareto front for multi-objective optimization.
        NOMAD::Pareto_Front  * _pareto_front;
        
        /// User search defined with NOMAD::Mads::set_user_search().
        NOMAD::Search        * _user_search;
        
        NOMAD::Search        * _model_search1;  ///< Model search #1.
        NOMAD::Search        * _model_search2;  ///< Model search #2.
        
        NOMAD::Search        * _VNS_search;     ///< VNS search.
        NOMAD::Search        * _cache_search;   ///< Cache search.
        NOMAD::L_Curve       * _L_curve;        ///< L-curve target.
        NOMAD::Extended_Poll * _extended_poll;  ///< Extended poll for categorical variables.
        bool                   _user_ext_poll;  ///< Flag for user-defined extended poll.
        
        // MADS flags (these are not parameters as users do not modify them):
        
        static bool _flag_check_bimads;   ///< Flag for the multi-objective test in \c run().
        static bool _flag_reset_mesh;     ///< Reset or not the mesh before a MADS run.
        static bool _flag_reset_barriers; ///< Reset or not the barriers before a MADS run.
        static bool _flag_p1_active;      ///< Flag equal to \c true if phase one is active.
        
        /*-----------------------------------------------------------------------------*/
        
        /// Initializations.
        void init ( void );
        
        /// Starting point evaluation.
        /**
         \param stop        Stop flag   -- \b IN/OUT.
         \param stop_reason Stop reason -- \b OUT.
         */
        void eval_x0 (  bool & stop , NOMAD::stop_type & stop_reason );
        
        /// One MADS iteration.
        /**
         \param stop           Stop flag                               -- \b IN/OUT.
         \param stop_reason    Stop reason                             -- \b OUT.
         \param success        Success for this iteration              -- \b OUT.
         \param new_feas_inc   Pointer to the new feasible incumbent   -- \b OUT.
         \param new_infeas_inc Pointer to the new infeasible incumbent -- \b OUT.
         */
        void iteration ( bool                     & stop           ,
                        NOMAD::stop_type         & stop_reason    ,
                        NOMAD::success_type      & success        ,
                        const NOMAD::Eval_Point *& new_feas_inc   ,
                        const NOMAD::Eval_Point *& new_infeas_inc   );
        
        /// The search step.
        /**
         \param stop           Stop flag                               -- \b IN/OUT.
         \param stop_reason    Stop reason                             -- \b OUT.
         \param success        Success for this step                   -- \b OUT.
         \param new_feas_inc   Pointer to the new feasible incumbent   -- \b OUT.
         \param new_infeas_inc Pointer to the new infeasible incumbent -- \b OUT.
         */
        void search ( bool                     & stop           ,
                     NOMAD::stop_type         & stop_reason    ,
                     NOMAD::success_type      & success        ,
                     const NOMAD::Eval_Point *& new_feas_inc   ,
                     const NOMAD::Eval_Point *& new_infeas_inc   );
        
        /// The poll step.
        /**
         \param stop           Stop flag                               -- \b IN/OUT.
         \param stop_reason    Stop reason                             -- \b OUT.
         \param success        Success for this step                   -- \b OUT.
         \param new_feas_inc   Pointer to the new feasible incumbent   -- \b OUT.
         \param new_infeas_inc Pointer to the new infeasible incumbent -- \b OUT.
         \param forbid_poll_size_stop Boolean used to check if the min poll
         size stopping criterion has to be
         disabled for integer variables -- \b OUT.
         */
        void poll ( bool					& stop,
                   NOMAD::stop_type		& stop_reason,
                   NOMAD::success_type		& success,
                   const NOMAD::Eval_Point *& new_feas_inc,
                   const NOMAD::Eval_Point *& new_infeas_inc,
                   bool					& forbid_poll_size_stop   );
        
        /// Sets the poll trial points from poll direction, poll center and mesh size
        /**
         \param  stop          Stop flag                                         -- \b IN/OUT.
         \param  stop_reason   Stop reason                                       -- \b OUT.
         \param  dirs          List of directions for the poll                   -- \b IN.
         \param  poll_center   the poll center (can be null)                     -- \b IN.
         \param  offset        Dir index offset for primary and sec. poll center -- \b IN.
         \param  sorting       If true than the points are for sorting           -- \b IN.
         */
        void set_poll_trial_points (  std::list<NOMAD::Direction> & dirs,
                                    size_t                           offset,
                                    const NOMAD::Eval_Point     &  poll_center,
                                    bool					    	& stop,
                                    NOMAD::stop_type				&stop_reason,
                                    bool							sorting);
        
        
        /// Compute a prospect point by optimization on quadratic models.
        /**
         \param  poll_center The poll center       -- \b IN.
         \param  dirs The directions that delimit the hypercube for optimization       -- \b IN.
         \param  prospect_point  The prospect point  -- \b OUT.
         \return A flag equal to \c true if the prospect direction has been computed.
         */
        bool optimize_quad_model ( const NOMAD::Eval_Point         & poll_center ,
                                  const std::list<NOMAD::Direction> & dirs    ,
                                  NOMAD::Point                    & prospect_point    )  ;
        
        
        /// Sets the poll directions from signature, poll center and mesh size
        /**
         \param dirs			List of directions for the poll			          -- \b OUT.
         \param i_pc			Poll type                                         -- \b IN.
         \param offset		Dir index offset for primary and sec. poll center -- \b IN.
         \param poll_center   The poll center                                   -- \b IN.
         \param stop 			Stop flag, true if cannot get direction   		  -- \b IN/OUT.
         \param stop_reason	Stop type										  -- \b OUT.
         */
        void set_poll_directions ( std::list<NOMAD::Direction> & dirs        ,
                                  NOMAD::poll_type              i_pc        ,
                                  size_t                        offset      ,
                                  const NOMAD::Eval_Point     & poll_center ,
                                  bool                        & stop        ,
                                  NOMAD::stop_type            & stop_reason   );
        
        /// Reduce the number of poll directions -> n
        /**
         \param dirs			List of directions for the poll			-- \b IN/OUT.
         \param  poll_center   the poll center                         -- \b IN.
         \return success for this step.
         */
        bool set_reduced_poll_to_n_directions(std::list<NOMAD::Direction>	& dirs,
                                              const NOMAD::Eval_Point		& poll_center);
        
        /// Compute the rank of a list of directions
        /**
         \param dirs List of directions for the poll -- \b IN/OUT.
         \return     rank>0 of the dirs if successfull or 0 if rank cannot be evaluated
         */
        int get_rank_from_dirs(const std::list<NOMAD::Direction> & dirs);
        
        
        /// Check the directions after the poll step.
        /**
         \param forbid_poll_size_stop Boolean equal to \c true if
         the \c MIN_POLL_SIZE parameter is valid for
         the last set of poll directions -- \b OUT.
         */
        void check_directions ( bool & forbid_poll_size_stop );
        
        /// Update of the success directions.
        /**
         - Occurs after the poll.
         \param new_inc    Pointer to the new incumbent -- \b IN (may be \c NULL).
         \param feasible   Flag equal to \c true if the incumbent is feasible -- \b IN.
         */
        void update_success_directions ( const NOMAD::Eval_Point         * new_inc    ,
                                        bool                              feasible  ) const;
        
        /// Launch a single-objective MADS run for multi-objective optimization.
        /**
         \param display_degree Display degree                                -- \b IN.
         \param mads_runs      Total number of MADS runs to execute          -- \b IN.
         \param overall_bbe    Global maximum number of blackbox evaluations -- \b IN.
         \param ev             Evaluator for multi-objective                 -- \b IN/OUT.
         \param stagnation_cnt Counter to detect a stagnation                -- \b IN/OUT.
         \param multi_stats    Stats for multi-objective                     -- \b IN/OUT.
         \param stop           Stop flag                                     -- \b IN/OUT.
         \param stop_reason    Stop reason                                   -- \b OUT.
         */
        void multi_launch_single_opt ( NOMAD::dd_type               display_degree ,
                                      int                          mads_runs      ,
                                      int                          overall_bbe    ,
                                      NOMAD::Multi_Obj_Evaluator & ev             ,
                                      int                        & stagnation_cnt ,
                                      NOMAD::Stats               & multi_stats    ,
                                      bool                       & stop           ,
                                      NOMAD::stop_type           & stop_reason      );
        
        /// Compute and set the minimal poll size for multi-objective optimization.
        /**
         \param lb        Lower bounds                        -- \b IN.
         \param ub        Upper bounds                        -- \b IN.
         \param delta_p_0 Initial poll size                   -- \b IN.
         \param delta_j   Delta criterion for multi-objective -- \b IN/OUT.
         */
        void multi_set_min_poll_size ( const NOMAD::Point & lb        ,
                                      const NOMAD::Point & ub        ,
                                      const NOMAD::Point & delta_p_0 ,
                                      NOMAD::Double        delta_j     );
        
        // Display mesh and poll sizes for a given signature.
        /**
         \param s The signature -- \b IN.
         */
        void display_deltas ( const NOMAD::Signature & s ) const;
        
        /// Displays at the beginning of an iteration.
        void display_iteration_begin ( void ) const;
        
        /// Displays at the end of an iteration.
        /**
         \param stop           Stop flag                               -- \b IN.
         \param stop_reason    Stop reason                             -- \b IN.
         \param success        Iteration success                       -- \b IN.
         \param new_feas_inc   Pointer to the new feasible incumbent   -- \b IN.
         \param new_infeas_inc Pointer to the new infeasible incumbent -- \b IN.
         */
        void display_iteration_end ( bool                      stop           ,
                                    NOMAD::stop_type          stop_reason    ,
                                    NOMAD::success_type       success        ,
                                    const NOMAD::Eval_Point * new_feas_inc   ,
                                    const NOMAD::Eval_Point * new_infeas_inc   ) const;
        
        
        /// Determine dynamic directions from a set of directions.
        /**
         - The computed opposite directions already include Delta^k_m.
         \param dirs          List of existing directions (no snap to bounds) -- \b IN.
         \param newDirs		New dynamic directions              -- \b OUT.
         \param poll_center   Poll center                         -- \b IN.
         \return true if new dynamic direction generated false otherwise
         */
        bool get_dynamic_directions (const std::list<NOMAD::Direction>	&	dirs,
                                     std::list<NOMAD::Direction>			&	newDirs,
                                     const NOMAD::Eval_Point				&	poll_center) ;
        
        
        
        /// Check if a set of directions include Ortho-MADS N+1 direction.
        /**
         \param dirs List of directions -- \b IN.
         \return A boolean equal to \c true if at
         least one direction in the set is
         of type Ortho-MADS N+1.
         */
        bool dirs_have_orthomads_np1 ( const std::list<NOMAD::Direction> & dirs );
        
        
        ///  Check if a dir needs to be obtained from model optimization
        /**
         \param dirs List of directions -- \b IN.
         \return A boolean equal to \c true if all directions are of type Ortho-MADS N+1 QUAD.
         */
        bool dir_from_model_opt( const std::list<NOMAD::Direction> & dirs);
        
        
        /// get a single direction using quad model optimization or sum of negatives
        /**
         \param dirs			Reduced poll directions	(no snap to bounds)	-- \b IN.
         \param poll_center	Poll center								    -- \b IN.
         \return new direction
         */
        NOMAD::Direction get_single_dynamic_direction (const std::list<NOMAD::Direction>	&	dirs,
                                                       const NOMAD::Eval_Point			&	poll_center) ;
        
        
        /*-----------------------------------------------------------------------------*/
        
    public:
        
        /// Constructor #1.
        /**
         - Basic version.
         \param p  Parameters                 -- \b IN.
         \param ev A pointer to the evaluator -- \b IN
         -- \b optional (default = \c NULL).
         */
        Mads ( NOMAD::Parameters & p , NOMAD::Evaluator * ev = NULL )
        : _p                   ( p                             ) ,
        _stats                 ( p.get_sgte_cost()             ) ,
        _ev_control            ( p , _stats , ev , NULL , NULL ) ,
        _ev_control_for_sorting( p , _stats , _ev_control.get_evaluator() , &(_ev_control.get_cache()) , &(_ev_control.get_sgte_cache()) ) ,
        _true_barrier          ( p , NOMAD::TRUTH              ) ,
        _sgte_barrier          ( p , NOMAD::SGTE               ) ,
        _mesh				   ( p.get_signature()->get_mesh() ) ,
        _pareto_front          ( NULL                          ) ,
        _user_search           ( NULL                          ) ,
        _model_search1         ( NULL                          ) ,
        _model_search2         ( NULL                          ) ,
        _VNS_search            ( NULL                          ) ,
        _cache_search          ( NULL                          ) ,
        _L_curve               ( NULL                          ) ,
        _extended_poll         ( NULL                          ) ,
        _user_ext_poll         ( false                         )   { init(); }
        
        /// Constructor #2.
        /**
         - Advanced version.
         \param p             Parameters                 -- \b IN.
         \param ev            A pointer to the evaluator -- \b IN (may be \c NULL).
         \param extended_poll A pointer to a NOMAD::Extended_Poll object
         -- \b IN (may be \c NULL).
         \param cache         A pointer to a cache       -- \b IN (may be \c NULL).
         \param sgte_cache    A pointer to a cache for surrogates
         -- \b IN (may be \c NULL).
         */
        Mads ( NOMAD::Parameters    & p             ,
              NOMAD::Evaluator     * ev            ,     // may be NULL
              NOMAD::Extended_Poll * extended_poll ,     // may be NULL
              NOMAD::Cache         * cache         ,     // may be NULL
              NOMAD::Cache         * sgte_cache      )   // may be NULL
        : _p                     ( p                                    ) ,
        _stats                 ( p.get_sgte_cost()                    ) ,
        _ev_control            ( p , _stats , ev , cache , sgte_cache ) ,
        _ev_control_for_sorting( p , _stats , _ev_control.get_evaluator() , cache , sgte_cache ) ,
        _true_barrier          ( p , NOMAD::TRUTH                     ) ,
        _sgte_barrier          ( p , NOMAD::SGTE                      ) ,
        _mesh                  ( p.get_signature()->get_mesh()		  ) ,
        _pareto_front          ( NULL                                 ) ,
        _user_search           ( NULL                                 ) ,
        _model_search1         ( NULL                                 ) ,
        _model_search2         ( NULL                                 ) ,
        _VNS_search            ( NULL                                 ) ,
        _cache_search          ( NULL                                 ) ,
        _L_curve               ( NULL                                 ) ,
        _extended_poll         ( extended_poll                        ) ,
        _user_ext_poll         ( (extended_poll!=NULL)                )   { init(); }
        
        /// Destructor.
        virtual ~Mads ( void );
        
        /// Algorithm execution for single-objective.
        /**
         \return Stop reason.
         */
        NOMAD::stop_type run ( void );
        
        /// Algorithm execution for multi-objective.
        /**
         \return Stop reason.
         */
        NOMAD::stop_type multi_run ( void );
        
        /// Force quit.
        /**
         Called by pressing Ctrl-C.
         \param signalValue Signal value -- \b IN.
         */
        static void force_quit ( int signalValue );
        
        /// Reset.
        /**
         - Also resets the user search.
         \param keep_barriers A boolean equal to \c true if NOMAD::Barrier objects
         have to be reseted
         -- \b IN -- \b optional (default = \c false).
         \param keep_stats    A boolean equal to \c true if the stats object
         has to be reseted
         -- \b IN -- \b optional (default = \c false).
         */
        void reset ( bool keep_barriers = false , bool keep_stats = false );
        
        /// Set user search.
        /**
         \param us A pointer to the user search -- \b IN.
         */
        void set_user_search  ( NOMAD::Search * us ) { _user_search  = us; }
        
        /// Set an extern Pareto front.
        /**
         \param pf A pointer to a Pareto front -- \b IN.
         */
        void set_pareto_front ( NOMAD::Pareto_Front * pf ) { _pareto_front = pf; }
        
        /// Set the flag for the multi-objective test.
        /**
         \param fcb The flag for the multi-objective test -- \b IN.
         */
        static void set_flag_check_bimads ( bool fcb ) { _flag_check_bimads = fcb; }
        
        /// Set the flag \c _flag_reset_mesh.
        /**
         \param frm The flag -- \b IN.
         */
        static void set_flag_reset_mesh ( bool frm ) { _flag_reset_mesh = frm; }
        
        /// Set the flag \c _flag_reset_barriers.
        /**
         \param frb The flag -- \b IN.
         */
        static void set_flag_reset_barriers ( bool frb ) { _flag_reset_barriers = frb; }
        
        /// Set the flag \c _flag_p1_active -- \b IN.
        /**
         \param fpa The flag.
         */
        static void set_flag_p1_active ( bool fpa ) { _flag_p1_active = fpa; }
        
        /// Access to the flags.
        /**
         \param flag_check_bimads   Multi-objective flag -- \b OUT.
         \param flag_reset_mesh     Mesh reset flag      -- \b OUT.
         \param flag_reset_barriers Reset barriers flag  -- \b OUT.
         \param flag_p1_active      Phase one flag       -- \b OUT.
         */
        static void get_flags ( bool & flag_check_bimads   ,
                               bool & flag_reset_mesh     ,
                               bool & flag_reset_barriers ,
                               bool & flag_p1_active        );
        
        /// Access to the stats.
        /**
         \return The stats.
         */
        NOMAD::Stats & get_stats ( void ) { return _stats; }
        
        /// Access to the evaluator control.
        /**
         \return The evaluator control.
         */
        NOMAD::Evaluator_Control & get_evaluator_control ( void ) { return _ev_control; }
        
        /// Access to the barrier for true function evaluations.
        /**
         \return The barrier for the true function evaluations.
         */
        NOMAD::Barrier & get_true_barrier ( void ) { return _true_barrier; }
        
        /// Access to the barrier for surrogate evaluations.
        /**
         \return The barrier for the surrogates.
         */
        NOMAD::Barrier & get_sgte_barrier ( void ) { return _sgte_barrier; }
        
        /// Access to the NOMAD::Extended_Poll object.
        /**
         \return A pointer to \c _extended_poll.
         */
        NOMAD::Extended_Poll * get_extended_poll ( void ) const { return _extended_poll; }
        
        
        
        /// Access to the Pareto front.
        /**
         \return A pointer to the Pareto front.
         */
        NOMAD::Pareto_Front * get_pareto_front ( void ) const { return _pareto_front; }
        
        /// Access to the active cache (truth or surrogate).
        /**
         \return The active cache.
         */
        const NOMAD::Cache & get_cache ( void ) const
        {
            return ( _p.get_opt_only_sgte() ) ?
            _ev_control.get_sgte_cache() : _ev_control.get_cache();
        }
        
        /// Access to the active barrier (truth or surrogate).
        /**
         \return The active barrier.
         */
        const NOMAD::Barrier & get_active_barrier ( void ) const
        {
            return ( _p.get_opt_only_sgte() ) ? _sgte_barrier : _true_barrier;
        }
        
        /// Access to the best feasible point.
        /**
         \return A pointer to the best feasible point;
         \return \c NULL if there is no feasible point.
         */
        const NOMAD::Eval_Point * get_best_feasible ( void ) const
        {
            return get_active_barrier().get_best_feasible();
        }
        
        /// Access to the best infeasible point.
        /**
         \return A pointer to the best infeasible point;
         \return \c NULL if there is no infeasible point.
         */
        const NOMAD::Eval_Point * get_best_infeasible( void ) const
        {
            return get_active_barrier().get_best_infeasible();
        }
        
        /// Access to the best infeasible point with minimun constraint violation.
        /**
         \return A pointer to the best infeasible point with min. viol.;
         \return \c NULL if there is no infeasible point.
         */
        const NOMAD::Eval_Point * get_best_infeasible_min_viol ( void ) const
        {
            return get_active_barrier().get_best_infeasible_min_viol();
        }
        
        /// Display model stats.
        /**
         \param out The NOMAD::Display object -- \b IN.
         */
        void display_model_stats ( const NOMAD::Display & out ) const;
        
        /// Display the Pareto front.
        /**
         Displays the front at the standard output or in a stats file.
         */
        void display_pareto_front ( void ) const;
        
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
    
    /// Display a NOMAD::Mads object.
    /**
     \param out The NOMAD::Display object -- \b IN.
     \param m   The NOMAD::Mads object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
     */
    inline const NOMAD::Display & operator << ( const NOMAD::Display & out ,
                                               const NOMAD::Mads    & m     )
    {
        m.display ( out );
        return out;
    }
}

#endif
