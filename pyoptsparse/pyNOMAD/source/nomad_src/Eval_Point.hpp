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
 \file   Eval_Point.hpp
 \brief  Evaluation point (headers)
 \author Sebastien Le Digabel
 \date   2010-04-14
 \see    Eval_Point.cpp
 */
#ifndef __EVAL_POINT__
#define __EVAL_POINT__

#include "Parameters.hpp"
#include "Cache_File_Point.hpp"
#include "Set_Element.hpp"

#ifdef WINDOWS
#ifndef isnan
inline bool isnan ( double x ) { return x != x; }
#endif
#endif

namespace NOMAD {
    
    /// Class for the representation of an evaluation point.
    /**
     An evaluation point gathers the point coordinates \c x, and the blackbox
     outputs at these coordinates \c f(x).
     */
    class Eval_Point : public NOMAD::Point {
        
    private:
        
        static int         _current_tag; ///< Current tag for all NOMAD::Eval_Point objects.
        int                _tag;         ///< Unique tag.
        
        static int			_current_bbe;		///< Current bbe for all NOMAD::Eval_Point objects
        static int			_current_sgte_bbe;  ///< Current sgte_bbe for all NOMAD::Eval_Point objects  .
        int					_bbe;				/// Recorded bbe for current NOMAD::Eval_Point object
        int					_sgte_bbe;		    /// Recorded sgte_bbe for current NOMAD::Eval_Point object
        int					_real_time;			/// Recorded time for current NOMAD::Eval_Point object
        
        NOMAD::Signature * _signature;   ///< Signature of the point.
        
        /**
         \c f is explicitely computed by a NOMAD::Evaluator object
         and is not saved in the cache file.
         */
        NOMAD::Double _f;
        
        /**
         \c h is explicitely computed by a NOMAD::Evaluator object
         and is not saved in the cache file.
         */
        NOMAD::Double _h;
        
        /// Flag equal to \c true if this eval point is in the cache.
        mutable bool _in_cache;
        
        /**
         Flag equal to \c true if the point has been evaluated
         during the current run.
         */
        mutable bool _current_run;
        
        /// Type of the evaluation (true function or surrogate).
        NOMAD::eval_type _eval_type;
        
        /// Direction from which the point has been constructed.
        /**
         May be \c NULL if no direction has been used.
         */
        NOMAD::Direction * _direction;
        
        
        /// Poll center from which the point has been constructed.
        /**
         May be \c NULL if no direction has been used.
         */
        const NOMAD::Eval_Point * _poll_center;
        
        /**
         Type of the poll center (feasible or not)
         from which the point has been constructed.
         */
        NOMAD::poll_center_type _poll_center_type;
        
        /// Evaluation status.
        /**
         Indicates if the evaluation failed, succeeded,
         or is in progress.
         */
        NOMAD::eval_status_type _eval_status;
        
        /// Flag equal to \c true if all EB constraints have been satisfied.
        bool _EB_ok;
        
        /// Blackbox outputs.
        NOMAD::Point _bb_outputs;
        
        /// User evaluation priority.
        /**
         - Decided by the user in his custom
         NOMAD::Evaluator::eval_x() function, via
         \c set_user_eval_priority().
         - Points with higher priorities
         will be evaluated first.
         */
        NOMAD::Double _user_eval_priority;
        
        /// Random evaluation priority.
        /**
         Same idea than \c _user_eval_priority
         for a random ordering of trial points.
         */
        NOMAD::Double _rand_eval_priority;
        
        /// Affectation operator.
        /**
         \param x The right-hand side object -- \b IN.
         \return \c *this as the result of the affectation.
         */
        Eval_Point & operator = ( const Eval_Point & x );
        
    public:
        
        /// Constructor #1.
        explicit Eval_Point ( void );
        
        /// Constructor #2.
        /**
         \param n Number of variables        -- \b IN.
         \param m Number of blackbox outputs -- \b IN.
         */
        Eval_Point ( int n , int m );
        
        /// Constructor #3.
        /**
         \param x  Coordinates of a point taken in a cache    -- \b IN.
         \param et Type of the evaluation (true or surrogate) -- \b IN.
         */
        explicit Eval_Point ( const NOMAD::Cache_File_Point & x  ,
                             NOMAD::eval_type                et   );
        
        
        ///Constructor 4
        /**
         \param x  Coordinates of a point                  -- \b IN.
         \param m Number of blackbox outputs -- \b IN.
         */
        explicit Eval_Point ( const NOMAD::Point & x, int m ) ;
        
        
        
        /// Copy constructor.
        /**
         \param x The copied object.
         */
        explicit Eval_Point ( const Eval_Point & x );
        
        /// Destructor.
        virtual ~Eval_Point ( void );
        
        /// Check the evaluation point.
        /**
         \param m  Number of blackbox outputs     -- \b IN.
         \param cf The reason for a check failure -- \b OUT.
         \return A boolean equal to \c true if the point is valid.
         */
        bool check ( int m , NOMAD::check_failed_type & cf ) const;
        
        /*---------------*/
        /*  GET methods  */
        /*---------------*/
        
        /// Size of the point in memory.
        /**
         \return Size of the point in memory, in bytes.
         */
        virtual int size_of ( void ) const;
        
        /// Access to the signature of the point.
        /**
         \return A pointer to the signature of the point.
         */
        NOMAD::Signature * get_signature ( void ) const;
        
        /// Access to the number of blackbox outputs.
        /**
         \return Number of blackbox outputs.
         */
        int get_m ( void ) const { return _bb_outputs.size(); }
        
        /// Access to the number of variables.
        /**
         \return Number of variables.
         */
        int get_n ( void ) const { return size(); }
        
        /// Check if the evaluation at this point is valid.
        /**
         \return A boolean equal to \c true if the evaluation is valid.
         */
        bool is_eval_ok ( void ) const { return (_eval_status == NOMAD::EVAL_OK); }
        
        /// Check if the evaluation at this point has been rejected by the user.
        /**
         \return A boolean equal to \c true if the evaluation has been rejected by the user.
         */
        bool check_rejected ( void ) const { return (_eval_status == NOMAD::EVAL_USER_REJECT); }
        
        
        /// Access to the tag of the point.
        /**
         \return The tag.
         */
        int get_tag ( void ) const { return _tag; }
        
        /// Access to the bbe of the point.
        /**
         \return The bbe.
         */
        int get_bbe ( void ) const { return _bbe; }
        
        /// Access to the bbe of the point.
        /**
         \return The sgte_bbe.
         */
        int get_sgte_bbe ( void ) const { return _sgte_bbe; }
        
        
        /// Access to the real time of the point evaluation.
        /**
         \return The real time.
         */
        int get_real_time( void ) const { return _real_time; }
        
        
        /// Access to the objective value \c f.
        /**
         \return The objective value.
         */
        const NOMAD::Double & get_f ( void ) const { return _f; }
        
        /// Access to the feasibility value \c h.
        /**
         \return The feasibility value.
         */
        const NOMAD::Double & get_h ( void ) const { return _h; }
        
        /// Access to the blackbox outputs.
        /**
         \return The \c m blackblack outputs.
         */
        const NOMAD::Point & get_bb_outputs ( void ) const { return _bb_outputs; }
        
        /// Access to the direction used to construct the point.
        /**
         \return The direction used to construct the point;
         may be \c NULL if no direction has been used.
         */
        const NOMAD::Direction * get_direction ( void ) const { return _direction; }
        
        /// Access to the poll_center used to construct the point.
        /**
         \return The poll center used to construct the point;
         may be \c NULL if no poll center has been used.
         */
        const NOMAD::Eval_Point * get_poll_center ( void ) const { return _poll_center ; }
        
        
        /// Check if the point has been generated during the current run.
        /**
         \return A boolean equal to \c true if the point has been
         generated during the current run.
         */
        bool get_current_run ( void ) const { return _current_run; }
        
        /// Check if the point is in cache.
        /**
         \return A boolean equal to \c true if the point is in cache.
         */
        bool is_in_cache ( void ) const { return _in_cache; }
        
        /// Check if the point respects the EB constraints.
        /**
         \return A boolean equal to \c true if the point respects the EB constraints.
         */
        bool is_EB_ok ( void ) const { return _EB_ok; }
        
        /// Access to the evaluation type.
        /**
         \return The evaluation type (true or surrogate).
         */
        NOMAD::eval_type get_eval_type ( void ) const { return _eval_type; }
        
        /// Access to the poll center type.
        /**
         \return The poll center type (feasible or not).
         */
        NOMAD::poll_center_type get_poll_center_type ( void ) const
        {
            return _poll_center_type;
        }
        
        /// Access to the evaluation status.
        /**
         \return The evaluation status
         (evaluation failed, succeeded, or is in progress).
         */
        NOMAD::eval_status_type get_eval_status ( void ) const { return _eval_status; }
        
        /// Access to the user evaluation priority.
        /**
         \return The user evaluation priority.
         */
        const NOMAD::Double & get_user_eval_priority ( void ) const
        {
            return _user_eval_priority;
        }
        
        /// Access to the random evaluation priority.
        /**
         \return The random evaluation priority.
         */
        const NOMAD::Double & get_rand_eval_priority ( void ) const
        {
            return _rand_eval_priority;
        }
        
        /// Check the point feasibility.
        /**
         The point is feasible if \c h \c <= \c h_min.
         \param h_min Feasibility threshold -- \b IN.
         \return A boolean equal to \c true if the point is feasible.
         */
        bool is_feasible ( const NOMAD::Double & h_min ) const
        {
            return ( _h.is_defined() && _h <= h_min );
        }
        
        /// Scaling.
        void scale ( void );
        
        /// Unscaling.
        void unscale ( void );
        
        /// Snap to bounds.
        /**
         \return A boolean equal to \c true if the snapping went well.
         */
        bool snap_to_bounds ( void );
        
        /// Treat the periodic variables.
        /**
         \param  new_dir A pointer to the modified direction used to generate
         this point from the poll center; may be \c NULL if no
         direction has been used -- \b OUT.
         \return A boolean equal to \c true if the treatment went well.
         */
        bool treat_periodic_variables ( NOMAD::Direction *& new_dir );
        
        /// Comparison operator.
        /**
         \param x Right-hand side object -- \b IN.
         \return A boolean equal to \c true if \c *this \c < \c x .
         */
        bool operator < ( const Eval_Point & x ) const;
        
        /*---------------*/
        /*  SET methods  */
        /*---------------*/
        
        /// Set the \c n and \c m.
        /**
         \param n Number of variables        -- \b IN.
         \param m Number of blackbox outputs -- \b IN.
         */
        void set ( int n , int m );
        
        /// Set the coordinates and \c m.
        /**
         \param x Coordinates of the point   -- \b IN.
         \param m Number of blackbox outputs -- \b IN.
         */
        void set ( const NOMAD::Point & x , int m );
        
        /// Set the tag.
        /**
         \param tag The tag -- \b IN.
         */
        void set_tag ( int tag );
        
        /// Set the real time.
        /**
         \param real_time The real time -- \b IN.
         */
        void set_real_time ( int real_time ){_real_time=real_time;}
        
        /// increment sgte_bbe.
        /**
         */
        void increment_sgte_bbe (void);
        
        /// increment bbe.
        /**
         */
        void increment_bbe (void);
        
        /// Set the blackbox evaluation number.
        /**
         \param bbe the evaluation number
         */
        void set_bbe (int bbe){_bbe=bbe;}
        
        /// Set the objective value \c f.
        /**
         \param f Objective value -- \b IN.
         */
        void set_f ( const NOMAD::Double & f ) { _f = f; }
        
        /// Set the feasibility value \c h.
        /**
         \param h Feasibility value -- \b IN.
         */
        void set_h ( const NOMAD::Double & h ) { _h = h; }
        
        /// Set the user evaluation priority.
        /**
         \param u User evaluation priority -- \b IN.
         */
        void set_user_eval_priority ( const NOMAD::Double & u ) { _user_eval_priority = u; }
        
        /// Set the random evaluation priority.
        /**
         \param r Random evaluation priority -- \b IN.
         */
        void set_rand_eval_priority ( const NOMAD::Double & r ) { _rand_eval_priority = r; }
        
        /// Set one blackbox output.
        /**
         \param i Index of the output to set -- \b IN.
         \param v Value of the output        -- \b IN.
         */
        void set_bb_output ( int i , const NOMAD::Double & v ) { _bb_outputs[i]= v; }
        
        /// Set all blackbox outputs.
        /**
         \param b The \c m blackbox outputs -- \b IN.
         */
        void set_bb_output ( const NOMAD::Point & b ) { _bb_outputs = b; }
        
        /// Set the evaluation status.
        /**
         \param e Evaluation status (failed, succeeded, or in progress)
         -- \b IN.
         */
        void set_eval_status ( const NOMAD::eval_status_type & e ) { _eval_status = e; }
        
        /// Set if the point respects the EB constraints.
        /**
         \param e A boolean equal to \c true if the point
         respects the EB constraints.
         */
        void set_EB_ok ( bool e ) { _EB_ok = e; }
        
        /// Set the evaluation type.
        /**
         \param e Evaluation type (true or surrogate) -- \b IN.
         */
        void set_eval_type ( NOMAD::eval_type e ) { _eval_type = e; }
        
        /// Set the type of the poll center.
        /**
         \param p Type of the poll center (feasible or not) -- \b IN.
         */
        void set_poll_center_type ( NOMAD::poll_center_type p ) { _poll_center_type = p; }
        
        /// Indicate if the point has been generated during the current run.
        /**
         \param c A boolean equal to \c true if the point
         has been generated during the current run -- \b IN.
         */
        void set_current_run ( bool c ) const { _current_run = c; }
        
        /// Indicate if the point is in cache.
        /**
         \param i A boolean equal to \c true if the point is in cache
         -- \b IN.
         */
        void set_in_cache ( bool i ) const { _in_cache = i; }
        
        
        /// Set the direction used to create the point.
        /**
         \param d A pointer to the direction; may be \c NULL -- \b IN.
         */
        void set_direction ( const NOMAD::Direction * d );
        
        
        /// Set the poll center used to create the point.
        /**
         \param d A pointer to the poll center; may be \c NULL -- \b IN.
         */
        void set_poll_center ( const NOMAD::Eval_Point * d );
        
        /// Set the signature.
        /**
         \param s A pointer to the signature -- \b IN.
         */
        void set_signature ( NOMAD::Signature * s );
        
        /// Check if there are nan's in the blackbox outputs:
        /**
         \return \c true if there is at least a nan.
         */
        bool check_nan ( void ) const;
        
        
        /// Reset the tags and bbes
        /**
         */		
        static void reset_tags_and_bbes ( void ) {_current_tag = 0;_current_bbe = 0;_current_sgte_bbe = 0;}
        
#ifdef MODEL_STATS
        
    private:
        
        mutable int           _mod_use; // 1: model search; 2: model eval sort
        mutable int           _nY;
        mutable NOMAD::Double _cond , _Yw , _mh , _mf;
        
    public:
        
        void set_mod_use ( int                 mod_use ) const { _mod_use = mod_use; }
        void set_nY       ( int                     nY ) const { _nY      = nY;      }
        void set_cond     ( const NOMAD::Double & cond ) const { _cond    = cond;    }
        void set_Yw       ( const NOMAD::Double &   Yw ) const { _Yw      = Yw;      }
        void set_mh       ( const NOMAD::Double &   mh ) const { _mh      = mh;      }
        void set_mf       ( const NOMAD::Double &   mf ) const { _mf      = mf;      }
        
        int                   get_mod_use ( void ) const    { return _mod_use; }
        int                   get_nY      ( void ) const    { return _nY;      }
        const NOMAD::Double & get_cond    ( void ) const    { return _cond;    }
        const NOMAD::Double & get_Yw      ( void ) const    { return _Yw;      }
        const NOMAD::Double & get_mh      ( void ) const    { return _mh;      }
        const NOMAD::Double & get_mf      ( void ) const    { return _mf;      }
        
        void set_model_data   ( const NOMAD::Eval_Point & x ) const;
        void clear_model_data ( void                        ) const;
        
#endif
        
        /// Display the tag of the point.
        /**
         \param out The NOMAD::Display object -- \b IN.
         */
        void display_tag ( const NOMAD::Display & out ) const;
        
        /// Display.
        /**
         \param out      The NOMAD::Display object -- \b IN.
         \param in_block If \c true, the point is displayed into an indented block
         -- \b IN -- \b optional (default = \c true ).
         */
        virtual void display_eval ( const NOMAD::Display & out , bool in_block = true ) const;
    };
    
    /// Display a NOMAD::Eval_Point object.
    /**
     \param out The NOMAD::Display object -- \b IN.
     \param x   The NOMAD::Eval_Point object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
     */
    inline const NOMAD::Display & operator << ( const NOMAD::Display    & out ,
                                               const NOMAD::Eval_Point & x     )
    {
        x.display_eval (out , true );
        return out;
    }
}

#endif
