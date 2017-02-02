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
 \file   Signature.hpp
 \brief  Evaluation point signature (headers)
 \author Sebastien Le Digabel
 \date   2010-04-22
 \see    Signature.cpp
 */
#ifndef __SIGNATURE__
#define __SIGNATURE__

#include "Variable_Group.hpp"

namespace NOMAD {
    
    /// Evaluation point signature.
    class Signature {
        
    private:
        
#ifdef MEMORY_DEBUG
        static int _cardinality;     ///< Number of NOMAD::Signature objects in memory.
        static int _max_cardinality; ///< Max number of NOMAD::Signature objects in memory.
#endif
        
        NOMAD::Point _lb;              ///< Lower bounds.
        NOMAD::Point _ub;              ///< Upper bounds.
        NOMAD::Point _scaling;         ///< Scaling.
        NOMAD::Point _fixed_variables; ///< Fixed variables.
        
        std::vector<NOMAD::bb_input_type> _input_types; ///< Input types.
        
        bool _all_continuous;  ///< Flag equal to \c true if all variables are continuous.
        bool _has_categorical; ///< Flag equal to \c true if there are categorical variables.
        
        static bool _warning_has_been_displayed;  ///< Flag equal to \c true if the warning has already been displayed
        
        std::vector<bool> _periodic_variables; ///< Periodic variables.
        
        /// Groups of variables.
        /**
         Include the directions/
         */
        std::list<NOMAD::Variable_Group*> _var_groups;
        
        /// Mesh associated to this signature.
        NOMAD::OrthogonalMesh *_mesh; ///< Orthogonal mesh
        
        /**
         Flag equal to \c true if the signature is standard
         (i.e created in Parameters class).
         */
        bool _std;
        
        // Feasible successful direction.
        /**
         Mutable since it is not used in
         \c operator \c < and have to be often changed.
         */
        mutable NOMAD::Direction _feas_success_dir;
        
        // Infeasible successful direction.
        /**
         Mutable since it is not used in
         \c operator \c < and have to be often changed.
         */
        mutable NOMAD::Direction _infeas_success_dir;
        
        
        // Direction group index
        int _dir_group_index;
        
        // Display
        NOMAD::Display _out;
        
        /*---------------------------------------------------------------------------*/
        
        /// Initializations.
        /**
         \param n                  Number of variables    -- \b IN.
         \param input_types        Types of the variables -- \b IN.
         \param lb                 Lower bounds           -- \b IN.
         \param ub                 Upper bounds           -- \b IN.
         \param scaling            Scaling                -- \b IN.
         \param fixed_variables    Fixed variables        -- \b IN.
         \param periodic_variables Periodic variables     -- \b IN.
         \param var_groups         Groups of variables    -- \b IN.
         */
        void init
        ( int                                                     n                  ,
         const std::vector<NOMAD::bb_input_type>               & input_types        ,
         const NOMAD::Point                                    & lb                 ,
         const NOMAD::Point                                    & ub                 ,
         const NOMAD::Point                                    & scaling            ,
         const NOMAD::Point                                    & fixed_variables    ,
         const std::vector<bool>                               & periodic_variables ,
         std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp> & var_groups           );
        
        /// Reset groups of variables.
        void reset_var_groups ( void );
        
        /// Clear.
        /**
         Private method called by the destructor and after an exception is thrown.
         */
        void clear ( void );
        
        /// Affectation operator.
        /**
         \param s The right-hand side object.
         */
        const Signature & operator = ( const Signature & s );
        
    public:
        
        /*---------------------------------------------------------------------------*/
        
        /// Exception class for an invalid signature.
        class Signature_Error : public NOMAD::Exception {
        public:
            /// Constructor.
            Signature_Error ( const std::string & file ,
                             int                 line ,
                             Signature         & s    ,
                             const std::string & msg    )
            : NOMAD::Exception ( file , line , msg ) { s.clear(); }
        };
        
        /*---------------------------------------------------------------------------*/
        
        /// Constructor #1.
        /**
         Advanced version.
         \param n							Number of variables			-- \b IN.
         \param input_types					Types of the variables		-- \b IN.
         \param use_smesh                   Type of mesh (SMesh, XMesh)	-- \b IN.
         \param anisotropic_mesh			Anisotropy of mesh (XMesh)	-- \b IN.
         \param initial_poll_size			Initial poll size			-- \b IN.
         \param min_poll_size				Minimum poll size			-- \b IN.
         \param min_mesh_size				Minimim mesh size			-- \b IN.
         \param mesh_update_basis			Mesh update basis			-- \b IN.
         \param poll_update_basis			Poll update basis			-- \b IN.
         \param mesh_coarsening_exponent	Mesh coarsening exponent	-- \b IN.
         \param mesh_refining_exponent		Mesh refining exponent		-- \b IN.
         \param initial_mesh_index			Initial mesh index			-- \b IN.
         \param lb							Lower bounds				-- \b IN.
         \param ub							Upper bounds				-- \b IN.
         \param scaling						Scaling						-- \b IN.
         \param fixed_variables				Fixed variables				-- \b IN.
         \param periodic_variables			Periodic variables			-- \b IN.
         \param var_groups					Groups of variables			-- \b IN.
         \param out							Display                     -- \b IN.
         */
        Signature
        ( int													n                  ,
         const std::vector<NOMAD::bb_input_type>				& input_types        ,
         const NOMAD::Point										& lb                 ,
         const NOMAD::Point										& ub                 ,
         bool													use_smesh	,
         bool													anisotropic_mesh	,
         const NOMAD::Point									    & initial_poll_size	,
         const NOMAD::Point										& min_poll_size		,
         const NOMAD::Point										& min_mesh_size		,
         NOMAD::Double											& mesh_update_basis	,
         NOMAD::Double											& poll_update_basis	,
         int													& mesh_coarsening_exponent,
         int													& mesh_refining_exponent,
         int													initial_mesh_index	,
         const NOMAD::Point										& scaling           ,
         const NOMAD::Point										& fixed_variables   ,
         const std::vector<bool>								& periodic_variables,
         std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp>		& var_groups         ,
         const NOMAD::Display									& out=NOMAD::Display()
         );
        
        /// Constructor #2.
        /**
         Basic version that will automatically construct groups of variables.
         \param n                  Number of variables                    -- \b IN.
         \param input_types        Types of the variables                 -- \b IN.
         \param initial_poll_size  Initial poll size                      -- \b IN.
         \param lb                 Lower bounds                           -- \b IN.
         \param ub                 Upper bounds                           -- \b IN.
         \param direction_types    Types of directions                    -- \b IN.
         \param sec_poll_dir_types Types of directions for secondary poll -- \b IN.
         \param out                Display                                -- \b IN.
         */
        Signature ( int                                       n                  ,
                   const std::vector<NOMAD::bb_input_type> & input_types        ,
                   const NOMAD::Point					   & initial_poll_size	,
                   const NOMAD::Point                      & lb                 ,
                   const NOMAD::Point                      & ub                 ,
                   const std::set<NOMAD::direction_type>   & direction_types    ,
                   const std::set<NOMAD::direction_type>   & sec_poll_dir_types ,
                   const NOMAD::Display                    & out = NOMAD::Display()
                   );
        
        /// Copy constructor.
        /**
         \param s The copied object.
         */
        Signature ( const Signature & s );
        
        /// Destructor.
        virtual ~Signature ( void );
        
        /// Reset.
        /**
         \param n                  Number of variables    -- \b IN.
         \param input_types        Types of the variables -- \b IN.
         \param lb                 Lower bounds           -- \b IN.
         \param ub                 Upper bounds           -- \b IN.
         \param scaling            Scaling                -- \b IN.
         \param fixed_variables    Fixed variables        -- \b IN.
         \param periodic_variables Periodic variables     -- \b IN.
         \param var_groups         Groups of variables    -- \b IN.
         */
        void reset
        ( int                                                     n                  ,
         const std::vector<NOMAD::bb_input_type>               & input_types        ,
         const NOMAD::Point                                    & lb                 ,
         const NOMAD::Point                                    & ub                 ,
         const NOMAD::Point                                    & scaling            ,
         const NOMAD::Point                                    & fixed_variables    ,
         const std::vector<bool>                               & periodic_variables ,
         std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp>	   & var_groups );
        
        /// Define a signature to be standard.
        void set_std ( void ) { _std = true; }
        
        /// Set a feasible successful direction.
        /**
         \param dir The direction -- \b IN.
         */
        void set_feas_success_dir ( const NOMAD::Direction & dir );
        
        /// Set an infeasible successful direction.
        /**
         \param dir The direction -- \b IN.
         */
        void set_infeas_success_dir ( const NOMAD::Direction & dir );
        
        /// Reset the feasible successful direction.
        void reset_feas_success_dir ( void ) const { _feas_success_dir.clear(); }
        
        /// Reset the infeasible successful direction.
        void reset_infeas_success_dir ( void ) const { _infeas_success_dir.clear(); }
        
        /// Scaling.
        /**
         Performed before an evaluation.
         \param x The scaled point -- \b IN/OUT.
         */
        void scale ( NOMAD::Point & x );
        
        /// Unscaling.
        /**
         Performed after an evaluation.
         \param x The unscaled point -- \b IN/OUT.
         */
        void unscale ( NOMAD::Point & x );
        
        /// Snap to bounds.
        /**
         - Supposes that \c this->treat_periodic_variables() has already been invoked.
         - If periodic variables have been treated, then bounds are
         satisfied and there is no need to snap anymore.
         \param x         The point to snap -- \b IN/OUT.
         \param direction A pointer to the direction associated to this point
         (may be \c NULL)  -- \b IN/OUT.
         \return A Boolean equal to \c true if \c x has been modified.
         */
        bool snap_to_bounds ( NOMAD::Point & x , NOMAD::Direction * direction );
        
        /// Treat the periodic variables.
        /**
         \param x The point to treat -- \b IN/OUT.
         \param old_dir A pointer to the direction associated to this point
         (before the treatment; may be \c NULL) -- \b IN.
         \param new_dir A pointer to the direction associated to this point
         (after the treatment; may be \c NULL) -- \b OUT.
         \return A boolean equal to \c true if \c x has been modified.
         */
        bool treat_periodic_variables ( NOMAD::Point            & x       ,
                                       const NOMAD::Direction *  old_dir ,
                                       NOMAD::Direction       *& new_dir   );
        
        /// Access to the lower bounds.
        /**
         \return The lower bounds.
         */
        const NOMAD::Point & get_lb ( void ) const { return _lb; }
        
        /// Access to the upper bounds.
        /**
         \return The upper bounds.
         */
        const NOMAD::Point & get_ub ( void ) const { return _ub; }
        
        /// Access to the mesh.
        /**
         \return The orthogonal mesh.
         */
        NOMAD::OrthogonalMesh * get_mesh ( void ) const { return _mesh; }
        
        
        /// Update the Mesh (poll and mesh sizes).
        /**
         \param success    Type of success of the iteration			-- \b IN.
         \param dir        Direction of the iteration (optional)	-- \b IN.
         */
        void update_mesh ( NOMAD::success_type success, const NOMAD::Direction * dir=NULL )
        {
            _mesh->update ( success, dir );
        }
        
        /// Complete reset of the mesh.
        /**
         */
        void reset_mesh ( void ) { _mesh->reset(); }
        
        /// Reset of the mesh indices.
        /**
         /param mesh_indices
         */
        void reset_mesh ( NOMAD::Point & mesh_indices  ) { _mesh->set_mesh_indices( mesh_indices); }
        
        
        /// Access to the scaling.
        /**
         \return The scaling.
         */
        const NOMAD::Point & get_scaling ( void ) const { return _scaling; }
        
        /// Access to the fixed variables.
        /**
         \return The fixed variables.
         */
        const NOMAD::Point & get_fixed_variables ( void ) const { return _fixed_variables; }
        
        
        /// Access to the number of fixed variables.
        /**
         \return The number of fixed variables.
         */
        int get_nb_fixed_variables ( void ) const;
        
        
        
        /// Access to the feasible successful direction.
        /**
         \return The feasible successful direction
         (may be undefined).
         */
        const NOMAD::Direction & get_feas_success_dir ( void ) const
        {
            return _feas_success_dir;
        }
        
        /// Access to the infeasible successful direction.
        /**
         \return The infeasible successful direction
         (may be undefined).
         */
        const NOMAD::Direction & get_infeas_success_dir ( void ) const
        {
            return _infeas_success_dir;
        }
        
        /// Access to the periodic variables.
        /**
         \return The periodic variables.
         */
        const std::vector<bool> & get_periodic_variables ( void ) const
        {
            return _periodic_variables;
        }
        
        /// Check if all variables are continuous.
        /**
         \return A boolean equal to \c true if all variables are continuous.
         */
        bool all_continuous ( void ) const { return _all_continuous; }
        
        /// Check if there are categorical variables.
        /**
         \return A boolean equal to \c true if there are categorical variables.
         */
        bool has_categorical ( void ) const { return _has_categorical; }
        
        /// Access to the number of categorical variables.
        /**
         \return Integer equal to the number of categorical variables.
         */
        int get_n_categorical ( void ) const ;
        
        /// Access to the number of variables.
        /**
         \return The number of variables.
         */
        int get_n ( void ) const { return static_cast<int> ( _input_types.size() ); }
        
        /// Access to the input types.
        /**
         \return The input types.
         */
        const std::vector<NOMAD::bb_input_type> & get_input_types ( void ) const
        {
            return _input_types;
        }
        
        /// Access to the input types.
        /**
         \return The input types.
         */
        const std::vector<NOMAD::bb_input_type> & get_input_type  ( void ) const
        {
            return _input_types;
        }
        
        /// Access to the groups of variables.
        /**
         \return The groups of variables.
         */
        const std::list<NOMAD::Variable_Group*> & get_var_groups ( void ) const
        {
            return _var_groups;
        }
        
        
        /// Check the compatibility of a point.
        /**
         - Only the number of variables is checked.
         - Other criteria (fixed variables, binary variables,...)
         are checked by NOMAD::Eval_Point::check().
         \param x The point to check -- \b IN.
         \return A boolean equal to \c true if the point is compatible
         with the signature.
         */
        bool is_compatible ( const NOMAD::Point & x ) const;
        
        /// Access to the directions.
        /**
         - The computed directions already include Delta^k_m.
         \param dirs              List of directions                      -- \b OUT.
         \param poll              Type of poll (primary or secondary)     -- \b IN.
         \param poll_center       Poll center                             -- \b IN.
         */
        void get_directions ( std::list<NOMAD::Direction>	& dirs			,
                             NOMAD::poll_type				poll			,
                             const NOMAD::Point				& poll_center	);
        
        
        /// Access to one direction for a given mesh.
        /**
         Used for example in the VNS search.
         \param dir          The direction  -- \b OUT.
         \param mesh_index   Mesh index ell -- \b IN.
         */
        void get_one_direction ( NOMAD::Direction & dir				,
                                int                mesh_index		) const;
        
        
        /// Comparison operator \c < .
        /**
         Successful directions are not considered.
         \param s The right-hand side object -- \b IN.
         \return A boolean equal to \c true if \c *this \c < \c s.
         */
        bool operator <  ( const Signature & s ) const;
        
        /// Comparison operator \c != .
        /**
         \param s The right-hand side object -- \b IN.
         \return A boolean equal to \c true if \c *this \c != \c s.
         */
        bool operator != ( const Signature & s ) const
        {
            return ( (*this < s) || (s < *this) );
        }
        
#ifdef MEMORY_DEBUG
        
        /// Access to the number of NOMAD::Signature objects in memory.
        /**
         \return Number of NOMAD::Signature objects in memory.
         */
        static int get_cardinality ( void ) { return Signature::_cardinality; }
        
        /// Access to the max number of NOMAD::Signature objects in memory.
        /**
         \return Max number of NOMAD::Signature objects in memory.
         */
        static int get_max_cardinality ( void ) { return Signature::_max_cardinality; }
#endif
        
        /// Display.
        /**
         \param out The NOMAD::Display object -- \b IN.
         */
        void display ( const NOMAD::Display & out ) const;
        
    };
    
    /// Display a NOMAD::Signature object.
    /**
     \param out The NOMAD::Display object                   -- \b IN.
     \param s   The NOMAD::Signature object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
     */
    inline const NOMAD::Display & operator << ( const NOMAD::Display   & out ,
                                               const NOMAD::Signature & s     ) {
        s.display ( out );
        return out;
    }
}

#endif
