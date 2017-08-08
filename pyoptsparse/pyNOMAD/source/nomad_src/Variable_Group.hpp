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
 \file   Variable_Group.hpp
 \brief  Group of variables (headers)
 \author Sebastien Le Digabel
 \date   2010-04-XX
 \see    Variable_Group.cpp
 */
#ifndef __VARIABLE_GROUP__
#define __VARIABLE_GROUP__

#include "Directions.hpp"

namespace NOMAD {
    
    /// Group of variables.
    /**
     A group can be composed of integer / continous / binary variables
     or only categorical variables.
     */
    class Variable_Group {
        
    private:
        
        std::set<int>          _var_indexes; ///< The variable indexes.
        NOMAD::Directions    * _directions;  ///< The directions.
        const NOMAD::Display & _out;         ///< Display.
        
        /// Affectation operator.
        /**
         \param vg The right-hand side object -- \b IN.
         */
        Variable_Group & operator = ( const Variable_Group & vg );
        
    public:
        
        /// Constructor.
        /**
         \param var_indexes        Variable indexes                   -- \b IN.
         \param direction_types    Direction types                    -- \b IN.
         \param sec_poll_dir_types Direction types for secondary poll -- \b IN.
         \param out                Display                            -- \b IN.
         */
        Variable_Group ( const std::set<int>                  & var_indexes        ,
                        const std::set<NOMAD::direction_type> & direction_types    ,
                        const std::set<NOMAD::direction_type> & sec_poll_dir_types ,
                        const NOMAD::Display                  & out                  )
        : _var_indexes ( var_indexes ),
        _directions  ( new Directions ( static_cast<int>(var_indexes.size()) ,
                                       direction_types                      ,
                                       sec_poll_dir_types                   ,
                                       out                                    ) ) ,
        _out         ( out )   {}
        
        /// Copy constructor.
        /**
         \param vg The copied object -- \b IN.
         */
        Variable_Group ( const Variable_Group & vg )
        : _var_indexes ( vg._var_indexes                    ) ,
        _directions  ( new Directions ( *vg._directions ) ) ,
        _out         ( vg._out                            )   {}
        
        /// Destructor.
        virtual ~Variable_Group ( void ) { delete _directions; }
        
        /// Comparison operator.
        /**
         \param vg The right-hand side object -- \b IN.
         \return A boolean equal to \c true if \c *this \c < \c vg.
         */
        bool operator < ( const Variable_Group & vg ) const;
        
        /// Check the group.
        /**
         Also removes fixed variables from the group.
         \param fixed_vars Fixed variables -- \b IN.
         \param bbit       Input types     -- \b IN.
         \param in_group   A pointer to a vector indicating if variables
         belong to a group (may be \c NULL) -- \b OUT.
         \param mod        True if the variable group has been modified    -- \b out
         \return A boolean equal to \c true if the verification is valid.
         */
        bool check ( const NOMAD::Point                      & fixed_vars ,
                    const std::vector<NOMAD::bb_input_type> & bbit       ,
                    std::vector<bool>                       * in_group   ,
                    bool 									 & mod            );
        
        /// Check if the directions are Ortho-MADS directions.
        /**
         \return A boolean equal to \c true if the directions are Ortho-MADS directions.
         */
        bool is_orthomads ( void ) const { return _directions->is_orthomads(); }
        
        
        /// Access to the variable indexes.
        /**
         \return The variable indexes.
         */
        const std::set<int> & get_var_indexes ( void ) const { return _var_indexes; }
        
        /// Access to the direction types.
        /**
         \return The direction types.
         */
        const std::set<NOMAD::direction_type> & get_direction_types ( void ) const
        {
            return _directions->get_direction_types();
        }
        
        /// Access to the direction types for the secondary poll.
        /**
         \return The direction types for the secondary poll.
         */
        const std::set<NOMAD::direction_type> & get_sec_poll_dir_types ( void ) const
        {
            return _directions->get_sec_poll_dir_types();
        }
        
        /// Access to the directions.
        /**
         \param dirs				List of directions                      -- \b OUT.
         \param poll				Type of poll (primary or secondary)     -- \b IN.
         \param mesh				Mesh									-- \b IN.
         */
        void get_directions ( std::list<NOMAD::Direction>	& dirs			,
                             NOMAD::poll_type				poll			,
                             const NOMAD::OrthogonalMesh	& mesh			);
        
        
        /// Access to one direction for a given mesh.
        /**
         Used for example in the VNS search.
         \param dir          The direction             -- \b OUT.
         */
        bool get_one_direction ( NOMAD::Direction & dir  )
        {
            return _directions->compute_dir_on_unit_sphere( dir );
        }
        
        
        /// Display.
        /**
         \param out The NOMAD::Display object -- \b IN.
         */
        void display ( const NOMAD::Display & out ) const;
        
        /// Display.
        /**
         Uses the \c this->_out member as NOMAD::Display object.
         */
        void display ( void ) const { display ( _out ); }
    };
    
    /// Display a NOMAD::Variable_Group object.
    /**
     \param out The NOMAD::Display object                        -- \b IN.
     \param vg  The NOMAD::Variable_Group object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
     */
    inline const NOMAD::Display & operator << ( const NOMAD::Display        & out ,
                                               const NOMAD::Variable_Group & vg    ) {
        vg.display ( out );
        return out;
    }
    
    /// Allow the comparison between two pointers to NOMAD::Variable_Group objects.
    /**
     Used with \c set<Variable_Group*,VG_Comp> in class NOMAD::Parameters.
     */
    struct VG_Comp {
        /// Comparison between two NOMAD::Variable_Group objects.
        bool operator() ( const Variable_Group * vg1 , const Variable_Group * vg2 ) const {
            return (*vg1 < *vg2);
        }
    };
}

#endif
