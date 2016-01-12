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
  \file   OrthogonalMesh.hpp
  \brief  Virtual class for the MADS orthogonal meshes (headers)
  \author Christophe Tribes
  \date   2010-04-06
  \see    Mesh.cpp XMesh.cpp
*/
#ifndef __ORTHOGONALMESH__
#define __ORTHOGONALMESH__

#include "Point.hpp"
#include "Direction.hpp"

namespace NOMAD {
	
	/// Virtual class for the MADS orthogonal meshes.
	/**
     - An orthogonal mesh in NOMAD is defined with the basic orthogonal directions the
	 mesh size parameter delta^k.
     - The poll size parameter Delta^k is not used to define the mesh but
	 to define the poll trial points. 
     - At each MADS iteration the mesh is updated.
     - Mesh and poll size parameters are stored as NOMAD::Point objects.
	 */
	class OrthogonalMesh {
		
		/*--------------------------------------------------------------*/
	private:

			
		/*--------------------------------------------------------------*/
		
		/// Private affectation operator.
		/**
		 \param m The right-hand side object -- \b IN.
		 */
		const OrthogonalMesh & operator = ( const OrthogonalMesh & m );
		
		
		/*--------------------------------------------------------------*/
	protected:
		
		
		NOMAD::Point	_delta_0;
		NOMAD::Point	_Delta_0;
		NOMAD::Point	_Delta_min;
		NOMAD::Point	_delta_min;
		
		NOMAD::Double	_update_basis;
		int				_coarsening_step;
		int				_refining_step  ;
		
		int				_n;
        int             _n_free_variables;
        
        int             _limit_mesh_index;   // Limit max or min of the mesh index for fine mesh (SMesh->max, XMesh->min)
		
		/// Constructor (called only by derived objects).
		/**
		 \param Delta_0						Initial poll size Delta_0						-- \b IN.
		 \param Delta_min					Minimal poll size Delta_min (may be undefined)	-- \b IN.
		 \param delta_min					Minimal mesh size delta_min (may be undefined)	-- \b IN.
		 \param fixed_variables				Fixed variables                                 -- \b IN.
		 \param update_basis				Mesh/poll update basis       (tau)				-- \b IN.
		 \param coarsening_step				Mesh/poll coarsening exponent (w+)				-- \b IN.
		 \param refining_step				Mesh/poll refining exponent   (w-)				-- \b IN.
         \param limit_mesh_index            Limit mesh index to trigger stopping criterion  -- \b IN.
		 */
		OrthogonalMesh (const NOMAD::Point	& Delta_0   ,
						const NOMAD::Point	& Delta_min ,
						const NOMAD::Point	& delta_min ,
                        const NOMAD::Point  & fixed_variables ,
						NOMAD::Double		update_basis,
						int					coarsening_step,
						int					refining_step,
                        int                 limit_mesh_index ) ;
		
		
		/// Copy constructor (called only by derived objects).
		/**
		 \param m The copied object -- \b IN.
		 */
		OrthogonalMesh ( const OrthogonalMesh & m )
		:	_delta_0			( m._delta_0			),
		_Delta_0				( m._Delta_0			),
		_Delta_min				( m._Delta_min			),
		_delta_min				( m._delta_min			),
		_update_basis			( m._update_basis		),
		_coarsening_step		( m._coarsening_step	),
		_refining_step			( m._refining_step		),
		_n						( m._n					),
        _n_free_variables       ( m._n_free_variables   ),
        _limit_mesh_index       ( m._limit_mesh_index){}
		
		
		/*--------------------------------------------------------------*/	  
	public:
		
        
		/// Destructor.
		virtual ~OrthogonalMesh ( void ){;}
		
		
		/// Update the Mesh (poll and mesh sizes).
		/**
		 \param success    Type of success of the iteration			-- \b IN.
		 \param dir        Direction of the iteration (optional)	-- \b IN.
		 */
		virtual void update ( NOMAD::success_type success, const NOMAD::Direction * dir=NULL) = 0;

		
		/// Update the provided mesh indices (the Mesh is unchanged).
		/**
		 \param success			Type of success of the iteration			-- \b IN.
		 \param mesh_indices	Provided mesh indices for update			-- \b IN/OUT.
		 \param dir				Direction of the iteration (optional)		-- \b IN.
		 */
		virtual void update ( NOMAD::success_type success, NOMAD::Point & mesh_indices, const NOMAD::Direction * dir=NULL ) const = 0;
		
		
		
		/// Reset the Mesh to its original sizes (poll and mesh sizes).
		virtual void reset ( void ) = 0;
		
		
		/// Access to the initial mesh size.
		/**
		 \return A NOMAD::Point for the initial mesh size.
		 */
		const NOMAD::Point & get_initial_mesh_size ( void ) const { return _delta_0; }
        
        /// Access to the initial poll size.
        /**
         \return A NOMAD::Point for the initial poll size.
         */
        const NOMAD::Point & get_initial_poll_size ( void ) const { return _Delta_0; }

				
		/// Access to the minimal mesh size.
		/**
		 \return A NOMAD::Point for the minimal mesh size.
		 */
		const NOMAD::Point & get_min_mesh_size ( void ) const { return _delta_min; }
		
        
		/// Access to the minimal poll size.
		/**
		 \return A NOMAD::Point for the minimal poll size.
		 */
		const NOMAD::Point & get_min_poll_size ( void ) const { return _Delta_min; }
		
		
		/// Test if mesh is finest so far.
		/**
		 \return True if mesh is the finest so far, False otherwise.
		 */	
		virtual bool is_finest(void) const = 0;
        

        /// Test if current mesh is finer than initial mesh (used by VNS search).
		/**
		 \return True if mesh size is smaller than initial mesh size for all components.
		 */
		bool is_finer_than_initial (void) const;

        
		/// Access to the mesh/poll update basis tau.
		/**
		 \return A double with the update basis tau.
		 */
		double get_update_basis ( void ) const { return _update_basis.value(); }
        
		
		/// Access to the mesh ratio after a success
		/**
		 \return A point with the ratio for each coordinate
		 */
		virtual NOMAD::Point get_mesh_ratio_if_success( void ) const = 0;
		
		
		/// Access to the number of variables.
		/**
		 \return An integer with the number of variables.
		 */
		int get_n ( void ) const { return _n; }
        

        /// Access to the number of free variables.
		/**
		 \return An integer with the number of free variables.
		 */
		int get_n_free_variables ( void ) const { return _n_free_variables; }

        
        /// Access to the mesh size parameter delta^k.
		/**
		 \return delta    The mesh size parameter delta^k -- \b OUT.
		 */
		NOMAD::Point get_delta ( void ) const
        {
            NOMAD::Point delta;
            get_delta(delta);
            return delta;
        }

		/// Access to the mesh size parameter delta^k.
		/**
		 \param delta    The mesh size parameter delta^k -- \b OUT.
		 \return A boolean equal to \c true if all values are
		 strictly inferior than the associated minimal
		 mesh size delta_min
		 (stopping criterion MIN_MESH_SIZE).
		 */
		virtual bool get_delta ( NOMAD::Point & delta ) const = 0;
		 
		/// Access to the largest mesh size.
		/**
		 \return  The largest mesh size  -- \b OUT.
		 */
		virtual NOMAD::Point get_delta_max ( void ) const = 0;
		
        
        /// Access to the poll size parameter Delta^k.
		/**
		 \return Delta    The poll size parameter Delta^k -- \b OUT.
		 */
		NOMAD::Point get_Delta ( void )
        {
            NOMAD::Point Delta;
            get_Delta(Delta);
            return Delta;
        }
        
		/// Access to the poll size parameter Delta^k.
		/**
		 \param Delta    The poll size parameter Delta^k -- \b OUT.
		 \return A boolean equal to \c true if all values are
		 strictly inferior than the associated minimal
		 mesh size delta_min
		 (stopping criterion MIN_POLL_SIZE).
		 */
		virtual bool get_Delta ( NOMAD::Point & Delta ) const = 0 ;
		
		
		/// Display.
		/**
		 \param out The NOMAD::Display object -- \b IN.
		 */
		virtual void display ( const NOMAD::Display & out ) const = 0;
		
		/// Check the stopping conditions on the minimal poll and mesh sizes.
		/**
		 \param stop           Stop flag                  -- \b IN/OUT.
		 \param stop_reason    Stop reason                -- \b OUT.
		 */
		virtual void check_min_mesh_sizes (	bool             & stop           ,
										   NOMAD::stop_type & stop_reason      ) const = 0;
		
		/// Access to the mesh indices per coordinate.
		/**
		 \return A point with the mesh index for each coordinate.
		 */
		virtual const NOMAD::Point get_mesh_indices ( void  ) const = 0;
		
		/// Manually set the mesh indices per coordinate (virtual).
		/**
		 \param r   The mesh index per coordinate -- \b IN.
		 */
		virtual void set_mesh_indices ( const NOMAD::Point & r ) =0 ;
        
 		/// Manually set the min mesh size per coordinate.
		/**
		 \param delta_min   The min mesh sizes per coordinate (can be undefined) -- \b IN.
		 */
		void set_min_mesh_sizes ( const NOMAD::Point & delta_min );
        
 		/// Manually set intial mesh size per coordinate.
		/**
		 \param d   The initial mesh sizes per coordinate -- \b IN.
		 */
        void set_delta_0 ( const NOMAD::Point & d );
        
  		/// Manually set intial poll size per coordinate.
		/**
		 \param d   The initial poll sizes per coordinate -- \b IN.
		 */
        void set_Delta_0 ( const NOMAD::Point & d );

		
		/// Access to the min mesh indices reached so far.
		/**
		 \return A point with the mesh index for each coordinate.
		 */
		virtual const NOMAD::Point get_min_mesh_indices ( void  ) const = 0;	
		
		/// Access to the max mesh indices reached so far.
		/**
		 \return A point with the mesh index for each coordinate.
		 */
		virtual const NOMAD::Point get_max_mesh_indices ( void  ) const = 0;	
		

        /// Access to the limit mesh index.
		/**
		 \return An integer with the limit mesh index.
		 */
		int get_limit_mesh_index ( void  ) const { return _limit_mesh_index;}
        
        
        /// Manually set the limit mesh index.
		/**
		 \param limit_mesh_index   The limit mesh index.
		 */
		virtual void set_limit_mesh_index ( int limit_mesh_index  ) = 0;

        
		/// Scale and project the ith component of a vector on the mesh
		/**
		 \param i	The vector component number			-- \b IN.
		 \param l	The vector component value			-- \b IN.
		 \return	The ith component of a vector after mesh scaling and projection
		 */
		virtual NOMAD::Double scale_and_project(int i, NOMAD::Double l) const = 0 ;
		
		
	};
	
	/// Display a NOMAD::OrthogonalMesh object.
	/**
     \param out The NOMAD::Display object -- \b IN.
     \param m   The NOMAD::OrthogonalMesh object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
	 */
	inline const NOMAD::Display & operator << ( const NOMAD::Display & out ,
											   const NOMAD::OrthogonalMesh    & m     )
	{
		m.display ( out );
		return out;
	}
}

#endif
