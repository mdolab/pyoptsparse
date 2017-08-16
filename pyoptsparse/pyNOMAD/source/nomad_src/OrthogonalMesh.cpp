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
 \brief  implementation
 \author Christophe Tribes
 \date   2014-06-19
 \see    SMesh.cpp XMesh.cpp
 */

#include "OrthogonalMesh.hpp"

/// Constructor (called only by derived objects).
NOMAD::OrthogonalMesh::OrthogonalMesh (const NOMAD::Point	& Delta_0   ,
                                       const NOMAD::Point	& Delta_min ,
                                       const NOMAD::Point	& delta_min ,
                                       const NOMAD::Point  & fixed_variables ,
                                       NOMAD::Double		update_basis,
                                       int					coarsening_step,
                                       int					refining_step,
                                       int                  limit_mesh_index ) :
_delta_0			( Delta_0 ),
_Delta_0			( Delta_0 ),
_Delta_min			( Delta_min ),
_delta_min			( delta_min ),
_update_basis		( update_basis ),
_coarsening_step	( coarsening_step ),
_refining_step		( refining_step ),
_limit_mesh_index   ( limit_mesh_index )
{
    bool chkMesh  = delta_min.is_defined();
    bool chkPoll  = Delta_min.is_defined();
    _n = Delta_0.size();
    
    _n_free_variables = _n - fixed_variables.nb_defined();
    
    // The delta_0 are obtained from the Delta_0
    _delta_0*=pow(_n_free_variables,-0.5);
    
    if ( !_Delta_0.is_complete() )
        throw NOMAD::Exception (  "OrthogonalMesh.hpp" , __LINE__ ,
                                "NOMAD::OrthogonalMesh::OrthogonalMesh(): delta_0 has undefined values" );
    
    if ( chkMesh && delta_min.size() != _n )
        throw NOMAD::Exception ( "OrthogonalMesh.hpp" , __LINE__ ,
                                "NOMAD::OrthogonalMesh::OrthogonalMesh(): delta_0 and delta_min have different sizes" );
    
    if ( chkPoll && Delta_min.size() != _n )
        throw NOMAD::Exception ( "OrthogonalMesh.hpp" , __LINE__ ,
                                "NOMAD::OrthogonalMesh::OrthogonalMesh(): Delta_0 and Delta_min have different sizes" );
    
    
    std::string error;
    for ( int k = 0 ; k < _n ; ++k )
    {
        // we check that Delta_min <= Delta_0 and that delta_min <= delta_0:
        if ( chkMesh &&
            _delta_min[k].is_defined()						&&
            _delta_0[k] < _delta_min[k]		)
        {
            error = "NOMAD::OrthogonalMesh::OrthogonalMesh(): delta_0 < delta_min";
            break;
        }
        if ( chkPoll &&
            _Delta_min[k].is_defined()						&&
            _Delta_0[k] < _Delta_min[k]     )
        {
            error = "NOMAD::OrthogonalMesh::OrthogonalMesh(): Delta_0 < Delta_min";
            break;
        }
        
    }
    
    if ( !error.empty() )
        throw NOMAD::Exception ( "OrthogonalMesh.hpp" , __LINE__ , error );
}



bool NOMAD::OrthogonalMesh::is_finer_than_initial (void) const
{
    NOMAD::Point delta;
    get_delta(delta);
    
    for (int i =0; i < _n ; ++i )
        if ( delta[i] >= _delta_0[i] )
            return false;
    
    return true;
}



/// Manually set the min mesh size per coordinate.
void NOMAD::OrthogonalMesh::set_min_mesh_sizes ( const NOMAD::Point & delta_min )
{
    
    // If delta_min undefined than _delta_min->undefined
    if ( ! delta_min.is_defined() )
    {
        _delta_min.clear();
        return;
    }
    
    // If delta_min defined test that everything is consistent
    if ( delta_min.size() != _n )
        throw NOMAD::Exception ( "XMesh.cpp" , __LINE__ ,
                                "NOMAD::OrthogonalMesh::set_delta_min() delta_min has dimension different than mesh dimension" );
    
    if ( !delta_min.is_complete() )
        throw NOMAD::Exception (  "OrthogonalMesh.hpp" , __LINE__ ,
                                "NOMAD::OrthogonalMesh::set_delta_min(): delta_min has some defined and undefined values" );
    
    std::string error;
    for ( int k = 0 ; k < _n ; ++k )
    {
        
        // we check that Delta_min <= Delta_0 and that delta_min <= delta_0:
        if ( delta_min[k].is_defined()						&&
            _delta_0[k] < delta_min[k]		)
        {
            error = "NOMAD::OrthogonalMesh::set_delta_min(): delta_0 < delta_min";
            break;
        }
        if ( delta_min[k].is_defined()						&&
            _Delta_0[k] < delta_min[k]     )
        {
            error = "NOMAD::OrthogonalMesh::set_delta_min(): Delta_0 < delta_min";
            break;
        }
        
    }
    
    if ( !error.empty() )
        throw NOMAD::Exception ( "OrthogonalMesh.hpp" , __LINE__ , error );
    
    _delta_min=delta_min;
    
}



/*-----------------------------------------------------------*/
/*             set delta_0                                   */
/*-----------------------------------------------------------*/
void NOMAD::OrthogonalMesh::set_delta_0 ( const NOMAD::Point & d )
{
    
    if ( d.size() != _delta_0.size() )
        throw NOMAD::Exception ( "OrthogonalMesh.hpp" , __LINE__ ,
                                "NOMAD::OrthogonalMesh::set_delta_0(): dimension of provided delta_0 must be consistent with their previous dimension" );
    
    _delta_0=d;
}

/*-----------------------------------------------------------*/
/*             set Delta_0                                   */
/*-----------------------------------------------------------*/
void NOMAD::OrthogonalMesh::set_Delta_0 ( const NOMAD::Point & d )
{
    
    if ( d.size() != _Delta_0.size() )
        throw NOMAD::Exception ( "XMesh.cpp" , __LINE__ ,
                                "NOMAD::XMesh::set_Delta_0(): dimension of provided Delta_0 must be consistent with their previous dimension" );
    
    _Delta_0=d;
}
