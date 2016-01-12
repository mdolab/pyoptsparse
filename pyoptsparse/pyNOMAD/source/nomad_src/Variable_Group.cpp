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
 \file   Variable_Group.cpp
 \brief  Group of variables (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-22
 \see    Variable_Group.hpp
 */
#include "Variable_Group.hpp"

/*-------------------------------------------------------------*/
/*     check (also removes fixed variables from the group)     */
/*-------------------------------------------------------------*/
bool NOMAD::Variable_Group::check
( const NOMAD::Point                      & fixed_vars ,
 const std::vector<NOMAD::bb_input_type> & bbit       ,
 std::vector<bool>                       * in_group   ,
 bool									 & mod           )
{
    
    // initial checks:
    if ( _var_indexes.empty() )
        return false;
    
    bool binary      = true;
    bool categorical = false;
    bool reset_dirs  = false;
    
    // other checks + filling of vector in_group:
    int                           n   = static_cast<int>(bbit.size());
    std::set<int>::const_iterator end = _var_indexes.end() ,
    it  = _var_indexes.begin();
    
    while ( it != end )
    {
        
        if ( *it < 0 || *it >= n )
            return false;
        
        // check if the variable is fixed:
        if ( fixed_vars[*it].is_defined() )
        {
            reset_dirs = true;
            _var_indexes.erase ( it++ );
            mod=true;
            continue;
        }
        
        if ( bbit[*it] == NOMAD::CATEGORICAL )
        {
            categorical = true;
            binary      = false;
        }
        else
        {
            if ( categorical )
                return false;
            if ( bbit[*it] != NOMAD::BINARY )
                binary = false;
        }
        
        if ( in_group )
            (*in_group)[*it] = true;
        
        ++it;
    }
    
    // modify the directions if a fixed variable has been removed:
    if ( reset_dirs )
    {
        
        std::set<NOMAD::direction_type> direction_types = _directions->get_direction_types();
        std::set<NOMAD::direction_type> sec_poll_dir_types = _directions->get_sec_poll_dir_types();
        
        delete _directions;
        
        _directions = new Directions ( static_cast<int>(_var_indexes.size()) ,
                                      direction_types                       ,
                                      sec_poll_dir_types                    ,
                                      _out                                    );
    }
    
    if ( binary )
        _directions->set_binary();
    
    else
    {
        
        // we check here that NOMAD::GPS_BINARY is not in
        // dir_types nor in sec_poll_dir_types:
        const std::set<NOMAD::direction_type> & direction_types
        = _directions->get_direction_types();
        const std::set<NOMAD::direction_type> & sec_poll_dir_types
        = _directions->get_sec_poll_dir_types();
        
        if ( direction_types.find    ( NOMAD::GPS_BINARY ) != direction_types.end   () ||
            sec_poll_dir_types.find ( NOMAD::GPS_BINARY ) != sec_poll_dir_types.end()    )
            return false;
        
        if ( categorical )
            _directions->set_categorical();
    }
    
    return true;
}

/*---------------------------------------------------------*/
/*                compute the directions                   */
/*---------------------------------------------------------*/
void NOMAD::Variable_Group::get_directions ( std::list<NOMAD::Direction>	& dirs		,
                                            NOMAD::poll_type				poll		,
                                            const NOMAD::OrthogonalMesh		& mesh		)
{
    
    _directions->compute ( dirs	,
                          poll	,
                          mesh );
    return;
}

/*---------------------------------------------------------*/
/*                     comparison operator                 */
/*---------------------------------------------------------*/
bool NOMAD::Variable_Group::operator < ( const NOMAD::Variable_Group & vg ) const
{
    // variable indexes:
    if ( _var_indexes.size() < vg._var_indexes.size() )
        return true;
    if ( _var_indexes.size() > vg._var_indexes.size() )
        return false;
    
    std::set<int>::const_iterator it1 ,
    it2 = vg._var_indexes.begin() ,
    end = _var_indexes.end();
    for ( it1 = _var_indexes.begin() ; it1 != end ; ++it1 , ++it2 ) {
        if ( *it1 < *it2 )
            return true;
        if ( *it1 > *it2 )
            return false;
    }
    
    // directions:
    return ( *_directions < *vg._directions );
}

/*---------------------------------------------------------*/
/*                          display                        */
/*---------------------------------------------------------*/
void NOMAD::Variable_Group::display ( const NOMAD::Display & out ) const
{
    out << "indexes: { ";
    std::set<int>::const_iterator end = _var_indexes.end();
    for ( std::set<int>::const_iterator it = _var_indexes.begin() ; it != end ; ++it )
        out << *it << " ";
    out << "}" << std::endl;
    if ( _directions->is_categorical() )
        out << "no directions (categorical variables)" << std::endl;
    else
        out << NOMAD::open_block ( "directions" )
        << *_directions
        << NOMAD::close_block();
}
