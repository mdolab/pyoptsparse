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
  \file   nomad.hpp
  \brief  NOMAD header file
  \author Sebastien Le Digabel
  \date   2010-04-12
*/
#ifndef __NOMAD__
#define __NOMAD__

#include "Mads.hpp"


namespace NOMAD {

  /// Display NOMAD information.
  /**
     \param out A NOMAD::Display object -- \b IN.
  */
  void display_info ( const NOMAD::Display & out );

  /// Display NOMAD version.
  /**
     \param out A NOMAD::Display object -- \b IN.
  */
  void display_version ( const NOMAD::Display & out );

  /// Display NOMAD usage.
  /**
     \param exeName Name of executable  -- \b IN.
     \param out A NOMAD::Display object -- \b IN.
  */
  void display_usage ( char* exeName, const NOMAD::Display & out );
  void display_usage ( const NOMAD::Display & out );	

#ifdef MEMORY_DEBUG
  /// Display NOMAD most important structures in memory.
  /**
     Is defined only in debug mode with flag MEMORY_DEBUG active.
     \param out A NOMAD::Display object -- \b IN.
  */
  void display_cardinalities ( const NOMAD::Display & out );
#endif
}

#endif
