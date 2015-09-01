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
  \file   Signature_Element.hpp
  \brief  Signature inside a set (headers)
  \author Sebastien Le Digabel
  \date   2010-04-12
*/
#ifndef __SIGNATURE_ELEMENT__
#define __SIGNATURE_ELEMENT__

#include "Set_Element.hpp"
#include "Signature.hpp"

namespace NOMAD {

  /// Signature inside a set of signatures.
  class Signature_Element : public NOMAD::Set_Element<NOMAD::Signature> {

  private:

    /// Affectation operator.
    /**
       \param se The right-hand side object -- \b IN.
    */
    Signature_Element & operator = ( const Signature_Element & se );

  public:

    /// Constructor.
    /**
       \param s A pointer to a NOMAD::Signature -- \b IN.
    */
    explicit Signature_Element ( const NOMAD::Signature * s )
      : NOMAD::Set_Element<NOMAD::Signature> ( s ) {}

    /// Copy constructor.
    /**
       \param se The copied object -- \b IN.
    */
    Signature_Element ( const Signature_Element & se )
      : NOMAD::Set_Element<NOMAD::Signature> ( se.get_element() ) {}
    
    /// Destructor.
    virtual ~Signature_Element ( void ) {}
    
    /// Comparison operator.
    /**
       \param se The right-hand side object -- \b IN.
    */
    virtual bool operator < ( const NOMAD::Set_Element<NOMAD::Signature> & se ) const
    {
      return ( *get_element() < *(se.get_element()) );
    }

    /// Access to the signature.
    /**
       \return A pointer to the signature.
    */
    NOMAD::Signature * get_signature ( void ) const
    {
      return const_cast<NOMAD::Signature *> ( get_element() );
    }
  };
}

#endif
