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
  \file   Exception.hpp
  \brief  Custom class for exceptions (headers)
  \author Sebastien Le Digabel
  \date   2010-03-29
  \see    Exception.cpp
*/
#ifndef __NOMAD_EXCEPTION__
#define __NOMAD_EXCEPTION__

#include <sstream>

namespace NOMAD {

  /// Custom class for exceptions.
  /**
     NOMAD uses this type of exceptions.
     It indicates the file and line number at which a throw is made.

     \b Example

     \code
     throw NOMAD::Exception ( __FILE__ , __LINE__ , "an error message" );
     \endcode
   */
  class Exception : public std::exception {

  private:

    mutable std::string _what;  ///< Error message.
    std::string         _file;  ///< File where the exception is thrown.
    int                 _line;  ///< Line number at which the exception is thrown.

  public:

    /// Constructor.
    /**
       \param file A string corresponding to the file where the
                     exception is thrown -- \b IN
       \param line An integer corresponding to the line number
                     at which the exception is thrown -- \b IN.
       \param msg  A string corresponding to the error message -- \b IN.
     */
    Exception ( const std::string & file , int line , const std::string & msg )
      : _what ( msg  ) ,
	_file ( file ) ,
	_line ( line )   {}

    /// Destructor.
    virtual ~Exception ( void ) throw() {}

    /// Access to the error message.
    /**
       \return A string with the error message.
    */
    const char * what ( void ) const throw();
  };
}

#endif
