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
  \file   Parameter_Entry.hpp
  \brief  Parameter entry (headers)
  \author Sebastien Le Digabel
  \date   2010-04-05
  \see    Parameter_Entry.cpp
*/
#ifndef __PARAMETER_ENTRY__
#define __PARAMETER_ENTRY__

#include "Display.hpp"
#include "Uncopyable.hpp"

namespace NOMAD {

  /// Parameter entry.
  /**
     - Describes the data relative to a parameter in a parameters file.
     - Objets of this class are stored in a NOMAD::Parameter_Entries object.
  */
  class Parameter_Entry : private NOMAD::Uncopyable {

  private:

    std::string            _name;   ///< Name of the parameter.
    std::list<std::string> _values; ///< List of values for the parameter.
    bool                   _ok;     ///< If the parameter is valid.
    bool                   _unique; ///< If the parameter is unique.
    Parameter_Entry      * _next;   ///< Acces to the next parameter.

    /// If the parameter has been interpreted.
    bool _has_been_interpreted; 

  public:

    /// Constructor.
    /**
       Ignores all entries after \c '#'.
       \param entry           A string describing the parameter entry -- \b IN.
       \param remove_comments A boolean equal to \c true if entries after
                              \c '#' are ignored -- \b IN
			      -- \b optional (default = \c true).
    */
    Parameter_Entry ( const std::string & entry , bool remove_comments = true );

    /// Destructor.
    virtual ~Parameter_Entry ( void ) {}

    /*---------------*/
    /*  GET methods  */
    /*---------------*/

    /// Access to the name of the parameter.
    /**
       \return The name.
    */
    const std::string & get_name ( void ) const { return _name; }

    /// Access to the parameter values.
    /**
       \return The parameter values as a list of strings.
    */
    const std::list<std::string> & get_values ( void ) const { return _values; }

    /// Access to the number of values of the parameter.
    /**
       \return The number of values.
    */
    int get_nb_values ( void ) const { return static_cast<int>(_values.size()); }

    /// Access to the \c _ok flag.
    /**
       This flag is equal to \c true if the parameter entry is well defined.
       \return A boolean equal to \c true if the parameter is valid.
    */
    bool is_ok ( void ) const { return _ok; }

    /// Access to the \c _unique flag.
    /**
       This flag is decided when a parameters file is read.
       \return A boolean equal to \c true if the parameter is unique
       in a parameters file.
    */
    bool is_unique ( void ) const { return _unique; }

    /// Access to another NOMAD::Parameter_Entry.
    /**
       NOMAD::Parameter_Entry objects are stored in a NOMAD::Parameter_Entries
       object. The link between elements is assumed by the \c _next member
       returned by this function.
       \return A pointer to the next entry.
    */
    Parameter_Entry * get_next ( void ) const { return _next; }
    
    /// Access to the \c _has_been_interpreted flag.
    /**
       \return A boolean equal to \c true if the parameter has already
               been interpreted.
    */
    bool has_been_interpreted ( void ) const { return _has_been_interpreted; }

    /*---------------*/
    /*  SET methods  */
    /*---------------*/

    /// Set the \c _next pointer.
    /**
       \param p A pointer to the next NOMAD::Parameter_Entry to be inserted -- \b IN.
    */
    void set_next   ( Parameter_Entry * p ) { _next = p; }

    /// Set the \c _unique flag.
    /**
       \param u Value of the flag -- \b IN.
    */
    void set_unique ( bool u ) { _unique = u; }

    /// Set the \c _has_been_interpreted flag. to \c true.
    void set_has_been_interpreted  ( void ) { _has_been_interpreted = true; }
    
    /// Comparison with another entry.
    /**
       The comparison is based on the parameter name.
       \param p The right-hand side object -- \b IN.
       \return A boolean equal to \c true if \c this->_name \c < \c p._name.
    */
    bool operator < ( const Parameter_Entry & p ) const { return _name < p._name; }
    
    /// Display.
    /**
       \param out The NOMAD::Display object -- \b IN.
    */
    void display ( const NOMAD::Display & out ) const;
  };

  /// Allows the comparison of two NOMAD::Parameter_Entry objects.
  struct Parameter_Entry_Comp {
    /// Comparison of two NOMAD::Parameter_Entry objects.
    /**
       \param  p1 Pointer to the first NOMAD::Parameter_Entry  -- \b IN.
       \param  p2 Pointer to the second NOMAD::Parameter_Entry -- \b IN.
       \return A boolean equal to \c true if \c *p1 \c < \c *p2.
    */
    bool operator() ( const Parameter_Entry * p1 , const Parameter_Entry * p2 ) const
    {
      return (*p1 < *p2);
    }
  };

  /// Display a NOMAD::Parameter_Entry object.
  /**
     \param out The NOMAD::Display object -- \b IN.
     \param e   The NOMAD::Parameter_Entry object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
  */
  inline const NOMAD::Display & operator << ( const NOMAD::Display         & out ,
					      const NOMAD::Parameter_Entry & e     ) 
  {
    e.display ( out );
    return out;
  }
}

#endif
