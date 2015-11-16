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
  \file   Double.hpp
  \brief  Custom class for double-precision reals (headers)
  \author Sebastien Le Digabel
  \date   2010-04-02
  \see    Double.cpp
*/
#ifndef __NOMAD_DOUBLE__
#define __NOMAD_DOUBLE__

#include "Exception.hpp"
#include "Display.hpp"

namespace NOMAD {

  /// Custom class for double-precision reals.
  /**
     - Allows comparisons on reals with custom precision.
     - Deals with undefined values.
     - Use \c value() to access the true \c double value.
  */
  class Double {
    
  private:
    
#ifdef MEMORY_DEBUG
    static int    _cardinality;      ///< Number of NOMAD::Double objects in memory.
    static int    _max_cardinality;  ///< Max number of NOMAD::Double objects in memory.
#endif

    double        _value;   ///< The value of the number.
    bool          _defined; ///< \c true if the number has a defined value.

    static double      _epsilon;    ///< Desired precision on comparisons.
    static std::string _inf_str;    ///< Infinity string.
    static std::string _undef_str;  ///< Undefined value string.

  public:
  
    /*-------------------------------------------------------------------*/

    /// Exception class for undefined NOMAD::Double objects.
    class Not_Defined : public NOMAD::Exception {
    public:
      /// Constructor.
      Not_Defined ( const std::string & file ,
		    int                 line ,
		    const std::string & msg    )
	: NOMAD::Exception ( file , line , msg ) {}
    };

    /// Exception class for divisons by zero with NOMAD::Double objects.
    class Invalid_Value : public NOMAD::Exception {
    public:
      /// Constructor.
      Invalid_Value ( const std::string & file ,
		      int                 line ,
		      const std::string & msg    )
	: NOMAD::Exception ( file , line , msg ) {}
    };

    /*-------------------------------------------------------------------*/
#ifdef MEMORY_DEBUG
    /// Access to the number of NOMAD::Double objects in memory.
    /**
       \return The number of NOMAD::Double objects in memory.
    */
    static int get_cardinality ( void ) { return NOMAD::Double::_cardinality; }

    /// Access to the max number of NOMAD::Double objects in memory.
    /**
       \return The max number of NOMAD::Double objects in memory.
    */
    static int get_max_cardinality ( void ) { return NOMAD::Double::_max_cardinality; }
#endif

    /// Constructor #1.
    explicit Double ( void );

    /// Constructor #2.
    /**
       From a \c double.
       \param v The \c double -- \b IN.
    */
    Double ( double v );

    /// Copy constructor.
    /**
       \param d The copied object -- \b IN.
    */
    Double ( const Double & d );

    /// Destructor.
    virtual ~Double ( void );

    /// Function \c atof (conversion from a string to a real).
    /**
       The value is determined by a string.
       \param s The string -- \b IN.
       \return  A boolean equal to \c true if the string was valid.
    */
    bool atof ( const std::string & s );

    /// Function \c atof (conversion from a string to a real) with a relative value.
    /**
       The value is determined by a string that may begin with \c 'r' to
       indicate a proportion (relative value).
       \param s The string -- \b IN.
       \param relative A boolean indicating if the conversion was relative -- \b OUT.
       \return A boolean equal to \c true if the string was valid.
    */
    bool relative_atof ( const std::string & s , bool & relative );

    /// Reset the Double.
    void clear ( void ) { _value = 0.0; _defined = false; }
    
    /// Reset the Double.
    void reset ( void ) { clear(); }

    /// Affectation operator #1.
    /**
       \param d The right-hand side object -- \b IN.
       \return \c *this as the result of the affectation.
    */
    Double & operator = ( const Double & d );

    /// Affectation operator #2.
    /**
       \param r The right-hand side \c double -- \b IN.
       \return \c *this as the result of the affectation.
    */
    Double & operator = ( double r );

    /// Access to the \c double value.
    /**
       \return The \c double value.
    */
    double value ( void ) const;

    /// Is the value defined ?
    /**
       \return A boolean equal to \c true if the NOMAD::Double is defined.
    */
    bool is_defined  ( void ) const { return _defined; }

    /// Is the value an integer ?
    /**
       \return A boolean equal to \c true if the NOMAD::Double is an integer.
    */
    bool is_integer ( void ) const;

    /// Is the value binary ?
    /**
       \return A boolean equal to \c true if the NOMAD::Double is in {0,1}.
    */
    bool is_binary ( void ) const;
    
    /// Sizeof.
    /**
       The size of an object in memory.
       \return An integer value for the size in bytes.
    */
    int size_of ( void ) const { return sizeof(_defined) + sizeof(_value); }

    /// Access to the NOMAD::Double precision.
    /**
       \return A \c double for the precision.
    */
    static double get_epsilon ( void )  { return Double::_epsilon; }

    /// Set the NOMAD::Double precision.
    /**
       \param eps The precision -- \b IN.
    */
    static void set_epsilon ( double eps );

    /// Access to the undefined value string.
    /**
       \return A \c std::string for the undefined value string.
    */
    static std::string get_undef_str ( void ) { return Double::_undef_str; }

    /// Set the NOMAD::Double undefined value string.
    /**
       \param undef_str The undefined value string -- \b IN.
    */
    static void set_undef_str ( const std::string & undef_str )
    {
      NOMAD::Double::_undef_str = undef_str;
    }

    /// Access to the infinity string.
    /**
       \return A \c std::string for the infinity string.
    */
    static std::string get_inf_str ( void ) { return Double::_inf_str; }

    /// Set the NOMAD::Double infinity string.
    /**
       \param inf_str The infinity string -- \b IN.
    */
    static void set_inf_str ( const std::string & inf_str )
    {
      NOMAD::Double::_inf_str = inf_str;
    }

    /// Rounding.
    /**
       \return An integer with the rounding result.
    */
    int round ( void ) const;
	  

	/// Ceil.
    /**
       \return A NOMAD::double with the rounding result.
    */
    const Double ceil ( void ) const;

	/// Floor.
    /**
       \return A NOMAD::double with the rounding result.
    */
    const Double floor ( void ) const;


    /// Absolute value.
    /**
       \return Another NOMAD::Double equal to max{\c -*this,\c *this}.
    */
    const Double abs ( void ) const;

    /// Square.
    /**
       \return Another NOMAD::Double equal to \c *this \c * \c *this.
    */
    const Double pow2 ( void ) const;

    /// Square root.
    /**
       \return Another NOMAD::Double equal to \c (*this)^0.5.
    */
    const Double sqrt ( void ) const;
    
    /// Relative error with another NOMAD::Double.
    /**
       \param  x The other real -- \b IN.
       \return A new object with the relative error value in \c [0;1].
    */
    const Double rel_err ( const Double & x ) const;

    /// Operator \c ++ (prefix position).
    /**
       Allows \c ++d;
       \return The NOMAD::Double plus 1.
    */
    Double & operator++ ( void );

    /// Operator \c ++ (suffix position).
    /**
       Allows \c d++;
       \return A copy of the NOMAD::Double before incrementation.
    */
    Double operator++ ( int n );

    /// Operator \c -- (prefix position).
    /**
       Allows \c --d;
       \return The NOMAD::Double minus 1.
    */    Double & operator-- ( void );

    /// Operator \c -- (suffix position).
    /**
       Allows \c d--;
       \return A copy of the NOMAD::Double before decrementation.
    */    Double operator-- ( int n );

    /// Operator \c +=.
    /**
       Allows \c d \c += \c d1.
       \param d1 The increment -- \b IN.
       \return The NOMAD::Double plus \c d1.
    */
    const Double & operator += ( const Double & d1 );

    /// Operator \c -=.
    /**
       Allows \c d \c -= \c d1.
       \param d1 The decrement -- \b IN.
       \return The NOMAD::Double minus \c d1.
    */
    const Double & operator -= ( const Double & d1 );

    /// Operator \c *=.
    /**
       Allows \c d \c *= \c d1.
       \param d1 The factor -- \b IN.
       \return The NOMAD::Double times \c d1.
    */
    const Double & operator *= ( const Double & d1 );

    /// Operator \c /=.
    /**
       Allows \c d \c /= \c d1. Throws a NOMAD::Exception::Invalid_Value if \c d1==0.
       \param d1 The denominator -- \b IN.
       \return The NOMAD::Double divided by \c d1.
    */
    const Double & operator /= ( const Double & d1 );

    /// Comparison <.
    /**
       The same as operator \c < but with consideration of undefined values.
       When comparing \c d1 and \c d2, \c (d1<d2) is equal to \c true for example
       if \c d1 is defined and if \c d2 is not.
       \param d The right-hand side of the comparison -- \b IN.
       \return A boolean equal to \c true if \c *this \c < \c d.
    */
    bool comp_with_undef ( const Double & d ) const;

    /// Projection to the mesh.
    /**
       Projection to the mesh of size delta
       ( \c *this \c = \c ref \c + \c k \c * \c delta ).
       \param ref    Reference for projection -- \b IN.
       \param delta  Mesh size parameter -- \b IN.
       \param lb     Lower bound -- \b IN -- \b optional (default=undefined real.)
       \param ub     Upper bound -- \b IN -- \b optional (default=undefined real.)
    */
    void project_to_mesh ( const Double & ref           ,
			   const Double & delta         ,
			   const Double & lb = Double() ,
			   const Double & ub = Double()   );

    /// Display.
    /**
       \param out The NOMAD::Display object -- \b IN.
    */
    void display ( const NOMAD::Display & out ) const;

    /// Display with format.
    /**
       \param out    The NOMAD::Display object -- \b IN.
       \param format The format as in C        -- \b IN.
    */
    void display ( const NOMAD::Display & out    ,
		   const std::string    & format   ) const;
  };

  /*---------------------------------------------------------------------------*/

  /// Display a NOMAD::Double object.
  /**
     \param out The NOMAD::Display object -- \b IN.
     \param d   The NOMAD::Double object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
  */
  inline const NOMAD::Display & operator << ( const NOMAD::Display & out ,
  					      const Double         & d     )
  {
    d.display ( out );
    return out;
  }
  
  /// Input.
  /**
     - Allows the input of NOMAD::Double objects with operator \c >>.
     - Can read undefined values (parameter \c UNDEF_STR with default \c "-".)
     - Example:
     \code
     NOMAD::Double d1 , d2;
     std::cout << "Enter d1 and d2: ";
     std::cin  >> d1 >> d2;
     std::cout << "d1 and d2 are equal to " << d1 << " and " << d2 << std::endl;
     \endcode
     \param in A \c std::istream object (can be a file) -- \b IN/OUT.
     \param d  The NOMAD::Double object to be read -- \b OUT.
     \return The modified \c std::istream object.
  */
  std::istream & operator >> ( std::istream & in , Double & d );

  /// Inverse operator.
  /**
     Allows operations such as \c d \c = \c -d.
     \param d The NOMAD::Double to be inversed -- \b IN.
     \return A new NOMAD::Double as the result of the operation.
  */
  inline const Double operator - ( const Double & d )
  {
    return NOMAD::Double (-d.value());
  }

  /// Operator \c + for two NOMAD::Double objects.
  /**
     Allows operations such as \c d \c = \c d1 \c + \c d2.
     \param d1 The first NOMAD::Double -- \b IN.
     \param d2 The second NOMAD::Double -- \b IN.
     \return   A third NOMAD::Double as the result of the operation.
  */
  inline const Double operator + ( const Double & d1 , const Double & d2 )
  {
    return NOMAD::Double ( d1.value() + d2.value() );
  }

  /// Operator \c - for two NOMAD::Double objects.
  /**
     Allows operations such as \c d \c = \c d1 \c - \c d2.
     \param d1 The first NOMAD::Double -- \b IN.
     \param d2 The second NOMAD::Double -- \b IN.
     \return   A third NOMAD::Double as the result of the operation.
  */
  inline const Double operator - ( const Double & d1 , const Double & d2 )
  {
    return NOMAD::Double (d1.value() - d2.value());
  }

  /// Operator \c * for two NOMAD::Double objects.
  /**
     Allows operations such as \c d \c = \c d1 \c * \c d2.
     \param d1 The first NOMAD::Double -- \b IN.
     \param d2 The second NOMAD::Double -- \b IN.
     \return   A third NOMAD::Double as the result of the operation.
  */
  inline const Double operator * ( const Double & d1 , const Double & d2 )
  {
    return NOMAD::Double ( d1.value() * d2.value() );
  }

  /// Operator \c / for two NOMAD::Double objects.
  /**
     Allows operations such as \c d \c = \c d1 \c / \c d2.
     \param d1 The first NOMAD::Double -- \b IN.
     \param d2 The second NOMAD::Double -- \b IN.
     \return   A third NOMAD::Double as the result of the operation.
  */
  const Double operator / ( const Double & d1 , const Double & d2 );

  /// Comparison operator \c ==.
  /**
     Allows the comparison \c d1 \c == \c d2.
     \param d1 The first NOMAD::Double -- \b IN.
     \param d2 The second NOMAD::Double -- \b IN.
     \return A boolean equal to true if \c d1 \c == \c d2.
  */
  inline bool operator == ( const Double & d1 , const Double & d2 )
  {
    return fabs ( d1.value() - d2.value() ) < NOMAD::Double::get_epsilon();
  }

  /// Comparison operator \c !=.
  /**
     Allows the comparison \c d1 \c != \c d2.
     \param d1 The first NOMAD::Double -- \b IN.
     \param d2 The second NOMAD::Double -- \b IN.
     \return A boolean equal to true if \c d1 \c != \c d2.
  */
  inline bool operator != ( const Double & d1 , const Double & d2 ) { return !(d1==d2); }

  /// Comparison operator \c <.
  /**
     Allows the comparison \c d1 \c < \c d2.
     \param d1 The first NOMAD::Double -- \b IN.
     \param d2 The second NOMAD::Double -- \b IN.
     \return A boolean equal to true if \c d1 \c < \c d2.
  */
  inline bool operator < ( const Double & d1 , const Double & d2 )
  {
    return d1.value() < d2.value() - NOMAD::Double::get_epsilon();
  }

  /// Comparison operator \c >.
  /**
     Allows the comparison \c d1 \c > \c d2.
     \param d1 The first NOMAD::Double -- \b IN.
     \param d2 The second NOMAD::Double -- \b IN.
     \return A boolean equal to true if \c d1 \c > \c d2.
  */
  inline bool operator > ( const Double & d1 , const Double & d2 )
  {
    return d1.value() > d2.value() + NOMAD::Double::get_epsilon();
  }

  /// Comparison operator \c <=.
  /**
     Allows the comparison \c d1 \c <= \c d2.
     \param d1 The first NOMAD::Double -- \b IN.
     \param d2 The second NOMAD::Double -- \b IN.
     \return A boolean equal to true if \c d1 \c <= \c d2.
  */
  inline bool operator <= ( const Double & d1 , const Double & d2 ) { return !(d1>d2); }

  /// Comparison operator \c >=.
  /**
     Allows the comparison \c d1 \c >= \c d2.
     \param d1 The first NOMAD::Double -- \b IN.
     \param d2 The second NOMAD::Double -- \b IN.
     \return A boolean equal to true if \c d1 \c >= \c d2.
  */  inline bool operator >= ( const Double & d1 , const Double & d2 )
  {
    return !(d1<d2);
  }
}
#endif
