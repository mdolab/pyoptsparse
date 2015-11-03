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
  \file   Point.hpp
  \brief  Custom class for points (headers)
  \author Sebastien Le Digabel
  \date   2010-04-05
  \see    Point.cpp
*/
#ifndef __POINT__
#define __POINT__

#include <numeric>
#include <functional>
#include "Double.hpp"

namespace NOMAD {

  /// \c accumulate operator for \c squared_norm().
  struct squared_norm_op {
    /** This function returns \c d1+d2^2 and is designed as an
	\c accumulate operator.
	\param d1 The NOMAD::Double \c d1 -- \b IN.
	\param d2 The NOMAD::Double \c d2 -- \b IN.
	\return A third NOMAD::Double equal to \c d1+d2^2.
    */
    NOMAD::Double operator () ( const NOMAD::Double & d1 , const NOMAD::Double & d2 ) {
      return d1 + d2.pow2();
    }
  };

  /// \c accumulate operator for \c size_of().
  struct size_of_op {
    /** This function returns \c i+size_of(d) and is designed as an
	\c accumulate operator.
	\param i The integer \c i -- \b IN.
	\param d The NOMAD::Double \c d -- \b IN.
	\return An integer equal to \c i+size_of(d).
    */
    int operator () ( int i , const NOMAD::Double & d ) {
      return i + d.size_of();
    }
  };

  /// Custom class for points.
  /**
     Coordinates are NOMAD::Double objects.
  */
  class Point {

  private:

#ifdef MEMORY_DEBUG
    static int _cardinality;     ///< Number of NOMAD::Point objects in memory.
    static int _max_cardinality; ///< Max number of NOMAD::Point objects in memory.
#endif

    /**
       Max number of coordinates to display.
       Default to 20, debug value at -1 (no limit).
    */
    static int _display_limit;
   
    int             _n;      ///< Dimension of the point.
    NOMAD::Double * _coords; ///< Coordinates of the point.

  public:

    /*-------------------------------------------------------------------*/

    /// Exception class for a bad access with NOMAD::Point objects.
    class Bad_Access : public NOMAD::Exception {
    public:
      /// Constructor.
      Bad_Access ( const std::string & file ,
		   int                 line ,
		   const std::string & msg    )
	: NOMAD::Exception ( file , line , msg ) {}
    };

    /// Exception class for NOMAD::Point objects that are not defined.
    class Not_Defined : public NOMAD::Exception {
    public:
      /// Constructor.
      Not_Defined ( const std::string & file ,
		    int                 line ,
		    const std::string & msg    )
	: NOMAD::Exception ( file , line , msg ) {}
    };

    /// Exception class for input errors with NOMAD::Point objects.
    class Bad_Input : public NOMAD::Exception {
    public:
      /// Constructor.
      Bad_Input ( const std::string & file ,
		  int                 line ,
		  const std::string & msg    )
	: NOMAD::Exception ( file , line , msg ) {}
    };
    
    /// Exception class for bad operations with NOMAD::Point objects.
    class Bad_Operation : public NOMAD::Exception {
    public:
      /// Constructor.
      Bad_Operation ( const std::string & file ,
		      int                 line ,
		      const std::string & msg    )
	: NOMAD::Exception ( file , line , msg ) {}
    };

    /*-------------------------------------------------------------------*/

#ifdef MEMORY_DEBUG
    /// Access to the number of NOMAD::Point objects in memory.
    /**
       \return The number of NOMAD::Point objects in memory.
    */
    static int get_cardinality ( void ) { return Point::_cardinality;     }

    /// Access to the max number of NOMAD::Point objects in memory.
    /**
       \return The max number of NOMAD::Point objects in memory.
    */
    static int get_max_cardinality ( void ) { return Point::_max_cardinality; }
#endif
    
    /// Constructor.
    /**
       \param n Dimension of the point -- \b IN --\b optional (default = 0).
       \param d Initial value for all coordinates
       -- \b IN -- \b optional (default = undefined real).
    */
    explicit Point ( int n = 0 , const NOMAD::Double & d = NOMAD::Double() );
    
    /// Copy constructor.
    /**
       \param p The copied object -- \b IN.
    */
    Point ( const Point & p );
    
    /// Affectation operator.
    /**
       \param p The right-hand side object -- \b IN.
       \return \c *this as the result of the affectation.
    */
    const Point & operator = ( const Point & p );

    /// Destructor.
    virtual ~Point ( void );

    /// Change the NOMAD::Point dimension.
    /**
       Sets also all coordinates to d.
       \param n New dimension -- \b IN --\b optional (default = 0).
       \param d Initial value for all coordinates
       -- \b IN -- \b optional (default = undefined real).
    */
    void reset ( int n = 0 , const NOMAD::Double & d = NOMAD::Double() );
    
    /// Clear the point.
    virtual void clear ( void ) { reset(); }
    
    /// Change the NOMAD::Point dimension.
    /**
       The values are kept.
       \param n New dimension of the point -- \b IN.
    */
    void resize ( int n );
    
    /// Set a new display limit.
    /**
       \param dl An integer for the new display limit -- \b IN.
    */
    static void set_display_limit ( int dl ) { Point::_display_limit = dl; }

    /// Access to the display limit.
    /**
       \return The display limit.
    */
    static int get_display_limit ( void ) { return Point::_display_limit; }
    
    /// Const operator \c [].
    /**
       \param i The index (0 for the first element) -- \b IN.
       \return The \c (i+1)th coordinate.
    */
    const NOMAD::Double & operator [] ( int i ) const;

    /// Non-const operator \c [].
    /**
       \param i The index (0 for the first element) -- \b IN.
       \return The \c (i+1)th coordinate.
    */
    NOMAD::Double & operator [] ( int i );
    
    /// Size of the point in memory.
    /**
       \return An integer for the size of the point in bytes.
    */
    virtual int size_of ( void ) const
    {
      return std::accumulate ( _coords                      ,
			       _coords+_n                   ,
			       static_cast<int>(sizeof(_n)) ,
			       size_of_op()                   );
    }

    /// Access to the dimension of the point.
    /**
       \return The dimension of the point.
    */
    int size ( void ) const { return _n; }

    /// Test if the point is empty (dimension equal to zero).
    /**
       \return A boolean equal to \c true if the point is empty.
    */
    bool empty ( void ) const { return _n==0; }
    
    /// Set all the coordinates to a specifi value.
    /**
       \param d The value for all coordinates -- \b IN.
    */
    void set ( const NOMAD::Double & d ) const
    { std::fill ( _coords , _coords+_n , d ); }

    /// Set the coordinates with an array of reals.
    /**
       \param n Dimension of the point -- \b IN.
       \param a Array of size \c n of reals -- \b IN.
    */
    void set ( int n , const NOMAD::Double * a );

    /// Check if all the coordinates are defined.
    /**
       \return A boolean equal to \c true if all the coordinates are defined.
    */
    bool is_complete ( void ) const;

    /// Check if at least one coordinate is defined.
    /**
       This virtual method is redefined in class NOMAD::Direction.
       \return A boolean equal to \c true if at least one coordinate is defined.
    */
    virtual bool is_defined ( void ) const;
    
    /// Count the number of defined values.
    /**
       \return The number of values that are defined.
    */
    int nb_defined ( void ) const;

    /// Squared norm of the point.
    /**
       \return A NOMAD::Double with the squared norm of the point.
    */
    const NOMAD::Double squared_norm ( void ) const
    {
      return std::accumulate ( _coords            ,
			       _coords+_n         ,
			       NOMAD::Double(0.0) ,
			       squared_norm_op()    );
    }

    /// Norm of the point.
    /**
       \return A NOMAD::Double with the norm of the point.
    */
    const NOMAD::Double norm ( void ) const { return squared_norm().sqrt(); }
    
    /// Dot product with another point \c x.
    /**
       \param x The other point \c x -- \b IN.
       \return The dot product \c *this \c . \c x.
    */
    const NOMAD::Double dot_product ( const Point & x ) const
    { return std::inner_product ( _coords , _coords+_n , x._coords , NOMAD::Double(0.0) ); }

    /// Angle with another point \c x.
    /**
       \param x The other point \c x -- \b IN.
       \return The angle between \c *this and \c x.
    */
    const NOMAD::Double get_angle ( const Point & x ) const;
    
    /// Mutiplication with a scalar.
    /**
       - This implements \c *this \c = \c d \c * \c *this.
       - The current object \c *this is modified.
       \param d The scalar -- \b IN.
       \return The point times \c d.
    */
    const Point & operator *= ( const NOMAD::Double & d );
    
    /// Multiplication with another point.
    /**
       - The multiplication is done coordinate by coordinate.
       - The current object \c *this is not modified.
       \param p The other point -- \b IN.
       \return A third point equal to \c *this \c .* \c p.
    */
    const Point operator * ( const Point & p ) const;

    /// Division with another point.
    /**
       - The division is done coordinate by coordinate.
       - The current object \c *this is not modified.
       \param p The other point -- \b IN.
       \return A third point equal to \c *this \c ./ \c p.
    */
    const Point operator / ( const Point & p ) const;

    /// Addition with another point.
    /**
       The current object \c *this is not modified.
       \param p The other point -- \b IN.
       \return A third point equal to \c *this \c + \c p.
    */
    const Point operator + ( const Point & p ) const;
    
    /// Substraction with another point.
    /**
       The current object \c *this is not modified.
       \param p The other point -- \b IN.
       \return A third point equal to \c *this \c - \c p.
    */
    const Point operator - ( const Point & p ) const;

    /// Negation.
    /**
       The current object \c *this is not modified.
       \return A new point equal to \c -*this.
    */
    const Point operator - ( void ) const;
    
    /// Comparison operator \c <.
    /**
       \param p The right-hand side object -- \b IN.
       \return A boolean equal to \c true if  \c *this \c < \c p.
    */
    virtual bool operator <  ( const Point & p ) const;


    /// Comparison operator \c ==.
    /**
       \param p The right-hand side object -- \b IN.
       \return A boolean equal to \c true if  \c *this \c == \c p.
    */
    bool operator == ( const Point & p ) const;

    /// Comparison operator \c !=.
    /**
       \param p The right-hand side object -- \b IN.
       \return A boolean equal to \c true if  \c *this \c != \c p.
    */
    bool operator != ( const Point & p ) const { return !(*this == p); }
    
    /// The same as operator \c < but with consideration of undefined values.
    /**
       \param p The right-hand side object -- \b IN.
       \return A boolean equal to \c true if \c *this \c < \c p.
    */
    bool comp_with_undef ( const Point & p ) const;
    
    /// Projection to the mesh.
    /**
       Projection to the mesh of size delta
       ( \c *this \c = \c ref \c + \c k \c * \c delta ).
       \param ref   Reference for projection -- \b IN.
       \param delta Mesh size parameter -- \b IN.
       \param lb    Lower bound -- \b IN -- \b optional
                    (default = undefined NOMAD::Point).
       \param ub    Upper bound -- \b IN -- \b optional
                    (default = undefined NOMAD::Point).
    */
    void project_to_mesh ( const Point & ref          ,
			   const Point & delta        ,
			   const Point & lb = Point() ,
			   const Point & ub = Point()   );

    /// Display.
    /**
       \param out The NOMAD::Display object -- \b IN.
       \param sep A string that is used as a separator between the coordinates
                  -- \b IN --\b optional (default = one space).
       \param w   An integer indicating a width for the display of
                  each coordinate -- \b IN -- \b optional
		  (default = -1, no limit).
       \param lim Max number of coordinates to display -- \b IN
              -- \b optional (default = -1, no limit).
    */
    virtual void display ( const NOMAD::Display & out       ,
			   const std::string    & sep = " " ,
			   int                    w   = -1  ,
			   int                    lim = -1    ) const;
	  	  
  };

  /*---------------------------------------------------------------------------*/

  /// Display a NOMAD::Point object.
  /**
     \param out The NOMAD::Display object -- \b IN.
     \param p   The NOMAD::Point object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
  */
  inline const NOMAD::Display & operator << ( const NOMAD::Display & out ,
					      const NOMAD::Point   & p     ) {
    p.display ( out , " " , 2 , NOMAD::Point::get_display_limit() );
    return out;
  }

  /// Input.
  /**
     - Allows the input of NOMAD::Double objects with operator \c >>.
     - Can read undefined coordinates (parameter \c UNDEF_STR with default \c "-".)
     - Example:
     \code
     NOMAD::Point x(3);
     std::cout << "Enter x (3 coordinates): ";
     std::cin  >> x;
     std::cout << "x is equal to " << x << std::endl;
     \endcode
     \param in A \c std::istream object (can be a file) -- \b IN/OUT.
     \param p  The NOMAD::Point object to be read -- \b OUT.
     \return The modified \c std::istream object.
  */
  std::istream & operator >> ( std::istream & in , Point & p );
}

#endif
