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
  \file   Display.hpp
  \brief  Custom class for display (headers)
  \author Sebastien Le Digabel
  \date   2010-03-30
  \see    Display.cpp
*/
#ifndef __DISPLAY__
#define __DISPLAY__

#include "utils.hpp"

namespace NOMAD {

  /// Custom display class.
  /**
     - This class is used instead of \c std::ostream ojects.
     - Use \c std::endl as new line character; \c '\\n' will ignore indentation.
     - Do not use \c << \c flush; : it would create a new indentation.
       Use method \c Display::flush() instead.
     
     \b Two \b examples \b for \b creating \b indented \b blocks:

     \code
     using namespace NOMAD;
     Display out ( std::cout );
     \endcode
     
     \b Example \b 1:

     \code
     out << "line #1" << std::endl
         << open_block()
         << "line #2" << std::endl << "line #3" << std::endl
         << open_block ( "begin of block 2" )
         << "line #4" << std::endl << "line #5" << std::endl
         << close_block ( "end of block 2" )
         << close_block()
         << std::endl;
     \endcode

     \b Example \b 2:

     \code
     out << "line #1" << std::endl;
     out.open_block();
     out << "line #2" << std::endl << "line #3" << std::endl;
     out.open_block ("begin of block 2");
     out << "line #4" << std::endl << "line #5" << std::endl;
     out.close_block ("end of block 2");
     out.close_block();
     out << std::endl;
     \endcode

     \b This \b displays \b twice:

     \verbatim
     line #1
     {
         line #2
	 line #3
         begin of block 2 {
	     line #4
	     line #5
	 } end of block 2
     }
     \endverbatim
  */
  class Display {

  private:

    std::ostream      & _out;           ///< Display.

    mutable std::string _indent_str;    ///< Indentation string (tabulations).
    mutable bool        _newline;       ///< Indent or not.
    
    std::string         _open_brace;    ///< Braces of the indentation blocks.
    std::string         _closed_brace;  ///< Defaults: "{" and "}".

    NOMAD::dd_type      _gen_dd;        ///< General display degree.
    NOMAD::dd_type      _search_dd;     ///< Search display degree.
    NOMAD::dd_type      _poll_dd;       ///< Poll display degree.
    NOMAD::dd_type      _iter_dd;       ///< Iterative display degree.

    /// Private affectation operator.
    /**
       \param out The right-hand side object -- \b IN.
    */
    const Display & operator = ( const Display & out );

  public:
    
    /// Constructor.
    /**
       \param out A \c std::ostream that will be used for all displays
                  (can be a \c std::ofstream)
		  -- \b IN -- \b optional (default = \c std::cout).
     */
    Display ( std::ostream & out = std::cout )
      : _out          ( out                   ) , // can be a std::ofstream
	_newline      ( true                  ) ,
	_open_brace   ( "{"                   ) ,
	_closed_brace ( "}"                   ) ,
	_gen_dd       ( NOMAD::NORMAL_DISPLAY ) ,
	_search_dd    ( NOMAD::NORMAL_DISPLAY ) ,
	_poll_dd      ( NOMAD::NORMAL_DISPLAY ) ,
	_iter_dd      ( NOMAD::NORMAL_DISPLAY )   {}

    /// Copy constructor.
    /**
       \param out The copied object -- \b IN.
    */
    Display ( const Display & out )
      : _out          ( out._out          ) ,
	_indent_str   ( out._indent_str   ) ,
	_newline      ( out._newline      ) ,
	_open_brace   ( out._open_brace   ) ,
	_closed_brace ( out._closed_brace ) ,
	_gen_dd       ( out._gen_dd       ) ,
	_search_dd    ( out._search_dd    ) ,
	_poll_dd      ( out._poll_dd      ) ,
	_iter_dd      ( out._iter_dd      )   {}

    /// Destructor.
    virtual ~Display ( void ) {}

    /// Flush.
    /**
       Must be used instead of \c out \c << \c std::flush.
    */
    void flush ( void ) const { _out << std::flush; }

    /*---------------*/
    /*  GET methods  */
    /*---------------*/

    /// Access to the indentation string.
    /**
      \return The indentation string.
    */
    const std::string get_indent_str ( void ) const { return _indent_str; }

    /// Access to the general display degree.
    /**
       \return _gen_dd.
    */
    NOMAD::dd_type get_gen_dd ( void ) const { return _gen_dd; }

    /// Access to the search display degree.
    /**
       \return _search_dd.
    */
    NOMAD::dd_type get_search_dd ( void ) const { return _search_dd; }

    /// Access to the poll display degree.
    /**
       \return _poll_dd.
    */
    NOMAD::dd_type get_poll_dd ( void ) const { return _poll_dd; }

    /// Access to the iterative display degree.
    /**
       \return _iter_dd.
    */
    NOMAD::dd_type get_iter_dd ( void ) const { return _iter_dd; }

    /// Get the display degree for a specific search type.
    /**
       \param search The search type.
       \return       The display degree.
    */
    NOMAD::dd_type get_display_degree ( NOMAD::search_type search ) const;

    /// Get the display degrees as a string of size 4.
    /**
       \param dd The string containing the display degrees -- \b OUT.
    */
    void get_display_degree ( std::string & dd ) const;

    /*---------------*/
    /*  SET methods  */
    /*---------------*/

    /// Set the indentation string.
    /**
      \param is The indentation string -- \b IN.
    */
    void set_indent_str ( const std::string & is ) { _indent_str = is; }

    /// Set the _open_brace string.
    /**
       \param ob The string -- \b IN.
    */
    void set_open_brace   ( const std::string & ob ) { _open_brace = ob; }

    /// Set the _closed_brace string.
    /**
       \param cb The string -- \b IN.
    */
    void set_closed_brace ( const std::string & cb ) { _closed_brace = cb; }

    /// Set the display degrees.
    /**
       \param gen_dd    General display degree   -- \b IN.
       \param search_dd Search display degree    -- \b IN.
       \param poll_dd   Poll display degree      -- \b IN.
       \param iter_dd   Iterative display degree -- \b IN.
    */
    void set_degrees ( NOMAD::dd_type gen_dd    ,
		       NOMAD::dd_type search_dd ,
		       NOMAD::dd_type poll_dd   ,
		       NOMAD::dd_type iter_dd     );

    /// Set all the display degrees to one given display degree.
    /**
       \param dd The 4 display degrees -- \b IN.
    */
    void set_degrees ( NOMAD::dd_type dd ) { set_degrees ( dd , dd , dd , dd ); }

    /// Open an indentation block.
    /**
       \param msg Message displayed as the block title
                  -- \b IN -- \b optional (default = empty string).
    */
    void open_block  ( const std::string & msg = "" ) const;
    
    /// Close an indentation block.
    /**
       \param msg Message displayed at the end of the block
                  -- \b IN -- \b optional (default = empty string).
    */
    void close_block ( const std::string & msg = "" ) const;

    /// Operator <<.
    template <class T>
    const Display & operator << ( const T & ) const;
    
    /// Defines the \c cout type.
    typedef std::basic_ostream<char, std::char_traits<char> > CoutType;

    /// Allows the use of \c out \c << \c endl (\c std::endl is used).
    /**
       \warning This considers also \c std::flush.
    */
    typedef CoutType& (*StandardEndLine)(CoutType&);

    /// Operator << for \c std::endl.
    const Display & operator << ( StandardEndLine ) const;  

    /// Set the display precision.
    /**
       \param p The display precision -- \b IN.
    */
    void precision ( int p ) const { _out.precision(p); }

    /// Get the current display precision.
    /**
       \return An integer for the current precision.
    */
    int precision ( void ) const { return static_cast<int>(_out.precision()); }

    /// Set the format flags (1/2).
    /**
       \param f The flags -- \b IN.
    */
    void flags ( std::ios_base::fmtflags f ) const { _out.flags(f); }

    /// Set the format flags (2/2).
    /**
       \param f The flags -- \b IN.
    */
    void setf ( std::ios_base::fmtflags f ) const { _out.setf(f); }

    /// Unset the format flags.
    /**
       \param f The flags -- \b IN.
    */
    void unsetf ( std::ios_base::fmtflags f ) const { _out.unsetf(f); }


    /// Get the current format flags.
    /**
       \return The format flags.
    */
    std::ios_base::fmtflags flags ( void ) const { return _out.flags(); }
    
    /*--------------------*/
    /*  type conversions  */
    /*--------------------*/

    /// Convert a NOMAD::dd_type to a character.
    /**
       \param dd The NOMAD::dd_type -- \b IN.
       \return   The character.
     */
    static char dd_to_char ( NOMAD::dd_type dd );

    /// Convert a NOMAD::dd_type to an integer.
    /**
       \param dd The NOMAD::dd_type -- \b IN.
       \return   The integer.
    */
    static int dd_to_int ( NOMAD::dd_type dd );

    /// Convert an integer to a NOMAD::dd_type.
    /**
       \param dd The integer -- \b IN.
       \return   The NOMAD::dd_type.
    */
    static NOMAD::dd_type int_to_dd ( int dd );

    /// Display a duration with a smart format.
    /**
       \param t Duration as an integer in seconds -- \b IN.
    */
    void display_time ( int t ) const;

    /// Display a boolean with format \c yes/no.
    /**
       \param b The boolean -- \b IN.
    */
    void display_yes_or_no ( bool b ) const { (*this) << ( (b) ? "yes" : "no" ); }

    /// Display a memory size.
    /**
       \param size The memory size.
    */
    void display_size_of ( float size ) const;

    /// Display an integer with a specific width.
    /**
       \param i     The integer to display -- \b IN.
       \param max_i Maximal value of \c i used to determine the display width
                    -- \b IN -- \b optional (default = \c -1).
    */
    void display_int_w ( int i , int max_i = -1 ) const;

    /// Get the keyword associated with a NOMAD::display_stats_type.
    /**
       \param dst The NOMAD::display_stats_type -- \b IN.
       \return    A string containing the keyword.
    */
    static std::string get_display_stats_keyword ( NOMAD::display_stats_type dst );

    /// Extract display format from a string.
    /**
       \param s      The string -- \b IN/OUT.
       \param format The format -- \b OUT.
    */
    static void extract_display_format ( std::string & s , std::string & format );

    /// Get the NOMAD::display_stats_type from a string.
    /**
       \param s The string -- \b IN.
       \return  The NOMAD::display_stats_type.
    */
    static NOMAD::display_stats_type get_display_stats_type ( const std::string & s );

  };

  /*-------------------------------------------------------------------------*/

  /// Open an indented block.
  /**
     Allows the use of \c out \c << \c open_block(msg).
  */
  class open_block {
  private:
    std::string _msg; ///< Message displayed as a block title.
  public:

    /// Constructor.
    /**
       Open an indented block.
       \param msg The block title
                  -- \b IN -- \b optional (default = empty string).
    */
    open_block ( const std::string & msg = "" ) : _msg ( msg ) {}

    /// Operator ().
    const Display & operator() ( const Display & out ) const {
      out.open_block ( _msg );
      return out;
    }
  };

  /*-------------------------------------------------------------------------*/

  /// Close an indented block.
  /**
     Allows the use of \c out \c << \c close_block(msg).
  */
  class close_block {
  private:
    std::string _msg; ///< Message displayed at the end of a block.
  public:

    /// Constructor.
    /**
       Close an indented block.
       \param msg Message displayed at the end of a block
                  -- \b IN -- \b optional (default = empty string).
    */
    close_block ( const std::string & msg = "" ) : _msg ( msg ) {}

    /// Operator ().
    const Display & operator() ( const Display & out ) const {
      out.close_block ( _msg );
      return out;
    }
  };

  /*-------------------------------------------------------------*/
  /*  display functions for enum types defined in 'defines.hpp'  */
  /*-------------------------------------------------------------*/

  /// Operator << for NOMAD::stop_type.
  std::ostream & operator << ( std::ostream & , NOMAD::stop_type );

  /// Operator << for NOMAD::dd_type.
  std::ostream & operator << ( std::ostream & , NOMAD::dd_type );
  
  /// Operator << for NOMAD::success_type.
  std::ostream & operator << ( std::ostream & , NOMAD::success_type );
  
  /// Operator << for NOMAD::bb_input_type.
  std::ostream & operator << ( std::ostream & , NOMAD::bb_input_type );

  /// Operator << for NOMAD::bb_output_type.
  std::ostream & operator << ( std::ostream & , NOMAD::bb_output_type );

  /// Operator << for NOMAD::interpolation_type.
  std::ostream & operator << ( std::ostream & , NOMAD::interpolation_type );

  /// Operator << for NOMAD::hnorm_type.
  std::ostream & operator << ( std::ostream & , NOMAD::hnorm_type );

  /// Operator << for NOMAD::search_type.
  std::ostream & operator << ( std::ostream & , NOMAD::search_type );

  /// Operator << for NOMAD::model_type.
  std::ostream & operator << ( std::ostream & , NOMAD::model_type );

  /// Operator << for NOMAD::TGP_mode_type.
  std::ostream & operator << ( std::ostream & , NOMAD::TGP_mode_type );

  /// Operator << for NOMAD::direction_type.
  std::ostream & operator << ( std::ostream & , NOMAD::direction_type );

  /// Operator << for NOMAD::check_failed_type.
  std::ostream & operator << ( std::ostream & , NOMAD::check_failed_type );

  /// Operator << for NOMAD::display_stats_type.
  std::ostream & operator << ( std::ostream & , NOMAD::display_stats_type );

  /// Operator << for NOMAD::eval_type.
  std::ostream & operator << ( std::ostream & , NOMAD::eval_type );

  /// Operator << for NOMAD::eval_status_type.
  std::ostream & operator << ( std::ostream & , NOMAD::eval_status_type );

  /// Operator << for NOMAD::multi_formulation_type.
  std::ostream & operator << ( std::ostream & , NOMAD::multi_formulation_type );

  /// Operator << for a vector of NOMAD::bb_intput_type.  
  std::ostream & operator << ( std::ostream                            & ,
			       const std::vector<NOMAD::bb_input_type> &   );

  /// Operator <<.
  template <class T>
  inline const NOMAD::Display & NOMAD::Display::operator << ( const T & t ) const
  {
    if ( _newline ) {
      _out << _indent_str;
      _newline = false;
    }
    _out << t;
    return *this;
  }

  /// Allows the use of \c out \c << \c endl.
  inline const NOMAD::Display & NOMAD::Display::operator << ( StandardEndLine m ) const
  {
    m ( _out ); // this could be a std::flush, so don't use it: instead use method flush()
    _newline = true;
    return *this;
  }

  /// Allows the use of \c out \c << \c open_block(msg).
  inline const NOMAD::Display & operator << ( const NOMAD::Display    & out ,
					      const NOMAD::open_block & ob    )
  {
    return ob ( out );
  }

  /// Allows the use of \c out \c << \c close_block(msg).
  inline const NOMAD::Display & operator << ( const NOMAD::Display     & out ,
					      const NOMAD::close_block & cb    )
  {
    return cb ( out );
  }

}

#endif
