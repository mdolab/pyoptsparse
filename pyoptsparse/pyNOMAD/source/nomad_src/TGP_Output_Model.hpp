/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search - version 3.7.0.beta        */
/*                                                                                     */
/*  Copyright (C) 2001-2014  Mark Abramson        - the Boeing Company, Seattle        */
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
  \file   TGP_Output_Model.hpp
  \brief  TGP (Bayesian treed Gaussian process) model for one output (headers)
  \author Sebastien Le Digabel
  \date   2011-02-07
  \see    TGP_Output_Model.cpp
*/
#ifdef USE_TGP

#ifndef __TGP_OUTPUT_MODEL__
#define __TGP_OUTPUT_MODEL__

#include "Eval_Point.hpp"
#include "tgp.h"

/*------------------------------*/
/*  TGP C functions prototypes  */
/*------------------------------*/
extern "C"
{
  unsigned long   three2lstate   ( int           * state );
  void          * newRNGstate    ( unsigned long   s     );
  void            deleteRNGstate ( void          * seed  );
  unsigned int *  GetImprovRank  ( int, int, double **, int, int, double * );
}

namespace NOMAD {

  // NOMAD-TGP callback function (called regularly by TGP).
  // void TGP_callback ( bool & TGP_interrupt ); // SLD -- 2014-09-04

  /// TGP models for one output.
  class TGP_Output_Model : private NOMAD::Uncopyable {

  private:

    const NOMAD::Display & _out; ///< Display.

    int           _p;            ///< Number of interpolation points.

    double      * _Z;            ///< Output vector (size = \c p).

    double        _Z_scaling[2]; ///< To scale/unscale \c Z.
    bool          _Z_is_scaled;  ///< If \c Z is scaled or not.

    bool          _is_binary;    ///< If the output has only two values.
    NOMAD::Point  _bin_values;   ///< Binary output values.

    bool          _is_fixed;     ///< Only one output value saved in \c bin_values[0].

    void        * _tgp_state;    ///< RNG (random number generator).

    Model       * _tgp_model;    ///< The TGP model.

    Temper      * _tgp_its;      ///< Importance tempering object.

    static double _ditemps[7];   ///< Importance tempering parameters.

    static bool   _force_quit;   ///< Flag equal to \c true if ctrl-c is pressed.

    /// Treat binary output.
    /**
       \param ZZ   Output vector  -- \b IN/OUT.
       \param nout Size of output -- \b IN.
    */
    void treat_binary_output ( double * ZZ , int nout ) const;

    /// Scale member \c _Z.
    void scale_Z ( void );

    /// Unscale member \c _Z.
    void unscale_Z ( void );

    /// Scale an output \c Z.
    /**
       \param Z Output to scale -- \b IN/OUT.
       \param n Size of output  -- \b IN.
    */
    void scale_Z ( double * Z , int n ) const;

    /// Unscale an output \c Z.
    /**
       \param Z Output to unscale -- \b IN/OUT.
       \param n Size of output    -- \b IN.
    */
    void unscale_Z ( double * Z , int n ) const;

  public:

    /// Constructor.
    /**
       \param X_pts     Interpolation points with output values          -- \b IN.
       \param bbo_index Output index                                     -- \b IN.
       \param seed      Random seed (\c 0 to have the same \c R results) -- \b IN.
       \param out       Display object                                   -- \b IN.
    */
    explicit TGP_Output_Model
    ( const std::list<const NOMAD::Eval_Point *> & X_pts     ,
      int                                          bbo_index ,
      int                                          seed      ,
      const NOMAD::Display                       & out         );

    /// Destructor.
    virtual ~TGP_Output_Model ( void );

    /// Force quit (called by pressing ctrl-c).
    static void force_quit ( void )
    {
      NOMAD::TGP_Output_Model::_force_quit = true;
    }

    /// Access to the \c force_quit flag.
    /**
       \return The \c force_quit flag.
    */
    static bool get_force_quit ( void ) { return NOMAD::TGP_Output_Model::_force_quit; }

    /// Compute the model.
    /**
       \param X              \c X matrix (\c X_pts \c x \c n): interpolation pts -- \b IN.
       \param XX             \c XX matrix (\c n_XX \c x \c n): prediction points -- \b IN.
       \param Xsplit         \x Xsplit matrix (\c X plus \c XX)                  -- \b IN.
       \param n              Dimension and number of columns of the matrices     -- \b IN.
       \param n_XX           Number of rows of \c XX                             -- \b IN.
       \param nsplit         Number of rows of \c Xsplit (\c X_pts \c + \c n_XX) -- \b IN.
       \param tgp_params     TGP parameters                                      -- \b IN.
       \param tgp_rect       TGP rectangle                                       -- \b IN.
       \param tgp_BTE        TGP \c B,\c T, and \c R parameters                  -- \b IN.
       \param tgp_linburn    TGP \c linburn parameter                            -- \b IN.
       \param tgp_verb       TGP \c verb parameter                               -- \b IN.
       \param ZZ             \c ZZ vector (size \c n_XX): prediction values      -- \b OUT.
       \param Ds2x           Expected reduction in predictive var (size \c n_XX) -- \b OUT.
       \param improv         Expected improvement of the obj. (size \c n_XX)     -- \b OUT.
    */
    void compute ( double ** X           ,
		   double ** XX          ,
		   double ** Xsplit      ,
		   int       n           ,
		   int       n_XX        ,
		   int       nsplit      ,
		   Params  * tgp_params  ,
		   double ** tgp_rect    ,
		   int     * tgp_BTE     ,
		   bool      tgp_linburn ,
		   bool      tgp_verb    ,
		   double  * ZZ          ,
		   double  * Ds2x        ,
		   int     * improv        );

    /// Prediction at one point.
    /**
       \param XX         \c XX matrix (\c n_XX \c x \c n): prediction points -- \b IN.
       \param n          Dimension and number of columns of the matrices     -- \b IN.
       \param ZZ         \c ZZ vector (size \c n_XX): prediction values      -- \b OUT.
       \param tgp_rect   TGP rectangle                                       -- \b IN.
       \return A boolean equal to \c true if the prediction went well.
    */
    bool predict ( double  * XX       ,
		   int       n        ,
		   double  & ZZ       ,
		   double ** tgp_rect   ) const;

    /// Access to \c Z.
    /**
       \return The \c Z vector.
    */
    const double * get_Z ( void ) const { return _Z; }

    /// Access to the \c fixed flag.
    /**
       \return The \c fixed flag.
    */
    bool is_fixed ( void ) const { return _is_fixed; }

    /// Access to the \c binary flag.
    /**
       \return The \c binary flag.
    */
    bool is_binary ( void ) const { return _is_binary; }

    /// Default display.
    void display ( void ) { display ( _out ); }

    /// Display.
    /**
       \param out Display -- \b IN.
    */
    void display ( const NOMAD::Display & out ) const;
  };

  /// Display a NOMAD::TGP_Output_Model object.
  /**
     \param out The NOMAD::Display object                          -- \b IN.
     \param s   The NOMAD::TGP_Output_Model object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
  */
  inline const NOMAD::Display & operator << ( const NOMAD::Display          & out ,
					      const NOMAD::TGP_Output_Model & s     ) {
    s.display ( out );
    return out;
  }

}

#endif
#endif
