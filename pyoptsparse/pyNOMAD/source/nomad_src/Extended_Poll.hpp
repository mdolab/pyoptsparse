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
  \file   Extended_Poll.hpp
  \brief  Extended poll for categorical variables (headers)
  \author Sebastien Le Digabel
  \date   2010-04-14
  \see    Extended_Poll.cpp
*/
#ifndef __EXTENDED_POLL__
#define __EXTENDED_POLL__

#include "Mads.hpp"
#include "Signature_Element.hpp"

namespace NOMAD {

  /// Generic class for the extended poll.
  /**
     This is an abstract class (it is not possible to create
     NOMAD::Extended_Poll objects).
  */
  class Extended_Poll : public NOMAD::Uncopyable {

  protected:

    /// Parameters (includes the standard signature).
    NOMAD::Parameters  & _p;
    
    /// Add an extended poll point.
    /**
       Get, check and register the extended point and its signature
       created by the user in \c construct_extended_points().
       \param ep Extended poll point  -- \b IN.
       \param s  Associated signature -- \b IN.
    */
    void add_extended_poll_point ( NOMAD::Point & ep , NOMAD::Signature & s );

    /*---------------------------------------------------------------------*/

  private:

    /// Executable for getting neighbors in batch mode.
    std::string _neighbors_exe;

    /// Set of signatures (does not include the standard signature).
    std::set<NOMAD::Signature_Element> _signatures;  

    /// Signatures used during one poll step.
    std::set<NOMAD::Signature_Element> _poll_signatures;

    /// Extended points for one poll step.
    std::list<NOMAD::Eval_Point *> _extended_points;

    /*---------------------------------------------------------------------*/

    /// Evaluation of an extended poll point.
    /**
       \param y              The extended poll point              -- \b IN/OUT.
       \param mads           NOMAD::Mads object invoking the poll -- \b IN/OUT.
       \param stop           Stop flag                            -- \b IN/OUT.
       \param stop_reason    Stop reason                          -- \b OUT.
       \param success        Type of success                      -- \b OUT.
       \param new_feas_inc   New feasible incumbent               -- \b IN/OUT.
       \param new_infeas_inc New infeasible incumbent             -- \b IN/OUT.
       \return A pointer to the evaluated point; may be \c NULL if
               the evaluation failed.
    */
    const NOMAD::Eval_Point * eval_epp
    ( NOMAD::Eval_Point        * y              ,
      Mads                     & mads           ,
      bool                     & stop           ,
      NOMAD::stop_type         & stop_reason    ,
      NOMAD::success_type      & success        ,
      const NOMAD::Eval_Point *& new_feas_inc   ,
      const NOMAD::Eval_Point *& new_infeas_inc   ) const;
    
    /// Sort the evaluated extended poll points.
    /**
       \param evaluated_pts The list of evaluated extended poll points -- \b IN/OUT.
    */
    void sort_epp ( const std::list<const NOMAD::Eval_Point *> & evaluated_pts );
    
    /// Check the extended poll trigger.
    /**
       \param old_bf A pointer to the old best feasible point   -- \b IN.
       \param old_bi A pointer to the old best infeasible point -- \b IN.
       \param y      A pointer to the extended poll center      -- \b IN.
       \return A boolean equal to \c true if the extended poll has to be performed.
    */
    bool check_trigger ( const NOMAD::Eval_Point * old_bf ,
			 const NOMAD::Eval_Point * old_bi ,
			 const NOMAD::Eval_Point * y        ) const;
    
    /// Check only the \c f values for the extended poll trigger.
    /**
       \param old_f Old \c f value -- \b IN.
       \param new_f New \c f value -- \b IN.
       \return A boolean equal to \c true if the extended poll has to be performed.
    */
    bool check_trigger_on_f ( const NOMAD::Double & old_f  ,
			      const NOMAD::Double & new_f    ) const;
    
    /// Create the descent parameters.
    /**
       \param y     Starting point         -- \b IN.
       \param stats Stats                  -- \b IN.
       \param descent_p Descent parameters -- \b IN/OUT.
    */
    void set_descent_parameters ( const NOMAD::Eval_Point  * y         ,
				  const NOMAD::Stats       & stats     ,
				  NOMAD::Parameters        & descent_p   ) const;
    
    /// Descent from the extended poll center.
    /**
       \param y               Extended poll center                 -- \b IN.
       \param mads            NOMAD::Mads object invoking the poll -- \b IN/OUT.
       \param nb_ext_poll_pts Number of extended poll points       -- \b IN/OUT.
       \param stop           Stop flag                                -- \b IN/OUT.
       \param stop_reason    Stop reason                              -- \b OUT.
       \param success        Type of success                          -- \b OUT.
       \param new_feas_inc   New feasible incumbent                   -- \b IN/OUT.
       \param new_infeas_inc New infeasible incumbent                 -- \b IN/OUT.
    */
    void descent ( const NOMAD::Eval_Point  * y               ,
		   Mads                     & mads            ,
		   int                      & nb_ext_poll_pts ,
		   bool                     & stop            ,
		   NOMAD::stop_type         & stop_reason     ,
		   NOMAD::success_type      & success         ,
		   const NOMAD::Eval_Point *& new_feas_inc    ,
		   const NOMAD::Eval_Point *& new_infeas_inc    );
    
    /*---------------------------------------------------------------------*/

  public:

    /// Constructor.
    /**
       \param p Parameters -- \b IN.
    */
    Extended_Poll ( NOMAD::Parameters & p ) : _p ( p ) {}
  
    /// Destructor.
    virtual ~Extended_Poll ( void );

    /// Construct the extended poll points.
    /**
       - Has to be implemented by every NOMAD::Extended_Poll subclass.
       - The extended poll points are the neighbors of \c xk where
         categorical variables have different values.
       - The default implementation of this method uses parameter NEIGHBORS_EXE.
       \param xk Poll center.
    */
    virtual void construct_extended_points ( const NOMAD::Eval_Point & xk );

    /// Set the neighbors executable name for the default implementation.
    /**
       \param   error_str A string containing a possible error message -- \b OUT.
       \return  \c true if no error.
    */
    bool set_neighbors_exe ( std::string & error_str );

    /// Reset.
    void reset ( void );

    /// Poll reset.
    /**
       Before the extended poll is launched.
    */
    void poll_reset ( void );

    /// Access to the poll signatures.
    /**
       \return The set of poll signatures.
    */
    const std::set<NOMAD::Signature_Element> & get_poll_signatures ( void ) const
    {
      return _poll_signatures;
    }

    /// Run the extended poll.
    /**
       \param mads            NOMAD::Mads object invoking this poll -- \b IN/OUT.
       \param nb_ext_poll_pts Number of extended poll points        -- \b OUT.
       \param stop            Stop flag                             -- \b IN/OUT.
       \param stop_reason     Stop reason                           -- \b OUT.
       \param success         Type of success                       -- \b OUT.
       \param new_feas_inc    New feasible incumbent                -- \b IN/OUT.
       \param new_infeas_inc  New infeasible incumbent              -- \b IN/OUT.
    */
    void run ( Mads                     & mads             ,
	       int                      & nb_ext_poll_pts  ,
	       bool                     & stop             ,
	       NOMAD::stop_type         & stop_reason      ,
	       NOMAD::success_type      & success          ,
	       const NOMAD::Eval_Point *& new_feas_inc     ,
	       const NOMAD::Eval_Point *& new_infeas_inc     );
  };
}

#endif
