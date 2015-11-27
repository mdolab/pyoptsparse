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
  \file   Slave.hpp
  \brief  Slave process (headers)
  \author Sebastien Le Digabel
  \date   2010-04-22
  \see    Slave.cpp
*/
#ifndef __SLAVE__
#define __SLAVE__

#include "Evaluator.hpp"

namespace NOMAD {

  /// Slave process for the parallel version.
  class Slave : private NOMAD::Uncopyable {

  private:

    static int  _rank;        /// Process rank.
    static int  _np;          /// Number of processes.

    static int  _data_sent;   /// Stats on the sent data.
    static int  _data_rcvd;   /// Stats on the received data.

    static bool _are_running; ///< \c true if the slaves are running.
    static bool _stop_ok;     ///< \c true if the slaves stopped without error.

    const NOMAD::Parameters * _p;  ///< Parameters (may be NULL).
    NOMAD::Evaluator        * _ev; ///< Evaluator  (may be NULL).

    /// Initializations.
    void init ( void ) const;

    /// Force quit.
    /**
       - Called by pressing Ctrl-C.
       - Does nothing: the slave will be stopped by the master.
       \param signalValue Signal value -- \b IN.
    */
    static void force_quit ( int signalValue ) {}

  public:

    /// Constructor.
    /**
       \param p  Parameters                -- \b IN.
       \param ev A pointer to an evaluator -- \b IN.
    */
    Slave ( const NOMAD::Parameters & p , NOMAD::Evaluator * ev )
      { _p=&p; _ev=ev; init(); }

    /// Destructor.
    virtual ~Slave ( void ) {}

    /// Run the slave code.
    void run ( void ) const;

    /// Access to the stat \c _data_sent.
    /**
       \return The stat \c _data_sent.
    */
    static int get_data_sent ( void ) { return Slave::_data_sent; }

    /// Access to the stat \c _data_rcvd.
    /**
       \return The stat \c _data_rcvd.
    */
    static int get_data_rcvd ( void ) { return Slave::_data_rcvd; }
    
    /// Access to the process rank.
    /**
       \return The process rank.
    */
    static int get_rank ( void );

    /// Access to the number of processes.
    /**
       \return The number of processes.
    */
    static int get_nb_processes ( void );

    /// Check if the slaves are running.
    /**
       \return A boolean equal to \c true if the slaves are running.
    */
    static bool are_running ( void ) { return _are_running; }
    
    /// Check if the current slave is the master.
    /**
       \return A boolean equal to \c true if the current slave is the master.
    */
    static bool is_master ( void ) { return ( Slave::get_rank() == 0 ); }

    /// Initialize all the slaves.
    /**
       \param out Display -- \b IN.
    */
    static void init_slaves ( const NOMAD::Display & out );

    /// Stop all the slaves.
    /**
       \param out Display -- \b IN.
    */
    static void stop_slaves ( const NOMAD::Display & out );

#ifdef USE_MPI

  private:
    
    /// Receive and evaluate a point from the master.
    /**
       \param count_eval Flag indicating if the evaluation has to be counted
                         or not -- \b OUT.
       \return A pointer to the point.
    */
    NOMAD::Eval_Point * eval_point ( bool & count_eval ) const;

    /// Send an evaluation result to the master.
    /**
       \param x          The evaluation point -- \b IN.
       \param count_eval Flag indicating if the evaluation has to be counted
                         or not -- \b IN.
    */
    void send_eval_result ( const NOMAD::Eval_Point * x , bool count_eval ) const;
    
    /// Send data.
    /**
       \param buf        Data to send  -- \b IN.
       \param count      Data quantity -- \b IN.
       \param datatype   Data type     -- \b IN.
       \param dest       Destination   -- \b IN.
       \param ready_send Flag equal to \c true if \c Rsend is used instead of \c Send
                                       -- \b IN.
    */
    static void send_data ( const void  * buf        ,
			    int           count      ,
			    MPI_Datatype  datatype   ,
			    int           dest       ,
			    bool          ready_send   );
    
    /// Receive data.
    /**
       \param buf      Data to receive                           -- \b OUT.
       \param count    Data quantity                             -- \b IN.
       \param datatype Data type                                 -- \b IN.
       \param source   Source (may be \c MPI_ANY_SOURCE)         -- \b IN.
       \param req      Pointer to a MPI request (may be \c NULL) -- \b IN/OUT.
       \return         The source.
    */
    static int receive_data ( void        * buf      ,
			      int           count    ,
			      MPI_Datatype  datatype ,
			      int           source   ,
			      MPI_Request * req        );

    /// Wait for a MPI request.
    /**
       \param req The request -- \b IN/OUT.
    */
    static void wait_request ( MPI_Request & req );
    
  public:

    /// Send an evaluation point to a slave.
    /**
       \param x          A pointer to the evaluation point to send -- \b IN.
       \param slave_rank Destination                               -- \b IN.
       \param h_max      Maximal feasibility value \c h_max        -- \b IN.
    */
    void send_eval_point ( const NOMAD::Eval_Point * x          ,
			   int                       slave_rank ,
			   const NOMAD::Double     & h_max        ) const;
    
    /// Receive an evaluation result from a slave.
    /**
       \param slave_rank Source                                      -- \b IN.
       \param x          A pointer to the evaluation point received  -- \b OUT.
       \param eval_ok    Flag indicating if the evaluation succeeded -- \b OUT.
       \param count_eval Flag indicating if the evaluation has to be counted
                         or not -- \b OUT.
    */
    void receive_eval_result ( int                 slave_rank ,
			       NOMAD::Eval_Point * x          ,
			       bool              & eval_ok    ,
			       bool              & count_eval    ) const;
    
    /// Receive a signal from a slave.
    /**
       \param signal The signal -- \b OUT.
       \return The source.
    */
    static int receive_signal ( char & signal )
    {
      return Slave::receive_data ( &signal , 1 , MPI_CHAR , MPI_ANY_SOURCE , NULL );
    }

    /// Send a signal to a slave.
    /**
       \param signal     The signal  -- \b IN.
       \param slave_rank Destination -- \b IN.
    */
    static void send_signal ( char signal , int slave_rank )
    {
      Slave::send_data ( &signal , 1 , MPI_CHAR , slave_rank , true );
    }

#endif
  };
}
#endif
