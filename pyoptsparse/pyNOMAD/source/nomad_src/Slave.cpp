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
 \file   Slave.cpp
 \brief  Slave process (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-22
 \see    Slave.hpp
 */
#include "Slave.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
int  NOMAD::Slave::_rank        = -1;
int  NOMAD::Slave::_np          = -1;
int  NOMAD::Slave::_data_sent   =  0;
int  NOMAD::Slave::_data_rcvd   =  0;
bool NOMAD::Slave::_are_running = false;
bool NOMAD::Slave::_stop_ok     = false;

/*----------------------------------------*/
/*        initializations (private)       */
/*----------------------------------------*/
void NOMAD::Slave::init ( void ) const
{
#ifdef USE_MPI
    MPI_Comm_rank ( MPI_COMM_WORLD, &NOMAD::Slave::_rank );
    MPI_Comm_size ( MPI_COMM_WORLD, &NOMAD::Slave::_np   );
#else
    NOMAD::Slave::_rank = 0;
    NOMAD::Slave::_np   = 1;
#endif
    
    // Slave::force_quit() will be called if ctrl-c is pressed:
    if ( !NOMAD::Slave::is_master() ) {
        
        NOMAD::Evaluator::force_quit();
        
        signal ( SIGTERM , NOMAD::Slave::force_quit );
        signal ( SIGINT  , NOMAD::Slave::force_quit );
#ifndef WINDOWS
        signal ( SIGPIPE , NOMAD::Slave::force_quit ); // (ctrl-c during a "| more")
#endif
    }
}

/*----------------------------------------*/
/*          get the process rank          */
/*               (static)                 */
/*----------------------------------------*/
int NOMAD::Slave::get_rank ( void )
{
    if ( NOMAD::Slave::_rank < 0 ) {
#ifdef USE_MPI
        MPI_Comm_rank ( MPI_COMM_WORLD, &NOMAD::Slave::_rank );
#else
        NOMAD::Slave::_rank = 0;
#endif
    }
    return NOMAD::Slave::_rank;
}

/*----------------------------------------*/
/*       get the number of processes      */
/*               (static)                 */
/*----------------------------------------*/
int NOMAD::Slave::get_nb_processes ( void )
{
    if ( NOMAD::Slave::_np < 0 ) {
#ifdef USE_MPI
        MPI_Comm_size ( MPI_COMM_WORLD, &NOMAD::Slave::_np );
#else
        NOMAD::Slave::_np = 1;
#endif
    }
    return NOMAD::Slave::_np;
}

/*----------------------*/
/*  run the slave code  */
/*----------------------*/
void NOMAD::Slave::run ( void ) const
{
#ifdef USE_MPI
    
    MPI_Request         req;
    char                signal     = 0;
    NOMAD::Eval_Point * x          = NULL;
    bool                count_eval = false;
    
    while ( true ) {
        
        // receive signal from master:
        // ---------------------------
        NOMAD::Slave::receive_data ( &signal , 1 , MPI_CHAR , 0 , &req );
        
        // slave is ready or not initialized:
        NOMAD::Slave::send_data ( &NOMAD::READY_SIGNAL , 1 , MPI_CHAR , 0 , false );
        
        NOMAD::Slave::wait_request ( req );
        
        // EVAL signal:
        // ------------
        if ( signal == NOMAD::EVAL_SIGNAL ) {
            
            // receive and evaluate the point:
            x = eval_point ( count_eval );
            
        }
        
        // RESULT signal:
        // --------------
        else if ( signal == NOMAD::RESULT_SIGNAL )
        {
            
            // send the evaluation result to the master:
            send_eval_result ( x , count_eval );
            
            delete x;
            x = NULL;
        }
        
        // STOP signal:
        // ------------
        else if ( signal == NOMAD::STOP_SIGNAL )
            break;
        
        // WAIT signal:
        // ------------
        // else if ( signal == NOMAD::WAIT_SIGNAL ) {
        // }
    }
    
    if ( x )
        delete x;
    
#endif
}

/*-----------------------------------------*/
/*        initialize all the slaves        */
/*                (static)                 */
/*-----------------------------------------*/
void NOMAD::Slave::init_slaves ( const NOMAD::Display & out )
{
#ifdef USE_MPI
    
    if ( !NOMAD::Slave::is_master() || NOMAD::Slave::_are_running )
        return;
    
    NOMAD::dd_type display_degree = out.get_gen_dd();
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out << std::endl << NOMAD::open_block ( "initializing slaves" );
    
    MPI_Status     status;
    MPI_Request ** req            = new MPI_Request * [ NOMAD::Slave::_np ];
    int            nb_initialized = 0;
    int            nb_slaves      = NOMAD::Slave::_np - 1;
    int            source;
    char           signal;
    NOMAD::Clock   clk;
    
    // 1. launch requests:
    for ( source = 1 ; source < NOMAD::Slave::_np ; ++source ) {
        req[source] = new MPI_Request;
        NOMAD::Slave::receive_data ( &signal , 1 , MPI_CHAR , source ,  req[source] );
        if ( display_degree == NOMAD::FULL_DISPLAY )
            out << "." << std::endl;
    }
    
    // 2. test requests (with a maximal delay of MAX_REQ_WAIT):
    int cnt = 0 , flag;
    while ( nb_initialized < nb_slaves && clk.get_real_time() < NOMAD::MAX_REQ_WAIT )
    {
        
        for ( source = 1 ; source < NOMAD::Slave::_np ; ++source )
        {
            
            if ( req[source] )
            {
                
                MPI_Test ( req[source] , &flag , &status );
                
                if ( flag )
                {
                    
                    MPI_Wait ( req[source] , &status );
                    
                    // send the WAIT signal:
                    NOMAD::Slave::send_data ( &NOMAD::WAIT_SIGNAL , 1 , MPI_CHAR , source , true );
                    
                    delete req[source];
                    req[source] = NULL;
                    ++nb_initialized;
                }
            }
        }
        // a constant is used in order to display only a few '.' :
        if ( display_degree == NOMAD::FULL_DISPLAY && cnt%1000000==0 )
            out << "." << std::endl;
        
        ++cnt;
    }
    
    // 3. delete requests:
    std::list<int> err_list;
    for ( source = 1 ; source < NOMAD::Slave::_np ; ++source )
    {
        if ( req[source] )
        {
            err_list.push_back ( source );
            MPI_Cancel ( req[source] );
            delete req[source];
        }
    }
    delete [] req;
    
    NOMAD::Slave::_are_running = true;
    NOMAD::Slave::_stop_ok     = false;
    
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out << NOMAD::close_block() << std::endl;
    
    if ( !err_list.empty() )
    {
        
        std::ostringstream oss;
        oss << "could not initialize slave";
        if ( err_list.size() > 1 )
        {
            oss << "s";
            std::list<int>::const_iterator it , end = err_list.end();
            for ( it = err_list.begin() ; it != end ; ++it )
                oss << " #" << *it;
        }
        else
            oss << " #" << *err_list.begin();
        
        throw NOMAD::Exception ( "Slave.cpp" , __LINE__ , oss.str() );
    }
    
#endif
}

/*-----------------------------------------*/
/*             stop the slaves             */
/*                (static)                 */
/*-----------------------------------------*/
void NOMAD::Slave::stop_slaves ( const NOMAD::Display & out )
{
#ifdef USE_MPI
    
    if ( !NOMAD::Slave::is_master() || NOMAD::Slave::_stop_ok )
        return;
    
    NOMAD::dd_type display_degree = out.get_gen_dd();
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out << std::endl << NOMAD::open_block ( "stopping slaves" );
    
    int  nb_stopped = 0;
    int  nb_slaves  = NOMAD::Slave::_np - 1;
    int  source;
    char signal;
    
    NOMAD::Clock clk;
    
    MPI_Status  status;
    MPI_Request ** req = new MPI_Request * [ NOMAD::Slave::_np ];
    
    // 1. launch requests:
    for ( source = 1 ; source < NOMAD::Slave::_np ; ++source ) {
        req[source] = new MPI_Request;
        NOMAD::Slave::receive_data ( &signal , 1 , MPI_CHAR , source ,  req[source] );
        if ( display_degree == NOMAD::FULL_DISPLAY )
            out << "." << std::endl;
    }
    
    // 2. test requests (with a maximal delay of MAX_REQ_WAIT):
    int cnt = 0 , flag;
    while ( nb_stopped < nb_slaves && clk.get_real_time() < NOMAD::MAX_REQ_WAIT ) {
        
        for ( source = 1 ; source < NOMAD::Slave::_np ; ++source ) {
            
            if ( req[source] ) {
                
                MPI_Test ( req[source] , &flag , &status );
                
                if ( flag ) {
                    
                    MPI_Wait ( req[source] , &status );
                    
                    // send the STOP signal:
                    NOMAD::Slave::send_data ( &NOMAD::STOP_SIGNAL , 1 , MPI_CHAR , source , true );
                    
                    delete req[source];
                    req[source] = NULL;
                    ++nb_stopped;
                }
            }
        }
        // a constant is used in order to display only a few '.' :
        if ( display_degree == NOMAD::FULL_DISPLAY && cnt%1000000==0 )
            out << "." << std::endl;
        ++cnt;
    }
    
    NOMAD::Slave::_are_running = false;
    NOMAD::Slave::_stop_ok     = true;
    
    // 3. delete requests:
    for ( source = 1 ; source < NOMAD::Slave::_np ; ++source ) {
        if ( req[source] ) {
            MPI_Cancel ( req[source] );
            delete req[source];
            NOMAD::Slave::_stop_ok = false;
        }
    }
    delete [] req;
    
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out << NOMAD::close_block() << std::endl;
    
#endif
}

#ifdef USE_MPI

/*------------------------------------------------------*/
/*               receive data (static, private)         */
/*------------------------------------------------------*/
int NOMAD::Slave::receive_data ( void        * buf      ,
                                int           count    ,
                                MPI_Datatype  datatype ,
                                int           source   ,  // may be MPI_ANY_SOURCE
                                MPI_Request * req        )
{
    int tag = ( NOMAD::Slave::is_master() ) ? source : NOMAD::Slave::get_rank();
    
    // immediate receive:
    if ( req ) {
        if ( source == MPI_ANY_SOURCE )
            throw NOMAD::Exception ( "Slave.cpp" , __LINE__ ,
                                    "Slave::receive_data(): immediate receive with no source" );
        MPI_Irecv ( buf , count , datatype , source , tag , MPI_COMM_WORLD , req );
    }
    
    // normal receive:
    else {
        MPI_Status status;
        if ( source == MPI_ANY_SOURCE )
            tag = MPI_ANY_TAG;
        MPI_Recv ( buf , count , datatype , source , tag , MPI_COMM_WORLD , &status );
        source = status.MPI_SOURCE;
    }
    
    // stats:
    int size;
    MPI_Type_size ( datatype , &size );
    NOMAD::Slave::_data_rcvd += count * size;
    
    return source;
}

/*------------------------------------------------------*/
/*              send data (static, private)             */
/*------------------------------------------------------*/
void NOMAD::Slave::send_data ( const void  * buf        ,
                              int           count      ,
                              MPI_Datatype  datatype   ,
                              int           dest       ,
                              bool          ready_send   )
{
    int tag = ( NOMAD::Slave::is_master() ) ? dest : NOMAD::Slave::get_rank();
    
    // ready send:
    if ( ready_send )
        MPI_Rsend ( const_cast<void*>(buf) , count , datatype ,
                   dest , tag , MPI_COMM_WORLD );
    
    // normal send:
    else
        MPI_Send ( const_cast<void*>(buf) , count , datatype ,
                  dest , tag , MPI_COMM_WORLD );
    
    // stats:
    int size;
    MPI_Type_size ( datatype , &size );
    NOMAD::Slave::_data_sent += count * size;
}

/*------------------------------------------------------*/
/*  receive and evaluate an Eval_Point from the master  */
/*  (private)                                           */
/*------------------------------------------------------*/
NOMAD::Eval_Point * NOMAD::Slave::eval_point ( bool & count_eval ) const
{
    // 1. receive the point:
    int         itab[3];
    MPI_Request req;
    NOMAD::Slave::receive_data ( itab                 , 3 , MPI_INT  , 0 , &req  );
    NOMAD::Slave::send_data    ( &NOMAD::READY_SIGNAL , 1 , MPI_CHAR , 0 , false );
    NOMAD::Slave::wait_request ( req );
    
    int      n    = itab[0];
    double * dtab = new double[n+1];
    
    NOMAD::Slave::receive_data ( dtab                 , n+1 , MPI_DOUBLE , 0 , &req  );
    NOMAD::Slave::send_data    ( &NOMAD::READY_SIGNAL , 1   , MPI_CHAR   , 0 , false );
    NOMAD::Slave::wait_request ( req );
    
    // 2. create the Eval_Point:
    int bb_nb_outputs=_p->get_bb_nb_outputs();
    NOMAD::Eval_Point * x = new NOMAD::Eval_Point ( n , bb_nb_outputs );
    for ( int i = 0 ; i < n ; ++i )
        (*x)[i] = dtab[i];
    NOMAD::Double h_max = dtab[n];
    
    x->set_tag       ( itab[2]                            );
    x->set_eval_type ( ( itab[1] > 0 ) ? NOMAD::SGTE : NOMAD::TRUTH );
    
    delete [] dtab;
    
    // 3. evaluate the point:
    bool eval_ok;
    try {
        eval_ok = _ev->eval_x ( *x , h_max , count_eval );
    }
    catch ( ... ) {
        eval_ok = false;
    }
    
    x->set_eval_status ( ( eval_ok ) ? NOMAD::EVAL_OK : NOMAD::EVAL_FAIL );
    
    return x;
}

/*-----------------------------------------------------*/
/*  send an evaluation result to the master (private)  */
/*-----------------------------------------------------*/
void NOMAD::Slave::send_eval_result ( const NOMAD::Eval_Point * x          ,
                                     bool                      count_eval   ) const
{
    // receive a signal from the master:
    char signal;
    NOMAD::Slave::receive_data ( &signal , 1 , MPI_CHAR , 0 , NULL );
    
    // send the evaluation result:
    int                  m    = _p->get_bb_nb_outputs();
    int                  s    = 2*m+2;
    double             * dtab = new double [s];
    const NOMAD::Point & bbo  = x->get_bb_outputs();
    
    // bb_outputs (m values):
    for ( int i = 0 ; i < m ; ++i ) {
        if ( bbo[i].is_defined() ) {
            dtab[i  ] = bbo[i].value();
            dtab[i+m] = 1.0;
        }
        else {
            dtab[i  ] = NOMAD::INF;
            dtab[i+m] = -1.0;
        }
    }
    
    // evaluation status:
    dtab[2*m] = ( x->get_eval_status() == NOMAD::EVAL_OK ) ? 1.0 : -1.0;
    
    // count_eval:
    dtab[s-1] = ( count_eval ) ? 1.0 : -1.0;
    
    // send the array:
    NOMAD::Slave::send_data ( dtab , s , MPI_DOUBLE , 0 , true );
    
    delete [] dtab;
}

/*---------------------------------------------*/
/*  receive an evaluation result from a slave  */
/*---------------------------------------------*/
void NOMAD::Slave::receive_eval_result ( int                 slave_rank ,
                                        NOMAD::Eval_Point * x          ,
                                        bool              & eval_ok    ,
                                        bool              & count_eval    ) const
{
    // send the RESULT signal to the slave:
    NOMAD::Slave::send_data ( &NOMAD::RESULT_SIGNAL , 1 , MPI_CHAR , slave_rank , true );
    
    // receive the evaluation result as a double array:
    int      m    =  _p->get_bb_nb_outputs();
    int      s    = 2*m+2;
    double * dtab = new double [s];
    
    MPI_Request req;
    NOMAD::Slave::receive_data ( dtab              , s , MPI_DOUBLE , slave_rank , &req  );
    NOMAD::Slave::send_data ( &NOMAD::READY_SIGNAL , 1 , MPI_CHAR   , slave_rank , false );
    NOMAD::Slave::wait_request ( req );
    
    // interpret the array:
    for ( int i = 0 ; i < m ; ++i )
        x->set_bb_output ( i , ( dtab[i+m] > 0.0 ) ? dtab[i] : NOMAD::Double() );
    
    eval_ok = ( dtab[2*m] > 0.0 );
    
    x->set_eval_status ( eval_ok ? NOMAD::EVAL_OK : NOMAD::EVAL_FAIL );
    
    count_eval = ( dtab[s-1] > 0.0 );
    
    delete [] dtab;
}

/*-----------------------------------------------------*/
/*            send an Eval_Point to a slave            */
/*-----------------------------------------------------*/
void NOMAD::Slave::send_eval_point ( const NOMAD::Eval_Point * x          ,
                                    int                       slave_rank ,
                                    const NOMAD::Double     & h_max        ) const
{
    char signal;
    int  itab[3];
    int  n = x->size();
    
    // n:
    itab[0] = n;
    
    // evaluation type (+1: sgte eval; -1: true eval):
    itab[1] = ( x->get_eval_type() == NOMAD::SGTE ) ? 1 : -1;
    
    // tag of the point:
    itab[2] = x->get_tag();
    
    // point coordinates:
    double * dtab = new double[n+1];
    for ( int i = 0 ; i < n ; ++i )
        dtab[i] = (*x)[i].value();
    dtab[n] = h_max.value();
    
    // wait for the slave signal:
    NOMAD::Slave::receive_data ( &signal , 1 , MPI_CHAR , slave_rank , NULL );
    
    // send n and evaluation type:
    NOMAD::Slave::send_data ( itab , 3 , MPI_INT , slave_rank , true );
    
    // wait for the slave signal:
    NOMAD::Slave::receive_data ( &signal , 1 , MPI_CHAR , slave_rank , NULL );
    
    // send the point coordinates:
    NOMAD::Slave::send_data ( dtab , n+1 , MPI_DOUBLE , slave_rank , true );
    
    delete [] dtab;
}

/*-----------------------------------------*/
/*          wait for a MPI request         */
/*                (static)                 */
/*-----------------------------------------*/
void NOMAD::Slave::wait_request ( MPI_Request & req )
{
    MPI_Status status;
    MPI_Wait ( &req , &status );
}

#endif
