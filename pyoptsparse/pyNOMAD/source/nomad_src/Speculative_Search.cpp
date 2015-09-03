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
 \file   Speculative_Search.cpp
 \brief  Speculative search (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-12
 \see    Speculative.hpp
 */
#include "Speculative_Search.hpp"

/*-------------------------------------------------------------*/
/*                     MADS speculative search                 */
/*-------------------------------------------------------------*/
/*  x_k = x_{k-1} + \Delta^m_{k-1} d                           */
/*  s_k = x_{k-1} +        \Delta^m_k d  if l_{k-1} > 0        */
/*    or  x_{k-1} + \tau * \Delta^m_k d  otherwise             */
/*-------------------------------------------------------------*/
/*  the directions that we use already contain \Delta^m:       */
/*     dir = \Delta^m_{k-1} d                                  */
/*-------------------------------------------------------------*/
/*  the test ' if (new_feas_inc || new_infeas_inc ) ' is true  */
/*  and has already been made in Mads.cpp                      */
/*-------------------------------------------------------------*/
void NOMAD::Speculative_Search::search ( NOMAD::Mads              & mads           ,
                                        int                      & nb_search_pts  ,
                                        bool                     & stop           ,
                                        NOMAD::stop_type         & stop_reason    ,
                                        NOMAD::success_type      & success        ,
                                        bool                     & count_search   ,
                                        const NOMAD::Eval_Point *& new_feas_inc   ,
                                        const NOMAD::Eval_Point *& new_infeas_inc   )
{
    // new_feas_inc and new_infeas_inc are used as inputs,
    // so do not initialize them to NULL here
    
    nb_search_pts = 0;
    success       = NOMAD::UNSUCCESSFUL;
    count_search  = !stop;
    
    if ( stop )
        return;
    
    const NOMAD::Display    & out = _p.out();
    NOMAD::dd_type display_degree = out.get_search_dd();
    if ( display_degree == NOMAD::FULL_DISPLAY ) {
        std::ostringstream oss;
        oss << NOMAD::SPEC_SEARCH;
        out << std::endl << NOMAD::open_block ( oss.str() ) << std::endl;
    }
    
    //int                       lkm1;  // l_{k-1}
    // int                       lk;    // l_k
    int                       n;
    NOMAD::Signature        * signature;
    NOMAD::Point              delta_m_k;
    NOMAD::Point              delta_m_km1;
    NOMAD::Point              factor;
    NOMAD::Point              xkm1;
    NOMAD::Eval_Point       * sk;
    const NOMAD::Eval_Point * x[2];
    x[0] = new_feas_inc;
    x[1] = new_infeas_inc;
    
    // Evaluator_Control:
    NOMAD::Evaluator_Control & ev_control = mads.get_evaluator_control();
    
    for ( int i = 0 ; i < 2 ; ++i )
    {
        if ( x[i] &&  x[i]->get_signature() )
        {
            
            const NOMAD::Direction * dir = x[i]->get_direction();
            if ( dir && ( dir->is_mads() || dir->get_type()==NOMAD::MODEL_SEARCH_DIR ) ) {
                
                // get the x_k's signature:
                signature = x[i]->get_signature();
                if ( !signature )
                    throw NOMAD::Exception ( "Speculative_Search.cpp" , __LINE__ ,
                                            "Speculative_Search::search(): could not get the signature" );
                
                xkm1 = *x[i] - *dir;
                
                factor = signature->get_mesh()->get_mesh_ratio_if_success();
                n=signature->get_n();
                for ( int k = 0 ; k < n ; ++k )
                {
                    if ( factor[k].is_defined() )
                    {
                        if ( factor[k] == 1 )
                        {
                            // factor determined based on mesh type ( default is 4 = mesh_update_basis = 2* poll_update_basis )
                            if ( dynamic_cast<NOMAD::XMesh *> (signature->get_mesh()) )
                                factor[k] = signature->get_mesh()->get_update_basis()*2.0;
                            else
                                factor[k] = signature->get_mesh()->get_update_basis();
                            
                        }
                        else if ( factor[k] == 0 )
                        {
                            if ( display_degree == NOMAD::FULL_DISPLAY )
                                out << "could not compute " << _type << " point: stop" << std::endl
                                << NOMAD::close_block ( "end of speculative search" );
                            stop        = true;
                            stop_reason = NOMAD::MESH_PREC_REACHED;
                            return;
                        }
                        
                    }
                    else
                        factor[k]=0;
                    
                }
                NOMAD::Point mesh_indices_k( x[i]->get_signature()->get_mesh()->get_mesh_indices() );
                signature->get_mesh()->update( NOMAD::FULL_SUCCESS, mesh_indices_k , dir );
                
                // speculative search point:
                NOMAD::Direction new_dir ( n , 0.0 , dir->get_type() );
                new_dir.Point::operator = ( factor * *dir );
                
                sk = new NOMAD::Eval_Point;
                sk->set ( n , _p.get_bb_nb_outputs() );
                sk->set_signature  ( signature );
                sk->set_direction  ( &new_dir );
                
                sk->Point::operator = ( xkm1 + new_dir );
                
                if ( display_degree == NOMAD::FULL_DISPLAY ) {
                    out << "trial point #" << sk->get_tag()
                    << ": ( ";
                    sk->Point::display ( out ," " , 2 , NOMAD::Point::get_display_limit() );
                    out << " )" << std::endl;
                }
                
                // add the new point to the list of search trial points:
                ev_control.add_eval_point ( sk                      ,
                                           display_degree          ,
                                           _p.get_snap_to_bounds() ,
                                           NOMAD::Double()         ,
                                           NOMAD::Double()         ,
                                           NOMAD::Double()         ,
                                           NOMAD::Double()           );
            }
        }
    }
    
    nb_search_pts = ev_control.get_nb_eval_points();
    
    // eval_list_of_points:s
    // --------------------
    new_feas_inc = new_infeas_inc = NULL;
    
    ev_control.eval_list_of_points ( _type                   ,
                                    mads.get_true_barrier() ,
                                    mads.get_sgte_barrier() ,
                                    mads.get_pareto_front() ,
                                    stop                    ,
                                    stop_reason             ,
                                    new_feas_inc            ,
                                    new_infeas_inc          ,
                                    success                   );
    
    if ( display_degree == NOMAD::FULL_DISPLAY ) {
        std::ostringstream oss;
        oss << "end of speculative search (" << success << ")";
        out << NOMAD::close_block ( oss.str() ) << std::endl;
    }
}
