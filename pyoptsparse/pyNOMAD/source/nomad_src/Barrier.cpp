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
 \file   Barrier.cpp
 \brief  Barrier for constraints handling (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-12
 \see    Barrier.hpp
 */
#include "Barrier.hpp"

/*---------------------------------------------------------*/
/*                    reset the barrier                    */
/*---------------------------------------------------------*/
void NOMAD::Barrier::reset ( void )
{
    _prefilter.clear();
    _filter.clear();
    _h_max           = _p.get_h_max_0();
    _best_feasible   = NULL;
    _ref             = NULL;
    _rho_leaps       = 0;
    _poll_center     = NULL;
    _sec_poll_center = NULL;
    
    if ( _peb_changes > 0 )
        _p.reset_PEB_changes();
    
    _peb_changes      = 0;
    _peb_filter_reset = 0;
    
    _peb_lop.clear();
    _all_inserted.clear();
    
    _one_eval_succ = _success = NOMAD::UNSUCCESSFUL;
}

/*---------------------------------------------------------*/
/*                         display                         */
/*---------------------------------------------------------*/
void NOMAD::Barrier::display ( const Display & out ) const
{
    if ( _eval_type == NOMAD::SGTE )
        out << "surrogate barrier" << std::endl;
    
    if ( _p.get_barrier_type() == NOMAD::EB )
        out << "extreme barrier (EB)" << std::endl;
    
    else {
        
        out << "type                       : "
        << ( (_p.get_barrier_type()==NOMAD::FILTER) ? "filter" : "progressive"  )
        << std::endl
        << "h_norm                     : " << _p.get_h_norm()   << std::endl
        << "h_min                      : " << _p.get_h_min()    << std::endl
        << "h_max                      : " << _h_max            << std::endl;
        if ( _p.get_barrier_type()==NOMAD::PB || _p.get_barrier_type()==NOMAD::PEB_P ) {
            out << "poll center  trigger rho   : " << _p.get_rho()      << std::endl
            << "number of trigger leaps    : " << _rho_leaps        << std::endl;
            if ( _p.get_barrier_type()==NOMAD::PEB_P )
                out << "number of PEB changes      : " << _peb_changes << std::endl
                << "number of PEB filter resets: " << _peb_filter_reset << std::endl;
        }
        if ( out.get_gen_dd() == NOMAD::FULL_DISPLAY )
            out << "number of pre-filter points: " << static_cast<int>(_prefilter.size())
            << std::endl;
        out << NOMAD::open_block ( "list of filter points ("
                                  + NOMAD::itos ( _filter.size() )
                                  + ")" )
        << std::endl;
        std::set<NOMAD::Filter_Point>::const_iterator end = _filter.end() , it;
        for ( it = _filter.begin() ; it != end ; ++it )
            out << *it->get_point() << std::endl;
        out.close_block();
    }
    
#ifdef DEBUG
    out << NOMAD::open_block ( "list of inserted points ("
                              + NOMAD::itos ( _all_inserted.size() )
                              + ")" )
    << std::endl;
    std::list<const NOMAD::Eval_Point *>::const_iterator it2 , end2 = _all_inserted.end();
    for ( it2 = _all_inserted.begin() ; it2 != end2 ; ++it2 )
        out << **it2 << std::endl;
    out.close_block();
#endif
}

/*------------------------------------------------------------*/
/*  barrier update: invoked by Evaluator_Control::eval_lop()  */
/*------------------------------------------------------------*/
void NOMAD::Barrier::update_and_reset_success ( void )
{
    if ( ( _p.get_barrier_type() == NOMAD::PB || _p.get_barrier_type() == NOMAD::PEB_P ) &&
        _success != NOMAD::UNSUCCESSFUL ) {
        
        if ( _success == NOMAD::PARTIAL_SUCCESS ) {
            
            if ( _filter.empty() )
                throw Barrier::Update_Error ( "Barrier.cpp" , __LINE__ ,
                                             "filter empty after a partial success" );
            
            std::set<NOMAD::Filter_Point>::const_iterator it = _filter.end();
            --it;
            
            while ( true ) {
                
                if ( it->get_point()->get_h().value() < _h_max.value() ) {
                    set_h_max ( it->get_point()->get_h() );
                    break;
                }
                
                if ( it == _filter.begin() )
                    throw Barrier::Update_Error ( "Barrier.cpp" , __LINE__ ,
                                                 "could not find a filter point with h < h_max after a partial success" );
                
                --it;
            }
        }
        
        _ref = get_best_infeasible();
        if ( _ref )
            set_h_max ( _ref->get_h() );
    }
    
    // reset success types:
    _one_eval_succ = _success = NOMAD::UNSUCCESSFUL;
}

/*---------------------------------------------------------*/
/*         insertion of an Eval_Point in the barrier       */
/*---------------------------------------------------------*/
void NOMAD::Barrier::insert ( const NOMAD::Eval_Point & x )
{
    // we compare the eval types (_SGTE_ or _TRUTH_) of x and *this:
    if ( x.get_eval_type() != _eval_type )
        throw Barrier::Insert_Error ( "Barrier.cpp" , __LINE__ ,
                                     "insertion of an Eval_Point into the bad Barrier object" );
    
    // basic check:
    if ( !x.is_eval_ok() ) {
        _one_eval_succ = NOMAD::UNSUCCESSFUL;
        return;
    }
    
    // pre-filter: if tag(x) is already in the pre-filter,
    // then return _UNSUCCESSFUL_:
    size_t size_before = _prefilter.size();
    _prefilter.insert ( x.get_tag() );
    if ( _prefilter.size() == size_before ) {
        _one_eval_succ = NOMAD::UNSUCCESSFUL;
        return;
    }
    
    // insertion in _all_inserted:
    _all_inserted.push_back ( &x );
    
    // other checks:
    const NOMAD::Double & h = x.get_h();
    if ( !x.is_EB_ok             () ||
        !x.get_f().is_defined   () ||
        !h.is_defined           () ||
        h.value() > _h_max.value()    ) {
        _one_eval_succ = NOMAD::UNSUCCESSFUL;
        return;
    }
    
    // insert_feasible or insert_infeasible:
    _one_eval_succ = ( x.is_feasible ( _p.get_h_min() ) ) ?
    insert_feasible ( x ) : insert_infeasible(x);
    
    if ( _one_eval_succ > _success )
        _success = _one_eval_succ;
}

/*---------------------------------------------------------*/
/*         update the barrier from another barrier         */
/*                (used by VNS search)                     */
/*---------------------------------------------------------*/
void NOMAD::Barrier::insert ( const Barrier & b )
{
    _one_eval_succ = _success = NOMAD::UNSUCCESSFUL;
    
    NOMAD::Eval_Point * modifiable_x;
    
    std::list<const NOMAD::Eval_Point *>::const_iterator it , end = b._all_inserted.end();
    for ( it = b._all_inserted.begin() ; it != end ; ++it ) {
        
        modifiable_x = &NOMAD::Cache::get_modifiable_point ( **it );
        
        modifiable_x->set_direction          ( NULL                              );
        modifiable_x->set_poll_center_type   ( NOMAD::UNDEFINED_POLL_CENTER_TYPE );
        modifiable_x->set_user_eval_priority ( NOMAD::Double()                   );
        modifiable_x->set_rand_eval_priority ( NOMAD::Double()                   );
        
        insert ( **it );
        
        if ( _one_eval_succ > _success )
            _success = _one_eval_succ;
    }
}

/*---------------------------------------------------------*/
/*          insertion of a feasible point (private)        */
/*---------------------------------------------------------*/
NOMAD::success_type NOMAD::Barrier::insert_feasible ( const NOMAD::Eval_Point & x )
{
    if ( !_best_feasible || ( x.get_f().value() < _best_feasible->get_f().value() ) ) {
        _best_feasible = &x;
        return NOMAD::FULL_SUCCESS;
    }
    return NOMAD::UNSUCCESSFUL;
}

/*---------------------------------------------------------*/
/*                 filter insertion (private)              */
/*---------------------------------------------------------*/
void NOMAD::Barrier::filter_insertion ( const NOMAD::Eval_Point & x , bool & insert )
{
    if ( _filter.empty() ) {
        _filter.insert (&x);
        insert = true;
    }
    else {
        
        insert = false;
        std::set<NOMAD::Filter_Point>::iterator it = _filter.begin();
        while ( it != _filter.end() ) {
            if ( x < *(it->get_point()) ) {
                _filter.erase(it++);
                insert = true;
                continue;
            }
            
            ++it;
        }
        
        if ( !insert ) {
            insert = true;
            std::set<NOMAD::Filter_Point>::iterator end = _filter.end();
            for ( it = _filter.begin() ; it != end ; ++it ) {
                if ( *(it->get_point()) < x ) {
                    insert = false;
                    break;
                }
            }
        }
        
        if ( insert )
            _filter.insert (&x);
    }
}

/*---------------------------------------------------------*/
/*         insertion of an infeasible point (private)      */
/*---------------------------------------------------------*/
NOMAD::success_type NOMAD::Barrier::insert_infeasible ( const NOMAD::Eval_Point & x )
{
    const NOMAD::Eval_Point * old_bi = get_best_infeasible();
    
    // filter insertion:
    // -----------------
    bool insert;
    filter_insertion ( x , insert );
    
    // filter:
    // -------
    if ( _p.get_barrier_type() == NOMAD::FILTER ) {
        const NOMAD::Eval_Point * bi = get_best_infeasible();
        if ( !bi )
            return NOMAD::UNSUCCESSFUL;
        if ( old_bi ) {
            if ( bi->get_h().value() < old_bi->get_h().value() )
                return NOMAD::FULL_SUCCESS;
            if ( insert )
                return NOMAD::PARTIAL_SUCCESS;
            return NOMAD::UNSUCCESSFUL;
        }
        return NOMAD::FULL_SUCCESS;
    }
    
    // progressive barrier:
    // --------------------
    
    // with PEB constraints, we remember all points that we tried to insert:
    if ( _p.get_barrier_type() == NOMAD::PEB_P )
        _peb_lop.push_back ( &x );
    
    // first infeasible successes are considered as partial successes
    // (improving iterations)
    if ( !_ref )
        return NOMAD::PARTIAL_SUCCESS;
    
    double hx = x.get_h().value();
    double fx = x.get_f().value();
    double hr = _ref->get_h().value();
    double fr = _ref->get_f().value();
    
    // failure:
    if ( hx > hr || ( hx == hr && fx >= fr ) )
        return NOMAD::UNSUCCESSFUL;
    
    // partial success:
    if ( fx > fr )
        return NOMAD::PARTIAL_SUCCESS;
    
    // full success:
    return NOMAD::FULL_SUCCESS;
}

/*---------------------------------------------------------*/
/*                    get_best_infeasible()                */
/*---------------------------------------------------------*/
const NOMAD::Eval_Point * NOMAD::Barrier::get_best_infeasible ( void ) const
{
    if ( _filter.empty() || _p.get_barrier_type() == NOMAD::EB )
        return NULL;
    
    if ( _p.get_barrier_type() == NOMAD::FILTER )
        return _filter.begin()->get_point();  // original
    
    return (--_filter.end())->get_point(); // original
}


/*---------------------------------------------------------*/
/*                    get_best_infeasible_min_viol()                */
/*---------------------------------------------------------*/
const NOMAD::Eval_Point * NOMAD::Barrier::get_best_infeasible_min_viol ( void ) const
{
    if ( _filter.empty() || _p.get_barrier_type() == NOMAD::EB )
        return NULL;
    
    if ( _p.get_barrier_type() == NOMAD::FILTER )
        return (--_filter.end())->get_point();
    
    return _filter.begin()->get_point();
}



/*---------------------------------------------------------*/
/*                  poll center selection                  */
/*---------------------------------------------------------*/
void NOMAD::Barrier::select_poll_center ( NOMAD::success_type last_it_success )
{
    const NOMAD::Eval_Point * best_infeasible = get_best_infeasible();
    
    _sec_poll_center = NULL;
    if ( !_best_feasible && !best_infeasible ) {
        _poll_center = NULL;
        return;
    }
    if ( !best_infeasible ) {
        _poll_center = _best_feasible;
        return;
    }
    if ( !_best_feasible ) {
        _poll_center = best_infeasible;
        return;
    }
    
    // filter:
    if ( _p.get_barrier_type() == NOMAD::FILTER ) {
        
        if ( !_poll_center ) {
            _poll_center = _best_feasible;
            return;
        }
        
        // switch:
        if ( last_it_success == NOMAD::UNSUCCESSFUL )
            _poll_center = ( _poll_center==best_infeasible ) ?
            _best_feasible : best_infeasible;
        return;
    }
    
    // progressive barrier:
    if ( _p.get_barrier_type() == NOMAD::PB || _p.get_barrier_type() == NOMAD::PEB_P ) {
        
        const NOMAD::Point * last_poll_center = _poll_center;
        
        if ( best_infeasible->get_f() < (_best_feasible->get_f() - _p.get_rho()) ) {
            _poll_center     = best_infeasible;
            _sec_poll_center = _best_feasible;
        }
        else {
            _poll_center     = _best_feasible;
            _sec_poll_center = best_infeasible;
        }
        if ( _poll_center != last_poll_center )
            ++_rho_leaps;
    }
}

/*---------------------------------------------------------*/
/*             change the value of _h_max (private)        */
/*---------------------------------------------------------*/
void NOMAD::Barrier::set_h_max ( const NOMAD::Double & h_max )
{
    // _h_max update:
    _h_max = h_max;
    
    // we remove all filter points x such that h(x) > h_max:
    if ( !_filter.empty() ) {
        
        if ( _filter.begin()->get_point()->get_h().value() > _h_max.value() ) {
            _filter.clear();
            return;
        }
        
        std::set<NOMAD::Filter_Point>::iterator it = _filter.end();
        do
            --it;
        while ( it != _filter.begin() &&
               it->get_point()->get_h().value() > _h_max.value() );
        ++it;
        _filter.erase ( it , _filter.end() );
    }
}

/*--------------------------------------------------------*/
/*  check the PEB constraints to change eventually their  */
/*  status from _PEB_P_ to _PEB_E_                        */
/*--------------------------------------------------------*/
void NOMAD::Barrier::check_PEB_constraints ( const NOMAD::Eval_Point & x , bool display )
{
    const NOMAD::Double                      & h_min = _p.get_h_min();
    const std::vector<NOMAD::bb_output_type> & bbot  = _p.get_bb_output_type();
    const NOMAD::Point                       & bbo   = x.get_bb_outputs();
    int                                        nb    = static_cast<int>(bbot.size());
    std::list<int>                             ks;
    
    for ( int k = 0 ; k < nb ; ++k ) {
        if ( bbot[k] == NOMAD::PEB_P && bbo[k] <= h_min ) {
            if ( display )
                _p.out() << std::endl
                << "change status of blackbox output " << k
                << " from progressive barrier constraint to extreme barrier constraint"
                << std::endl;
            ++_peb_changes;
            _p.change_PEB_constraint_status (k);
            ks.push_back(k);
        }
    }
    
    // at least one constraint changed status, so we have to update the filter
    // and remove all points that have their h value changed to infinity
    // (it can add new dominant points from the list _peb_lop):
    if ( !ks.empty() ) {
        
        std::list<int>::const_iterator it_k , end_k = ks.end() , begin_k = ks.begin();
        
        // we inspect the filter points if some have to be removed; if so,
        // all filter candidates (stored in _peb_lop) will be re-inserted
        // into the filter:
        bool reset_filter = false;
        std::set<NOMAD::Filter_Point>::const_iterator end = _filter.end() , it;
        
        for ( it = _filter.begin() ; it != end ; ++it ) {
            
            const NOMAD::Point & bbo_cur = it->get_point()->get_bb_outputs();
            for ( it_k = begin_k ; it_k != end_k ; ++it_k )
                if ( bbo_cur[*it_k] > h_min ) {
                    reset_filter = true;
                    break;
                }
            if ( reset_filter )
                break;
        }
        
        if ( reset_filter ) {
            
            if ( display )
                _p.out() << std::endl << "PEB change of status: filter reset" << std::endl;
            
            ++_peb_filter_reset;
            
            _filter.clear();
            bool insert;
            
            std::list<const NOMAD::Eval_Point *>::const_iterator end2 = _peb_lop.end  ();
            std::list<const NOMAD::Eval_Point *>::iterator       it2  = _peb_lop.begin();
            
            while ( it2 != end2 ) {
                
                insert = true;
                const NOMAD::Point & bbo_cur = (*it2)->get_bb_outputs();
                
                for ( it_k = begin_k ; it_k != end_k ; ++it_k )
                    if ( bbo_cur[*it_k] > h_min ) {
                        insert = false;
                        break;
                    }
                
                // if insert==true: this point is potentially a new filter point:
                if ( insert ) {
                    filter_insertion ( **it2 , insert );
                    ++it2;
                }
                
                // if insert==false: it means that the current filter point
                // has to be removed from filter and from _peb_lop, and
                // in addition, its h is put to INF:
                else {
                    NOMAD::Cache::get_modifiable_point ( **it2 ).set_h ( NOMAD::Double() );
                    _peb_lop.erase(it2++);
                }
            }
        }
    }
}
