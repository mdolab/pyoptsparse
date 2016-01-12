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
  \file   Stats.cpp
  \brief  Algorithm stats (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-22
  \see    Stats.hpp
*/
#include "Stats.hpp"

/*---------------------------------------------------------*/
/*                     affectation operator                */
/*---------------------------------------------------------*/
NOMAD::Stats & NOMAD::Stats::operator = ( const NOMAD::Stats & s )
{
  _eval              = s._eval;
  _sim_bb_eval       = s._sim_bb_eval;
  _sgte_eval         = s._sgte_eval;
  _bb_eval           = s._bb_eval;
  _block_eval         = s._block_eval;	
  _failed_eval       = s._failed_eval;
  _cache_hits        = s._cache_hits;
  _interrupted_eval  = s._interrupted_eval;
  _iterations        = s._iterations;
  _nb_poll_searches  = s._nb_poll_searches;
  _poll_pts          = s._poll_pts;
  _poll_success      = s._poll_success;
  _nb_success_dyn_dir= s._nb_success_dyn_dir;
  _nb_ext_polls      = s._nb_ext_polls;
  _ext_poll_pts      = s._ext_poll_pts;
  _ext_poll_succ     = s._ext_poll_succ;
  _ext_poll_bb_eval  = s._ext_poll_bb_eval;
  _ext_poll_descents = s._ext_poll_descents;
  _nb_spec_searches  = s._nb_spec_searches;
  _spec_pts          = s._spec_pts;
  _spec_success      = s._spec_success;
  _nb_LH_searches    = s._nb_LH_searches;
  _LH_pts            = s._LH_pts;
  _LH_success        = s._LH_success;
  _nb_cache_searches = s._nb_cache_searches;
  _CS_pts            = s._CS_pts;
  _CS_success        = s._CS_success;
  _model_stats       = s._model_stats;
  _nb_VNS_searches   = s._nb_VNS_searches;
  _VNS_pts           = s._VNS_pts;
  _VNS_success       = s._VNS_success;
  _VNS_bb_eval       = s._VNS_bb_eval;
  _VNS_sgte_eval     = s._VNS_sgte_eval;
  _nb_usr_searches   = s._nb_usr_searches;
  _usr_srch_pts      = s._usr_srch_pts;
  _usr_srch_success  = s._usr_srch_success;
  _p1_iterations     = s._p1_iterations;
  _p1_bbe            = s._p1_bbe;
  _mads_runs         = s._mads_runs;
  _clock             = s._clock;
  _stat_sum          = s._stat_sum;
  _stat_avg          = s._stat_avg;
  _cnt_avg           = s._cnt_avg;

#ifdef USE_MPI
  _asynchronous_success = s._asynchronous_success;
  _MPI_data_size        = s._MPI_data_size;
#endif
  
  return *this;
}

/*---------------------------------------------------------*/
/*                      reset the stats                    */
/*---------------------------------------------------------*/
void NOMAD::Stats::reset ( void )
{
  _eval                =
    _sim_bb_eval       =
    _sgte_eval         =
    _bb_eval           =
	_block_eval         =
    _failed_eval       =
    _cache_hits        =
    _interrupted_eval  =
    _iterations        =
    _nb_poll_searches  =
    _poll_pts          =
    _poll_success      =
    _nb_ext_polls      =
    _ext_poll_pts      =
    _ext_poll_succ     =
    _ext_poll_bb_eval  =
    _ext_poll_descents =
    _nb_spec_searches  =
    _spec_pts          =
    _spec_success      =
    _nb_LH_searches    =
    _LH_pts            =
    _LH_success        =
    _nb_cache_searches =
    _CS_pts            =
    _CS_success        =
    _nb_VNS_searches   =
    _VNS_pts           =
    _VNS_success       =
    _VNS_bb_eval       =
    _VNS_sgte_eval     =
    _nb_usr_searches   =
    _usr_srch_pts      =
    _usr_srch_success  =
    _p1_iterations     =
    _p1_bbe            =
	_nb_success_dyn_dir =
    _mads_runs         = 0;

  _model_stats.reset();

  _stat_sum.clear();
  _stat_avg.clear();
  _cnt_avg = 0;

#ifdef USE_MPI
  _asynchronous_success = 0;
  _MPI_data_size        = -1;
#endif

  _clock.reset();
}

/*---------------------------------------------------------*/
/*          update stats from another Stats object         */
/*---------------------------------------------------------*/
/*      for_search==true means that this method has been   */
/*      invoked from a search                              */
/*---------------------------------------------------------*/
void NOMAD::Stats::update ( const NOMAD::Stats & s , bool for_search )
{
  _eval              += s._eval;
  _sim_bb_eval       += s._sim_bb_eval;
  _sgte_eval         += s._sgte_eval;
  _bb_eval           += s._bb_eval;
  _block_eval        += s._block_eval;	
  _failed_eval       += s._failed_eval;
  _cache_hits        += s._cache_hits;
  _interrupted_eval  += s._interrupted_eval;
  _nb_ext_polls      += s._nb_ext_polls;
  _ext_poll_pts      += s._ext_poll_pts;
  _ext_poll_succ     += s._ext_poll_succ;
  _ext_poll_bb_eval  += s._ext_poll_bb_eval;
  _ext_poll_descents += s._ext_poll_descents;
  _nb_LH_searches    += s._nb_LH_searches;
  _LH_pts            += s._LH_pts;
  _LH_success        += s._LH_success;
  _nb_cache_searches += s._nb_cache_searches;
  _CS_pts            += s._CS_pts;
  _CS_success        += s._CS_success;
  _nb_usr_searches   += s._nb_usr_searches;
  _usr_srch_pts      += s._usr_srch_pts;
  _usr_srch_success  += s._usr_srch_success;
  _nb_success_dyn_dir += s._nb_success_dyn_dir;

#ifdef USE_MPI
  _asynchronous_success += s._asynchronous_success;
  _MPI_data_size         = s._MPI_data_size;
#endif

  // _stat_sum and _stat_avg:
  int tmp = _cnt_avg + s._cnt_avg;
  update_stat_sum ( s._stat_sum );
  update_stat_avg ( s._stat_avg );
  _cnt_avg = tmp;

  // specific updates when for_search==false:
  if ( !for_search ) {
    _nb_poll_searches  += s._nb_poll_searches;
    _poll_pts          += s._poll_pts;
    _poll_success      += s._poll_success;
    _nb_spec_searches  += s._nb_spec_searches;
    _spec_pts          += s._spec_pts;
    _spec_success      += s._spec_success;
    _nb_VNS_searches   += s._nb_VNS_searches;
    _VNS_pts           += s._VNS_pts;
    _VNS_success       += s._VNS_success;
    _VNS_bb_eval       += s._VNS_bb_eval;
    _VNS_sgte_eval     += s._VNS_sgte_eval;
    _p1_iterations     += s._p1_iterations;
    _p1_bbe            += s._p1_bbe;
    _iterations        += s._iterations;
  }
}

/*---------------------------------------------------------*/
/*                      update _stat_sum                   */
/*---------------------------------------------------------*/
void NOMAD::Stats::update_stat_sum ( const NOMAD::Double & d )
{
  if ( !d.is_defined() )
    return;
  if ( _stat_sum.is_defined() )
    _stat_sum += d;
  else
    _stat_sum  = d;
}

/*---------------------------------------------------------*/
/*                      update _stat_avg                   */
/*---------------------------------------------------------*/
void NOMAD::Stats::update_stat_avg ( const NOMAD::Double & d )
{
  if ( !d.is_defined() )
    return;
  if ( _stat_avg.is_defined() )
    _stat_avg += d;
  else
    _stat_avg  = d;
  ++_cnt_avg;
}

/*---------------------------------------------------------*/
/*                          display                        */
/*---------------------------------------------------------*/
void NOMAD::Stats::display ( const NOMAD::Display & out ) const
{
  out << "MADS iterations                 : " << _iterations;
  if ( _p1_iterations > 0 )
    out << " (phase one: " << _p1_iterations << ")";
  out << std::endl;
  if ( _sgte_cost > 0 )
    out << "bb evaluations (with sgte cost) : " << get_bb_eval();
  else
    out << "blackbox evaluations            : " << _bb_eval;
  out << std::endl;
  if ( _block_eval > 0)
	  out << "Block of evaluations            : " << _block_eval;
  if ( _p1_bbe > 0 )
    out << " (phase one: " << _p1_bbe << ")";
  out << std::endl;
  if ( _sim_bb_eval != _bb_eval )
    out << "simulated blackbox evaluations  : " << _sim_bb_eval         << std::endl;
  out << "evaluations                     : "   << _eval                << std::endl;
  if ( _sgte_cost > 0 || _sgte_eval > 0 )
    out << "surrogate evaluations           : " << _sgte_eval           << std::endl;
  out << "failed evaluations              : "   << _failed_eval;
  if ( _failed_eval > 0 && _failed_eval==_bb_eval+_sgte_eval )
    out << " (all evaluations failed)";
  out << std::endl;

  out << "interrupted sequences of eval.  : "   << _interrupted_eval    << std::endl;

  if ( _mads_runs > 1 )
    out << "number of MADS runs             : " << _mads_runs           << std::endl;

  out << "cache hits                      : "   << _cache_hits          << std::endl

      << "number of poll searches         : "   << _nb_poll_searches    << std::endl;
  if ( _nb_poll_searches > 0 )
    out << "poll successes                  : " << _poll_success        << std::endl
	  << "poll points                     : " << _poll_pts            << std::endl;
	  
  out << "dyn. direction successes        : " << _nb_success_dyn_dir <<  std::endl;
	  
  if ( _nb_ext_polls > 0 )
    out << "number of extended polls        : " << _nb_ext_polls        << std::endl
	<< "number of ext. poll descents    : " << _ext_poll_descents   << std::endl
	<< "number of ext. poll successes   : " << _ext_poll_succ       << std::endl
	<< "number of ext. poll points      : " << _ext_poll_pts        << std::endl
	<< "number of ext. poll bb eval     : " << _ext_poll_bb_eval    << std::endl;  
  out << "number of speculative searches  : "   << _nb_spec_searches    << std::endl;
  if ( _nb_spec_searches > 0 )
    out << "speculative search successes    : " << _spec_success        << std::endl
	<< "speculative search points       : " << _spec_pts            << std::endl;
  out << "number of user searches         : "   << _nb_usr_searches     << std::endl;
  if ( _nb_usr_searches > 0 )
    out << "user search successes           : " << _usr_srch_success    << std::endl
	<< "user search points              : " << _usr_srch_pts        << std::endl;
  out << "number of LH searches           : "   << _nb_LH_searches      << std::endl;
  if ( _nb_LH_searches > 0 )
    out << "LH search successes             : " << _LH_success          << std::endl
	<< "LH search points                : " << _LH_pts              << std::endl;
  out << "number of cache searches        : "   << _nb_cache_searches   << std::endl;
  if ( _nb_cache_searches > 0 )
    out << "cache search successes          : " << _CS_success          << std::endl
	<< "cache search points             : " << _CS_pts              << std::endl;	
  out << "number of VNS searches          : "   << _nb_VNS_searches     << std::endl;
  if ( _nb_VNS_searches > 0 ) {
    out << "VNS search successes            : " << _VNS_success         << std::endl
	<< "VNS search points               : " << _VNS_pts             << std::endl
	<< "VNS blackbox evaluations        : " << _VNS_bb_eval         << std::endl;
    if ( _VNS_sgte_eval > 0 )
      out << "VNS surrogate evaluations       : " << _VNS_sgte_eval     << std::endl;
  }
  if ( _model_stats.get_nb_models() > 0 ) {
#ifdef DEBUG
    out << NOMAD::open_block ( "model stats" )
	<< _model_stats
	<< NOMAD::close_block();
#else
    out << "number of models built          : "
	<< _model_stats.get_nb_models() << std::endl
	<< "number of model searches        : "
	<< _model_stats.get_MS_nb_searches() << std::endl;
    if ( _model_stats.get_MS_nb_searches() > 0 ) {
      out << "model search successes          : "
	  << _model_stats.get_MS_success() << std::endl
	  << "model search points             : "
	  << _model_stats.get_MS_pts() << std::endl
	  << "model search blackbox eval.     : "
	  << _model_stats.get_MS_bb_eval() << std::endl;
      if ( _model_stats.get_MS_sgte_eval() > 0 )
	out << "model search sgte evaluations   : "
	    << _model_stats.get_MS_sgte_eval() << std::endl;
    }
#endif
  }
  else
    out << "no model has been constructed" << std::endl;
#ifdef USE_MPI
  out << "number of asynchronous successes: " << _asynchronous_success  << std::endl;
  out << "total size of MPI communications: ";
  if ( _MPI_data_size < 0 )
    out << "-";
  else
    out.display_size_of ( _MPI_data_size );
  out << std::endl;
#endif

  out << "wall-clock time                 : ";
  out.display_time ( _clock.get_real_time() );
  out << std::endl;
  
  if ( _stat_sum.is_defined() )
    out << "stat sum                        : " << _stat_sum            << std::endl;
  if ( _stat_avg.is_defined() )
    out << "stat avg                        : " << get_stat_avg()       << std::endl;
  }
