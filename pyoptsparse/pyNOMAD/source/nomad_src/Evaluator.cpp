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
 \file   Evaluator.cpp
 \brief  Evaluation of blackbox functions (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-14
 \see    Evaluator.hpp
 */
#include "Evaluator.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
bool NOMAD::Evaluator::_force_quit = false;

/*-----------------------------------------------------------------*/
/*                            constructor                          */
/*-----------------------------------------------------------------*/
/*  Parameters::_bb_exe is used to construct _bb_exe and _bb_nbo:  */
/*    Parameters::_bb_exe can have similar blackbox names, while   */
/*    . _bb_exe will have unique blackbox names                    */
/*    . _bb_exe includes the blackbox path                         */
/*-----------------------------------------------------------------*/
NOMAD::Evaluator::Evaluator ( const NOMAD::Parameters & p )
: _p					( p     ) ,
_is_multi_obj		( false ) ,
_is_model_evaluator	( false )
{
    NOMAD::Evaluator::_force_quit = false;
    
    if ( _p.get_bb_exe().empty() )
        return;
    
    // _bbe_exe and _bb_nbo construction:
    std::list<std::string>::const_iterator it = _p.get_bb_exe().begin();
    _bb_exe.push_back(*it);
    _bb_nbo.push_back(1);
    ++it;
    
    std::list<std::string>::const_iterator end = _p.get_bb_exe().end();
    while ( it != end )
    {
        if ( *it != _bb_exe[_bb_exe.size()-1] )
        {
            _bb_exe.push_back(*it);
            _bb_nbo.push_back(1);
        }
        else
            ++_bb_nbo[_bb_exe.size()-1];
        ++it;
    }
    
    // we check that _bb_exe contains unique names and we add the problem path:
    size_t k , l , n = _bb_exe.size() , nm1 = n-1;
    for ( k = 0 ; k < nm1 ; ++k )
    {
        for ( l = k+1 ; l < n ; ++l )
            if ( _bb_exe[k] == _bb_exe[l] )
                throw NOMAD::Exception ( "Evaluator.cpp" , __LINE__ ,
                                        "problem with executable names" );
    }
    
    // construction of _sgte_exe:
    bool        has_sgte_exe = _p.has_sgte_exe();
    std::string err;
    if ( has_sgte_exe )
    {
        for ( k = 0 ; k < n ; ++k )
        {
            
            _sgte_exe.push_back ( _p.get_sgte_exe(_bb_exe[k]) );
            
            if ( _sgte_exe[_sgte_exe.size()-1].empty() )
            {
                err = "blackbox executable \'" + _bb_exe[k] + "\' has no surrogate";
                throw NOMAD::Exception ( "Evaluator.cpp" , __LINE__ , err );
            }
        }
    }
    
    // process blakc-box executables (check and add problem path):
    for ( k = 0 ; k < n ; ++k )
    {
        process_bb_exe_name ( _bb_exe[k] );
        if ( has_sgte_exe )
            process_bb_exe_name ( _sgte_exe[k] );
    }
    
    // blackbox names and indexes display:
#ifdef DEBUG
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank);
    if ( rank == 0 )
    {
#else
        {
#endif
            const NOMAD::Display & out = _p.out();
            if ( !_bb_exe.empty() )
            {
                out << std::endl
                << NOMAD::open_block ( "blackbox executables" );
                for ( k = 0 ; k < n ; ++k )
                {
                    out << NOMAD::open_block ( "bb #" + NOMAD::itos(k) )
                    << _bb_exe[k] << std::endl
                    << "number of outputs=" << _bb_nbo[k] << std::endl
                    << NOMAD::close_block();
                }
                out.close_block();
            }
            if ( !_sgte_exe.empty() )
            {
                out << std::endl
                << NOMAD::open_block ( "surrogate executables" );
                for ( k = 0 ; k < n ; ++k )
                    out << "sgte #" << static_cast<int>(k) << ": "
                    << _sgte_exe[k] << std::endl;
                out.close_block();
            }
        }
#endif
    }
    
    /*----------------------------------------------------------------*/
    /*            process a blackbox executable name (private)        */
    /*----------------------------------------------------------------*/
    void NOMAD::Evaluator::process_bb_exe_name ( std::string & bb_exe ) const
    {
        std::string            err;
        std::list<std::string> bb_exe_words;
        
        NOMAD::get_words ( bb_exe , bb_exe_words );
        
        if ( bb_exe_words.empty() )
        {
            err = "problem with executable \'" + bb_exe + "\'";
            throw NOMAD::Exception ( "Evaluator.cpp" , __LINE__ , err );
        }
        
        std::string problem_dir = _p.get_problem_dir();
        
        // bb_exe is composed of several words (it is a command):
        if ( bb_exe_words.size() > 1 )
        {
            
            bb_exe.clear();
            
            std::list<std::string>::const_iterator it  = bb_exe_words.begin() ,
            end = bb_exe_words.end();
            while (true)
            {
                
                if ( (*it)[0] != '$' )
                {
                    bb_exe += "\"" + problem_dir;
                    bb_exe += *it + "\"";
                }
                else
                    bb_exe += it->substr ( 1 , it->size()-1 );
                
                ++it;
                
                if ( it == end )
                    break;
                
                bb_exe += " ";
            }
        }
        
        // bb_exe is just composed of one name (it is an executable):
        else {
            if ( bb_exe[0] != '$' )
                bb_exe = problem_dir + bb_exe;
            else
                bb_exe = bb_exe.substr ( 1 , bb_exe.size()-1 );
            if ( !NOMAD::check_exe_file ( bb_exe ) )
            {
                err = "\'" + bb_exe + "\' is not a valid executable file";
                throw NOMAD::Exception ( "Evaluator.cpp" , __LINE__ , err );
            }
            if ( bb_exe[0] != '$' )
                bb_exe = "\"" + bb_exe + "\"";
        }
    }
    
    /*-----------------------------------------------------------------------*/
    /*  check the constraints to decide if an evaluation have to be stopped  */
    /*-----------------------------------------------------------------------*/
    /*   . checked when h > h_max or if a 'EB' constraint is violated        */
    /*   . private method                                                    */
    /*-----------------------------------------------------------------------*/
    bool NOMAD::Evaluator::interrupt_evaluations ( const NOMAD::Eval_Point & x     ,
                                                  const NOMAD::Double     & h_max   ) const
    {
        int                                        nbo     = _p.get_bb_nb_outputs();
        const NOMAD::Point                       & bbo     = x.get_bb_outputs();
        const std::vector<NOMAD::bb_output_type> & bbot    = _p.get_bb_output_type();
        NOMAD::Double                              h       = 0.0;
        bool                                       check_h = h_max.is_defined();
        
        for ( int i = 0 ; i < nbo ; ++i ) {
            
            if ( bbo[i].is_defined()                                 &&
                ( bbot[i] == NOMAD::EB || bbot[i] == NOMAD::PEB_E ) &&
                bbo[i] > _p.get_h_min() )
                return true;
            
            if ( check_h && bbo[i].is_defined() &&
                (bbot[i] == NOMAD::FILTER ||
                 bbot[i] == NOMAD::PB     ||
                 bbot[i] == NOMAD::PEB_P     ) ) {
                    
                    if ( bbo[i] > _p.get_h_min() ) {
                        switch ( _p.get_h_norm() ) {
                            case NOMAD::L1:
                                h += bbo[i];
                                break;
                            case NOMAD::L2:
                                h += bbo[i].pow2();
                                break;
                            case NOMAD::LINF:
                                if ( bbo[i] > h )
                                    h = bbo[i];
                                break;
                        }
                        
                        if ( _p.get_h_norm() == NOMAD::L2 ) {
                            if ( h > h_max.pow2() )
                                return true;
                        }
                        else if ( h > h_max )
                            return true;
                    }
                }
        }
        return false;
    }
    
    /*--------------------------------------------------------*/
    /*    compute f(x) from the blackbox outputs of a point   */
    /*      (define a Multi_Obj_Evaluator to treat more than  */
    /*       one objective)                                   */
    /*--------------------------------------------------------*/
    void NOMAD::Evaluator::compute_f ( NOMAD::Eval_Point & x ) const
    {
        if ( x.get_bb_outputs().size() != _p.get_bb_nb_outputs() ) {
            std::ostringstream err;
            err << "Evaluator::compute_f(x): x has a wrong number of blackbox outputs ("
            << x.get_bb_outputs().size() << " != "
            << _p.get_bb_nb_outputs() << ")";
            throw NOMAD::Exception ( "Evaluator.cpp" , __LINE__ , err.str() );
        }
        
        x.set_f ( x.get_bb_outputs()[*(_p.get_index_obj().begin())] );
    }
    
    /*--------------------------------------------------------*/
    /*    compute h(x) from the blackbox outputs of a point   */
    /*    set also the flag 'EB_ok' of the point              */
    /*--------------------------------------------------------*/
    void NOMAD::Evaluator::compute_h ( NOMAD::Eval_Point & x ) const
    {
        if ( x.get_bb_outputs().size() != _p.get_bb_nb_outputs() )
        {
            std::ostringstream err;
            err << "Evaluator::compute_h(x): x has a wrong number of blackbox outputs ("
            << x.get_bb_outputs().size() << " != "
            << _p.get_bb_nb_outputs() << ")";
            throw NOMAD::Exception ( "Evaluator.cpp" , __LINE__ , err.str() );
        }
        
        int                                        nbo  = _p.get_bb_nb_outputs();
        const std::vector<NOMAD::bb_output_type> & bbot = _p.get_bb_output_type();
        const NOMAD::Point                       & bbo  = x.get_bb_outputs();
        NOMAD::Double                              h    = 0.0 , bboi;
        
        x.set_EB_ok ( true );
        
        for ( int i = 0 ; i < nbo ; ++i )
        {
            
            bboi = bbo[i];
            
            if ( bboi.is_defined()                                  &&
                (bbot[i] == NOMAD::EB || bbot[i] == NOMAD::PEB_E ) &&
                bboi > _p.get_h_min() )
            {
                h.clear();
                x.set_h     ( h     );
                x.set_EB_ok ( false );
                return;
            }
            
            if ( bboi.is_defined() &&
                ( bbot[i] == NOMAD::FILTER ||
                 bbot[i] == NOMAD::PB     ||
                 bbot[i] == NOMAD::PEB_P     ) )
            {
                if ( bboi > _p.get_h_min() )
                {
                    switch ( _p.get_h_norm() )
                    {
                        case NOMAD::L1:
                            h += bboi;
                            break;
                        case NOMAD::L2:
                            h += bboi * bboi;
                            break;
                        case NOMAD::LINF:
                            if ( bboi > h )
                                h = bboi;
                            break;
                    }
                }
            }
        }
        
        if ( _p.get_h_norm() == NOMAD::L2 )
            h = h.sqrt();
        
        x.set_h ( h );
    }
    
    /*-------------------------------------------------------------------*/
    /*      . evaluate the black boxes at a given Eval_Point             */
    /*      . the function returns true if the evaluation did not fail   */
    /*      . set count_eval=true to count the evaluation                */
    /*        (unless the output value CNT_EVAL is defined and set to 1  */
    /*         by the blackbox)                                          */
    /*-------------------------------------------------------------------*/
    bool NOMAD::Evaluator::eval_x ( NOMAD::Eval_Point	& x          ,
                                   const NOMAD::Double	& h_max      ,
                                   bool					& count_eval   ) const
    {
        count_eval = false;
        
        if ( _bb_exe.empty() || !x.is_complete() )
            throw NOMAD::Exception ( "Evaluator.cpp" , __LINE__ ,
                                    "Evaluator: no BB_EXE is defined (blackbox executable names)" );
        
        bool sgte = x.get_eval_type() == NOMAD::SGTE;
        if ( sgte && _sgte_exe.empty() )
            throw NOMAD::Exception ( "Evaluator.cpp" , __LINE__ ,
                                    "Evaluator: no SGTE_EXE is defined (surrogate executable names)" );
        
        int         pid     = NOMAD::get_pid();
        int         seed    = _p.get_seed();
        std::string tmp_dir = _p.get_tmp_dir();
        
        std::ostringstream oss;
        oss << "." << seed;
        if ( pid != seed )
            oss << "." << pid;
        oss << "." << x.get_tag() << ".";
        const std::string & sint = oss.str();
        
        // for the parallel version: no need to include the process rank in the names
        // as the point tags are unique for all the processes: each process creates
        // its own points and uses Eval_Point::set_tag()
        
        // blackbox input file writing:
        // ----------------------------
        std::string bb_input_file_name =
        tmp_dir + NOMAD::BLACKBOX_INPUT_FILE_PREFIX
        + sint + NOMAD::BLACKBOX_INPUT_FILE_EXT;
        
        std::string bb_output_file_name =
        tmp_dir + NOMAD::BLACKBOX_OUTPUT_FILE_PREFIX
        + sint + NOMAD::BLACKBOX_OUTPUT_FILE_EXT;
        
        std::ofstream fout ( bb_input_file_name.c_str() );
        if ( fout.fail() )
        {
            std::string err = "could not create file blackbox input file " + bb_input_file_name + ". \n \n #### Please check that write permission are granted for the working directory. #### ";
            throw NOMAD::Exception ( "Evaluator.cpp" , __LINE__ , err );
       }
        
        // include seed:
        if ( _p.get_bb_input_include_seed() )
            fout << seed << " ";
        
        // include tag:
        if ( _p.get_bb_input_include_tag() )
            fout << x.get_tag() << " ";
        
        fout.setf ( std::ios::fixed );
        fout.precision ( NOMAD::DISPLAY_PRECISION_BB );
        x.Point::display ( fout , " " , -1 , -1 );
        fout << std::endl;
        
        fout.close();
        
        if ( fout.fail() )
            return false;
        
        x.set_eval_status ( NOMAD::EVAL_IN_PROGRESS );
        
        std::string   cmd , bb_exe;
        std::ifstream fin;
        bool          failed;
        NOMAD::Double d;
        int           j , nbbok;
        int           ibbo = 0;
        
        // system call to evaluate the blackbox:
        // -------------------------------------
        size_t bn = _bb_exe.size();
        for ( size_t k = 0 ; k < bn ; ++k )
        {
            
            // executable name:
            bb_exe = ( sgte ) ? _sgte_exe[k] : _bb_exe[k];
            
            // system command:
            cmd = bb_exe + " " + bb_input_file_name;
            
            // redirection ? if no, the blackbox has to create
            // the output file 'bb_output_file_name':
            if ( _p.get_bb_redirection() )
                cmd += " > " + bb_output_file_name;
            
#ifdef DEBUG
#ifdef USE_MPI
            int rank;
            MPI_Comm_rank ( MPI_COMM_WORLD, &rank);
            _p.out() << "command(rank=" << rank
            << ") = \'" << cmd << "\'" << std::endl;
#else
            _p.out() << "command=\'" << cmd << "\'" << std::endl;
#endif
#endif
            
            // the evaluation:
            {
                int signal = system ( cmd.c_str() );
                
                // catch the ctrl-c signal:
                if ( signal == SIGINT )
                    raise ( SIGINT );
                
                // other evaluation error:
                failed = ( signal != 0 );
                count_eval = true;
            }
            
            // the evaluation failed (we stop the evaluations):
            if ( failed )
            {
                x.set_eval_status ( NOMAD::EVAL_FAIL );
                break;
            }
            
            // reading of the blackbox output file:
            // ------------------------------------
            else
            {
                
                // bb-output file reading:
                fin.open ( bb_output_file_name.c_str() );
                
                failed          = false;
                bool is_defined = true;
                bool is_inf     = false;
                
                // loop on the number of outputs for this blackbox:
                nbbok = _bb_nbo[k];
                for ( j = 0 ; j < nbbok ; ++j )
                {
                    
                    fin >> d;
                    
                    if ( !d.is_defined() )
                    {
                        is_defined = false;
                        break;
                    }
                    
                    if ( fin.fail() )
                    {
                        failed = true;
                        break;
                    }
                    
                    if ( d.value() >= NOMAD::INF )
                    {
                        is_inf = true;
                        break;
                    }
                    
                    x.set_bb_output ( ibbo++ , d );
                }
                
                fin.close();
                
                // the evaluation failed:
                if ( failed || !is_defined || is_inf )
                {
                    x.set_eval_status ( NOMAD::EVAL_FAIL );
                    break;
                }
                
                // stop the evaluations if h > h_max or if a 'EB' constraint is violated:
                if ( k < _bb_exe.size() - 1 && interrupt_evaluations ( x , h_max ) )
                    break;
            }
        }
        
        if ( x.get_eval_status() == NOMAD::EVAL_IN_PROGRESS )
            x.set_eval_status ( NOMAD::EVAL_OK );
        
        // delete the blackbox input and output files:
        // -------------------------------------------
        remove ( bb_input_file_name.c_str () );
        remove ( bb_output_file_name.c_str() );
        
        // check the CNT_EVAL output:
        // --------------------------
        int index_cnt_eval = _p.get_index_cnt_eval();
        if ( index_cnt_eval >= 0 && x.get_bb_outputs()[index_cnt_eval] == 0.0 )
            count_eval = false;
        
        return x.is_eval_ok();
    }
    
    
    /*-------------------------------------------------------------------*/
    /*      . evaluate the black boxes at a list of given Eval_Points  */
    /*-------------------------------------------------------------------*/
    bool NOMAD::Evaluator::eval_x ( std::list<NOMAD::Eval_Point *>	& list_eval,
                                   const NOMAD::Double				& h_max ,
                                   std::list<bool>					& list_count_eval) const
    {
        
        std::list<NOMAD::Eval_Point *>::iterator it;
        std::list<NOMAD::Eval_Point *>::iterator it_begin=list_eval.begin();
        std::list<NOMAD::Eval_Point *>::iterator it_end=list_eval.end();
        std::list<bool>::iterator it_count=list_count_eval.begin();
        
        if ( list_eval.size() !=list_count_eval.size())
            throw NOMAD::Exception ( "Evaluator.cpp" , __LINE__ ,
                                    "Evaluator: inconsistent size of list" );
        
        if ( _bb_exe.empty())
            throw NOMAD::Exception ( "Evaluator.cpp" , __LINE__ ,
                                    "Evaluator: no BB_EXE is defined (blackbox executable names)" );
        
        bool sgte = ((*it_begin)->get_eval_type() == NOMAD::SGTE);
        if ( sgte && _sgte_exe.empty() )
            throw NOMAD::Exception ( "Evaluator.cpp" , __LINE__ ,
                                    "Evaluator: no SGTE_EXE is defined (surrogate executable names)" );
        
        
        
        int         pid     = NOMAD::get_pid();
        int         seed    = _p.get_seed();
        std::string tmp_dir = _p.get_tmp_dir();
        
        std::ostringstream oss;
        oss << "." << seed;
        if ( pid != seed )
            oss << "." << pid;
        
        for (it=it_begin;it!=it_end;++it)
        {
            if (!(*it)->is_complete() )
                throw NOMAD::Exception ( "Evaluator.cpp" , __LINE__ ,
                                        "Evaluator: points provided for evaluations are incomplete " );
            
        }
        // add the tag of the first point
        oss << "." << (*it_begin)->get_tag();
        
        oss << "." ;
        const std::string & sint = oss.str();
        
        // for the parallel version: no need to include the process rank in the names
        // as the point tags are unique for all the processes: each process creates
        // its own points and uses Eval_Point::set_tag()
        
        // blackbox input file writing:
        // ----------------------------
        std::string bb_input_file_name =
        tmp_dir + NOMAD::BLACKBOX_INPUT_FILE_PREFIX
        + sint + NOMAD::BLACKBOX_INPUT_FILE_EXT;
        
        std::string bb_output_file_name =
        tmp_dir + NOMAD::BLACKBOX_OUTPUT_FILE_PREFIX
        + sint + NOMAD::BLACKBOX_OUTPUT_FILE_EXT;
        
        std::ofstream fout ( bb_input_file_name.c_str() );
        if ( fout.fail() ) {
            std::string err = "could not open file blackbox input file " + bb_input_file_name;
            throw NOMAD::Exception ( "Evaluator.cpp" , __LINE__ , err );
        }
        
        
        for (it=it_begin;it!=it_end;++it)
        {
            // include seed:
            if ( _p.get_bb_input_include_seed() )
                fout << seed << " ";
            
            // include tag:
            if ( _p.get_bb_input_include_tag() )
                fout << (*it)->get_tag() << " ";
            
            fout.setf ( std::ios::fixed );
            fout.precision ( NOMAD::DISPLAY_PRECISION_BB );
            (*it)->Point::display ( fout , " " , -1 , -1 );
            fout << std::endl;
        }
        
        fout.close();
        
        if ( fout.fail() )
            return false;
        
        for (it=it_begin;it!=it_end;++it)
            (*it)->set_eval_status ( NOMAD::EVAL_IN_PROGRESS );
        
        
        std::string   cmd , bb_exe;
        std::ifstream fin;
        bool          failed;
        NOMAD::Double d;
        int           j , nbbok;
        int           ibbo = 0;
        
        // system call to evaluate the blackboxes:
        // -------------------------------------
        size_t bn = _bb_exe.size();
        for ( size_t k = 0 ; k < bn ; ++k )
        {
            
            // executable name:
            bb_exe = ( sgte ) ? _sgte_exe[k] : _bb_exe[k];
            
            // system command:
            cmd = bb_exe + " " + bb_input_file_name;
            
            // redirection ? if no, the blackbox has to create
            // the output file 'bb_output_file_name':
            if ( _p.get_bb_redirection() )
                cmd += " > " + bb_output_file_name;
            
            
            // the evaluation:
            {
                int signal = system ( cmd.c_str() );
                
                // catch the ctrl-c signal:
                if ( signal == SIGINT )
                    raise ( SIGINT );
                
                // other evaluation error:
                failed = ( signal != 0 );
            }
            
            // the evaluation failed (we stop the evaluations):
            if ( failed )
            {
                it_count=list_count_eval.begin();
                for (it=it_begin;it!=it_end;++it,++it_count)
                {
                    (*it)->set_eval_status ( NOMAD::EVAL_FAIL );
                    (*it_count)=true;    //
                }
                break;
            }
            
            // reading of the blackbox output file:
            // ------------------------------------
            else
            {
                
                // bb-output file reading:
                fin.open ( bb_output_file_name.c_str() );
                
                string s;
                bool is_defined,is_inf;
                
                bool list_all_failed_eval=true;
                bool list_all_interrupt=true;
                
                // loop on the points
                it_count=list_count_eval.begin();
                for (it=it_begin;it!=it_end;++it,++it_count)
                {
                    failed		= false;
                    is_defined	= true;
                    is_inf		= false;
                    
                    // loop on the number of outputs for this blackbox:
                    nbbok = _bb_nbo[k];
                    ibbo=0;
                    for ( j = 0 ; j < nbbok ; ++j )
                    {
                        
                        fin >> s;
                        
                        if ( fin.fail() )
                        {
                            failed = true;
                            break;
                        }
                        
                        toupper(s);
                        if (s.compare("REJECT")==0)
                        {
                            *it_count=false;   // Rejected points are not counted
                            (*it)->set_eval_status(NOMAD::EVAL_USER_REJECT);
                            break;
                        }
                        else
                        {
                            d.atof(s);
                            (*it_count)=true;
                        }
                        //
                        
                        if  (s.compare("FAIL")==0)
                        {
                            failed = true;
                            break;
                        }
                        
                        
                        if ( !d.is_defined() )
                        {
                            is_defined = false;
                            break;
                        }
                        
                        
                        if ( d.value() >= NOMAD::INF ) {
                            is_inf = true;
                            break;
                        }
                        
                        (*it)->set_bb_output ( ibbo++ , d );
                    }
                    
                    
                    // the evaluation failed:
                    if ( failed || !is_defined || is_inf )
                    {
                        (*it)->set_eval_status ( NOMAD::EVAL_FAIL );
                        
                    } else
                        list_all_failed_eval=false;
                    
                    
                    // stop the evaluations if h > h_max or if a 'EB' constraint is violated:
                    if ( !( k < _bb_exe.size() - 1 && interrupt_evaluations ( *(*it) , h_max ) && list_all_interrupt ))
                        list_all_interrupt=false;
                }
                
                fin.close();
                
                if (list_all_failed_eval || list_all_interrupt)
                    break;
                
            }
        }
        
        
        // delete the blackbox input and output files:
        // -------------------------------------------
        remove ( bb_input_file_name.c_str () );
        remove ( bb_output_file_name.c_str() );
        
        bool at_least_one_eval_ok=false;
        int index_cnt_eval = _p.get_index_cnt_eval();
        
        
        // update eval status and check that at least one was ok
        it_count=list_count_eval.begin();
        for (it=it_begin;it!=it_end;++it,++it_count)
        {
            if ( (*it)->get_eval_status() == NOMAD::EVAL_IN_PROGRESS )
                (*it)->set_eval_status ( NOMAD::EVAL_OK );
            
            if (!at_least_one_eval_ok && (*it)->is_eval_ok())
                at_least_one_eval_ok=true;
            
            // count_eval from bb_outputs:
            // --------------------------
            if ( index_cnt_eval >= 0 && (*it)->get_bb_outputs()[index_cnt_eval]==0)
                *it_count=false;
        }
        
        return at_least_one_eval_ok;
    }
