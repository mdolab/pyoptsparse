/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search - version 3.7.1      */
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
 \file   nomad.cpp
 \brief  NOMAD main file
 \author Sebastien Le Digabel
 \date   2010-04-12
 */
#include "nomad.hpp"

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv )
{
    
    // display:
    NOMAD::Display out ( std::cout );
    out.precision ( NOMAD::DISPLAY_PRECISION_STD );
    
    
    std::string error;
    {
        // NOMAD initializations:
        NOMAD::begin ( argc , argv );
        
        // usage:
        if ( argc < 2 )
        {
            NOMAD::display_usage ( argv[0],std::cerr );
            NOMAD::end();
            return EXIT_FAILURE;
        }
        
        // parameters file:
        std::string param_file_name = argv[1];
        std::string opt             = param_file_name;
        NOMAD::toupper ( opt );
        
        // display version if option '-v' has been specified:
        if ( opt == "-U" )
        {
            NOMAD::display_usage ( argv[0], out );
            NOMAD::end();
            return EXIT_SUCCESS;
        }
        
        
        // display version if option '-v' has been specified:
        if ( opt == "-V" || opt =="-VERSION")
        {
            NOMAD::display_version ( out );
            NOMAD::end();
            return EXIT_SUCCESS;
        }
        
        // display info if option '-i' has been specified:
        if ( opt == "-I" || opt == "-INFO" )
        {
            NOMAD::display_info  ( out );
            NOMAD::display_usage ( argv[0], out );
            NOMAD::end();
            return EXIT_SUCCESS;
        }
        
        // parameters creation:
        NOMAD::Parameters p ( out );
        
        // display help on parameters if option '-h' has been specified:
        if ( opt == "-H" || opt == "-HELP" )
        {
            p.help ( argc , argv );
            NOMAD::end();
            return EXIT_SUCCESS;
        }
        
        // display developer help on parameters if option '-d' has been specified:
        if ( opt == "-D" )
        {
            p.help ( argc , argv,true );
            NOMAD::end();
            return EXIT_SUCCESS;
        }
        
        
        // check the number of processess:
#ifdef USE_MPI
        if ( NOMAD::Slave::get_nb_processes() < 2 )
        {
            std::cerr << "ERROR: Incorrect command to run with MPI." << std::endl;
            NOMAD::display_usage ( argv[0], std::cerr );
            NOMAD::end();
            return EXIT_FAILURE;
        }
#endif
        
        try {
            
            
            // read parameters file:
            p.read ( param_file_name );
            
            // parameters check:
            p.check();
            
            // display NOMAD info:
            if ( p.get_display_degree() > NOMAD::MINIMAL_DISPLAY)
                NOMAD::display_info ( out );
            
            // parameters display:
            if ( NOMAD::Slave::is_master() &&
                p.get_display_degree() == NOMAD::FULL_DISPLAY )
                out << std::endl
                << NOMAD::open_block ( "parameters" ) << std::endl
                << p
                << NOMAD::close_block();
            
            // algorithm creation and execution:
            NOMAD::Mads mads ( p , NULL );
            if ( p.get_nb_obj() == 1 )
                mads.run();
            else
                mads.multi_run();
            
#ifdef MODEL_STATS
            mads.display_model_stats ( out );
#endif
            
        }
        catch ( std::exception & e )
        {
            if ( NOMAD::Slave::is_master() )
            {
                error = std::string ( "NOMAD has been interrupted: " ) + e.what();
                std::cerr << std::endl << error << std::endl << std::endl;
            }
        }
        
        NOMAD::Slave::stop_slaves ( out );
        NOMAD::end();
    }
    
#ifdef MEMORY_DEBUG
    NOMAD::display_cardinalities ( out );
#endif
    
    return ( error.empty() ) ? EXIT_SUCCESS : EXIT_FAILURE;
}

/*-----------------------------------------------------*/
/*  display NOMAD most important structures in memory  */
/*-----------------------------------------------------*/
#ifdef MEMORY_DEBUG
void NOMAD::display_cardinalities ( const NOMAD::Display & out )
{
    
#ifdef USE_MPI
    if ( !NOMAD::Slave::is_master() )
        return;
#endif
    
    // compute the biggest int value for appropriate display width:
    int max = (NOMAD::Double::get_max_cardinality() > NOMAD::Point::get_max_cardinality())
    ? NOMAD::Double::get_max_cardinality() : NOMAD::Point::get_max_cardinality();
    if ( NOMAD::Direction::get_max_cardinality() > max )
        max = NOMAD::Direction::get_max_cardinality();
    if ( NOMAD::Set_Element<NOMAD::Eval_Point>::get_max_cardinality() > max )
        max = NOMAD::Set_Element<NOMAD::Eval_Point>::get_max_cardinality();
    if ( NOMAD::Set_Element<NOMAD::Signature>::get_max_cardinality() > max )
        max = NOMAD::Set_Element<NOMAD::Signature>::get_max_cardinality();
    if ( NOMAD::Cache_File_Point::get_max_cardinality() > max )
        max = NOMAD::Cache_File_Point::get_max_cardinality();
    
    // cardinalities display:
    // ----------------------
    out << std::endl
    << NOMAD::open_block ( "important objects in memory" );
    
    // NOMAD::Signature:
    out << "Signature              : ";
    out.display_int_w ( NOMAD::Signature::get_cardinality() , max );
    out << " (max=";
    out.display_int_w ( NOMAD::Signature::get_max_cardinality() , max );
    out << ")" << std::endl;
    
    // NOMAD::Double:
    out << "Double                 : ";
    out.display_int_w ( NOMAD::Double::get_cardinality() , max );
    out << " (max=";
    out.display_int_w ( NOMAD::Double::get_max_cardinality() , max );
    out << ")" << std::endl;
    
    // NOMAD::Point:
    out << "Point                  : ";
    out.display_int_w ( NOMAD::Point::get_cardinality() , max );
    out << " (max=";
    out.display_int_w ( NOMAD::Point::get_max_cardinality() , max );
    out << ")" << std::endl;
    
    // NOMAD::Direction:
    out << "Direction              : ";
    out.display_int_w ( NOMAD::Direction::get_cardinality() , max );
    out << " (max=";
    out.display_int_w ( NOMAD::Direction::get_max_cardinality() , max );
    out << ")" << std::endl;
    
    // Set_Element<Eval_Point>:
    out << "Set_Element<Eval_Point>: ";
    out.display_int_w (NOMAD::Set_Element<NOMAD::Eval_Point>::get_cardinality(), max);
    out << " (max=";
    out.display_int_w (NOMAD::Set_Element<NOMAD::Eval_Point>::get_max_cardinality(), max);
    out << ")" << std::endl;
    
    // Set_Element<NOMAD::Signature>:
    out << "Set_Element<Signature> : ";
    out.display_int_w (NOMAD::Set_Element<NOMAD::Signature>::get_cardinality(), max);
    out << " (max=";
    out.display_int_w (NOMAD::Set_Element<NOMAD::Signature>::get_max_cardinality(), max);
    out << ")" << std::endl;
    
    // NOMAD::Cache_File_Point:
    out << "Cache_File_Point       : ";
    out.display_int_w ( NOMAD::Cache_File_Point::get_cardinality() , max );
    out << " (max=";
    out.display_int_w ( NOMAD::Cache_File_Point::get_max_cardinality() , max );
    out << ")" << std::endl;
    
    out << NOMAD::close_block();
}
#endif

/*------------------------------------------*/
/*            display NOMAD version         */
/*------------------------------------------*/
void NOMAD::display_version ( const NOMAD::Display & out )
{
#ifdef USE_MPI
    if ( !NOMAD::Slave::is_master() )
        return;
#endif
    out << std::endl << "NOMAD - version "
    << NOMAD::VERSION << " - www.gerad.ca/nomad"
    << std::endl << std::endl;
}

/*------------------------------------------*/
/*          display NOMAD information       */
/*------------------------------------------*/
void NOMAD::display_info ( const NOMAD::Display & out )
{
#ifdef USE_MPI
    if ( !NOMAD::Slave::is_master() )
        return;
#endif
    NOMAD::display_version ( out );
    out << NOMAD::open_block ( "Copyright (C) 2001-2015" )
    << "Mark A. Abramson     - The Boeing Company"              << std::endl
    << "Charles Audet        - Ecole Polytechnique de Montreal" << std::endl
    << "Gilles Couture       - Ecole Polytechnique de Montreal" << std::endl
    << "John E. Dennis, Jr.  - Rice University"                 << std::endl
    << "Sebastien Le Digabel - Ecole Polytechnique de Montreal" << std::endl
    << "Christophe Tribes    - Ecole Polytechnique de Montreal" << std::endl
    << NOMAD::close_block()
    << std::endl
    << "Funded in part by AFOSR and Exxon Mobil."               << std::endl
    << std::endl
    << "License   : \'" << NOMAD::LGPL_FILE       << "\'" << std::endl
    << "User guide: \'" << NOMAD::USER_GUIDE_FILE << "\'" << std::endl
    << "Examples  : \'" << NOMAD::EXAMPLES_DIR    << "\'" << std::endl
    << "Tools     : \'" << NOMAD::TOOLS_DIR       << "\'" << std::endl
    << std::endl
    << "Please report bugs to nomad@gerad.ca"
    << std::endl;
}

/*------------------------------------------*/
/*             display NOMAD usage          */
/*------------------------------------------*/
void NOMAD::display_usage ( char* exeName, const NOMAD::Display & out )
{
#ifdef USE_MPI
    if ( !NOMAD::Slave::is_master() )
        return;
    out << std::endl
    << "Run NOMAD.MPI  : mpirun -np p " << exeName << " parameters_file" << std::endl
    << "Info           : " << exeName << " -i"                           << std::endl
    << "Help           : " << exeName << " -h keyword(s) (or 'all')"     << std::endl
    << "Developer help : " << exeName << " -d keyword(s) (or 'all')"     << std::endl
    << "Version        : " << exeName << " -v"                           << std::endl
    << "Usage          : " << exeName << " -u"                          << std::endl
    << std::endl;  
#else
    out << std::endl
    << "Run NOMAD      : " << exeName << " parameters_file"          << std::endl
    << "Info           : " << exeName << " -i"                       << std::endl
    << "Help           : " << exeName << " -h keyword(s) (or 'all')" << std::endl
    << "Developer help : " << exeName << " -d keyword(s) (or 'all')" << std::endl
    << "Version        : " << exeName << " -v"                       << std::endl
    << "Usage          : " << exeName << " -u"                       << std::endl
    << std::endl; 
#endif
}
