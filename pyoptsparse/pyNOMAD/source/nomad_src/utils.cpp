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
  \file   utils.cpp
  \brief  Utility functions
  \author Sebastien Le Digabel
  \date   2010-03-23
  \see    utils.hpp
*/
#include "utils.hpp"

/*---------------------------------------------------------------*/
/*               construct the first n prime numbers             */
/*---------------------------------------------------------------*/
void NOMAD::construct_primes ( int n , int * primes )
{
  bool   is_prime;
  double kk;
  int    i = 0 , k = 2 , j;
  while ( true ) {
    is_prime = true;
    kk = sqrt ( static_cast<double>(k) );
    for ( j = 2 ; j <= kk ; ++j )
      if ( k%j==0 ) {
	is_prime = false;
	break;
      }
    if ( is_prime ) {
      primes[i++] = k;
      if ( i==n )
	break;
    }
    ++k;
  }
}

/*--------------------------------------------*/
/*  decompose a string (sentence) into words  */
/*--------------------------------------------*/
void NOMAD::get_words ( const std::string & sentence , std::list<std::string> & words )
{
  std::string        s;
  std::istringstream in ( sentence );
   while ( true ) {
    in >> s;
    if ( in.fail() )
      break;
    words.push_back ( s );
  }
}

/*---------------------------------------------------------------*/
/*  get the pid (process id): used as random seed or unique tag  */
/*---------------------------------------------------------------*/
int NOMAD::get_pid ( void )
{
#ifdef _MSC_VER
  return _getpid();
#else
  return getpid();
#endif
}

/*---------------------------------------*/
/*    called at the beginning of NOMAD   */
/*---------------------------------------*/
void NOMAD::begin ( int argc , char ** argv )
{
#ifdef USE_MPI
  MPI_Init ( &argc , &argv );
#endif
}

/*---------------------------------------*/
/*      called at the end of NOMAD       */
/*---------------------------------------*/
void NOMAD::end ( void )
{
#ifdef USE_MPI
  MPI_Finalize();
#endif
}

/*-----------------------------------------------------------------*/
/*              check if a file exists and is executable           */
/*-----------------------------------------------------------------*/
bool NOMAD::check_exe_file ( const std::string & file_name )
{
#ifdef WINDOWS
  // don't check on Windows:
  return true;
#else
  return ( access ( file_name.c_str() , X_OK ) == 0 );
#endif
}

/*-----------------------------------------------------------------*/
/*              check if a file exists and is readable             */
/*-----------------------------------------------------------------*/
bool NOMAD::check_read_file ( const std::string & file_name )
{
#ifdef _MSC_VER
  return ( _access ( file_name.c_str() , 4 ) == 0 );
#else
  return ( access ( file_name.c_str() , R_OK ) == 0 );
#endif
}

/*-----------------------------------------------------------------*/
/*                         NOMAD::itos                             */
/*-----------------------------------------------------------------*/
std::string NOMAD::itos ( int i )
{
  std::ostringstream oss;
  oss << i;
  return oss.str();
}

/*-----------------------------------------------------------------*/
/*                         NOMAD::itos                             */
/*-----------------------------------------------------------------*/
std::string NOMAD::itos ( size_t i )
{
    std::ostringstream oss;
    oss << i;
    return oss.str();
}


/*-----------------------------------------------------------------*/
/*                         NOMAD::toupper - 1/2                    */
/*-----------------------------------------------------------------*/
void NOMAD::toupper ( std::string & s )
{
  size_t ns = s.size();
  for ( size_t i = 0 ; i < ns ; ++i )
    s[i] = std::toupper(s[i]);
}

/*-----------------------------------------------------------------*/
/*                         NOMAD::toupper - 2/2                    */
/*-----------------------------------------------------------------*/
void NOMAD::toupper ( std::list<std::string> & ls )
{
  std::list<std::string>::iterator       it;
  std::list<std::string>::const_iterator end = ls.end();
  for ( it = ls.begin() ; it != end ; ++it )
    NOMAD::toupper ( *it );
}

/*-----------------------------------------------------------------*/
/*                             NOMAD::atoi                         */
/*-----------------------------------------------------------------*/
bool NOMAD::atoi ( const std::string & s , int & i )
{
  i = -1;
  if ( s.empty() )
    return false;

  size_t n = s.size();

  if ( s[0] == '-' ) {
    if ( n > 1 && s[1] == '-' )
      return false;
    std::string ss = s;
    ss.erase(ss.begin());
    if ( NOMAD::atoi ( ss , i ) ) {
      i = -i;
      return true;
    }
    return false;
  }

  for ( size_t k = 0 ; k < n ; ++k )
    if ( !isdigit(s[k]) )
      return false;
  i = std::atoi(s.c_str());
  return true;
}

bool NOMAD::atoi ( char c , int & i ) {
  std::string s = "-";
  s[0] = c;
  return NOMAD::atoi(s,i);
}

/*-------------------------------------------------------------------*/
/*  if a NOMAD::bb_output_type variable corresponds to a constraint  */
/*-------------------------------------------------------------------*/
bool NOMAD::bbot_is_constraint ( NOMAD::bb_output_type bbot )
{
  return ( bbot == NOMAD::EB     ||
	   bbot == NOMAD::PB     ||
	   bbot == NOMAD::PEB_P  ||
	   bbot == NOMAD::PEB_E  ||
	   bbot == NOMAD::FILTER    );
}

/*-----------------------------------------------------------------------*/
/*  if a NOMAD::direction_type variable corresponds to a MADS direction  */
/*-----------------------------------------------------------------------*/
bool NOMAD::dir_is_mads ( NOMAD::direction_type dt )
{
  return (dt == NOMAD::ORTHO_1        ||
		  dt == NOMAD::ORTHO_2        ||
		  dt == NOMAD::ORTHO_NP1_QUAD ||
		  dt == NOMAD::ORTHO_NP1_NEG  ||
		  dt == NOMAD::DYN_ADDED      ||
		  dt == NOMAD::ORTHO_2N       ||
		  dt == NOMAD::LT_1           ||
		  dt == NOMAD::LT_2           ||
		  dt == NOMAD::LT_2N          ||
		  dt == NOMAD::LT_NP1      );
}

/*----------------------------------------------------------------------*/
/*  if a NOMAD::direction_type variable corresponds to a GPS direction  */
/*----------------------------------------------------------------------*/
bool NOMAD::dir_is_gps ( NOMAD::direction_type dt )
{
  return ( dt == NOMAD::GPS_BINARY             ||
	   dt == NOMAD::GPS_2N_STATIC          ||
	   dt == NOMAD::GPS_2N_RAND            ||
	   dt == NOMAD::GPS_NP1_STATIC_UNIFORM ||
	   dt == NOMAD::GPS_NP1_STATIC         ||
	   dt == NOMAD::GPS_NP1_RAND_UNIFORM   ||
	   dt == NOMAD::GPS_NP1_RAND              );
}

/*--------------------------------------------------------------------------*/
/*  if a NOMAD::direction_type variable corresponds to a LT-MADS direction  */
/*--------------------------------------------------------------------------*/
bool NOMAD::dir_is_ltmads ( NOMAD::direction_type dt )
{
  return ( dt == NOMAD::LT_1     ||
	   dt == NOMAD::LT_2     ||
	   dt == NOMAD::LT_2N    ||
	   dt == NOMAD::LT_NP1      );
}

/*-------------------------------------------------------------------------*/
/*  if a NOMAD::direction_type variable corresponds to a random direction  */
/*-------------------------------------------------------------------------*/
bool NOMAD::dir_is_random ( NOMAD::direction_type dt )
{
  return ( dt == NOMAD::GPS_NP1_RAND         ||
	   dt == NOMAD::GPS_NP1_RAND_UNIFORM ||
	   dt == NOMAD::GPS_2N_RAND          ||
	   dt == NOMAD::LT_1                 ||
	   dt == NOMAD::LT_2                 ||
	   dt == NOMAD::LT_2N                ||
	   dt == NOMAD::LT_NP1                  );
}

/*-----------------------------------------------------------------------------*/
/*  if a NOMAD::direction_type variable corresponds to a Ortho-MADS direction  */
/*-----------------------------------------------------------------------------*/
bool NOMAD::dir_is_orthomads ( NOMAD::direction_type dt )
{
  return (dt == NOMAD::ORTHO_1        ||
		  dt == NOMAD::ORTHO_2        ||
		  dt == NOMAD::ORTHO_NP1_QUAD ||
		  dt == NOMAD::ORTHO_NP1_NEG  ||
		  dt == NOMAD::ORTHO_2N         );
}

/*---------------------------------------------------------------------*/
/*  check if a set of directions includes one Ortho-MADS direction     */
/*  (true if at least one direction in the set is of type Ortho-MADS)  */
/*---------------------------------------------------------------------*/
bool NOMAD::dirs_have_orthomads ( const std::set<NOMAD::direction_type> & dir_types )
{
  std::set<NOMAD::direction_type>::const_iterator it , end = dir_types.end();
  for ( it = dir_types.begin() ; it != end ; ++it )
    if ( NOMAD::dir_is_orthomads (*it) )
      return true;
  return false;
}

/*---------------------------------------------------------------------*/
/*  check if a set of directions includes one Ortho-MADS N+1 direction */
/*  (true if at least one direction in the set is of type Ortho-MADS)  */
/*---------------------------------------------------------------------*/
bool NOMAD::dirs_have_orthomads_np1( const std::set<NOMAD::direction_type> & dir_types )
{
	std::set<NOMAD::direction_type>::const_iterator it , end = dir_types.end();
	for ( it = dir_types.begin() ; it != end ; ++it )
		if ( (*it)==NOMAD::ORTHO_NP1_QUAD || (*it)==NOMAD::ORTHO_NP1_NEG)
			return true;
	return false;
}

/*---------------------------------------------------*/
/*  returns true if one of the string of ls is in s  */
/*---------------------------------------------------*/
bool NOMAD::string_find ( const std::string & s , const std::list<std::string> & ls )
{
  std::list<std::string>::const_iterator it , end = ls.end();
  for ( it = ls.begin() ; it != end ; ++it )
    if ( NOMAD::string_find ( s , *it ) )
      return true;
  return false;
}

/*---------------------------------------------------*/
/*        search a string into another string        */
/*---------------------------------------------------*/
bool NOMAD::string_find ( const std::string & s1 , const std::string & s2 )
{
  return ( s1.find(s2) < s1.size() );
}

/*--------------------------------------------------------------------------*/
/*        returns true if one of the string in ls is an element  of s       */
/*--------------------------------------------------------------------------*/
bool NOMAD::string_match ( const std::string & s , const std::list<std::string> & ls )
{
	std::list<std::string>::const_iterator it , end = ls.end();
	for ( it = ls.begin() ; it != end ; ++it )
		if ( s.compare(*it) == 0 )
			return true;
	return false;
}


/*-----------------------------------------------------------------*/
/*         interpret a list of strings as a direction type         */
/*-----------------------------------------------------------------*/
bool NOMAD::strings_to_direction_type ( const std::list<std::string> & ls ,
					NOMAD::direction_type        & dt   )
{
	
	dt = NOMAD::UNDEFINED_DIRECTION;
	
	if ( ls.empty() || ls.size() > 4 )
		return false;
	
	std::list<std::string>::const_iterator it = ls.begin() , end = ls.end();
	std::string                            s  = *it;
	NOMAD::toupper ( s );
	
	// no direction:
	if ( s == "NONE" ) {
		dt = NOMAD::NO_DIRECTION;
		return true;
	}
	
	// Ortho-MADS with 1, 2, n+1 (plus QUAD or NEG), or 2n directions:
	if ( s == "ORTHO" ) 
	{
		++it;
		if ( it == end )
		{
			dt = NOMAD::ORTHO_NP1_QUAD;  // Default for ORTHO
			return true;
		}
		if ( *it == "1" )
		{
			dt = NOMAD::ORTHO_1;
			return true;
		}
		if ( *it == "2" )
		{
			dt = NOMAD::ORTHO_2;
			return true;
		}
		s = *it;
		NOMAD::toupper ( s );
		if ( s == "2N" )
		{
			dt = NOMAD::ORTHO_2N;
			return true;
		}
		if ( s == "N+1" )
		{
			++it;
			if (it==end)
			{
				dt = NOMAD::ORTHO_NP1_QUAD;   // Default for ORTHO N+1
				return true;
			}
			s = *it;
			NOMAD::toupper ( s );
			if ( s=="QUAD" )
			{
				dt= NOMAD::ORTHO_NP1_QUAD;
				return true;
			}
			if ( s=="NEG" )
			{
				dt=NOMAD::ORTHO_NP1_NEG;
				return true;
			}
			
		}
		
		return false;
	}
	
	// LT-MADS with 1, 2 or 2n directions:
	if ( s == "LT" ) {
		++it;
		if ( it == end ) {
			dt = NOMAD::LT_2N;
			return true;
		}
		if ( *it == "1" ) {
			dt = NOMAD::LT_1;
			return true;
		}
		if ( *it == "2" ) {
			dt = NOMAD::LT_2;
			return true;
		}
		s = *it;
		NOMAD::toupper ( s );
		if ( s == "N+1" ) {
			dt = NOMAD::LT_NP1;
			return true;
		}
		if ( s == "2N" ) {
			dt = NOMAD::LT_2N;
			return true;
		}
		return false;
	}
	
	// GPS:
	if ( s == "GPS" ) {
		++it;
		if ( it == end ) {
			dt = NOMAD::GPS_2N_STATIC;
			return true;
		}
		s = *it;
		NOMAD::toupper ( s );
		
		// GPS for binary variables:
		if ( s == "BINARY" || s == "BIN" ) {
			dt = NOMAD::GPS_BINARY;
			return true;
		}
		
		// GPS, n+1 directions:
		if ( s == "N+1" ) {
			++it;
			if ( it == end ) {
				dt = NOMAD::GPS_NP1_STATIC;
				return true;
			}
			s = *it;
			NOMAD::toupper ( s );
			
			// GPS, n+1, static:
			if ( s == "STATIC" ) {
				++it;
				if ( it == end ) {
					dt = NOMAD::GPS_NP1_STATIC;
					return true;
				}
				s = *it;
				NOMAD::toupper ( s );
				if ( s == "UNIFORM" ) {
					dt = NOMAD::GPS_NP1_STATIC_UNIFORM;
					return true;
				}
				return false;
			}
			
			// GPS, n+1, random:
			if ( s == "RAND" || s == "RANDOM" ) {
				++it;
				if ( it == end ) {
					dt = NOMAD::GPS_NP1_RAND;
					return true;
				}
				s = *it;
				NOMAD::toupper ( s );
				if ( s == "UNIFORM" ) {
					dt = NOMAD::GPS_NP1_RAND_UNIFORM;
					return true;
				}
				return false;
			}
			return false;
		}
		
		// 2n directions:
		if ( s == "2N" ) {
			++it;
			if ( it == end ) {
				dt = NOMAD::GPS_2N_STATIC;
				return true;
			}
			s = *it;
			NOMAD::toupper ( s );
			if ( s == "STATIC" ) {
				dt = NOMAD::GPS_2N_STATIC;
				return true;
			}
			if ( s == "RAND" || s == "RANDOM" ) {
				dt = NOMAD::GPS_2N_RAND;
				return true;
			}
			return false;
		}
		return false;
	}
	return false;
}

/*---------------------------------------*/
/*   convert a string into a hnorm_type  */
/*---------------------------------------*/
bool NOMAD::string_to_hnorm_type ( const std::string & s , NOMAD::hnorm_type & hn )
{
  std::string ss = s;
  NOMAD::toupper(ss);
  if ( ss == "L1" ) {
    hn = NOMAD::L1;
    return true;
  }
  if ( ss == "L2" ) {
    hn = NOMAD::L2;
    return true;
  }
  if ( ss == "LINF" ) {
    hn = NOMAD::LINF;
    return true;
  }
  return false;
}

/*-----------------------------------------*/
/*  convert a string into a TGP_mode_type  */
/*-----------------------------------------*/
bool NOMAD::string_to_TGP_mode_type ( const std::string & s , NOMAD::TGP_mode_type & m )
{
  std::string ss = s;
  NOMAD::toupper(ss);
  if ( ss == "FAST" ) {
    m = NOMAD::TGP_FAST;
    return true;
  }
  if ( ss == "PRECISE" ) {
    m = NOMAD::TGP_PRECISE;
    return true;
  }
  if ( ss == "USER" ) {
    m = NOMAD::TGP_USER;
    return true;
  }
  return false;
}

/*--------------------------------------------------*/
/*  convert a string into a multi_formulation_type  */
/*--------------------------------------------------*/
bool NOMAD::string_to_multi_formulation_type ( const std::string             & s   ,
					       NOMAD::multi_formulation_type & mft   )
{
  std::string ss = s;
  NOMAD::toupper(ss);
  if ( ss == "NORMALIZED" ) {
    mft = NOMAD::NORMALIZED;
    return true;
  }
  if ( ss == "PRODUCT" ) {
    mft = NOMAD::PRODUCT;
    return true;
  }
  if ( ss == "DIST_L1" ) {
    mft = NOMAD::DIST_L1;
    return true;
  }
  if ( ss == "DIST_L2" ) {
    mft = NOMAD::DIST_L2;
    return true;
  }
  if ( ss == "DIST_LINF" ) {
    mft = NOMAD::DIST_LINF;
    return true;
  }
  return false;
}

/*-------------------------------------------*/
/*   convert a string into a bb_output_type  */
/*-------------------------------------------*/
bool NOMAD::string_to_bb_output_type ( const std::string     & s    ,
				       NOMAD::bb_output_type & bbot   )
{
  std::string ss = s;
  NOMAD::toupper(ss);
  
  if ( ss == "OBJ" ) {
    bbot = NOMAD::OBJ;
    return true;
  }
  if ( ss == "EB" ) {
    bbot = NOMAD::EB;
    return true;
  }
  if ( ss == "PB" || ss == "CSTR" ) {
    bbot = NOMAD::PB;
    return true;
  }
  if ( ss == "PEB" ) {
    bbot = NOMAD::PEB_P;
    return true;
  }
  if ( ss == "F" ) {
    bbot = NOMAD::FILTER;
    return true;
  }
  if ( ss == "STAT_AVG" ) {
    bbot = NOMAD::STAT_AVG;
    return true;
  }
  if ( ss == "STAT_SUM" ) {
    bbot = NOMAD::STAT_SUM;
    return true;
  }
  if ( ss == "CNT_EVAL" ) {
    bbot = NOMAD::CNT_EVAL;
    return true;
  }
  if ( ss == "NOTHING" || ss == "-" ) {
    bbot = NOMAD::UNDEFINED_BBO;
    return true;
  }
  return false;
}

/*-----------------------------------------------------------------*/
/*                convert a string into a bb_input_type            */
/*-----------------------------------------------------------------*/
bool NOMAD::string_to_bb_input_type ( const std::string    & s    ,
				      NOMAD::bb_input_type & bbit   )
{
  std::string ss = s;
  NOMAD::toupper ( ss );
  if ( ss=="R" || ss=="REAL" ) {
    bbit = NOMAD::CONTINUOUS;
    return true;
  }
  if ( ss=="C" || ss=="CAT" ) {
    bbit = NOMAD::CATEGORICAL;
    return true;
  }
  if ( ss=="B" || ss=="BIN" ) {
    bbit = NOMAD::BINARY;
    return true;
  }
  if ( ss=="I" || ss=="INT" ) {
    bbit = NOMAD::INTEGER;
    return true;
  }
  return false;
}

/*-----------------------------------------------------------------*/
/*                 convert a string into a model_type              */
/*-----------------------------------------------------------------*/
bool NOMAD::string_to_model_type ( const std::string & s  ,
				   NOMAD::model_type & mt   )
{
  std::string ss = s;
  NOMAD::toupper ( ss );
  if ( ss=="TGP" || ss=="TGP_MODEL" ) {
    mt = NOMAD::TGP_MODEL;
    return true;
  }
  if ( ss=="QUADRATIC" || ss=="QUADRATIC_MODEL" ) {
    mt = NOMAD::QUADRATIC_MODEL;
    return true;
  }

  mt = NOMAD::NO_MODEL;
  return false;
}

/*----------------------------------------------------------------------*/
/*         convert a string in {"YES","NO","Y","N"} to a bool           */
/*         value of return: -1: error                                   */
/*                           0: bool=false                              */
/*                           1: bool=true                               */
/*----------------------------------------------------------------------*/
int NOMAD::string_to_bool ( const std::string & ss )
{
  std::string s = ss;
  NOMAD::toupper ( s );
  if ( s=="Y" || s=="YES" || s=="1" || s=="TRUE" )
    return 1;
  if ( s=="N" || s=="NO" || s=="0" || s=="FALSE" )
    return 0;
  return -1;
}

/*---------------------------------------------------*/
/*  convert a string 'i-j' to the integers i and j   */
/*---------------------------------------------------*/
bool NOMAD::string_to_index_range ( const std::string & s           ,
				    int               & i           ,
				    int               & j           ,
				    int               * n           ,
				    bool                check_order   )
{
  if ( s.empty() )
    return false;
  
  if ( s == "*" ) {
    if ( !n )
      return false;
    i = 0;
    j = *n-1;
    return true;
  }
  
  if ( s[0] == '-' ) {
      
    size_t ns = s.size();
    if ( ns > 1 && s[1] == '-' )
      return false;
    
    std::string ss = s;
    ss.erase ( ss.begin() );
      
    if ( NOMAD::string_to_index_range ( ss , i , j , n , false ) ) {
      i = -i;
      return true;
    }
    return false;
  }

  std::istringstream in (s);
  std::string        s1;

  getline ( in , s1 , '-' );

  if (in.fail())
    return false;

  size_t k , n1 = s1.size();

  if ( n1 >= s.size() - 1 ) {
    for ( k = 0 ; k < n1 ; ++k )
      if (!isdigit(s1[k]))
	return false;
    if ( !NOMAD::atoi ( s1 , i ) )
      return false;
    if ( n1 == s.size() ) {
      j = i;
      return true;
    }
    if (n) {
      j = *n-1;
      return true;
    }
    return false;
  }

  std::string s2;
  getline (in, s2);

  if (in.fail())
    return false;

  size_t n2 = s2.size();
  for ( k = 0 ; k < n2 ; ++k )
    if ( !isdigit(s2[k]) )
      return false;

  if ( !NOMAD::atoi ( s1, i ) || !NOMAD::atoi ( s2 , j ) )
    return false;

  return !check_order || i <= j;
}


int NOMAD::get_rank(double ** M,
					size_t m,
					size_t n)
{
	double  * W = new double  [n];
	double ** V = new double *[n];
	for (size_t i = 0 ; i < n ; ++i )
	{
		V[i]=new double [n];
	}
	
	std::string error_msg;
	NOMAD::SVD_decomposition ( error_msg , M , W , V , static_cast<int>(m) , static_cast<int>(n) );

	for (size_t i=0;i<n;++i)
		delete [] V[i];
	delete [] V;

	
 	if (! error_msg.empty())
	{
		delete [] W;
		return -1;
	}
	
	int rank=0;
	for (size_t i=0;i<n;i++)
	{
		if (fabs(W[i])>NOMAD::SVD_EPS)
			rank++;
	}
	
	delete [] W;
	return rank;
	
}


/*--------------------------------------------------------------*/
/*                        SVD decomposition                     */
/*  inspired and recoded from an old numerical recipes version  */
/*--------------------------------------------------------------*/
/*                                                              */
/*           M = U . W . V'                                     */
/*                                                              */
/*           M ( m x n )                                        */
/*           U ( m x n )                                        */
/*           W ( n x n )                                        */
/*           V ( n x n )                                        */
/*                                                              */
/*           U.U' = U'.U = I if m = n                           */
/*           U'.U = I        if m > m                           */
/*                                                              */
/*           V.V' = V'.V = I                                    */
/*                                                              */
/*           W diagonal                                         */
/*                                                              */
/*           M is given as first argument and becomes U         */
/*           W is given as a size-n vector                      */
/*           V is given, not V'                                 */
/*                                                              */
/*--------------------------------------------------------------*/
bool NOMAD::SVD_decomposition ( std::string & error_msg ,
				double     ** M         ,
				double      * W         ,
				double     ** V         ,
				int           m         ,
				int           n         ,
				int           max_mpn     ) // default=1500
{
	error_msg.clear();
	
	if ( max_mpn > 0 && m+n > max_mpn ) {
		error_msg = "SVD_decomposition() error: m+n > " + NOMAD::itos ( max_mpn );
		return false;
	}
	
	double * rv1   = new double[n];
	double   scale = 0.0;
	double   g     = 0.0;
	double   norm  = 0.0;
	
	int      nm1   = n - 1;
	
	bool   flag;
	int    i , j , k , l , its , jj , nm = 0;
	double s , f , h , tmp , c , x , y , z , absf , absg , absh;
	
	const int NITER = 30;
	
	// Initialization W and V
	for (i=0; i < n; ++i)
	{
		W[i]=0.0;
		for (j=0; j < n ; ++j)
			V[i][j]=0.0;
	}
	
	// Householder reduction to bidiagonal form:
	for ( i = 0 ; i < n ; ++i )
	{
		l      = i + 1;
		rv1[i] = scale * g;
		g      = s = scale = 0.0;
		if ( i < m )
		{
			for ( k = i ; k < m ; ++k )
				scale += fabs ( M[k][i] );
			if ( scale ) {
				for ( k = i ; k < m ; ++k ) 
				{
					M[k][i] /= scale;
					s += M[k][i] * M[k][i];
				}
				f       = M[i][i];
				g       = ( f >= 0.0 ) ? -fabs(sqrt(s)) : fabs(sqrt(s));
				h       = f * g - s;
				M[i][i] = f - g;
				for ( j = l ; j < n ; ++j ) {
					for ( s = 0.0 , k = i ; k < m ; ++k )
						s += M[k][i] * M[k][j];
					f = s / h;
					for ( k = i ; k < m ; ++k )
						M[k][j] += f * M[k][i];
				}
				for ( k = i ; k < m ; ++k )
					M[k][i] *= scale;
			}
		}
		W[i] = scale * g;
		g    = s = scale = 0.0;
		if ( i < m && i != nm1 ) 
		{
			for ( k = l ; k < n ; ++k )
				scale += fabs ( M[i][k] );
			if ( scale ) 
			{
				for ( k = l ; k < n ; ++k )
				{
					M[i][k] /= scale;
					s       += M[i][k] * M[i][k];
				}
				f       = M[i][l];
				g       = ( f >= 0.0 ) ? -fabs(sqrt(s)) : fabs(sqrt(s));
				h       = f * g - s;
				M[i][l] = f - g;
				for ( k = l ; k < n ; ++k )
					rv1[k] = M[i][k] / h;
				for ( j = l ; j < m ; ++j ) {
					for ( s=0.0,k=l ; k < n ; ++k )
						s += M[j][k] * M[i][k];
					for ( k=l ; k < n ; ++k )
						M[j][k] += s * rv1[k];
				}
				for ( k = l ; k < n ; ++k )
					M[i][k] *= scale;
			}
		}
		tmp  = fabs ( W[i] ) + fabs ( rv1[i] );
		norm = ( norm > tmp ) ? norm : tmp;
	}
	
	// accumulation of right-hand transformations:
	for ( i = nm1 ; i >= 0 ; --i )
	{
		if ( i < nm1 ) 
		{
			if ( g ) 
			{
				for ( j = l ; j < n ; ++j )
					V[j][i] = ( M[i][j] / M[i][l] ) / g;
				for ( j = l ; j < n ; ++j ) {
					for ( s = 0.0 , k = l ; k < n ; ++k )
						s += M[i][k] * V[k][j];
					for ( k = l ; k < n ; ++k )
						V[k][j] += s * V[k][i];
				}
			}
			for ( j = l ; j < n ; ++j )
				V[i][j] = V[j][i] = 0.0;
		}
		V[i][i] = 1.0;
		g       = rv1[i];
		l       = i;
	}
	
	// accumulation of left-hand transformations:
	for ( i = ( ( m < n ) ? m : n ) - 1 ; i >= 0 ; --i ) 
	{
		l = i + 1;
		g = W[i];
		for ( j = l ; j < n ; ++j )
			M[i][j] = 0.0;
		if ( g ) {
			g = 1.0 / g;
			for ( j = l ; j < n ; ++j ) 
			{
				for ( s = 0.0 , k = l ; k < m ; ++k )
					s += M[k][i] * M[k][j];
				f = ( s / M[i][i] ) * g;
				for ( k = i ; k < m ; ++k )
					M[k][j] += f * M[k][i];
			}
			for ( j = i ; j < m ; ++j )
				M[j][i] *= g;
		}
		else
			for ( j = i ; j < m ; ++j )
				M[j][i] = 0.0;
		++M[i][i];
	}
	
	// diagonalization of the bidiagonal form:
	for ( k = nm1 ; k >= 0 ; --k )
	{
		for ( its = 1 ; its <= NITER ; its++ ) 
		{
			flag = true;
			for ( l = k ; l >= 0 ; l-- ) 
			{
				nm = l - 1;
				if ( nm < 0 || fabs ( rv1[l]) + norm == norm )
				{
					flag = false;
					break;
				}
				if ( fabs ( W[nm] ) + norm == norm )
					break;
			}
			if ( flag )
			{
				c = 0.0;
				s = 1.0;
				for ( i = l ; i <= k ; i++ )
				{
					f      = s * rv1[i];
					rv1[i] = c * rv1[i];
					if ( fabs(f) + norm == norm )
						break;
					g = W[i];
					
					absf = fabs(f);
					absg = fabs(g);
					h    = ( absf > absg ) ?
					absf * sqrt ( 1.0 + pow ( absg/absf , 2.0 ) ) :
					( ( absg==0 ) ? 0.0 : absg * sqrt ( 1.0 + pow ( absf/absg , 2.0 ) ) );
					
					W[i] =  h;
					h    =  1.0 / h;
					c    =  g * h;
					s    = -f * h;
					for ( j = 0 ; j < m ; ++j ) {
						y = M[j][nm];
						z = M[j][ i];
						M[j][nm] = y * c + z * s;
						M[j][ i] = z * c - y * s;
					}
				}
			}
			z = W[k];
			if ( l == k) {
				if ( z < 0.0 ) {
					W[k] = -z;
					for ( j = 0 ; j < n ; j++ )
						V[j][k] = -V[j][k];
				}
				break;  // this 'break' is always active if k==0
			}
			if ( its == NITER ) 
			{
				error_msg = "SVD_decomposition() error: no convergence in " +
	            NOMAD::itos ( NITER ) + " iterations";
				delete [] rv1;
				return false;
			}
			x  = W[l];
			nm = k - 1;
			y  = W[nm];
			g  = rv1[nm];
			h  = rv1[k];
			f  = ( (y-z) * (y+z) + (g-h) * (g+h) ) / ( 2.0 * h * y );
			
			absf = fabs(f);
			g    = ( absf > 1.0 ) ?
			absf * sqrt ( 1.0 + pow ( 1.0/absf , 2.0 ) ) :
			sqrt ( 1.0 + pow ( absf , 2.0 ) );
			
			f = ( (x-z) * (x+z) +
				 h * ( ( y / ( f + ( (f >= 0)? fabs(g) : -fabs(g) ) ) ) - h ) ) / x;
			c = s = 1.0;
			
			for ( j = l ; j <= nm ; ++j ) {
				i = j + 1;
				g = rv1[i];
				y = W[i];
				h = s * g;
				g = c * g;
				
				absf = fabs(f);
				absh = fabs(h);
				z    = ( absf > absh ) ?
				absf * sqrt ( 1.0 + pow ( absh/absf , 2.0 ) ) :
				( ( absh==0 ) ? 0.0 : absh * sqrt ( 1.0 + pow ( absf/absh , 2.0 ) ) );
				
				rv1[j] = z;
				c      = f / z;
				s      = h / z;
				f      = x * c + g * s;
				g      = g * c - x * s;
				h      = y * s;
				y     *= c;
				for ( jj = 0 ; jj < n ; ++jj )
				{
					x = V[jj][j];
					z = V[jj][i];
					V[jj][j] = x * c + z * s;
					V[jj][i] = z * c - x * s;
				}
				
				absf = fabs(f);
				absh = fabs(h);
				z    = ( absf > absh ) ?
				absf * sqrt ( 1.0 + pow ( absh/absf , 2.0 ) ) :
				( ( absh==0 ) ? 0.0 : absh * sqrt ( 1.0 + pow ( absf/absh , 2.0 ) ) );
				
				W[j] = z;
				
				if ( z ) 
				{
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for ( jj = 0 ; jj < m ; ++jj )
				{
					y = M[jj][j];
					z = M[jj][i];
					M[jj][j] = y * c + z * s;
					M[jj][i] = z * c - y * s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			W  [k] = x;
		}
	}
	
	delete [] rv1;
	return true;
}
