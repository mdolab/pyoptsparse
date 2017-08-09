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
  \file   Pareto_Front.cpp
  \brief  Pareto front (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-09
  \see    Pareto_Front.hpp
*/
#include "Pareto_Front.hpp"

/*------------------------------------------------------*/
/*                 insertion of a point                 */
/*   (returns true if the point is a new Pareto point)  */
/*------------------------------------------------------*/
bool NOMAD::Pareto_Front::insert ( const NOMAD::Eval_Point & x )
{

  NOMAD::Pareto_Point pp ( &x );
    
  if ( _pareto_pts.empty() ) {
    _pareto_pts.insert (pp);
    return true;
  }

  bool insert = false;

  std::set<NOMAD::Pareto_Point>::iterator it = _pareto_pts.begin();
  while ( it != _pareto_pts.end() ) {
    if ( pp.dominates (*it) ) {
      _pareto_pts.erase(it++);
      insert = true;

      continue;
    }
    ++it;
  }
    
  if ( !insert ) {
    insert = true;
    std::set<NOMAD::Pareto_Point>::iterator end = _pareto_pts.end();
    for ( it = _pareto_pts.begin() ; it != end ; ++it ) {
      if ( it->dominates (pp) ) {
	insert = false;
	break;
      }
    }
  }
      
  if ( insert ) {
    _pareto_pts.insert ( pp );
    return true;
  }
  return false;
}

/*-------------------------------------------------------*/
/*            get the point that minimizes f2(x)         */
/*  (this is simply the last point of the Pareto front)  */
/*-------------------------------------------------------*/
const NOMAD::Eval_Point * NOMAD::Pareto_Front::get_best_f2 ( void ) const
{
  if ( _pareto_pts.empty() )
    return NULL;

  std::set<NOMAD::Pareto_Point>::const_iterator it = _pareto_pts.end();
  --it;

  return it->get_element();
}

/*------------------------------------------------------*/
/*                get the reference point               */
/*------------------------------------------------------*/
NOMAD::Point * NOMAD::Pareto_Front::get_ref ( const NOMAD::Pareto_Point *& xj ,
					      NOMAD::Double & delta_j ) const
{
  xj = NULL;
  delta_j.clear();

  int p = static_cast<int>(_pareto_pts.size());
    
  // no points in the front:
  if ( p == 0 )
    return NULL;

  // just one point in the front:
  if ( p == 1 ) {
    xj      = &(*_pareto_pts.begin());
    delta_j = 1.0 / ( xj->get_w() + 1 );  // delta=1.0
    return NULL;
  }
    
  std::set<NOMAD::Pareto_Point>::const_iterator it = _pareto_pts.begin();
  NOMAD::Point * ref = new NOMAD::Point ( 2 );

  NOMAD::Double f1xm1;  // f_1 ( x_{j-1} )
  NOMAD::Double f1x;    // f_1 ( x_j     )
  NOMAD::Double f1xp1;  // f_1 ( x_{j+1} )

  NOMAD::Double f2xm1;  // f_2 ( x_{j-1} )
  NOMAD::Double f2x;    // f_2 ( x_j     )
  NOMAD::Double f2xp1;  // f_2 ( x_{j+1} )
    
  // two points in the front:
  if ( p == 2 ) {

    f1xm1 = it->get_f1();
    f2xm1 = it->get_f2();

    ++it;
    xj = &(*it);
    
    f1x = xj->get_f1();
    f2x = xj->get_f2();

    delta_j = ( (f1x-f1xm1).pow2() + (f2x-f2xm1).pow2() ) / ( xj->get_w() + 1.0 );
    
    const_cast<NOMAD::Pareto_Point *>(xj)->update_w();
    
    (*ref)[0] = f1x;
    (*ref)[1] = f2xm1;

    return ref;
  }

  // more than two points in the front:
  std::set<NOMAD::Pareto_Point>::const_iterator end = _pareto_pts.end();
    
  const NOMAD::Pareto_Point * prev , * cur , * next;

  NOMAD::Double delta;

  prev = &(*it);
  ++it;
    
  while ( true ) {

    cur = &(*it);

    ++it;
    if ( it == end )
      break;
    
    next = &(*it);
   
    f1xm1 = prev->get_f1();
    f2xm1 = prev->get_f2();
    
    f1x   = cur->get_f1();
    f2x   = cur->get_f2();

    f1xp1 = next->get_f1();
    f2xp1 = next->get_f2();

    delta = ( (f1x-f1xm1).pow2() + (f2x-f2xm1).pow2() +
	      (f1x-f1xp1).pow2() + (f2x-f2xp1).pow2()   ) / ( cur->get_w() + 1.0 );
      
    if ( !delta_j.is_defined() || delta > delta_j ) {
      xj      = cur;
      delta_j = delta;
      (*ref)[0] = f1xp1;
      (*ref)[1] = f2xm1;
    }
    
    prev = cur;
  }
  
  const_cast<Pareto_Point *>(xj)->update_w();
  
  return ref;
}

/*------------------------------------------------------*/
/*                  compute delta and surf              */
/*------------------------------------------------------*/
void NOMAD::Pareto_Front::get_delta_surf ( NOMAD::Double      & delta_j  ,
					   NOMAD::Double      & surf     ,
					   const NOMAD::Point & f_bounds   ) const
{
  bool def = f_bounds.is_complete();
  NOMAD::Double f1_min , f1_max , f2_min , f2_max;

  if ( def ) {

    if ( f_bounds.size() == 4 ) {

      f1_min = f_bounds[0];
      f1_max = f_bounds[1];
      f2_min = f_bounds[2];
      f2_max = f_bounds[3];

      if ( f1_min >= f1_max || f2_min >= f2_max ) {
	f1_min.clear();
	f1_max.clear();
	f2_min.clear();
	f2_max.clear();
	def = false;
      }
    }
    else
      def = false;
  }

  delta_j.clear();
  surf.clear();

  int p = static_cast<int> ( _pareto_pts.size() );

  // no point in the front:
  if ( p == 0 ) {
    if ( def )
      surf = 1.0;
    return;
  }

  const NOMAD::Pareto_Point * xj;
  NOMAD::Double               f1x;  // f_1 ( x_j )
  NOMAD::Double              f2x;  // f_2 ( x_j )

  NOMAD::Double surf_frame = (def) ?
    ( f2_max - f2_min ) * ( f1_max - f1_min )
    : NOMAD::Double();

  // just one point in the front:
  if ( p == 1 ) {
    xj      = &(*_pareto_pts.begin());
    delta_j = 1.0 / ( xj->get_w() + 1 );  // delta=1.0
    f1x     = xj->get_f1();
    f2x     = xj->get_f2();

    if ( !def || f1x > f1_max || f1x < f1_min || f2x > f2_max || f2x < f2_min )
      return;

    surf = (   ( f2_max - f2_min ) * (    f1x - f1_min )
	     + (    f2x - f2_min ) * ( f1_max - f1x    ) ) / surf_frame;

    return;
  }

  std::set<NOMAD::Pareto_Point>::const_iterator it = _pareto_pts.begin();

  NOMAD::Double f1xm1;  // f_1 ( x_{j-1} )
  NOMAD::Double f1xp1;  // f_1 ( x_{j+1} )

  NOMAD::Double f2xm1;  // f_2 ( x_{j-1} )
  NOMAD::Double f2xp1;  // f_2 ( x_{j+1} )
    
  // two points in the front:
  if ( p == 2 ) {

    f1xm1 = it->get_f1();
    f2xm1 = it->get_f2();

    if ( def && ( f1xm1 < f1_min ||
		  f1xm1 > f1_max ||
		  f2xm1 < f2_min ||
		  f2xm1 > f2_max    ) )
      def = false;

    ++it;
    xj = &(*it);
    
    f1x = xj->get_f1();
    f2x = xj->get_f2();

    if ( def && ( f1x < f1_min ||
		  f1x > f1_max ||
		  f2x < f2_min ||
		  f2x > f2_max    ) )
      def = false;
    
    delta_j = ( (f1x-f1xm1).pow2() + (f2x-f2xm1).pow2() ) / ( xj->get_w() + 1.0 );

    if ( def )
      surf  = (    ( f2xm1  - f2_min ) * ( f1x    - f1xm1  )
		 + ( f2_max - f2_min ) * ( f1xm1  - f1_min )
	         + ( f2x    - f2_min ) * ( f1_max - f1x    )  ) / surf_frame;

    return;
  }

  // more than two points in the front:
  std::set<NOMAD::Pareto_Point>::const_iterator end = _pareto_pts.end();
    
  const NOMAD::Pareto_Point * prev , * cur , * next;

  NOMAD::Double delta;

  prev  = &(*it);
  f1xm1 = prev->get_f1();
  f2xm1 = prev->get_f2();

  ++it;

  cur  = &(*it);
  f1x  = cur->get_f1();

  if ( def && ( f1xm1 < f1_min ||
		f1xm1 > f1_max ||
		f2xm1 < f2_min ||
		f2xm1 > f2_max ||
	        f1x   < f1_min ||
		f1x   > f1_max    ) )
    def = false;

  if ( def )
    surf =  ( f2xm1  - f2_min ) * ( f1x   - f1xm1  )
          + ( f2_max - f2_min ) * ( f1xm1 - f1_min );

  while ( true ) {

    cur = &(*it);

    ++it;
    if ( it == end )
      break;
    
    next = &(*it);
   
    f1xm1 = prev->get_f1();
    f2xm1 = prev->get_f2();
    
    f1x   = cur->get_f1();
    f2x   = cur->get_f2();

    f1xp1 = next->get_f1();
    f2xp1 = next->get_f2();


    if ( def &&
	 ( f1xm1 < f1_min || f1xm1 > f1_max || f2xm1 < f2_min || f2xm1 > f2_max ||
	   f1x   < f1_min || f1x   > f1_max || f2x   < f2_min || f2x   > f2_max ||
	   f1xp1 < f1_min || f1xp1 > f1_max || f2xp1 < f2_min || f2xp1 > f2_max    ) )
      def = false;

    delta = ( (f1x-f1xm1).pow2() + (f2x-f2xm1).pow2() +
 	      (f1x-f1xp1).pow2() + (f2x-f2xp1).pow2()   ) / ( cur->get_w() + 1.0 );
      
    if ( !delta_j.is_defined() || delta > delta_j )
      delta_j = delta;

    if ( def )
      surf += ( f2x - f2_min ) * ( f1xp1 - f1x );
    
    prev = cur;
  }

  if ( def ) {
    surf += ( f2xp1 - f2_min ) * ( f1_max - f1xp1 );
    surf  = surf / surf_frame;
  }
  else
    surf.clear();
}

/*------------------------------------------------------*/
/*               display the Pareto points              */
/*------------------------------------------------------*/
void NOMAD::Pareto_Front::display ( const NOMAD::Display & out ) const
{
  size_t nb  = _pareto_pts.size();
  int    cnt = 0;
  std::set<NOMAD::Pareto_Point>::const_iterator it , end = _pareto_pts.end();
  for ( it = _pareto_pts.begin() ; it != end ; ++it ) {
    out << "#";
    out.display_int_w ( cnt++ , static_cast<int>(nb) );
    out << " ";
    it->display ( out );
    out << std::endl;
  }
}

/*------------------------------------------------------------------*/
/*  . begin() and next() methods, to browse the front               */
/*                                                                  */
/*  . example:                                                      */
/*                                                                  */
/*        const Eval_Point * cur = pareto_front.begin();            */
/*        while (cur) {                                             */
/*          ...                                                     */
/*          cur = pareto_front.next();                              */
/*        }                                                         */
/*------------------------------------------------------------------*/

// begin():
// --------
const NOMAD::Eval_Point * NOMAD::Pareto_Front::begin ( void ) const
{
  if ( _pareto_pts.empty() )
    return NULL;
  _it = _pareto_pts.begin();
  return _it->get_element();
}

// next(): (supposes that begin() has been called)
// -------
const NOMAD::Eval_Point * NOMAD::Pareto_Front::next ( void ) const
{
  if ( _pareto_pts.empty() )
    return NULL;
  ++_it;
  if ( _it == _pareto_pts.end() )
    return NULL;
  return _it->get_element();
}
