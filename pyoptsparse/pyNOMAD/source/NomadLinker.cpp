#include "NomadLinker.h"
#include "Callback.h"

#include <iostream>
#include "nomad.hpp"
using namespace NOMAD;

NomadLinker::NomadLinker() : callback_(NULL) {}

NomadLinker::~NomadLinker()
{
}

void NomadLinker::setCallback(Callback &callback)
{
    callback_ = &callback;
}

std::vector<double> NomadLinker::call( int n, int m, double* x, int xdim, double* lb, int
		lbdim, double* ub, int ubdim, double min_poll_size, double min_mesh_size, int max_bbe, int display_degree, int print_file )
{
    if( !callback_ )
    {
        std::cerr << "No callback is set.\n";
    }
    else
    {
        /*----------------------------------------*/
        /*  The problem: the evaluation is made   */
        /*  by calling a FORTRAN routine          */
        /*----------------------------------------*/
        class My_Evaluator : public Evaluator 
        {
            private:
                NomadLinker* nlobject;
            public:
                My_Evaluator( const Parameters &p , NomadLinker &nltransfer ) : Evaluator( p )
                {
                    nlobject = &nltransfer;
                }

                ~My_Evaluator( void ) {}

                bool eval_x( Eval_Point &x, const Double &h_max, bool &count_eval ) const
                {
                    int n = x.size();
                    int m = x.get_bb_outputs().size();
                    int i;

                    // call the FORTRAN routine:
                    std::vector<double> xvec( n );
                    for( i=0; i<n; i++ )
                    {
                        xvec[i] = x[i].value();
                    }
                    nlobject->passcallback()->call( *nlobject, xvec );
                    std::vector<double> fxvec = nlobject->getSolution();
                    
                    for ( i=0; i<m; i++ )
                        x.set_bb_output( i, fxvec[i] );

                    count_eval = true; // count a black-box evaluation

                    return true;       // the evaluation succeeded
                }
        };
		
		/*----------------------------------------*/
        /*  Problem Setup                         */
        /*----------------------------------------*/
        
	    // display:
        Display out( std::cout );
        out.precision( DISPLAY_PRECISION_STD );

        try {
            int i;

            // parameters creation:
            Parameters p( out );

            p.set_DIMENSION( n );             // number of variables

            vector<bb_output_type> bbot( m ); // definition of output types:
            bbot[0] = CNT_EVAL;
            bbot[1] = OBJ;                    // first output : objective value to minimize
            for( i=2; i<m; i++ )      // other outputs: constraints cj <= 0
                bbot[i] = EB;
            p.set_BB_OUTPUT_TYPE( bbot );

            // starting point and bounds:
            Point px0( n );
            Point plb( n );
            Point pub( n );
            for( i=0; i<n; i++ )
            {
                px0[i] = x[i];
                if( lb[i] > -1e20 )
                    plb[i] = lb[i];
                if( ub[i] < 1e20 )
                    pub[i] = ub[i];
            }
            p.set_X0( px0 );
            p.set_LOWER_BOUND( plb );
            p.set_UPPER_BOUND( pub );
            
            // Stopping criteria
            if(min_poll_size)
                p.set_MIN_POLL_SIZE( min_poll_size);
            if(min_mesh_size)
                p.set_MIN_MESH_SIZE( min_mesh_size);

            // maximum number of black-box evaluations:
            if( max_bbe > 0 )
                p.set_MAX_BB_EVAL( max_bbe );

            // display degree:
            p.set_DISPLAY_DEGREE( display_degree );
            // History file
            if( print_file > 0 )
                p.set_HISTORY_FILE( "NOMAD.out" );
            // parameters validation:
            p.check();

            // custom evaluator creation:
            My_Evaluator ev( p, *this);

            // algorithm creation and execution:
            Mads mads( p, &ev );
            mads.run();

            // get the solution:
            const Eval_Point * bf = mads.get_best_feasible();
            if( bf )
                for( i=0; i<n; i++ )
	                x[i] = (*bf)[i].value();
        }
        catch( exception &e )
        {
            cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
        }
		
        return solution;
    }
}

void NomadLinker::setSolution( std::vector<double> transfer)
{
    solution = transfer;
}

std::vector<double> NomadLinker::getSolution()
{
    return solution;
}

