/**
 * FOR CODE GENERATION, BOTH RUNNING AND TERMINAL COST FUNCSTIONS MUST BE SET
 */

#include <acado_code_generation.hpp>
#include <acado_toolkit.hpp>
#include <acado_gnuplot.hpp>

USING_NAMESPACE_ACADO

int main( ){

	AlgebraicState  tau;    /*      Thruster effort         */
	Control         f;      /*      Thruster command        */
//	Control         s;      /*      Slack variable          */
	DMatrix TCM(6,6);
	TCM(0,0) = 1;           TCM(0,1) =  1;          TCM(0,2) = 0.0008;      TCM(0,3) = 0.0008;      TCM(0,4) = 0;           TCM(0,5) = 0;
	TCM(1,0) = 0;           TCM(1,1) =  0;          TCM(1,2) = -1;          TCM(1,3) = -1;          TCM(1,4) = -0.0008;     TCM(1,5) = -0.0008;
	TCM(2,0) = 0;           TCM(2,1) =  0;          TCM(2,2) = 0;           TCM(2,3) = 0;           TCM(2,4) = -1;          TCM(2,5) = -1;
	TCM(3,0) = 0;           TCM(3,1) =  0;          TCM(3,2) = 0/*0.0325*/; TCM(3,3) = 0/*0.0325*/; TCM(3,4) = 0;           TCM(3,5) = 0;
	TCM(4,0) = -0.0316;     TCM(4,1) =  -0.316;     TCM(4,2) = 0;           TCM(4,3) = 0;           TCM(4,4) = 0.4235;      TCM(4,5) = -0.556;
	TCM(5,0) = -0.4210;     TCM(5,1) =  0.4200;     TCM(5,2) = -0.5735;     TCM(5,3) = 0.936;       TCM(5,4) = 0;           TCM(5,5) = 0;


}
