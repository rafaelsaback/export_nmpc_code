/**
 * FOR CODE GENERATION, BOTH RUNNING AND TERMINAL COST FUNCSTIONS MUST BE SET
 */

#include <acado_code_generation.hpp>
#include <acado_toolkit.hpp>
#include <acado_gnuplot.hpp>
#include <acado/matrix_vector/matrix_vector.hpp>

USING_NAMESPACE_ACADO

void SetAbsV(Expression &absV, DifferentialState &v);
void SetAbsU(Expression &absU, Control &u_thruster);
void SetGravityBuoyancy(Expression &gravityBuoyancy, DifferentialState &n,
        DVector &rg, DVector &rb, double &weight, double &buoyancy);
void SetJv(Expression &Jv, DifferentialState &v, DifferentialState &n);
DMatrix readVariable(const std::string &folder, const std::string file);

int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        std::cerr << "You should provide a system as argument, either gazebo or flat_fish." << std::endl;
        return 1;
    }

    //========================================================================================
    //                                  SETTING VARIABLES
    //========================================================================================

    std::string system = argv[1];

    double weight;
    double buoyancy;
    std::string dataFolder;

    if(system == "gazebo")
    {
        weight = 4600;
        buoyancy = 4618;
        dataFolder = "./config/gazebo/";
    }
    else if(system == "flatfish" || system == "flat_fish")
    {
        weight = 2800;
        buoyancy = 2812;
        dataFolder = "./config/flatfish/";
    }
    else
    {
        std::cerr << "Unknown system " << system << "." << std::endl;
        return 1;
    }

    double uncertainty = 0;

     /* true  - Export the code to the specified folder
      * false - Simulates the controller inside Acado's environment */
    std::string rockFolder = "/home/rafaelsaback/rock/1_dev/control/uwv_model_pred_control/src";

    /***********************
     * CONTROLLER SETTINGS
     ***********************/

    // Prediction and control horizon
    double horizon         = 2;
    double sampleTime      = 0.1;
    int    iteractions     = horizon / sampleTime;
    bool   positionControl = false;

    //========================================================================================
    //                                  DEFINING VARIABLES
    //========================================================================================

    DifferentialEquation f;                     // Variable that will carry 12 ODEs (6 for position and 6 for velocity)

    DifferentialState v("Velocity", 6, 1);      // Velocity States
    DifferentialState n("Pose", 6, 1);          // Pose States
    DifferentialState tau("Efforts", 6, 1);     // Effort
    Control tauDot("Efforts rate", 6, 1);       // Effort rate of change

    DVector rg(3);                              // Center of gravity
    DVector rb(3);                              // Center of buoyancy

    DMatrix M(6, 6);                            // Inertia matrix
    DMatrix Minv(6, 6);                         // Inverse of inertia matrix
    DMatrix Dl(6, 6);                           // Linear damping matrix
    DMatrix Dq(6, 6);                           // Quadratic damping matrix

    // H-representation of the output domain's polyhedron
    DMatrix AHRep(12,6);
    DMatrix BHRep(12,1);

    Expression linearDamping(6, 6);             // Dl*v
    Expression quadraticDamping(6, 6);          // Dq*|v|*v
    Expression gravityBuoyancy(6);              // g(n)
    Expression absV(6, 6);                      // |v|
    Expression Jv(6);                           // J(n)*v
    Expression vDot(6);                         // AUV Fossen's equation

    Function J;                                 // Cost function
    Function Jn;                                // Terminal cost function

    AHRep = readVariable("./config/", "a_h_rep.dat");
    BHRep = readVariable("./config/", "b_h_rep.dat");
    M = readVariable(dataFolder, "m_matrix.dat");
    Dl = readVariable(dataFolder, "dl_matrix.dat");
    Dq = readVariable(dataFolder, "dq_matrix.dat");

    /***********************
     * LEAST-SQUARES PROBLEM
     ***********************/

    // Cost Functions
    if (positionControl)
    {
        J  << n;
        Jn << n;
    }
    else
    {
        J  << v;
        Jn << v;
    }

    J << tauDot;

    // Weighting Matrices for simulational controller
    DMatrix Q  = eye<double>(J.getDim());
    DMatrix QN = eye<double>(Jn.getDim());

    Q = readVariable("./config/", "q_matrix.dat");

    //========================================================================================
    //                                  MOTION MODEL EQUATIONS
    //========================================================================================
    rg(0) = 0;  rg(1) = 0;  rg(2) = 0;
    rb(0) = 0;  rb(1) = 0;  rb(2) = 0;

    M  *= (2 -(1 + uncertainty));
    Dl *= 1 + uncertainty;
    Dq *= 1 + uncertainty;

    SetAbsV(absV, v);

    Minv = M.inverse();
    linearDamping = Dl * v;
    quadraticDamping = Dq * absV * v;
    SetGravityBuoyancy(gravityBuoyancy, n, rg, rb, weight, buoyancy);
    SetJv(Jv, v, n);

    // Dynamic Equation
    vDot = Minv * (tau - linearDamping - quadraticDamping - gravityBuoyancy);

    // Differential Equations
    f << dot(v) == vDot;
    f << dot(n) == Jv;
    f << dot(tau) == tauDot;

    //========================================================================================
    //                              OPTIMAL CONTROL PROBLEM (OCP)
    //========================================================================================

    OCP ocp(0, horizon, iteractions);

    // Weighting Matrices for exported controller
    BMatrix Qexp = eye<bool>(J.getDim());
    BMatrix QNexp = eye<bool>(Jn.getDim());

    ocp.minimizeLSQ(Qexp, J);
    ocp.minimizeLSQEndTerm(QNexp, Jn);
    ocp.subjectTo(f);

    // Individual degrees of freedom constraints
    ocp.subjectTo(tau(3) == 0);
    ocp.subjectTo(tau(4) == 0);

    /* Note: If the constraint for Roll is not defined the MPC does not
     * work since there's no possible actuation in that DOF.
     * Always consider Roll equal to zero. */

    // Polyhedron set of inequalities
    ocp.subjectTo(AHRep*tau - BHRep <= 0);


    //========================================================================================
    //                                  CODE GENERATION
    //========================================================================================

    // Export code to ROCK
    OCPexport mpc(ocp);

    int Ni = 1;

    mpc.set(HESSIAN_APPROXIMATION, GAUSS_NEWTON);
    mpc.set(DISCRETIZATION_TYPE, SINGLE_SHOOTING);
    mpc.set(INTEGRATOR_TYPE, INT_RK4);
    mpc.set(NUM_INTEGRATOR_STEPS, iteractions * Ni);
    mpc.set(HOTSTART_QP, YES);
    mpc.set(QP_SOLVER, QP_QPOASES);
    mpc.set(GENERATE_TEST_FILE, NO);
    mpc.set(GENERATE_MAKE_FILE, NO);
    mpc.set(GENERATE_MATLAB_INTERFACE, NO);
    mpc.set(GENERATE_SIMULINK_INTERFACE, NO);
    mpc.set(CG_USE_VARIABLE_WEIGHTING_MATRIX, YES);
    mpc.set(CG_HARDCODE_CONSTRAINT_VALUES, YES);

    if (mpc.exportCode(rockFolder.c_str()) != SUCCESSFUL_RETURN)
        exit( EXIT_FAILURE);

    mpc.printDimensionsQP();

    return EXIT_SUCCESS;
}

void SetAbsV(Expression &absV, DifferentialState &v)
{
    for (uint i = 0; i < absV.getNumCols(); i++)
        absV(i, i) = (1 / (1 + exp(-100 * v(i))) - 1 / (1 + exp(100 * v(i))))
        * v(i);
}

void SetAbsU(Expression &absU, Control &u)
{
    for (uint i = 0; i < absU.getNumCols(); i++)
        absU(i, i) = (1 / (1 + exp(-100 * u(i))) - 1 / (1 + exp(100 * u(i))))
        * u(i);
}

void SetGravityBuoyancy(Expression &gravityBuoyancy, DifferentialState &n,
        DVector &rg, DVector &rb, double &weight, double &buoyancy)
{
    gravityBuoyancy(0) = -(weight - buoyancy) * sin(n(4));
    gravityBuoyancy(1) = (weight - buoyancy) * (cos(n(4)) * sin(n(3)));
    gravityBuoyancy(2) = (weight - buoyancy) * (cos(n(4)) * cos(n(3)));
    gravityBuoyancy(3) = (rg(1) * weight - rb(1) * buoyancy) * cos(n(4)) * cos(n(3))
                    - (rg(2) * weight - rb(2) * buoyancy) * cos(n(4)) * sin(n(3));
    gravityBuoyancy(4) =  -(rg(2) * weight - rb(2) * buoyancy) * sin(n(4))
                    - (rg(0) * weight - rb(0) * buoyancy) * cos(n(4)) * cos(n(3));
    gravityBuoyancy(5) = (rg(0) * weight - rb(0) * buoyancy) * cos(n(4))
                    * sin(n(3)) + (rg(1) * weight - rb(1) * buoyancy) * sin(n(4));
}

void SetJv(Expression &Jv, DifferentialState &v, DifferentialState &n)
{
    Expression J(6, 6);

    J(0, 3) = 0; J(0, 4) = 0; J(0, 5) = 0;
    J(1, 3) = 0; J(1, 4) = 0; J(1, 5) = 0;
    J(2, 3) = 0; J(2, 4) = 0; J(2, 5) = 0;
    J(3, 0) = 0; J(3, 1) = 0; J(3, 2) = 0;
    J(4, 0) = 0; J(4, 1) = 0; J(4, 2) = 0;
    J(5, 0) = 0; J(5, 1) = 0; J(5, 2) = 0;

    J(0, 0) = cos(n(5)) * cos(n(4));
    J(0, 1) = -sin(n(5)) * cos(n(3)) + cos(n(5)) * sin(n(4)) * sin(n(3));
    J(0, 2) = sin(n(5)) * sin(n(3)) + cos(n(5)) * cos(n(3)) * sin(n(4));
    J(1, 0) = sin(n(5)) * cos(n(4));
    J(1, 1) = cos(n(5)) * cos(n(3)) + sin(n(3)) * sin(n(4)) * sin(n(5));
    J(1, 2) = -cos(n(5)) * sin(n(3)) + sin(n(4)) * sin(n(5)) * cos(n(3));
    J(2, 0) = -sin(n(4));
    J(2, 1) = cos(n(4)) * sin(n(3));
    J(2, 2) = cos(n(4)) * cos(n(3));

    J(3, 3) = 1;
    J(3, 4) = sin(n(3)) * tan(n(4));
    J(3, 5) = cos(n(3)) * tan(n(4));
    J(4, 3) = 0;
    J(4, 4) = cos(n(3));
    J(4, 5) = -sin(n(3));
    J(5, 3) = 0;
    J(5, 4) = sin(n(3)) / cos(n(4));
    J(5, 5) = cos(n(3)) / cos(n(4));

    Jv = J * v;
}

DMatrix readVariable(const std::string &folder, const std::string filename)
{
    std::ifstream file;
    DMatrix variable;

    file.open(folder + filename);

    if(!file)
        throw std::runtime_error("Couldn't open the file " + folder + filename + ".");

    variable.read(file);
    file.close();
    return variable;
}
