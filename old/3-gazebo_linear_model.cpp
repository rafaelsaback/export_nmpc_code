/**
 * FOR CODE GENERATION, BOTH RUNNING AND TERMINAL COST FUNCSTIONS MUST BE SET
 */

#include <acado_code_generation.hpp>
#include <acado_toolkit.hpp>
#include <acado_gnuplot.hpp>

USING_NAMESPACE_ACADO

void SetAbsV(Expression &absV, DifferentialState &v);
void SetAbsU(Expression &absU, Control &u_thruster);
void SetGravityBuoyancy(Expression &gravityBuoyancy, DifferentialState &n,
        DVector &rg, DVector &rb, double &weight, double &buoyancy);
void SetJv(Expression &Jv, DifferentialState &v, DifferentialState &n);

int main()
{

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

    Expression linearDamping(6, 6);             // Dl*v
    Expression quadraticDamping(6, 6);          // Dq*|v|*v
    Expression gravityBuoyancy(6);              // g(n)
    Expression absV(6, 6);                      // |v|
    Expression Jv(6);                           // J(n)*v
    Expression vDot(6);                         // AUV Fossen's equation

    Function J;                                 // Cost function
    Function Jn;                                // Terminal cost function

    // Initiliazing the matrices as null
    M.setZero();
    Minv.setZero();
    Dl.setZero();
    Dq.setZero();

    //========================================================================================
    //                                  SETTING VARIABLES
    //========================================================================================

    /**
     * true  - Export the code to the specified folder
     * false - Simulates the controller inside Acado's environment
     */
    bool exportCode = true;
    std::string rockFolder =
            "/home/rafaelsaback/flatfish/bir_rock/control/uwv_model_pred_control/src";

    /***********************
     * SIMULATION PARAMETERS
     ***********************/
    double simulationTime = 10;
    DVector reference     = zeros<double>(18);
    reference(4) = -0.1;
//    reference(8) = -5;
//    reference(2)          = -2;

    /***********************
     * CONTROLLER SETTINGS
     ***********************/

    // Prediction and control horizon
    double horizon         = 2;
    double sampleTime      = 0.1;
    int    iteractions     = horizon / sampleTime;
    bool   positionControl = false;

    /***********************
     *       PLOTS
     ***********************/
    bool plotPosition = true;
    bool plotVelocity = true;
    bool plotTau      = true;

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

    if(positionControl)
    {
        // Output Error Weights
        Q(0, 0) = 1; // SURGE
        Q(1, 1) = 1; // SWAY
        Q(2, 2) = 1; // HEAVE
        Q(3, 3) = 1e-6; // ROLL
        Q(4, 4) = 1; // PITCH
        Q(5, 5) = 1; // YAW

        // Effort Weights
        Q(6, 6) = 1e-4;   // SURGE
        Q(7, 7) = 1e-4;   // SWAY
        Q(8, 8) = 1e-5;   // HEAVE
        Q(9, 9) = 1e+6;      // ROLL
        Q(10, 10) = 1e-13; // PITCH
        Q(11, 11) = 1e-4; // YAW

        // Effort Rate Weights
        Q(12, 12) = 1e-4;   // SURGE
        Q(13, 13) = 1e-4;   // SWAY
        Q(14, 14) = 1e-5;   // HEAVE
        Q(15, 15) = 1e+6;      // ROLL
        Q(16, 16) = 1e-8; // PITCH
        Q(17, 17) = 1e-4; // YAW
    }
    else
    {
        // Output Error Weights
        Q(0, 0) = 1; // SURGE
        Q(1, 1) = 1; // SWAY
        Q(2, 2) = 1; // HEAVE
        Q(3, 3) = 1e-6; // ROLL
        Q(4, 4) = 1; // PITCH
        Q(5, 5) = 1; // YAW

        // Effort Weights
        Q(6, 6) = 1e-7;   // SURGE
        Q(7, 7) = 1e-7;   // SWAY
        Q(8, 8) = 1e-7;   // HEAVE
        Q(9, 9) = 1e+6;      // ROLL
        Q(10, 10) = 1e-7; // PITCH
        Q(11, 11) = 1e-7; // YAW
    }

    /***********************
     * 	  CONSTRAINTS
     ***********************/

    // Variables for allowing the constraints
    bool constrainTau    = true;

    // Constraints' values
    DVector tauLimit    = 100 * ones<double>(6); // Maximum/minimun effort allowed
    tauLimit[4] = 48;   // Pitch
    tauLimit[5] = 140;  // Yaw

    /***********************
     * 	VEHICLE PARAMETERS
     ***********************/

    double weight = 4600;
    double buoyancy = 4605;

    rg(0) = 0;  rg(1) = 0;  rg(2) = 0;
//    rb(0) = -0.00698;  rb(1) = 0;  rb(2) = 0.10186;
    rb(0) = -0.0054289;  rb(1) = 0;  rb(2) = 0.079226;

    M(0, 0) = 460;  M(1, 1) = 460;  M(2, 2) = 460;
    M(3, 3) = 35;   M(4, 4) = 35;   M(5, 5) = 35;
    Dl(0, 0) = 50;  Dl(1, 1) = 50;  Dl(2, 2) = 50;
    Dl(3, 3) = 45;  Dl(4, 4) = 45;  Dl(5, 5) = 45;
    Dq(0, 0) = 40;  Dq(1, 1) = 40;  Dq(2, 2) = 40;
    Dq(3, 3) = 35;  Dq(4, 4) = 35;  Dq(5, 5) = 35;

//    double uncertainty = 1.3;
//    M *= (2-uncertainty);
//    Dl *= uncertainty;
//    Dq *= uncertainty;

    //========================================================================================
    //                                  MOTION MODEL EQUATIONS
    //========================================================================================

    SetAbsV(absV, v);

    Minv = M.inverse();
    linearDamping = Dl * v;
    quadraticDamping = Dq * absV * v;
    SetGravityBuoyancy(gravityBuoyancy, n, rg, rb, weight, buoyancy);
    SetJv(Jv, v, n);

    // Dynamic Equation
    //vDot = Minv * (tau - linearDamping - quadraticDamping - gravityBuoyancy);
    vDot = Minv * (tau - linearDamping - gravityBuoyancy);

    // Differential Equations
    f << dot(v) == vDot;
    f << dot(n) == Jv;
    f << dot(tau) == tauDot;

    //========================================================================================
    //                              OPTIMAL CONTROL PROBLEM (OCP)
    //========================================================================================

    OCP ocp(0, horizon, iteractions);

    if (exportCode)
    {
        // Weighting Matrices for exported controller
        BMatrix Qexp = eye<bool>(J.getDim());
        BMatrix QNexp = eye<bool>(Jn.getDim());

        ocp.minimizeLSQ(Qexp, J);
        ocp.minimizeLSQEndTerm(QNexp, Jn);
    }
    else
    {
        DVector R = zeros<double>(J.getDim(), 1);
        for(uint i = 0; i < reference.size(); i++)
        {
            R[i] = reference[i];
        }
        ocp.minimizeLSQ(Q, J, R);
        ocp.minimizeLSQEndTerm(QN, Jn);
    }

    ocp.subjectTo(f);

    //	CONTROL CONSTRAINTS

    if (constrainTau)
    {
        ocp.subjectTo(-tauLimit(0) <= tau(0) <= tauLimit(0));
        ocp.subjectTo(-tauLimit(1) <= tau(1) <= tauLimit(1));
        ocp.subjectTo(-tauLimit(2) <= tau(2) <= tauLimit(2));
        ocp.subjectTo(tau(3) == 0);
        ocp.subjectTo(-tauLimit(4) <= tau(4) <= tauLimit(4));
        ocp.subjectTo(-tauLimit(5) <= tau(5) <= tauLimit(5));
    }


    //========================================================================================
    //                          CODE GENERATION / CONTROLLER SIMULATION
    //========================================================================================

    //	EXPORTING THE CODE TO ROCK
    if (exportCode)
    {
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
        mpc.set(CG_HARDCODE_CONSTRAINT_VALUES, NO);

        if (mpc.exportCode(rockFolder.c_str()) != SUCCESSFUL_RETURN)
            exit( EXIT_FAILURE);

        mpc.printDimensionsQP();

        return EXIT_SUCCESS;
    }
    //	SIMULATING THE CONTROLLER INSIDE ACADO
    else
    {
        //========================================================================================
        //                               SETTING UP THE PROCESS
        //========================================================================================

        OutputFcn identity;
        DynamicSystem dynamicSystem(f, identity);

        Process process(dynamicSystem, INT_RK45);
        process.set(ABSOLUTE_TOLERANCE, 1.0e-2);
        process.set(INTEGRATOR_TOLERANCE, 1.0e-2);

        //========================================================================================
        //                                      MPC
        //========================================================================================

        RealTimeAlgorithm alg(ocp, sampleTime);
        alg.set(MAX_NUM_ITERATIONS, 3);

        VariablesGrid refGrid(dynamicSystem.getNumDynamicEquations(), 0,
                simulationTime, simulationTime / sampleTime);
        refGrid.setZero();

        for (uint i = 0; i < refGrid.getNumPoints(); i++)
        {
            for (uint j = 0; j < reference.size(); j++)
                if(i <= 20)
                    refGrid(i, j) = i*reference(j)/20;
                else
                    refGrid(i, j) = reference(j);
        }

        StaticReferenceTrajectory trajectory(refGrid);
        Controller controller(alg, trajectory);

        //========================================================================================
        //                              SIMULATION ENVIRONMENT
        //========================================================================================

        SimulationEnvironment sim(0, simulationTime, process, controller);

        DVector x0(dynamicSystem.getNumDynamicEquations());
        x0.setZero();
        x0(0) = 0;

        if (sim.init(x0) != SUCCESSFUL_RETURN)
            exit( EXIT_FAILURE);
        if (sim.run() != SUCCESSFUL_RETURN)
            exit( EXIT_FAILURE);


        // ...AND PLOT THE RESULTS
        // ----------------------------------------------------------
        VariablesGrid sampledProcessOutput;
        sim.getSampledProcessOutput(sampledProcessOutput);

        VariablesGrid feedbackControl;
        sim.getFeedbackControl(feedbackControl);

        if (plotPosition)
        {
            GnuplotWindow positionWindow;
            positionWindow.addSubplot(sampledProcessOutput(6),
                    "Surge Position [m]");
//            positionWindow.addSubplot(sampledProcessOutput(7),
//                    "Sway Position [m]");
            positionWindow.addSubplot(sampledProcessOutput(8),
                    "Heave Position [m]");
//            positionWindow.addSubplot( sampledProcessOutput(9),
//                    "Roll Position [rad]" );
            positionWindow.addSubplot( sampledProcessOutput(10),
                    "Pitch Position [rad]" );
            positionWindow.addSubplot(sampledProcessOutput(11),
                    "Yaw Position [rad]");
            positionWindow.plot();
        }

        if (plotVelocity)
        {
            GnuplotWindow velocityWindow;
            velocityWindow.addSubplot(sampledProcessOutput(0),
                    "Surge Velocity [m/s]");
//            velocityWindow.addSubplot(sampledProcessOutput(1),
//                    "Sway Velocity [m/s]");
            velocityWindow.addSubplot(sampledProcessOutput(2),
                    "Heave Velocity [m/s]");
//            velocityWindow.addSubplot( sampledProcessOutput(3),
//                    "Roll Velocity [rad/s]" );
            velocityWindow.addSubplot( sampledProcessOutput(4),
                    "Pitch Velocity [rad/s]" );
            velocityWindow.addSubplot(sampledProcessOutput(5),
                    "Yaw Velocity [rad/s]");
            velocityWindow.plot();
        }
            GnuplotWindow effortWindow;
            effortWindow.addSubplot(sampledProcessOutput(12),
                    "Surge Effort [N]");
//            effortWindow.addSubplot(sampledProcessOutput(13),
//                    "Sway Effort [N]");
            effortWindow.addSubplot(sampledProcessOutput(14),
                    "Heave Effort [N]");
//            effortWindow.addSubplot( sampledProcessOutput(15),
//                    "Roll Effort [N]" );
            effortWindow.addSubplot( sampledProcessOutput(16),
                    "Pitch Effort [N]" );
            effortWindow.addSubplot(sampledProcessOutput(17),
                    "Yaw Effort [N]");
            effortWindow.plot();

        if (plotTau)
        {
            GnuplotWindow tauDotWindow;
            tauDotWindow.addSubplot(feedbackControl(0), "Surge Effort Increment");
//            tauDotWindow.addSubplot(feedbackControl(1), "Sway Effort Increment");
            tauDotWindow.addSubplot(feedbackControl(2), "Heave Effort Increment");
//            tauDotWindow.addSubplot(feedbackControl(3), "Roll Effort Increment" );
            tauDotWindow.addSubplot(feedbackControl(4),"Pitch Effort Increment" );
            tauDotWindow.addSubplot(feedbackControl(5), "Yaw Effort Increment");
            tauDotWindow.plot();
        }
    }
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
//    gravityBuoyancy(0) = -(weight - buoyancy) * sin(n(4));
//    gravityBuoyancy(1) = (weight - buoyancy) * (cos(n(4)) * sin(n(3)));
//    gravityBuoyancy(2) = (weight - buoyancy) * (cos(n(4)) * cos(n(3)));
//    gravityBuoyancy(3) = (rg(1) * weight - rb(1) * buoyancy) * cos(n(4)) * cos(n(3))
//                    - (rg(2) * weight - rb(2) * buoyancy) * cos(n(4)) * sin(n(3));
//    gravityBuoyancy(4) =  -(rg(2) * weight - rb(2) * buoyancy) * sin(n(4))
//                    - (rg(0) * weight - rb(0) * buoyancy) * cos(n(4)) * cos(n(3));
//    gravityBuoyancy(5) = (rg(0) * weight - rb(0) * buoyancy) * cos(n(4))
//                    * sin(n(3)) + (rg(1) * weight - rb(1) * buoyancy) * sin(n(4));

    gravityBuoyancy(0) = 0;
    gravityBuoyancy(1) = 0;
    gravityBuoyancy(2) = (weight - buoyancy);
    gravityBuoyancy(3) = 0;
    gravityBuoyancy(4) = 0;
    gravityBuoyancy(5) = 0;
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
