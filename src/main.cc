/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   main.cc                                           Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/04/10 15:22:36 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "unsteady_solver.hh"

//
void loadConfig(double &i_x0, double &i_x1, double &i_y0, double &i_y1,
                int &i_n_x, int &i_n_y, double &i_tStop, double &i_cfl,
                double &i_rho_0, double &i_temp_0, double &i_c_p, double &i_mu,
                double &i_k, double &i_a, double &i_tempDiff,
                std::string &i_repoDir, std::string &i_simulationName, int &i_plottingFactor)
{
    std::ifstream in("./config");

    std::string parameter;

    // dummie disposable variables
    double double_value;
    int int_value;
    std::string string_value;

    while (!in.eof())
    {
        // if no read string fits the required parameter, do nothing
        in >> parameter;

        if (parameter == "x0")
        {
            in >> double_value;
            i_x0 = double_value;
        }
        else if (parameter == "x1")
        {
            in >> double_value;
            i_x1 = double_value;
        }
        else if (parameter == "y0")
        {
            in >> double_value;
            i_y0 = double_value;
        }
        else if (parameter == "y1")
        {
            in >> double_value;
            i_y1 = double_value;
        }
        else if (parameter == "n_x")
        {
            in >> int_value;
            i_n_x = int_value;
        }
        else if (parameter == "n_y")
        {
            in >> int_value;
            i_n_y = int_value;
        }
        else if (parameter == "tStop")
        {
            in >> double_value;
            i_tStop = double_value;
        }
        else if (parameter == "cfl")
        {
            in >> double_value;
            i_cfl = double_value;
        }
        else if (parameter == "rho_0")
        {
            in >> double_value;
            i_rho_0 = double_value;
        }
        else if (parameter == "temp_0")
        {
            in >> double_value;
            i_temp_0 = double_value;
        }
        else if (parameter == "c_p")
        {
            in >> double_value;
            i_c_p = double_value;
        }
        else if (parameter == "mu")
        {
            in >> double_value;
            i_mu = double_value;
        }
        else if (parameter == "k")
        {
            in >> double_value;
            i_k = double_value;
        }
        else if (parameter == "a")
        {
            in >> double_value;
            i_a = double_value;
        }
        else if (parameter == "tempDiff")
        {
            in >> double_value;
            i_tempDiff = double_value;
        }
        else if (parameter == "repoDir")
        {
            in >> string_value;
            i_repoDir = string_value;
        }
        else if (parameter == "simulationName")
        {
            in >> string_value;
            i_simulationName = string_value;
        }
        else if (parameter == "plottingFactor")
        {
            in >> int_value;
            i_plottingFactor = int_value;
        }
    }
    in.close();
}
//

//
int main()
{
    //
    // #############################################
    // PARAMETERS & VARIABLE FIELDS #############
    // space and time discretisation parameters
    double x0;
    double x1;
    double y0;
    double y1;
    int n_x;
    int n_y;
    double tStop;
    double cfl;

    // physical parameters of the fluid
    double rho_0;
    double temp_0;
    double c_p;
    double mu;
    double k;

    // hydro and buoyancy dynamics parameters
    double a;
    double tempDiff;

    // miscellaneous
    std::string repoDir;
    std::string simulationName;
    int plottingFactor;
    loadConfig(x0, x1, y0, y1,
               n_x, n_y, tStop, cfl,
               rho_0, temp_0, c_p, mu,
               k, a, tempDiff,
               repoDir, simulationName, plottingFactor);

    // initially all zero w.r.t. the gauge pressure
    Eigen::VectorXd uVecInitialCond = Eigen::VectorXd::Zero((n_x - 1) * n_y);
    Eigen::VectorXd vVecInitialCond = Eigen::VectorXd::Zero(n_x * (n_y - 1));
    Eigen::VectorXd pVecInitialCond = Eigen::VectorXd::Zero(n_x * n_y);
    Eigen::VectorXd tempVecInitialCond = temp_0 * Eigen::VectorXd::Ones(n_x * n_y); // room temp
    // END OF DECLARATION #######################
    // #############################################
    //

    //
    // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    // SOLVER CLASS & SOLVE $$$$$$$$$$$$$$$$$$$$$
    // initialisation and declaration
    UnsteadySolver testSolver(x0, x1, y0, y1, n_x, n_y, tStop, cfl);
    testSolver.setPhysicalParameters(rho_0, temp_0, tempDiff, c_p, mu, k, a);
    testSolver.setCompDomainVec(uVecInitialCond, vVecInitialCond, pVecInitialCond, tempVecInitialCond);
    // testSolver.pressureShift();
    testSolver.setOutputDataAttributes(repoDir, simulationName);

    // linear solvers
    testSolver.constructMatrixA_u();
    testSolver.constructMatrixA_v();
    testSolver.constructMatrixA_p();
    testSolver.constructMatrixA_temp();
    testSolver.constructMatrixA_streamFunc();

    double t = 0.0;
    int frame = 0;
    testSolver.writeDataToFiles(t); // initial state
    while (testSolver.okToContinue())
    {
        // load vectors on current field status:
        testSolver.constructLoadVecTemp();
        testSolver.constructLoadVecU();
        testSolver.constructLoadVecV();
        testSolver.constructLoadVecP();

        // pre-projection decoupled equation:
        testSolver.solveForTemp_Next(); // affect v

        // step 1: predictor
        testSolver.solveForU_Star();
        testSolver.solveForV_Star();

        // step 2: corrector
        testSolver.solveForP_Next();

        // step 3: constrain projection
        testSolver.solveForU_Next();
        testSolver.solveForV_Next();

        // post-projection decoupled equations:
        testSolver.constructVorticityVec();
        testSolver.solveForStreamFunc();

        testSolver.checkIfContinue(t);

        // progress checking and data logging
        frame++;
        t += testSolver.dt();
        if (frame % plottingFactor == 0)
        {
            std::cout << "t " << t << "/" << testSolver.tStop() << std::endl;
            testSolver.writeDataToFiles(t);
        }
    }
    // END OF SOLVING $$$$$$$$$$$$$$$$$$$$$$$$$$$
    // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    //

    return 0;
}