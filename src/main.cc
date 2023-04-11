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

void loadConfig(int &i_N, double &i_tStop, double &i_cfl, double &i_re, double &i_a, std::string &i_repoDir, std::string &i_simulationName, int &i_plottingFactor)
{
    std::ifstream in("./config.txt");

    std::string parameter;

    // dummie disposable variables
    double double_value;
    int int_value;
    std::string string_value;

    while (!in.eof())
    {
        // if no read string fits the required parameter, do nothing
        in >> parameter;

        if (parameter == "N")
        {
            in >> int_value;
            i_N = int_value;
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
        else if (parameter == "re")
        {
            in >> double_value;
            i_re = double_value;
        }
        else if (parameter == "a")
        {
            in >> double_value;
            i_a = double_value;
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

int main()
{
    int N;
    double tStop;
    double cfl;
    double re;
    double a;
    std::string repoName;
    std::string simulationName;
    int plottingFactor;
    loadConfig(N, tStop, cfl, re, a, repoName, simulationName, plottingFactor);

    // initially all zero w.r.t. the gauge pressure
    Eigen::VectorXd uVecInitialCond = Eigen::VectorXd::Zero(N * (N - 1));
    Eigen::VectorXd vVecInitialCond = Eigen::VectorXd::Zero(N * (N - 1));
    Eigen::VectorXd pVecInitialCond = Eigen::VectorXd::Zero(N * N);

    UnsteadySolver testSolver(N, tStop, cfl, re, a);
    testSolver.setCompDomainVec(uVecInitialCond, vVecInitialCond, pVecInitialCond);
    testSolver.setOutputDataAttributes(repoName, simulationName);
    testSolver.constructMatrixA_uv();
    testSolver.constructMatrixA_p();

    double t = 0.0;
    int frame = 0;
    testSolver.writeDataToFiles(t);
    while (!testSolver.reachedSteady())
    {
        testSolver.constructLoadVecU();
        testSolver.constructLoadVecV();
        testSolver.constructLoadVecP();

        testSolver.solveForU_Star();
        testSolver.solveForV_Star();

        testSolver.solveForP_Next();

        testSolver.solveForU_Next();
        testSolver.solveForV_Next();

        testSolver.checkIfSteady();

        frame++;
        t += testSolver.dt();
        if (frame % plottingFactor == 0)
        {
            testSolver.writeDataToFiles(t);
        }
    }

    return 0;
}