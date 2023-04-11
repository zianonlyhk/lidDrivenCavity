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

int main()
{
    int N = 20;
    double tStop = 0.25;
    double cfl = 0.8;
    double re = 0.1;
    double a = 1.0;

    // initially all zero w.r.t. the gauge pressure
    Eigen::VectorXd uVecInitialCond = Eigen::VectorXd::Zero(N * (N - 1));
    Eigen::VectorXd vVecInitialCond = Eigen::VectorXd::Zero(N * (N - 1));
    Eigen::VectorXd pVecInitialCond = Eigen::VectorXd::Zero(N * N);

    UnsteadySolver testSolver(N, tStop, cfl, re, a);
    testSolver.setCompDomainVec(uVecInitialCond, vVecInitialCond, pVecInitialCond);
    testSolver.setOutputDataAttributes((std::string) "/Users/zianhuang/Room214N/dev/mphil/lidDrivenCavity/", (std::string) "testSolver");
    testSolver.constructMatrixA_uv();
    testSolver.constructMatrixA_p();

    double t = 0.0;
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

        t += testSolver.dt();
        testSolver.writeDataToFiles(t);
    }

    return 0;
}