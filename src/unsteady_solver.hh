/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   unsteady_solver.hh                                Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/04/10 15:22:45 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#ifndef UNSTEADY_SOLVER_HH
#define UNSTEADY_SOLVER_HH

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <Eigen/Dense>
#include <Eigen/Sparse>

class UnsteadySolver
{
public:
    // SETTING AND INITIALISING

    UnsteadySolver();
    // initial conditions part 1:
    UnsteadySolver(int N, double tStop, double cfl, double re, double a);
    // initial conditions part 2:
    void setCompDomainVec(const Eigen::VectorXd &uVec, const Eigen::VectorXd &vVec, const Eigen::VectorXd &pVec);
    // other attributes
    void setOutputDataAttributes(std::string simulationName, std::string repoDir);

    // PUBLIC METHODS

    int getVecIdxU(int i, int j);
    int getVecIdxV(int i, int j);
    int getVecIdxP(int i, int j);
    // solving the linear system at each time step
    void constructMatrixA_uv();
    void constructMatrixA_p();
    void constructMatrixA_streamFunc();

    void constructLoadVecU();
    void constructLoadVecV();
    void constructLoadVecP();

    void constructVorticityVec();

    void solveForU_Star(); // step 1 update fieldVecU
    void solveForV_Star(); // step 1 update fieldVecV
    void solveForP_Next(); // step 2 update fieldVecP
    void solveForU_Next(); // step 3 update fieldVecU again
    void solveForV_Next(); // step 3 update fieldVecV again

    void solveForStreamFunc();

    void checkIfBreak(double time);

    void writeDataToFiles(double time);

    // ACCESSING PRIVATE MEMBERS
    const double tStop();
    const double dt();
    const double reachedSteady();

private:
    // CONSTANT ATTRIBUTES
    // -----------------------------

    double m_tol = 1e-5;
    double m_re;
    double m_cfl;
    double m_a;
    // discretisation parameters
    int m_N;
    // spatial parameters
    double m_x0 = 0.0;
    double m_x1 = 1.0;
    double m_y0 = m_x0;
    double m_y1 = m_x1;
    double m_dx;
    // temporal parameters
    double m_tStart = 0.0;
    double m_tStop;
    bool m_reachedSteady = false;

    std::string m_repoDir;
    std::string m_simulationName;
    std::ofstream m_uResults;
    std::ofstream m_vResults;
    std::ofstream m_pResults;
    std::ofstream m_vecResults;
    std::ofstream m_streamFuncResults;

    // DYNAMIC ATTRIBUTES
    // -----------------------------

    double m_dt;

    // COMPUTATIONAL DOMAIN ATTRIBUTES
    // -----------------------------

    // sparse matrices at each time step
    Eigen::SparseMatrix<double> m_sMatrixA_uv;
    Eigen::SparseMatrix<double> m_sMatrixA_p;
    Eigen::SparseMatrix<double> m_sMatrixA_streamFunc;
    // sparse matrix solvers
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> m_solverUV;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> m_solverP;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> m_solverStreamFunc;

    // field vectors at each step
    Eigen::VectorXd m_uVec;
    Eigen::VectorXd m_vVec;
    Eigen::VectorXd m_pVec;
    Eigen::VectorXd m_vortVec;
    Eigen::VectorXd m_streamFuncVec;

    // previous time step field vectors
    Eigen::VectorXd m_uVecPrev;
    Eigen::VectorXd m_vVecPrev;
    // field vector difference
    Eigen::VectorXd m_uVecDiff;
    Eigen::VectorXd m_vVecDiff;

    // load vectors at each step
    Eigen::VectorXd m_uLoadVec;
    Eigen::VectorXd m_vLoadVec;
    Eigen::VectorXd m_pLoadVec;
};

#endif