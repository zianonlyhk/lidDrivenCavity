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
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Sparse>

class UnsteadySolver
{
public:
    // SETTING AND INITIALISING

    UnsteadySolver();
    // initial conditions part 1:
    UnsteadySolver(int N, double tStop, double cfl, double re);
    // initial conditions part 2:
    void setCompDomainVec(const Eigen::VectorXd &uVec, const Eigen::VectorXd &vVec, const Eigen::VectorXd &pVec);

    // PUBLIC METHODS

    void updateDt();

    // solving the linear system at each time step
    void constructMatrixA_uv(double delTime);
    void constructMatrixA_p(); // only once
    void constructLoadVecU();
    void constructBcVecU();
    void constructLoadVecV();
    void constructBcVecV();

    void solveForU_Star(); // step 1 update fieldVecU
    void solveForV_Star(); // step 1 update fieldVecV
    void solveForP_Next(); // step 2 update fieldVecP
    void solveForU_Next(); // step 3 update fieldVecU again
    void solveForV_Next(); // step 3 update fieldVecV again

    // ACCESSING PRIVATE MEMBERS
    const double tStop();
    const double dt();

private:
    // CONSTANT ATTRIBUTES
    // -----------------------------

    double m_tol = 1e-5;
    double m_re;
    double m_cfl;
    // discretisation parameters
    int m_N;
    // spatial parameters
    double m_x0 = 0.0;
    double m_x1 = 1.0;
    double m_y0 = 0.0;
    double m_y1 = 1.0;
    double m_dx;
    // temporal parameters
    double m_tStart = 0.0;
    double m_tStop;

    // DYNAMIC ATTRIBUTES
    // -----------------------------

    double m_dt;

    // COMPUTATIONAL DOMAIN ATTRIBUTES
    // -----------------------------

    // sparse matrices at each time step
    Eigen::SparseMatrix<double> m_sMatrixA_uv;
    Eigen::SparseMatrix<double> m_sMatrixA_p;
    // sparse matrix solvers
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> m_solverUV;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> m_solverP;

    // field vectors at each step
    Eigen::VectorXd m_uVec;
    Eigen::VectorXd m_vVec;
    Eigen::VectorXd m_pVec;
    // previous time step field vectors
    Eigen::VectorXd m_uVecPrev;
    Eigen::VectorXd m_vVecPrev;
    Eigen::VectorXd m_pVecPrev;
    // field vector difference
    Eigen::VectorXd m_uVecDiff;
    Eigen::VectorXd m_vVecDiff;
    Eigen::VectorXd m_pVecDiff;

    // load vectors at each step
    Eigen::VectorXd m_uLoadVec;
    Eigen::VectorXd m_vLoadVec;
    Eigen::VectorXd m_pLoadVec;
    // boundary condition vectors at each step
    Eigen::VectorXd m_uBcVec;
    Eigen::VectorXd m_vBcVec;
    Eigen::VectorXd m_pBcVec;
};

#endif