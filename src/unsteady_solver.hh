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
#include <fstream>
#include <iostream>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>

class UnsteadySolver
{
public:
    //
    // #############################################
    // SETTING AND INITIALISING #################
    UnsteadySolver();
    // initial conditions part 1:
    UnsteadySolver(double x0, double x1, double y0, double y1, int n_x, int n_y, double tStop, double cfl);
    void setPhysicalParameters(double rho_0, double temp_0, double tempDiff, double c_p, double mu, double k, double a);

    // initial conditions part 2:
    void setCompDomainVec(const Eigen::VectorXd &uVec, const Eigen::VectorXd &vVec, const Eigen::VectorXd &pVec, const Eigen::VectorXd &tempVec);
    void pressureShift();

    // other attributes
    void setOutputDataAttributes(std::string simulationName, std::string repoDir);
    // END OF INITIALISATION ####################
    // #############################################
    //

    //
    // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    // METHODS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    // indexing public methods
    int getVecIdxU(int i, int j);
    int getVecIdxV(int i, int j);
    int getVecIdxP(int i, int j);
    // solving the linear system at each time step
    void constructMatrixA_u();
    void constructMatrixA_v();
    void constructMatrixA_p();
    void constructMatrixA_temp();
    void constructMatrixA_streamFunc();

    void constructLoadVecU();
    void constructLoadVecV();
    void constructLoadVecP();
    void constructLoadVecTemp();

    void constructVorticityVec();

    void solveForTemp_Next();
    void solveForU_Star(); // step 1 update fieldVecU
    void solveForV_Star(); // step 1 update fieldVecV
    void solveForP_Next(); // step 2 update fieldVecP
    void solveForU_Next(); // step 3 update fieldVecU again
    void solveForV_Next(); // step 3 update fieldVecV again

    void solveForStreamFunc();

    void checkIfContinue(double time);

    void writeDataToFiles(double time);

    // ACCESSING PRIVATE MEMBERS
    const double tStop();
    const double dt();
    const double okToContinue();
    // END OF METHODS $$$$$$$$$$$$$$$$$$$$$$$$$$$
    // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    //

private:
    //
    // ---------------------------------------------
    // CONSTANT ATTRIBUTES ----------------------
    double m_tol = 5e-6;
    double m_a; // lid flow speed
    // space discretisation parameters
    int m_Nx;
    int m_Ny;
    double m_x0;
    double m_x1;
    double m_y0;
    double m_y1;
    double m_Lx;
    double m_Ly;
    double m_dx;
    double m_dy;
    // temporal parameters
    double m_dt;
    double m_cfl;
    double m_tStart = 0.0;
    double m_tStop;
    bool m_okToContinue = true;
    // parameters introduced in heat transfer
    double m_rho_0;
    double m_temp_0;
    double m_tempDiff;
    double m_Cp;
    // linear approximation coefficients
    double m_mu;
    double m_k;

    std::string m_repoDir;
    std::string m_simulationName;
    std::ofstream m_uResults;
    std::ofstream m_vResults;
    std::ofstream m_pResults;
    std::ofstream m_tempResults;
    std::ofstream m_vecResults;
    std::ofstream m_streamFuncResults;
    // END CONSTANT ATTRIBUTES ------------------
    // ---------------------------------------------
    //

    //
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // COMPUTATIONAL DOMAIN %%%%%%%%%%%%%%%%%%%%%
    // sparse matrices at each time step
    Eigen::SparseMatrix<double> m_sMatrixA_u;
    Eigen::SparseMatrix<double> m_sMatrixA_v;
    Eigen::SparseMatrix<double> m_sMatrixA_p;
    Eigen::SparseMatrix<double> m_sMatrixA_temp;
    Eigen::SparseMatrix<double> m_sMatrixA_streamFunc;
    // sparse matrix solvers
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > m_solverU;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > m_solverV;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > m_solverP;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > m_solverTemp;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > m_solverStreamFunc;

    // field vectors at each step
    Eigen::VectorXd m_uVec;
    Eigen::VectorXd m_vVec;
    Eigen::VectorXd m_pVec;
    Eigen::VectorXd m_tempVec;
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
    Eigen::VectorXd m_tempLoadVec;
    // END OF COMPUTATIONAL DOMAIN %%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //
};

#endif