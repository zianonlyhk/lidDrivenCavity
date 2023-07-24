/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   unsteady_solver.cc                                Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/04/10 15:22:53 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "unsteady_solver.hh"

// #define GRAVITY -9.81
#define GRAVITY -0.981 // for smaller Lx and Ly

UnsteadySolver::UnsteadySolver() {}
UnsteadySolver::UnsteadySolver(double i_x0, double i_x1, double i_y0, double i_y1, int i_n_x, int i_n_y, double i_tStop, double i_cfl)
{
    m_x0 = i_x0;
    m_x1 = i_x1;
    m_y0 = i_y0;
    m_y1 = i_y1;
    m_Nx = i_n_x;
    m_Ny = i_n_y;
    m_tStop = i_tStop;
    m_cfl = i_cfl;

    m_Lx = (m_x1 - m_x0);
    m_Ly = (m_y1 - m_y0);
    m_dx = m_Lx / m_Nx;
    m_dy = m_Ly / m_Ny;

    m_uLoadVec = Eigen::VectorXd::Zero((m_Nx - 1) * m_Ny);
    m_vLoadVec = Eigen::VectorXd::Zero(m_Nx * (m_Ny - 1));
    m_pLoadVec = Eigen::VectorXd::Zero(m_Nx * m_Ny);
    m_tempLoadVec = Eigen::VectorXd::Zero(m_Nx * m_Ny);

    m_vortVec = Eigen::VectorXd::Zero((m_Nx - 1) * (m_Ny - 1));
    m_streamFuncVec = Eigen::VectorXd::Zero((m_Nx - 1) * (m_Ny - 1));
}

void UnsteadySolver::setPhysicalParameters(double i_rho_0, double i_temp_0, double i_tempDiff, double i_c_p, double i_mu, double i_k, double i_a)
{
    m_rho_0 = i_rho_0;
    m_temp_0 = i_temp_0;
    m_tempDiff = i_tempDiff;
    m_Cp = i_c_p;
    m_mu = i_mu;
    m_k = i_k;
    m_a = i_a;

    // arbitrary choice of delta time
    // m_dt = 100 * std::min(m_dx, m_dy) * m_mu / (m_rho_0 * m_a * std::max(m_Lx, m_Ly));
    m_dt = 0.01;
}

void UnsteadySolver::setCompDomainVec(const Eigen::VectorXd &uVec, const Eigen::VectorXd &vVec, const Eigen::VectorXd &pVec, const Eigen::VectorXd &tempVec)
{
    m_uVec = uVec;
    m_uVecPrev = uVec;

    m_vVec = vVec;
    m_vVecPrev = vVec;

    m_pVec = pVec;

    m_tempVec = tempVec;
}

void UnsteadySolver::pressureShift()
{
    for (int i = 0; i < m_Nx; ++i)
    {
        for (int j = 0; j < m_Ny; ++j)
        {
            // comsol page on "The Boussinesq Approximation"
            m_pVec(getVecIdxP(i, j)) = m_pVec(getVecIdxP(i, j)) + m_rho_0 * GRAVITY * j * m_dy;
        }
    }
}

void UnsteadySolver::setOutputDataAttributes(std::string simulationName, std::string repoDir)
{
    m_simulationName = simulationName;
    m_repoDir = repoDir;

    m_uResults.open(simulationName + (std::string) "data/" + repoDir + (std::string) "_uResults.dat");
    m_vResults.open(simulationName + (std::string) "data/" + repoDir + (std::string) "_vResults.dat");
    m_pResults.open(simulationName + (std::string) "data/" + repoDir + (std::string) "_pResults.dat");
    m_tempResults.open(simulationName + (std::string) "data/" + repoDir + (std::string) "_tempResults.dat");
    m_vecResults.open(simulationName + (std::string) "data/" + repoDir + (std::string) "_vecResults.dat");
    m_streamFuncResults.open(simulationName + (std::string) "data/" + repoDir + (std::string) "_streamFuncResults.dat");
}

int UnsteadySolver::getVecIdxU(int i, int j)
{
    return i + j * (m_Nx - 1);
}

int UnsteadySolver::getVecIdxV(int i, int j)
{
    return i * (m_Ny - 1) + j;
}

int UnsteadySolver::getVecIdxP(int i, int j)
{
    return i + j * m_Nx;
}

void UnsteadySolver::constructMatrixA_u()
{
    double alpha = m_mu / m_rho_0 * m_dt / m_dx / m_dy;

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m_Nx - 1, m_Nx - 1);
    B.diagonal() = -4 * Eigen::VectorXd::Ones(m_Nx - 1);
    B.diagonal(1) = Eigen::VectorXd::Ones(m_Nx - 2);
    B.diagonal(-1) = Eigen::VectorXd::Ones(m_Nx - 2);

    Eigen::MatrixXd A_diag = Eigen::MatrixXd::Zero(m_Nx - 1, m_Nx - 1);
    A_diag.diagonal() = Eigen::VectorXd::Ones(m_Nx - 1);
    A_diag = A_diag - alpha * B;

    Eigen::MatrixXd A_viceDiag = Eigen::MatrixXd::Zero(m_Nx - 1, m_Nx - 1);
    A_viceDiag.diagonal() = -alpha * Eigen::VectorXd::Ones(m_Nx - 1);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero((m_Nx - 1) * m_Ny, (m_Nx - 1) * m_Ny);

    for (int i = 0; i < m_Ny; i++)
    {
        A.block(i * (m_Nx - 1), i * (m_Nx - 1), (m_Nx - 1), (m_Nx - 1)) = A_diag;
        if (i > 0)
            A.block((i - 1) * (m_Nx - 1), i * (m_Nx - 1), (m_Nx - 1), (m_Nx - 1)) = A_viceDiag;
        if (i < m_Ny - 1)
            A.block((i + 1) * (m_Nx - 1), i * (m_Nx - 1), (m_Nx - 1), (m_Nx - 1)) = A_viceDiag;
    }

    m_sMatrixA_u = A.sparseView();
    m_solverU.compute(m_sMatrixA_u);

    // DEBUG
    // bool db = true;
    bool db = false;
    if (db)
    {
        std::cout << "alpha is " << alpha << std::endl;
        std::cout << A << std::endl;
        std::cout << std::endl;
        std::cout << m_sMatrixA_u << std::endl;
    }
}

void UnsteadySolver::constructMatrixA_v()
{
    double alpha = m_mu / m_rho_0 * m_dt / m_dx / m_dy;

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m_Ny - 1, m_Ny - 1);
    B.diagonal() = -4 * Eigen::VectorXd::Ones(m_Ny - 1);
    B.diagonal(1) = Eigen::VectorXd::Ones(m_Ny - 2);
    B.diagonal(-1) = Eigen::VectorXd::Ones(m_Ny - 2);

    Eigen::MatrixXd A_diag = Eigen::MatrixXd::Zero(m_Ny - 1, m_Ny - 1);
    A_diag.diagonal() = Eigen::VectorXd::Ones(m_Ny - 1);
    A_diag = A_diag - alpha * B;

    Eigen::MatrixXd A_viceDiag = Eigen::MatrixXd::Zero(m_Ny - 1, m_Ny - 1);
    A_viceDiag.diagonal() = -alpha * Eigen::VectorXd::Ones(m_Ny - 1);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m_Nx * (m_Ny - 1), m_Nx * (m_Ny - 1));

    for (int i = 0; i < m_Nx; i++)
    {
        A.block(i * (m_Ny - 1), i * (m_Ny - 1), (m_Ny - 1), (m_Ny - 1)) = A_diag;
        if (i > 0)
            A.block((i - 1) * (m_Ny - 1), i * (m_Ny - 1), (m_Ny - 1), (m_Ny - 1)) = A_viceDiag;
        if (i < m_Nx - 1)
            A.block((i + 1) * (m_Ny - 1), i * (m_Ny - 1), (m_Ny - 1), (m_Ny - 1)) = A_viceDiag;
    }

    m_sMatrixA_v = A.sparseView();
    m_solverV.compute(m_sMatrixA_v);

    // DEBUG
    // bool db = true;
    bool db = false;
    if (db)
    {
        std::cout << "alpha is " << alpha << std::endl;
        std::cout << A << std::endl;
        std::cout << std::endl;
        std::cout << m_sMatrixA_v << std::endl;
    }
}

void UnsteadySolver::constructMatrixA_p()
{
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m_Nx, m_Nx);
    B.diagonal() = -4 * Eigen::VectorXd::Ones(m_Nx);
    // repeating boundary conditions so possible, for LR
    B(0, 0) = -3;
    B(m_Nx - 1, m_Nx - 1) = -3;
    B.diagonal(1) = Eigen::VectorXd::Ones(m_Nx - 1);
    B.diagonal(-1) = Eigen::VectorXd::Ones(m_Nx - 1);

    Eigen::MatrixXd A_viceDiag = Eigen::MatrixXd::Zero(m_Nx, m_Nx);
    A_viceDiag.diagonal() = Eigen::VectorXd::Ones(m_Nx);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m_Nx * m_Ny, m_Nx * m_Ny);

    // repeating bc for the top and bottom included here:
    for (int i = 0; i < m_Ny; i++)
    {
        A.block(i * m_Nx, i * m_Nx, m_Nx, m_Nx) = B;
        if (i > 0)
        {
            A.block((i - 1) * m_Nx, i * m_Nx, m_Nx, m_Nx) = A_viceDiag;
        }
        else
        {
            A.block(i * m_Nx, i * m_Nx, m_Nx, m_Nx) += A_viceDiag;
        }
        if (i < m_Ny - 1)
        {
            A.block((i + 1) * m_Nx, i * m_Nx, m_Nx, m_Nx) = A_viceDiag;
        }
        else
        {
            A.block(i * m_Nx, i * m_Nx, m_Nx, m_Nx) += A_viceDiag;
        }
    }

    m_sMatrixA_p = A.sparseView();
    m_solverP.compute(m_sMatrixA_p);

    // DEBUG
    // bool db = true;
    bool db = false;
    if (db)
    {
        std::cout << A << std::endl;
        std::cout << std::endl;
        std::cout << m_sMatrixA_p << std::endl;
    }
}

void UnsteadySolver::constructMatrixA_temp()
{
    // look at the comsol page "Heat Transfer: Conservation of Energy"
    double alpha = m_dt / m_dx / m_dy * m_k / m_rho_0 / m_Cp;

    // notes here: built-in bc on the left right boundary absent here
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m_Nx, m_Nx);
    B.diagonal() = -4 * Eigen::VectorXd::Ones(m_Nx);
    B.diagonal(1) = Eigen::VectorXd::Ones(m_Nx - 1);
    B.diagonal(-1) = Eigen::VectorXd::Ones(m_Nx - 1);

    Eigen::MatrixXd A_diag = Eigen::MatrixXd::Zero(m_Nx, m_Nx);
    A_diag.diagonal() = Eigen::VectorXd::Ones(m_Nx);
    A_diag = A_diag - alpha * B;

    Eigen::MatrixXd A_viceDiag = Eigen::MatrixXd::Zero(m_Nx, m_Nx);
    A_viceDiag.diagonal() = -alpha * Eigen::VectorXd::Ones(m_Nx);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m_Nx * m_Ny, m_Nx * m_Ny);

    // notes here: transmissive bc on the top and bottom
    for (int i = 0; i < m_Ny; i++)
    {
        A.block(i * m_Nx, i * m_Nx, m_Nx, m_Nx) = A_diag;
        if (i > 0)
        {
            A.block((i - 1) * m_Nx, i * m_Nx, m_Nx, m_Nx) = A_viceDiag;
        }
        else
        {
            A.block(i * m_Nx, i * m_Nx, m_Nx, m_Nx) += A_viceDiag;
        }
        if (i < m_Ny - 1)
        {
            A.block((i + 1) * m_Nx, i * m_Nx, m_Nx, m_Nx) = A_viceDiag;
        }
        else
        {
            A.block(i * m_Nx, i * m_Nx, m_Nx, m_Nx) += A_viceDiag;
        }
    }

    m_sMatrixA_temp = A.sparseView();
    m_solverTemp.compute(m_sMatrixA_temp);

    // DEBUG
    // bool db = true;
    bool db = false;
    if (db)
    {
        std::cout << "alpha is " << alpha << std::endl;
        std::cout << A << std::endl;
        std::cout << std::endl;
        std::cout << m_sMatrixA_temp << std::endl;
    }
}

void UnsteadySolver::constructMatrixA_streamFunc()
{
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m_Nx - 1, m_Nx - 1);
    B.diagonal() = -4 * Eigen::VectorXd::Ones(m_Nx - 1);
    B.diagonal(1) = Eigen::VectorXd::Ones(m_Nx - 2);
    B.diagonal(-1) = Eigen::VectorXd::Ones(m_Nx - 2);

    Eigen::MatrixXd A_viceDiag = Eigen::MatrixXd::Zero(m_Nx - 1, m_Nx - 1);
    A_viceDiag.diagonal() = Eigen::VectorXd::Ones(m_Nx - 1);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero((m_Nx - 1) * (m_Ny - 1), (m_Nx - 1) * (m_Ny - 1));

    for (int i = 0; i < m_Ny - 1; i++)
    {
        A.block(i * (m_Nx - 1), i * (m_Nx - 1), (m_Nx - 1), (m_Nx - 1)) = B;
        if (i > 0)
        {
            A.block((i - 1) * (m_Nx - 1), i * (m_Nx - 1), (m_Nx - 1), (m_Nx - 1)) = A_viceDiag;
        }
        if (i < m_Ny - 2)
        {
            A.block((i + 1) * (m_Nx - 1), i * (m_Nx - 1), (m_Nx - 1), (m_Nx - 1)) = A_viceDiag;
        }
    }

    m_sMatrixA_streamFunc = A.sparseView();
    m_solverStreamFunc.compute(m_sMatrixA_streamFunc);
}

void UnsteadySolver::constructLoadVecTemp()
{
    // define here to add the missing bc terms in its A matrix
    double alpha = m_dt / m_dx / m_dy * m_k / m_rho_0 / m_Cp;

    double centre, north, south, east, west;
    double uLeft, uRight, vUp, vDown;

    for (int i = 0; i < m_Nx; ++i)
    {
        for (int j = 0; j < m_Ny; ++j)
        {
            centre = m_tempVec(getVecIdxP(i, j));

            // Direction ###
            if (j < m_Ny - 1)
            {
                north = m_tempVec(getVecIdxP(i, j + 1));
            }
            else // top boundary conditions
            {
                north = centre;
            }

            if (j > 0)
            {
                south = m_tempVec(getVecIdxP(i, j - 1));
            }
            else // bottom
            {
                south = centre;
            }

            if (i < m_Nx - 1)
            {
                east = m_tempVec(getVecIdxP(i + 1, j));
            }
            else // right
            {
                east = m_temp_0 - m_tempDiff; // right bc
            }

            if (i > 0)
            {
                west = m_tempVec(getVecIdxP(i - 1, j));
            }
            else // left
            {
                west = m_temp_0 + m_tempDiff; // left bc
            }

            // u components ###
            if (i < m_Nx - 1 && i > 0) // in the middle
            {
                uLeft = m_uVec(getVecIdxU(i - 1, j));
                uRight = m_uVec(getVecIdxU(i, j));
            }
            else if (i == m_Nx - 1) // right
            {
                uLeft = m_uVec(getVecIdxU(i - 1, j));
                uRight = 0.0;
            }
            else // left
            {
                uLeft = 0.0;
                uRight = m_uVec(getVecIdxU(i, j));
            }

            // v components ###
            if (j < m_Ny - 1 && j > 0) // in the middle
            {
                vDown = m_vVec(getVecIdxV(i, j - 1));
                vUp = m_vVec(getVecIdxV(i, j));
            }
            else if (j == m_Ny - 1) // up
            {
                vDown = m_vVec(getVecIdxV(i, j - 1));
                vUp = 0.0;
            }
            else // down
            {
                vDown = 0.0;
                vUp = m_vVec(getVecIdxV(i, j));
            }

            // compromise with the average
            m_tempLoadVec(getVecIdxP(i, j)) = centre - m_dt / ((m_dx + m_dy) / 2) / 2 * (0.5 * (uLeft + uRight) * (east - west) + 0.5 * (vUp + vDown) * (north - south));

            if (i == 0) // left
            {
                m_tempLoadVec(getVecIdxP(i, j)) += alpha * (m_temp_0 + m_tempDiff);
            }
            if (i == m_Nx - 1) // right
            {
                m_tempLoadVec(getVecIdxP(i, j)) += alpha * (m_temp_0 - m_tempDiff);
            }
        }
    }

    // DEBUG
    // bool db = true;
    bool db = false;
    if (db)
    {
        std::cout << m_tempLoadVec << std::endl;
    }
}

void UnsteadySolver::constructLoadVecU()
{
    double centre, north, south, east, west;
    double nokia1, nokia3, nokia7, nokia9;

    for (int j = 0; j < m_Ny; ++j)
    {
        for (int i = 0; i < m_Nx - 1; ++i)
        {
            centre = m_uVec(getVecIdxU(i, j));

            // Direction ###
            if (j < m_Ny - 1)
            {
                north = m_uVec(getVecIdxU(i, j + 1));
            }
            else // top
            {
                north = 2 * m_a - centre;
            }

            if (j > 0)
            {
                south = m_uVec(getVecIdxU(i, j - 1));
            }
            else // bottom
            {
                south = -centre;
            }

            if (i < m_Nx - 2)
            {
                east = m_uVec(getVecIdxU(i + 1, j));
            }
            else // right
            {
                east = 0.0;
            }

            if (i > 0)
            {
                west = m_uVec(getVecIdxU(i - 1, j));
            }
            else // left
            {
                west = 0.0;
            }

            // Nokia ###
            if (j < m_Ny - 1 && j > 0) // middle
            {
                nokia1 = m_vVec(getVecIdxV(i, j));
                nokia3 = m_vVec(getVecIdxV(i + 1, j));
                nokia7 = m_vVec(getVecIdxV(i, j - 1));
                nokia9 = m_vVec(getVecIdxV(i + 1, j - 1));
            }
            else if (j == m_Ny - 1) // top
            {
                nokia1 = 0.0;
                nokia3 = 0.0;
                nokia7 = m_vVec(getVecIdxV(i, j - 1));
                nokia9 = m_vVec(getVecIdxV(i + 1, j - 1));
            }
            else // bottom
            {
                nokia1 = m_vVec(getVecIdxV(i, j));
                nokia3 = m_vVec(getVecIdxV(i + 1, j));
                nokia7 = 0.0;
                nokia9 = 0.0;
            }

            m_uLoadVec(getVecIdxU(i, j)) = centre - m_dt / ((m_dx + m_dy) / 2) / 2 * (centre * (east - west) + (north - south) / 4 * (nokia1 + nokia3 + nokia7 + nokia9));

            // left right omitted since they are zero due to staggered grid
            if (j == m_Ny - 1) // top
            {
                m_uLoadVec(getVecIdxU(i, j)) += m_dt * m_mu / m_rho_0 / m_dx / m_dy * north;
            }
            else if (j == 0) // bottom
            {
                m_uLoadVec(getVecIdxU(i, j)) += m_dt * m_mu / m_rho_0 / m_dx / m_dy * south;
            }
        }
    }

    // DEBUG
    // bool db = true;
    bool db = false;
    if (db)
    {
        std::cout << m_uLoadVec << std::endl;
    }
}

void UnsteadySolver::constructLoadVecV()
{
    double centre, north, south, east, west;
    double nokia1, nokia3, nokia7, nokia9;

    for (int i = 0; i < m_Nx; ++i)
    {
        for (int j = 0; j < m_Ny - 1; ++j)
        {
            centre = m_vVec(getVecIdxV(i, j));

            // Direction ###
            if (j < m_Ny - 2)
            {
                north = m_vVec(getVecIdxV(i, j + 1));
            }
            else // top
            {
                north = 0.0;
            }

            if (j > 0)
            {
                south = m_vVec(getVecIdxV(i, j - 1));
            }
            else // bottom
            {
                south = 0.0;
            }

            if (i < m_Nx - 1)
            {
                east = m_vVec(getVecIdxV(i + 1, j));
            }
            else // right
            {
                east = -centre;
            }

            if (i > 0)
            {
                west = m_vVec(getVecIdxV(i - 1, j));
            }
            else // left
            {
                west = -centre;
            }

            // Nokia ###
            if (i < m_Nx - 1 && i > 0) // middle
            {
                nokia1 = m_uVec(getVecIdxU(i - 1, j + 1));
                nokia3 = m_uVec(getVecIdxU(i, j + 1));
                nokia7 = m_uVec(getVecIdxU(i - 1, j));
                nokia9 = m_uVec(getVecIdxU(i, j));
            }
            else if (i == m_Nx - 1) // right
            {
                nokia1 = m_uVec(getVecIdxU(i - 1, j + 1));
                nokia3 = 0.0;
                nokia7 = m_uVec(getVecIdxU(i - 1, j));
                nokia9 = 0.0;
            }
            else // left
            {
                nokia1 = 0.0;
                nokia3 = m_uVec(getVecIdxU(i, j + 1));
                nokia7 = 0.0;
                nokia9 = m_uVec(getVecIdxU(i, j));
            }

            // m_vLoadVec(getVecIdxV(i, j)) = centre - m_dt / ((m_dx + m_dy) / 2) / 2 * (centre * (north - south) + (east - west) / 4 * (nokia1 + nokia3 + nokia7 + nokia9)) - m_dt * (m_tempVec(getVecIdxP(i, j)) - m_temp_0) / m_temp_0 * GRAVITY;
            // m_vLoadVec(getVecIdxV(i, j)) = centre - m_dt / ((m_dx + m_dy) / 2) / 2 * (centre * (north - south) + (east - west) / 4 * (nokia1 + nokia3 + nokia7 + nokia9)) + m_dt / m_dx / m_dy * (GRAVITY - (0.5 * (m_tempVec(getVecIdxP(i, j)) + m_tempVec(getVecIdxP(i, j + 1))) - m_temp_0) / m_temp_0 * GRAVITY);
            m_vLoadVec(getVecIdxV(i, j)) = centre - m_dt / ((m_dx + m_dy) / 2) / 2 * (centre * (north - south) + (east - west) / 4 * (nokia1 + nokia3 + nokia7 + nokia9));

            // top bottom omitted since they are zero due to staggered grid
            if (i == m_Nx - 1) // right
            {
                m_vLoadVec(getVecIdxV(i, j)) += m_dt * m_mu / m_rho_0 / m_dx / m_dy * east;
            }
            else if (i == 0) // left
            {
                m_vLoadVec(getVecIdxV(i, j)) += m_dt * m_mu / m_rho_0 / m_dx / m_dy * west;
            }
        }
    }

    // DEBUG
    // bool db = true;
    bool db = false;
    if (db)
    {
        std::cout << m_vLoadVec << std::endl;
    }
}

void UnsteadySolver::solveForTemp_Next()
{
    m_tempVec = m_solverTemp.solve(m_tempLoadVec);

    // DEBUG
    // bool db = true;
    bool db = false;
    if (db)
    {
        std::cout << m_tempVec << std::endl;
    }
}

void UnsteadySolver::solveForU_Star()
{
    m_uVecPrev = m_uVec;
    m_uVec = m_solverU.solve(m_uLoadVec);

    // DEBUG
    // bool db = true;
    bool db = false;
    if (db)
    {
        std::cout << m_uVec << std::endl;
    }
}

void UnsteadySolver::solveForV_Star()
{
    m_vVecPrev = m_vVec;
    m_vVec = m_solverV.solve(m_vLoadVec);

    // DEBUG
    // bool db = true;
    bool db = false;
    if (db)
    {
        std::cout << m_vVec << std::endl;
    }
}

void UnsteadySolver::constructLoadVecP()
{
    double north, south, east, west;

    for (int i = 0; i < m_Nx; ++i)
    {
        for (int j = 0; j < m_Ny; ++j)
        {
            // Direction ###
            if (j < m_Ny - 1)
            {
                north = m_vVec(getVecIdxV(i, j));
            }
            else // top
            {
                north = 0.0;
            }

            if (j > 0)
            {
                south = m_vVec(getVecIdxV(i, j - 1));
            }
            else // bottom
            {
                south = 0.0;
            }

            if (i < m_Nx - 1)
            {
                east = m_uVec(getVecIdxU(i, j));
            }
            else // right
            {
                east = 0.0;
            }

            if (i > 0)
            {
                west = m_uVec(getVecIdxU(i - 1, j));
            }
            else // left
            {
                west = 0.0;
            }

            m_pLoadVec(getVecIdxP(i, j)) = ((m_dx + m_dy) / 2) / m_dt * (east - west + north - south);
        }
    }

    // DEBUG
    // bool db = true;
    bool db = false;
    if (db)
    {
        std::cout << m_pLoadVec << std::endl;
    }
}

void UnsteadySolver::solveForP_Next()
{
    m_pVec = m_solverP.solve(m_pLoadVec);

    // DEBUG
    // bool db = true;
    bool db = false;
    if (db)
    {
        std::cout << m_pVec << std::endl;
    }
}

void UnsteadySolver::solveForU_Next()
{
    double pGradX;

    for (int j = 0; j < m_Ny; ++j)
    {
        for (int i = 0; i < m_Nx - 1; ++i)
        {
            // Direction ###
            if (i < m_Nx - 2 && i > 0)
            {
                pGradX = (m_pVec(getVecIdxP(i + 1, j)) - m_pVec(getVecIdxP(i, j))) / m_dx;
            }
            else // left and right
            {
                pGradX = 0.0;
            }

            m_uVec(getVecIdxU(i, j)) -= m_dt * pGradX;
        }
    }

    // DEBUG
    // bool db = true;
    bool db = false;
    if (db)
    {
        std::cout << m_uVec << std::endl;
    }
}

void UnsteadySolver::solveForV_Next()
{
    double pGradY;

    for (int i = 0; i < m_Nx; ++i)
    {
        for (int j = 0; j < m_Ny - 1; ++j)
        {
            // Direction ###
            if (j < m_Ny - 2 && j > 0)
            {
                pGradY = (m_pVec(getVecIdxP(i, j + 1)) - m_pVec(getVecIdxP(i, j))) / m_dy;
            }
            else // top and bottom
            {
                pGradY = 0.0;
            }

            m_vVec(getVecIdxV(i, j)) -= m_dt * pGradY;
        }
    }

    // DEBUG
    // bool db = true;
    bool db = false;
    if (db)
    {
        std::cout << m_vVec << std::endl;
    }
}

void UnsteadySolver::constructVorticityVec()
{
    for (int i = 0; i < m_Nx - 1; ++i)
    {
        for (int j = 0; j < m_Ny - 1; ++j)
        {
            m_vortVec((m_Nx - 1) * j + i) = -(m_vVec(getVecIdxV(i + 1, j)) - m_vVec(getVecIdxV(i, j)) - m_uVec(getVecIdxU(i, j + 1)) + m_uVec(getVecIdxU(i, j)));
        }
    }
}

void UnsteadySolver::solveForStreamFunc()
{
    m_streamFuncVec = m_solverStreamFunc.solve(m_vortVec);
}

void UnsteadySolver::checkIfContinue(double time)
{
    m_uVecDiff = m_uVec - m_uVecPrev;
    m_vVecDiff = m_vVec - m_vVecPrev;
    double uContribution = abs(m_uVecDiff.lpNorm<Eigen::Infinity>() / m_uVec.lpNorm<Eigen::Infinity>());
    double vContribution = abs(m_vVecDiff.lpNorm<Eigen::Infinity>() / m_vVec.lpNorm<Eigen::Infinity>());

    if (std::isnan(uContribution + vContribution) && time != 0.0) // initial nan at t=0.0
    {
        std::cout << uContribution + vContribution << " encountered!" << std::endl;
    }

    // if (uContribution + vContribution < m_tol)
    // {
    //     std::cout << "Steady state reached." << std::endl;
    //     m_okToContinue = false;
    // }

    if (time > m_tStop)
    {
        std::cout << "Stop time reached." << std::endl;
        m_okToContinue = false;
    }
}

void UnsteadySolver::writeDataToFiles(double time)
{
    // writing x velocity
    for (int i = 0; i < m_Nx - 1; ++i)
    {
        for (int j = 0; j < m_Ny; ++j)
        {
            m_uResults << time << ' ' << m_x0 + m_dx / 2 + i * m_dx << ' ' << m_y0 + j * m_dx << ' ' << m_uVec(getVecIdxU(i, j)) << std::endl;
        }
        m_uResults << std::endl;
    }
    m_uResults << std::endl;

    // writing y velocity
    for (int j = 0; j < m_Ny - 1; ++j)
    {
        for (int i = 0; i < m_Nx; ++i)
        {
            m_vResults << time << ' ' << m_x0 + i * m_dx << ' ' << m_y0 + m_dx / 2 + j * m_dx << ' ' << m_vVec(getVecIdxV(i, j)) << std::endl;
        }
        m_vResults << std::endl;
    }
    m_vResults << std::endl;

    // writing p
    for (int j = 0; j < m_Ny; ++j)
    {
        for (int i = 0; i < m_Nx; ++i)
        {
            m_pResults << time << ' ' << m_x0 + i * m_dx << ' ' << m_y0 + j * m_dx << ' ' << m_pVec(getVecIdxP(i, j)) << std::endl;
        }
        m_pResults << std::endl;
    }
    m_pResults << std::endl;

    // writing temp
    for (int j = 0; j < m_Ny; ++j)
    {
        for (int i = 0; i < m_Nx; ++i)
        {
            m_tempResults << time << ' ' << m_x0 + i * m_dx << ' ' << m_y0 + j * m_dx << ' ' << m_tempVec(getVecIdxP(i, j)) << std::endl;
        }
        m_tempResults << std::endl;
    }
    m_tempResults << std::endl;

    // writing arrows
    double currX, currY, xVel, yVel, velMag;
    for (int j = 0; j < m_Ny - 1; j += 2)
    {
        for (int i = 0; i < m_Nx - 1; i += 2)
        {
            currX = m_x0 + i * m_dx;
            currY = m_y0 + j * m_dy;
            xVel = m_uVec(getVecIdxU(i, j));
            yVel = m_vVec(getVecIdxV(i, j));
            velMag = sqrt(xVel * xVel + yVel * yVel);
            m_vecResults << time << ' ' << currX << ' ' << currY << ' ' << 1.5 * m_dx * xVel / velMag << ' ' << 1.5 * m_dy * yVel / velMag << std::endl;
        }
        m_vecResults << std::endl;
    }
    m_vecResults << std::endl;

    // writing streamFunc
    for (int j = 0; j < m_Ny - 1; ++j)
    {
        for (int i = 0; i < m_Nx - 1; ++i)
        {
            currX = m_x0 + m_dx / 2 + i * m_dx;
            currY = m_y0 + m_dx / 2 + j * m_dx;
            m_streamFuncResults << time << ' ' << currX << ' ' << currY << ' ' << m_streamFuncVec((m_Nx - 1) * j + i) << std::endl;
        }
        m_streamFuncResults << std::endl;
    }
    m_streamFuncResults << std::endl;

    if (!m_okToContinue)
    {
        m_uResults.close();
        m_vResults.close();
        m_pResults.close();
        m_tempResults.close();
        m_vecResults.close();
        m_streamFuncResults.close();
    }
}

const double UnsteadySolver::tStop() { return m_tStop; }
const double UnsteadySolver::dt() { return m_dt; }
const double UnsteadySolver::okToContinue() { return m_okToContinue; }