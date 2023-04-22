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

UnsteadySolver::UnsteadySolver() {}
UnsteadySolver::UnsteadySolver(int N, double tStop, double cfl, double re, double a)
{
    m_N = N;
    m_tStop = tStop;
    m_cfl = cfl;
    m_re = re;
    m_a = a;

    m_dx = (m_x1 - m_x0) / m_N;

    m_dt = 100 * m_dx / m_re;
    // m_dt = m_re * m_dx * m_dx;

    m_vortVec = Eigen::VectorXd::Zero((m_N - 1) * (m_N - 1));
    m_streamFuncVec = Eigen::VectorXd::Zero((m_N - 1) * (m_N - 1));

    m_uLoadVec = Eigen::VectorXd::Zero(m_N * (m_N - 1));
    m_vLoadVec = Eigen::VectorXd::Zero(m_N * (m_N - 1));
    m_pLoadVec = Eigen::VectorXd::Zero(m_N * m_N);
    m_tempLoadVec = Eigen::VectorXd::Zero(m_N * m_N);
}
void UnsteadySolver::setTemp(double initTemp, double tempDiff)
{
    m_initTemp = initTemp;
    m_tempDiff = tempDiff;
}

void UnsteadySolver::setCompDomainVec(const Eigen::VectorXd &uVec, const Eigen::VectorXd &vVec, const Eigen::VectorXd &pVec)
{
    m_uVec = uVec;
    m_uVecPrev = uVec;

    m_vVec = vVec;
    m_vVecPrev = vVec;

    m_pVec = pVec;
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
    return j * (m_N - 1) + i;
}

int UnsteadySolver::getVecIdxV(int i, int j)
{
    return i * (m_N - 1) + j;
}

int UnsteadySolver::getVecIdxP(int i, int j)
{
    return j * m_N + i;
}

void UnsteadySolver::constructMatrixA_uv()
{
    double alpha = m_dt / m_re / m_dx / m_dx;

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m_N - 1, m_N - 1);
    B.diagonal() = -4 * Eigen::VectorXd::Ones(m_N - 1);
    B.diagonal(1) = Eigen::VectorXd::Ones(m_N - 2);
    B.diagonal(-1) = Eigen::VectorXd::Ones(m_N - 2);

    Eigen::MatrixXd A_diag = Eigen::MatrixXd::Zero(m_N - 1, m_N - 1);
    A_diag.diagonal() = Eigen::VectorXd::Ones(m_N - 1);
    A_diag = A_diag - alpha * B;

    Eigen::MatrixXd A_viceDiag = Eigen::MatrixXd::Zero(m_N - 1, m_N - 1);
    A_viceDiag.diagonal() = -alpha * Eigen::VectorXd::Ones(m_N - 1);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m_N * (m_N - 1), m_N * (m_N - 1));

    for (int i = 0; i < m_N; i++)
    {
        A.block(i * (m_N - 1), i * (m_N - 1), (m_N - 1), (m_N - 1)) = A_diag;
        if (i > 0)
            A.block((i - 1) * (m_N - 1), i * (m_N - 1), (m_N - 1), (m_N - 1)) = A_viceDiag;
        if (i < m_N - 1)
            A.block((i + 1) * (m_N - 1), i * (m_N - 1), (m_N - 1), (m_N - 1)) = A_viceDiag;
    }

    m_sMatrixA_uv = A.sparseView();
    m_solverUV.compute(m_sMatrixA_uv);

    // DEBUG
    // bool db = true;
    bool db = false;
    if (db)
    {
        std::cout << "alpha is " << alpha << std::endl;
        std::cout << A << std::endl;
        std::cout << std::endl;
        std::cout << m_sMatrixA_uv << std::endl;
    }
}

void UnsteadySolver::constructMatrixA_p()
{
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m_N, m_N);
    B.diagonal() = -4 * Eigen::VectorXd::Ones(m_N);
    B(0, 0) = -3;
    B(m_N - 1, m_N - 1) = -3;
    B.diagonal(1) = Eigen::VectorXd::Ones(m_N - 1);
    B.diagonal(-1) = Eigen::VectorXd::Ones(m_N - 1);

    Eigen::MatrixXd A_viceDiag = Eigen::MatrixXd::Zero(m_N, m_N);
    A_viceDiag.diagonal() = Eigen::VectorXd::Ones(m_N);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m_N * m_N, m_N * m_N);

    for (int i = 0; i < m_N; i++)
    {
        A.block(i * m_N, i * m_N, m_N, m_N) = B;
        if (i > 0)
        {
            A.block((i - 1) * m_N, i * m_N, m_N, m_N) = A_viceDiag;
        }
        else
        {
            A.block(i * m_N, i * m_N, m_N, m_N) += A_viceDiag;
        }
        if (i < m_N - 1)
        {
            A.block((i + 1) * m_N, i * m_N, m_N, m_N) = A_viceDiag;
        }
        else
        {
            A.block(i * m_N, i * m_N, m_N, m_N) += A_viceDiag;
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
    // new coefficient
    double alpha = (m_rho * m_cp * m_lambda / m_a / m_l) * m_dt / m_dx / m_dx;

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m_N, m_N);
    B.diagonal() = -4 * Eigen::VectorXd::Ones(m_N);
    B.diagonal(1) = Eigen::VectorXd::Ones(m_N - 1);
    B.diagonal(-1) = Eigen::VectorXd::Ones(m_N - 1);

    Eigen::MatrixXd A_diag = Eigen::MatrixXd::Zero(m_N, m_N);
    A_diag.diagonal() = Eigen::VectorXd::Ones(m_N);
    A_diag = A_diag - alpha * B;

    Eigen::MatrixXd A_viceDiag = Eigen::MatrixXd::Zero(m_N, m_N);
    A_viceDiag.diagonal() = -alpha * Eigen::VectorXd::Ones(m_N);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m_N * m_N, m_N * m_N);

    for (int i = 0; i < m_N; i++)
    {
        A.block(i * m_N, i * m_N, m_N, m_N) = A_diag;
        if (i > 0)
        {
            A.block((i - 1) * m_N, i * m_N, m_N, m_N) = A_viceDiag;
        }
        else
        {
            A.block(i * m_N, i * m_N, m_N, m_N) += A_viceDiag;
        }
        if (i < m_N - 1)
        {
            A.block((i + 1) * m_N, i * m_N, m_N, m_N) = A_viceDiag;
        }
        else
        {
            A.block(i * m_N, i * m_N, m_N, m_N) += A_viceDiag;
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
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m_N - 1, m_N - 1);
    B.diagonal() = -4 * Eigen::VectorXd::Ones(m_N - 1);
    B.diagonal(1) = Eigen::VectorXd::Ones(m_N - 2);
    B.diagonal(-1) = Eigen::VectorXd::Ones(m_N - 2);

    Eigen::MatrixXd A_viceDiag = Eigen::MatrixXd::Zero(m_N - 1, m_N - 1);
    A_viceDiag.diagonal() = Eigen::VectorXd::Ones(m_N - 1);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero((m_N - 1) * (m_N - 1), (m_N - 1) * (m_N - 1));

    for (int i = 0; i < m_N - 1; i++)
    {
        A.block(i * (m_N - 1), i * (m_N - 1), (m_N - 1), (m_N - 1)) = B;
        if (i > 0)
        {
            A.block((i - 1) * (m_N - 1), i * (m_N - 1), (m_N - 1), (m_N - 1)) = A_viceDiag;
        }
        if (i < m_N - 2)
        {
            A.block((i + 1) * (m_N - 1), i * (m_N - 1), (m_N - 1), (m_N - 1)) = A_viceDiag;
        }
    }

    m_sMatrixA_streamFunc = A.sparseView();
    m_solverStreamFunc.compute(m_sMatrixA_streamFunc);
}

void UnsteadySolver::constructLoadVecTemp()
{
    double alpha = (m_rho * m_cp * m_lambda / m_a / m_l) * m_dt / m_dx / m_dx;

    double centre, north, south, east, west;
    double uLeft, uRight, vUp, vDown;

    for (int i = 0; i < m_N; ++i)
    {
        for (int j = 0; j < m_N; ++j)
        {
            centre = m_tempVec(getVecIdxP(i, j));

            // Direction ###
            if (j < m_N - 1) // if there is anything up
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

            if (i < m_N - 1)
            {
                east = m_tempVec(getVecIdxP(i + 1, j));
            }
            else // right
            {
                east = m_initTemp - m_tempDiff;
            }

            if (i > 0)
            {
                west = m_tempVec(getVecIdxP(i - 1, j));
            }
            else // left
            {
                west = m_initTemp + m_tempDiff;
            }

            // u components ###
            if (i < m_N - 1 && i > 0) // in the middle
            {
                uLeft = m_uVec(getVecIdxU(i - 1, j));
                uRight = m_uVec(getVecIdxU(i, j));
            }
            else if (i == m_N - 1) // right
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
            if (j < m_N - 1 && j > 0) // in the middle
            {
                vDown = m_vVec(getVecIdxV(i, j - 1));
                vUp = m_vVec(getVecIdxV(i, j));
            }
            else if (j == m_N - 1) // up
            {
                vDown = m_vVec(getVecIdxV(i, j - 1));
                vUp = 0.0;
            }
            else // down
            {
                vDown = 0.0;
                vUp = m_vVec(getVecIdxV(i, j));
            }

            m_tempLoadVec(getVecIdxP(i, j)) = centre - m_dt / m_dx / 2 * (0.5 * (uLeft + uRight) * (east - west) + 0.5 * (vUp + vDown) * (north - south));

            if (i == 0) // left
            {
                m_tempLoadVec(getVecIdxP(i, j)) += alpha * (m_initTemp + m_tempDiff);
            }
            if (i == m_N - 1) // right
            {
                m_tempLoadVec(getVecIdxP(i, j)) += alpha * (m_initTemp - m_tempDiff);
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

    for (int j = 0; j < m_N; ++j)
    {
        for (int i = 0; i < m_N - 1; ++i)
        {
            centre = m_uVec(getVecIdxU(i, j));

            // Direction ###
            if (j < m_N - 1)
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

            if (i < m_N - 2)
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
            if (j < m_N - 1 && j > 0) // middle
            {
                nokia1 = m_vVec(getVecIdxV(i, j));
                nokia3 = m_vVec(getVecIdxV(i + 1, j));
                nokia7 = m_vVec(getVecIdxV(i, j - 1));
                nokia9 = m_vVec(getVecIdxV(i + 1, j - 1));
            }
            else if (j == m_N - 1) // top
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

            m_uLoadVec(getVecIdxU(i, j)) = centre - m_dt / m_dx / 2 * (centre * (east - west) + (north - south) / 4 * (nokia1 + nokia3 + nokia7 + nokia9));

            if (j == m_N - 1) // top
            {
                m_uLoadVec(getVecIdxU(i, j)) += m_dt / m_re / m_dx / m_dx * north;
            }
            else if (j == 0) // bottom
            {
                m_uLoadVec(getVecIdxU(i, j)) += m_dt / m_re / m_dx / m_dx * south;
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

    for (int i = 0; i < m_N; ++i)
    {
        for (int j = 0; j < m_N - 1; ++j)
        {
            centre = m_vVec(getVecIdxV(i, j));

            // Direction ###
            if (j < m_N - 2)
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

            if (i < m_N - 1)
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
            if (i < m_N - 1 && i > 0) // middle
            {
                nokia1 = m_uVec(getVecIdxU(i - 1, j + 1));
                nokia3 = m_uVec(getVecIdxU(i, j + 1));
                nokia7 = m_uVec(getVecIdxU(i - 1, j));
                nokia9 = m_uVec(getVecIdxU(i, j));
            }
            else if (i == m_N - 1) // right
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

            m_vLoadVec(getVecIdxV(i, j)) = centre - m_dt / m_dx / 2 * (centre * (north - south) + (east - west) / 4 * (nokia1 + nokia3 + nokia7 + nokia9));

            if (i == m_N - 1) // right
            {
                m_vLoadVec(getVecIdxV(i, j)) += m_dt / m_re / m_dx / m_dx * east;
            }
            else if (i == 0) // left
            {
                m_vLoadVec(getVecIdxV(i, j)) += m_dt / m_re / m_dx / m_dx * west;
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
    m_uVec = m_solverUV.solve(m_uLoadVec);

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
    m_vVec = m_solverUV.solve(m_vLoadVec);

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

    double tempScaling = m_lambda / m_cp / m_rho / m_a / m_l;
    double tempC, tempN, tempS, tempR, tempL;

    for (int i = 0; i < m_N; ++i)
    {
        for (int j = 0; j < m_N; ++j)
        {
            tempC = m_tempVec(getVecIdxP(i, j));

            // Direction ###
            if (j < m_N - 1)
            {
                north = m_vVec(getVecIdxV(i, j));
                tempN = m_tempVec(getVecIdxP(i, j + 1));
            }
            else // top
            {
                north = 0.0;
                tempN = tempC;
            }

            if (j > 0)
            {
                south = m_vVec(getVecIdxV(i, j - 1));
                tempS = m_tempVec(getVecIdxP(i, j - 1));
            }
            else // bottom
            {
                south = 0.0;
                tempS = tempC;
            }

            if (i < m_N - 1)
            {
                east = m_uVec(getVecIdxU(i, j));
                tempR = m_tempVec(getVecIdxP(i + 1, j));
            }
            else // right
            {
                east = 0.0;
                tempR = m_initTemp - m_tempDiff;
            }

            if (i > 0)
            {
                west = m_uVec(getVecIdxU(i - 1, j));
                tempL = m_tempVec(getVecIdxP(i - 1, j));
            }
            else // left
            {
                west = 0.0;
                tempL = m_initTemp + m_tempDiff;
            }

            m_pLoadVec(getVecIdxP(i, j)) = m_dx / m_dt * (east - west + north - south) - 1 / m_dt * tempScaling / tempC * (tempN + tempS + tempR + tempL - 4 * tempC);
            // m_pLoadVec(getVecIdxP(i, j)) = m_dx / m_dt * (east - west + north - south);
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

    for (int j = 0; j < m_N; ++j)
    {
        for (int i = 0; i < m_N - 1; ++i)
        {
            // Direction ###
            if (i < m_N - 2 && i > 0)
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

    for (int i = 0; i < m_N; ++i)
    {
        for (int j = 0; j < m_N - 1; ++j)
        {
            // Direction ###
            if (j < m_N - 2 && j > 0)
            {
                pGradY = (m_pVec(getVecIdxP(i, j + 1)) - m_pVec(getVecIdxP(i, j))) / m_dx;
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
    for (int i = 0; i < m_N - 1; ++i)
    {
        for (int j = 0; j < m_N - 1; ++j)
        {
            m_vortVec((m_N - 1) * j + i) = -(m_vVec(getVecIdxV(i + 1, j)) - m_vVec(getVecIdxV(i, j)) - m_uVec(getVecIdxU(i, j + 1)) + m_uVec(getVecIdxU(i, j)));
        }
    }
}

void UnsteadySolver::solveForStreamFunc()
{
    m_streamFuncVec = m_solverStreamFunc.solve(m_vortVec);
}

void UnsteadySolver::checkIfBreak(double time)
{
    m_uVecDiff = m_uVec - m_uVecPrev;
    m_vVecDiff = m_vVec - m_vVecPrev;
    double uContribution = abs(m_uVecDiff.lpNorm<Eigen::Infinity>() / m_uVec.lpNorm<Eigen::Infinity>());
    double vContribution = abs(m_vVecDiff.lpNorm<Eigen::Infinity>() / m_vVec.lpNorm<Eigen::Infinity>());

    if (std::isnan(uContribution + vContribution))
    {
        std::cout << uContribution + vContribution << " encountered!" << std::endl;
    }

    if (uContribution + vContribution < m_tol)
    {
        m_reachedSteady = true;
    }

    if (time > m_tStop)
    {
        m_reachedSteady = true;
    }
}

void UnsteadySolver::writeDataToFiles(double time)
{
    // writing x velocity
    for (int j = 0; j < m_N; ++j)
    {
        for (int i = 0; i < m_N - 1; ++i)
        {
            m_uResults << time << ' ' << m_x0 + m_dx / 2 + i * m_dx << ' ' << m_y0 + j * m_dx << ' ' << m_uVec(getVecIdxU(i, j)) << std::endl;
        }
        m_uResults << std::endl;
    }
    m_uResults << std::endl;

    // writing y velocity
    for (int j = 0; j < m_N - 1; ++j)
    {
        for (int i = 0; i < m_N; ++i)
        {
            m_vResults << time << ' ' << m_x0 + i * m_dx << ' ' << m_y0 + m_dx / 2 + j * m_dx << ' ' << m_vVec(getVecIdxV(i, j)) << std::endl;
        }
        m_vResults << std::endl;
    }
    m_vResults << std::endl;

    // writing p
    for (int j = 0; j < m_N; ++j)
    {
        for (int i = 0; i < m_N; ++i)
        {
            m_pResults << time << ' ' << m_x0 + i * m_dx << ' ' << m_y0 + j * m_dx << ' ' << m_pVec(getVecIdxP(i, j)) << std::endl;
        }
        m_pResults << std::endl;
    }
    m_pResults << std::endl;

    // writing temp
    for (int j = 0; j < m_N; ++j)
    {
        for (int i = 0; i < m_N; ++i)
        {
            m_tempResults << time << ' ' << m_x0 + i * m_dx << ' ' << m_y0 + j * m_dx << ' ' << m_tempVec(getVecIdxP(i, j)) << std::endl;
        }
        m_tempResults << std::endl;
    }
    m_tempResults << std::endl;

    // writing arrows
    double currX, currY, xVel, yVel, velMag;
    for (int j = 0; j < m_N - 1; j += 2)
    {
        for (int i = 0; i < m_N - 1; i += 2)
        {
            currX = m_x0 + 3 * m_dx / 4 + i * m_dx;
            currY = m_y0 + 3 * m_dx / 4 + j * m_dx;
            xVel = m_uVec(getVecIdxU(i, j));
            yVel = m_vVec(getVecIdxV(i, j));
            velMag = sqrt(xVel * xVel + yVel * yVel);
            m_vecResults << time << ' ' << currX << ' ' << currY << ' ' << 0.015 * xVel / velMag << ' ' << 0.015 * yVel / velMag << std::endl;
        }
        m_vecResults << std::endl;
    }
    m_vecResults << std::endl;

    // writing streamFunc
    for (int j = 0; j < m_N - 1; ++j)
    {
        for (int i = 0; i < m_N - 1; ++i)
        {
            currX = m_x0 + 3 * m_dx / 4 + i * m_dx;
            currY = m_y0 + 3 * m_dx / 4 + j * m_dx;
            m_streamFuncResults << time << ' ' << currX << ' ' << currY << ' ' << m_streamFuncVec((m_N - 1) * j + i) << std::endl;
        }
        m_streamFuncResults << std::endl;
    }
    m_streamFuncResults << std::endl;

    if (m_reachedSteady)
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
const double UnsteadySolver::reachedSteady() { return m_reachedSteady; }