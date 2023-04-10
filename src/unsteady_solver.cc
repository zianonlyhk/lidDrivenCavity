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
}

void UnsteadySolver::setCompDomainVec(const Eigen::VectorXd &uVec, const Eigen::VectorXd &vVec, const Eigen::VectorXd &pVec)
{
    m_uVec = uVec;
    m_uVecPrev = uVec;

    m_vVec = vVec;
    m_vVecPrev = vVec;

    m_pVec = pVec;
    m_pVecPrev = pVec;
}

void UnsteadySolver::updateDt()
{
    double maxSpeed = std::max(m_uVec.lpNorm<Eigen::Infinity>(), m_vVec.lpNorm<Eigen::Infinity>());
    double diffusionPart = m_dx / maxSpeed;

    double convectionPart = m_re * m_dx * m_dx;

    m_dt = m_cfl * std::min(diffusionPart, convectionPart);
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

void UnsteadySolver::constructMatrixA_uv(double delTime)
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
    B.diagonal(1) = Eigen::VectorXd::Ones(m_N - 1);
    B.diagonal(-1) = Eigen::VectorXd::Ones(m_N - 1);

    Eigen::MatrixXd A_viceDiag = Eigen::MatrixXd::Zero(m_N, m_N);
    A_viceDiag.diagonal() = Eigen::VectorXd::Ones(m_N);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m_N * m_N, m_N * m_N);

    for (int i = 0; i < m_N; i++)
    {
        A.block(i * m_N, i * m_N, m_N, m_N) = B;
        if (i > 0)
            A.block((i - 1) * m_N, i * m_N, m_N, m_N) = A_viceDiag;
        if (i < m_N - 1)
            A.block((i + 1) * m_N, i * m_N, m_N, m_N) = A_viceDiag;
    }

    m_sMatrixA_p = A.sparseView();

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

void UnsteadySolver::constructLoadVecU()
{
}

void UnsteadySolver::constructBcVecU()
{
}

void UnsteadySolver::constructLoadVecV()
{
}

void UnsteadySolver::constructBcVecV()
{
}

const double UnsteadySolver::tStop() { return m_tStop; }
const double UnsteadySolver::dt() { return m_dt; }