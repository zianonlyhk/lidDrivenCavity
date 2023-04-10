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
UnsteadySolver::UnsteadySolver(int N, double tStop, double cfl, double re)
{
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
}

void UnsteadySolver::constructMatrixA_uv(double delTime)
{
}
