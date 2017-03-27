#pragma once

#include <iostream>
#include <map>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <fstream>

namespace NuTo
{

void WriteSimulationParameters(const std::map<std::string, std::string>& parameters, const std::string& fileName)
{
    std::cout << "Writing parameters to file: " << fileName + "parameters.txt" << std::endl;
    std::ofstream file(fileName + std::string("parameters.txt"));

    if (not file.is_open())
    {
        std::cout << "cant open file: " << fileName << std::endl;
    }

    for (const auto& pair : parameters)
        file << pair.first << "\n" << pair.second << "\n\n";

    file.close();
}


bool BiCGStab(const Eigen::MatrixXd& mat, const Eigen::VectorXd& rhs, Eigen::VectorXd& x,
              Eigen::FullPivLU<Eigen::MatrixXd>& precond, double tol_error, int& iters)
{

    using std::sqrt;
    using std::abs;
    using RealScalar = double;
    using Scalar     = double;
    using VectorType = Eigen::VectorXd;
    RealScalar tol   = tol_error;
    using Index      = int;
    Index maxIters   = iters;

    Index n       = mat.cols();
    VectorType r  = rhs - mat * x;
    VectorType r0 = r;

    RealScalar r0_sqnorm  = r0.squaredNorm();
    RealScalar rhs_sqnorm = rhs.squaredNorm();
    if (rhs_sqnorm == 0)
    {
        x.setZero();
        return true;
    }
    Scalar rho   = 1;
    Scalar alpha = 1;
    Scalar w     = 1;

    VectorType v = VectorType::Zero(n), p = VectorType::Zero(n);
    VectorType y(n), z(n);
    VectorType kt(n), ks(n);

    VectorType s(n), t(n);

    RealScalar tol2 = tol * tol * rhs_sqnorm;
    RealScalar eps2 = Eigen::NumTraits<Scalar>::epsilon() * Eigen::NumTraits<Scalar>::epsilon();
    Index i         = 0;
    Index restarts  = 0;

    while (r.squaredNorm() > tol2 && i < maxIters)
    {
        Scalar rho_old = rho;

        rho = r0.dot(r);
        if (abs(rho) < eps2 * r0_sqnorm)
        {
            // The new residual vector became too orthogonal to the arbitrarily chosen direction r0
            // Let's restart with a new r0:
            r   = rhs - mat * x;
            r0  = r;
            rho = r0_sqnorm = r.squaredNorm();
            if (restarts++ == 0)
                i = 0;
        }
        Scalar beta = (rho / rho_old) * (alpha / w);
        p           = r + beta * (p - w * v);

        y = precond.solve(p);

        v.noalias() = mat * y;

        alpha = rho / r0.dot(v);
        s     = r - alpha * v;

        z           = precond.solve(s);
        t.noalias() = mat * z;

        RealScalar tmp = t.squaredNorm();
        if (tmp > RealScalar(0))
            w = t.dot(s) / tmp;
        else
            w = Scalar(0);
        x += alpha * y + w * z;
        r = s - w * t;
        ++i;
    }
    tol_error = sqrt(r.squaredNorm() / rhs_sqnorm);
    iters     = i;
    return true;
}



} // namespace NuTo
