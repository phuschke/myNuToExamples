//
// Created by phuschke on 3/16/17.
//

#include "../AuxFunctions.h"
#include <eigen3/Eigen/Sparse>
int main()
{
    std::cout << "Start PreconditionedBiCGStab" << std::endl;


    Eigen::MatrixXd mat(2, 2);
    mat << 1, 0, 0, 1;

    Eigen::VectorXd rhs(2);
    Eigen::VectorXd x(2);
    rhs << 1, 1;

    Eigen::FullPivLU<Eigen::MatrixXd> precond(200*mat);

    int iters = 10;
    NuTo::BiCGStab(mat, rhs, x, precond, 1.e-6, iters);

    std::cout << "x \n" << x << "\n";
    std::cout << "iterations \n" << iters << "\n";

    Eigen::SparseMatrix<double> sparseMat(2,2);
    sparseMat.insert(0,0) = 1;
    sparseMat.insert(1,0) = 2;
    sparseMat.insert(1,1) = 3;
    sparseMat.insert(0,1) = 4;
    std::cout << sparseMat << std::endl;


    sparseMat.conservativeResize(5,5);
    std::cout << sparseMat << std::endl;
}