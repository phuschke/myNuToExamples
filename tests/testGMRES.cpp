//
// Created by phuschke on 3/16/17.
//

#include <iostream>

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/QR>






int main() {

    int dim = 3;
    Eigen::MatrixXd A(dim, dim);
    A.setRandom();

    Eigen::DiagonalPreconditioner<double> precond(A);



}
