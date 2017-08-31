//
// Created by phuschke on 3/16/17.
//

#include <iostream>
//! @brief Print information about the section


#include <eigen3/Eigen/Sparse>

int main()
{

    std::vector<int> vec01(2, 33);
    std::vector<int> vec02(4, 22);

    vec01.insert(vec01.end(), vec02.begin(), vec02.end());


    Eigen::SparseMatrix<double> A(2, 2);
    A.insert(0, 0) = 1;
    A.insert(1, 0) = 1;
    A.insert(0, 1) = 1;
    A.insert(1, 1) = 1;

    Eigen::SparseMatrix<double> B = 2 * A;
    Eigen::SparseMatrix<double> C = 3 * A;
    Eigen::SparseMatrix<double> D = 4 * A;

    auto E = A;
    E.resize(4, 4);

    std::cout << "A \n" << A << std::endl;
    std::cout << "B \n" << B << std::endl;
    std::cout << "C \n" << C << std::endl;
    std::cout << "D \n" << D << std::endl;
    std::cout << "E \n" << E << std::endl;


    int startRowId = 0;
    int startColId = 0;

    for (int rowId = 0; rowId < A.rows(); ++rowId)
        for (int colId = 0; colId < A.cols(); ++colId)
            E.insert(startRowId + rowId, startColId + colId) = A.coeff(rowId, colId);

    startRowId = 0;
    startColId = A.cols();

    for (int rowId = 0; rowId < B.rows(); ++rowId)
        for (int colId = 0; colId < B.cols(); ++colId)
            E.insert(startRowId + rowId, startColId + colId) = B.coeff(rowId, colId);

    startRowId = A.rows();
    startColId = 0;

    for (int rowId = 0; rowId < C.rows(); ++rowId)
        for (int colId = 0; colId < C.cols(); ++colId)
            E.insert(startRowId + rowId, startColId + colId) = C.coeff(rowId, colId);

    startRowId = A.cols();
    startColId = A.rows();

    for (int rowId = 0; rowId < D.rows(); ++rowId)
        for (int colId = 0; colId < D.cols(); ++colId)
            E.insert(startRowId + rowId, startColId + colId) = D.coeff(rowId, colId);

    std::cout << "E \n" << E << std::endl;

    E.resize(4, 4);
    std::cout << "E \n" << E << std::endl;

    C.resize(2, 2);
    C.insert(1, 0) = 443;
    C.insert(1, 1) = 21;
    std::cout << "C.nonZeros() \n" << C.nonZeros() << std::endl;
    std::cout << "C.innerVec \n" << C.innerSize() << std::endl;
    std::cout << "C.innerVec \n" << C.outerSize() << std::endl;
    std::cout << "C.innerVec \n" << C.outerIndexPtr()[0] << std::endl;
    std::cout << "C.innerVec \n" << C.outerIndexPtr()[1] << std::endl;
    std::cout << "C.innerVec \n" << C.outerIndexPtr()[2] << std::endl;

    std::cout << "C.innerVec \n" << C.outerIndexPtr()[2] << std::endl;


    Eigen::SparseMatrix<double> mat(2, 2);
    mat = C;
    for (int k = 0; k < mat.outerSize(); ++k)
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
        {
            long int row = it.row();
            long int col = it.col();
            E.insert(row, col) = it.value();
            std::cout << "it.value()" << it.value() << std::endl;
            std::cout << "it.row()" << it.row() << std::endl;
            std::cout << "it.col()" << it.col() << std::endl;
        }


    std::cout << "E \n" << E << std::endl;
}
