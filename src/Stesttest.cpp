//
// Created by haitan on 09/03/2020.
//

#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"


int main () {

    std::complex<double> c1(1.0, 1.0);
    std::complex<double> c2(0.0, 1.0);
    std::complex<double> c3(1.0, 0.0);
    std::complex<double> c4(0.5, 0.5);

//    Eigen::VectorXd b(2);
//    b(0) = 1.0; b(1) = 1.0;

    Eigen::MatrixXcd A(2,2);
    A(0,0) = c1;
    A(0,1) = c2;
    A(1,0) = c3;
    A(1,1) = c4;

    Eigen::MatrixXcd b(2,2);
    b(0,0) = 2.0*c3;
    b(0,1) = 0.5*c1;
    b(1,0) = c1;
    b(1,1) = 0.3*c2;


    Eigen::MatrixXcd x = A.colPivHouseholderQr().solve(b);

    std::cout << "Here's matrix A:\n" <<std::endl;
    std::cout << A << std::endl;
    std::cout << "Here's matrix b:\n" <<std::endl;
    std::cout << b << std::endl;
    std::cout << "Here's the solution:\n" <<std::endl;
    std::cout << x << std::endl;
//    Eigen::VectorXcd r = A*b;
//
//    std::cout << "Complex matrix: " << std::endl;
//    std::cout << c << std::endl;
//    std::cout << "Real vector: " << std::endl;
//    std::cout << b << std::endl;
//
//    std::cout << "Matrix complex per vector real multiplication: " << std::endl;
//    std::cout << r << std::endl;

    return 0;
}




