//
// Created by haitan on 09/03/2020.
//

#include <iostream>
#include <vector>
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

    std::vector<std::string> myvector = {"POD","DMD","RDMD"};
    std::vector<std::string>::iterator it;

    it = find (myvector.begin(), myvector.end(), "RDMD");
    if (it != myvector.end())
        std::cout << "Element found in myvector: " << *it << '\n';
    else
        std::cout << "Element not found in myvector\n";

    Eigen::VectorXd bb(3);
    bb(0)=3.2;
    bb(1)=3.2;
    bb(2)=3.2;
    int index;
    double minval = bb.minCoeff(&index);

    std::cout << "Min val in bb " << minval << " at " << index << std::endl;

    return 0;
}




