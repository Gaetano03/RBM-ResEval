//
// Created by haitan on 05/03/2020.
//

#ifndef PRE_PROCESS_HPP
#define PRE_PROCESS_HPP

#include "read_Inputs.hpp"
#include <cmath>

Eigen::VectorXd IC ( Eigen::MatrixXd &sn_set, prob_settings settings, int nC, int Nr, std::string flag = "NO" );

Eigen::VectorXi Inverse_POS (const Eigen::MatrixXd &sn_set, int Nsamples);

#endif //MODES_PRE_PROCESS_HPP
