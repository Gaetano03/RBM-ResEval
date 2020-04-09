//
// Created by haitan on 05/03/2020.
//

#ifndef PRE_PROCESS_HPP
#define PRE_PROCESS_HPP

#include "read_Inputs.hpp"
#include <cmath>

Eigen::VectorXd IC ( Eigen::MatrixXd &sn_set, prob_settings settings, int nC, int Nr);

void get_MinMax_ConsVar (const Eigen::MatrixXd sn_set, const prob_settings &settings, const int nC, double &rho_max,
                         double &rho_min, double &rhoU_max, double &rhoU_min, double &rhoV_max, double &rhoV_min,
                         double &rhoW_max, double &rhoW_min, double &rhoE_max, double &rhoE_min, double &tke_min,
                         double &tke_max, double &omega_min, double &omega_max, double &nuTilde_min, double &nuTilde_max);

void Direct_Normalization(Eigen::MatrixXd sn_set, const prob_settings &settings, const int nC);

Eigen::VectorXi Inverse_POS (const Eigen::MatrixXd &sn_set, int Nsamples);

#endif //MODES_PRE_PROCESS_HPP
