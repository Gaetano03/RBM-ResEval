#include "Pre-Process.hpp"


Eigen::VectorXd IC( prob_settings settings, int nC, int Nr) {
    double M = settings.Mach;
    double Re = settings.Re;
    double T = settings.T;
    double alpha = M_PI*settings.alpha/double(180);
    double beta =  M_PI*settings.beta/double(180);
    if ( settings.ndim == 2 ) beta = 0.0;

    double length = 1.0;
    double R = 287.058;
    double gamma = 1.4;
    double mu_ref = 1.716E-5;
    double T_ref = 273.15;
    double S = 110.4;

    double mu = mu_ref*std::pow(T/T_ref,1.5)*(T_ref + S)/(T + S);
    double V_magn = M*std::sqrt(gamma*R*T);
    double rho = Re*mu/(V_magn*length);
    double rhoU = rho*V_magn*std::cos(alpha)*std::cos(beta);
    double rhoV = rho*V_magn*std::sin(alpha);
    double rhoW = rho*V_magn*std::cos(alpha)*std::sin(beta);
    double rhoE = rho*(R/(gamma-1)*T + 0.5*V_magn*V_magn);

    //That only matters for turbulent calculation
    //since these values are used as default values in SU2, they are not present in config file but hard coded
    double n_turb = 0.05;
    double mu_turb2lam_ratio = 10.0;
    double nuTilde=0.1*mu/rho;

    double tke = 1.5*n_turb*n_turb*V_magn*V_magn;
    double omega = rho*tke/(std::max(mu*mu_turb2lam_ratio,1.e-25));
    // double rhotke = rho*tke;
    // double rhoomega = rho*omega;
    double rhotke = tke;
    double rhoomega = omega;

    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC*Nr);

    if ( nC == 4 ) //Laminar 2D Navier-Stokes
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
    } else if (nC== 5) // Turbolent 2D Spalart Allmaras
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(4*Nr, Nr) = nuTilde*Eigen::MatrixXd::Ones(Nr,1);

    } else if ( nC == 6 ) //Turbulent 2D Navier-Stokes
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(4*Nr, Nr) = rhotke*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(5*Nr, Nr) = rhoomega*Eigen::MatrixXd::Ones(Nr,1);
    } else if ( nC == 7 ) //Turbulent 3D Navier-Stokes
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoW*Eigen::MatrixXd::Ones(Nr,1); //no sideslip angle
        Ic.segment(4*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(5*Nr, Nr) = rhotke*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(6*Nr, Nr) = rhoomega*Eigen::MatrixXd::Ones(Nr,1);
    } else {
        std::cout << "Set well numbe rof conservative Variables" << std::endl;
        exit(EXIT_FAILURE);
    }

    return Ic;
}

