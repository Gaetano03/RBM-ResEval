#include "Pre-Process.hpp"

Eigen::VectorXd IC ( Eigen::MatrixXd &sn_set, prob_settings settings, int nC, int Nr, std::string flag ) {

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

    double rho_max, rho_min, rhoU_max, rhoU_min, rhoV_max, rhoV_min, rhoW_max, rhoW_min,
        rhoE_max, rhoE_min, tke_min, tke_max, omega_min, omega_max, nu_min, nu_max; //add turbulence

    if ( nC == 2 )
    {
        Ic.head(Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;

    } else if ( nC == 3 )
    {
        Ic.head(Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoW*Eigen::MatrixXd::Ones(Nr,1);
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;

    } else if ( nC == 4 ) //Laminar 2D Navier-Stokes
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;

        if ( flag == "YES" ) {
            //Introduce an if on the number of conservative variables

            rho_max = sn_set.middleRows(0, Nr).maxCoeff();
            rho_min = sn_set.middleRows(0, Nr).minCoeff();
            sn_set.middleRows(0, Nr) = (sn_set.middleRows(0, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rho_min )/(rho_max - rho_min);

            rhoU_max = sn_set.middleRows(Nr, Nr).maxCoeff();
            rhoU_min = sn_set.middleRows(Nr, Nr).minCoeff();
            sn_set.middleRows(Nr, Nr) = (sn_set.middleRows(Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoU_min)/(rhoU_max - rhoU_min);

            rhoV_max = sn_set.middleRows(2*Nr, Nr).maxCoeff();
            rhoV_min = sn_set.middleRows(2*Nr, Nr).minCoeff();
            sn_set.middleRows(2*Nr, Nr) = (sn_set.middleRows(2*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoV_min)/(rhoV_max - rhoV_min);

            rhoE_max = sn_set.middleRows(3*Nr, Nr).maxCoeff();
            rhoE_min = sn_set.middleRows(3*Nr, Nr).minCoeff();
            sn_set.middleRows(3*Nr, Nr) = (sn_set.middleRows(3*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoE_min)/(rhoE_max - rhoE_min);

        }


    } else if ( nC== 5 ) // Turbolent 2D Spalart Allmaras
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(4*Nr, Nr) = nuTilde*Eigen::MatrixXd::Ones(Nr,1);
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;

    } else if ( nC == 6 && settings.ndim == 2 ) //Turbulent 2D Navier-Stokes
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(4*Nr, Nr) = rhotke*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(5*Nr, Nr) = rhoomega*Eigen::MatrixXd::Ones(Nr,1);
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;

        if ( flag == "YES" ) {
            //Introduce an if on the number of conservative variables

            rho_max = sn_set.middleRows(0, Nr).maxCoeff();
            rho_min = sn_set.middleRows(0, Nr).minCoeff();
            sn_set.middleRows(0, Nr) = (sn_set.middleRows(0, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rho_min )/(rho_max - rho_min);

            rhoU_max = sn_set.middleRows(Nr, Nr).maxCoeff();
            rhoU_min = sn_set.middleRows(Nr, Nr).minCoeff();
            sn_set.middleRows(Nr, Nr) = (sn_set.middleRows(Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoU_min)/(rhoU_max - rhoU_min);

            rhoV_max = sn_set.middleRows(2*Nr, Nr).maxCoeff();
            rhoV_min = sn_set.middleRows(2*Nr, Nr).minCoeff();
            sn_set.middleRows(2*Nr, Nr) = (sn_set.middleRows(2*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoV_min)/(rhoV_max - rhoV_min);

            rhoE_max = sn_set.middleRows(3*Nr, Nr).maxCoeff();
            rhoE_min = sn_set.middleRows(3*Nr, Nr).minCoeff();
            sn_set.middleRows(3*Nr, Nr) = (sn_set.middleRows(3*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoE_min)/(rhoE_max - rhoE_min);

            tke_max = sn_set.middleRows(4*Nr, Nr).maxCoeff();
            tke_min = sn_set.middleRows(4*Nr, Nr).minCoeff();
            sn_set.middleRows(4*Nr, Nr) = (sn_set.middleRows(4*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*tke_min)/(tke_max - tke_min);

            omega_max = sn_set.middleRows(5*Nr, Nr).maxCoeff();
            omega_min = sn_set.middleRows(5*Nr, Nr).minCoeff();
            sn_set.middleRows(5*Nr, Nr) = (sn_set.middleRows(5*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*omega_min)/(omega_max - omega_min);

        }

    } else if ( nC == 7 ) //Turbulent 3D Navier-Stokes
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoW*Eigen::MatrixXd::Ones(Nr,1); //no sideslip angle
        Ic.segment(4*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(5*Nr, Nr) = rhotke*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(6*Nr, Nr) = rhoomega*Eigen::MatrixXd::Ones(Nr,1);
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;

        if ( flag == "YES" ) {
            //Introduce an if on the number of conservative variables

            rho_max = sn_set.middleRows(0, Nr).maxCoeff();
            rho_min = sn_set.middleRows(0, Nr).minCoeff();
            sn_set.middleRows(0, Nr) = (sn_set.middleRows(0, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rho_min )/(rho_max - rho_min);

            rhoU_max = sn_set.middleRows(Nr, Nr).maxCoeff();
            rhoU_min = sn_set.middleRows(Nr, Nr).minCoeff();
            sn_set.middleRows(Nr, Nr) = (sn_set.middleRows(Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoU_min)/(rhoU_max - rhoU_min);

            rhoV_max = sn_set.middleRows(2*Nr, Nr).maxCoeff();
            rhoV_min = sn_set.middleRows(2*Nr, Nr).minCoeff();
            sn_set.middleRows(2*Nr, Nr) = (sn_set.middleRows(2*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoV_min)/(rhoV_max - rhoV_min);

            rhoW_max = sn_set.middleRows(3*Nr, Nr).maxCoeff();
            rhoW_min = sn_set.middleRows(3*Nr, Nr).minCoeff();
            sn_set.middleRows(3*Nr, Nr) = (sn_set.middleRows(3*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoW_min)/(rhoW_max - rhoW_min);

            rhoE_max = sn_set.middleRows(4*Nr, Nr).maxCoeff();
            rhoE_min = sn_set.middleRows(4*Nr, Nr).minCoeff();
            sn_set.middleRows(4*Nr, Nr) = (sn_set.middleRows(4*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoE_min)/(rhoE_max - rhoE_min);

            tke_max = sn_set.middleRows(5*Nr, Nr).maxCoeff();
            tke_min = sn_set.middleRows(5*Nr, Nr).minCoeff();
            sn_set.middleRows(5*Nr, Nr) = (sn_set.middleRows(5*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*tke_min)/(tke_max - tke_min);

            omega_max = sn_set.middleRows(6*Nr, Nr).maxCoeff();
            omega_min = sn_set.middleRows(6*Nr, Nr).minCoeff();
            sn_set.middleRows(6*Nr, Nr) = (sn_set.middleRows(6*Nr, Nr) - Eigen::MatrixXd::Ones(Nr, settings.Ns)*omega_min)/(omega_max - omega_min);

        }

    } else {
        std::cout << "Set well number of Variables for subtracting initial condition" << std::endl;
        exit(EXIT_FAILURE);
    }

    return Ic;
}


Eigen::VectorXi Inverse_POS (const Eigen::MatrixXd &sn_set, int Nsamples) {

    int Ns = sn_set.cols();
    int Ndof = sn_set.rows();
    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(Ns);
    for ( int i = 0; i < Ns; i++ )
        norm_sn_set(i) = sn_set.col(i).norm()/(double)Ndof;

    Eigen::VectorXd ysamples = Eigen::VectorXd::LinSpaced(Nsamples, norm_sn_set.minCoeff(), norm_sn_set.maxCoeff());
    std::vector<int> I_POS = {};

//    std::cout << "Samples on y : "<< ysamples.transpose() << std::endl;

    for (int i = 0; i < Nsamples; i++) {
        Eigen::VectorXd ftime = norm_sn_set - Eigen::VectorXd::Ones(Ns)*ysamples(i);
//        std::cout << "ftime : "<< ftime.transpose() << std::endl;
        for ( int j = 0; j < Ns-1; j++) {
            if ( (ftime(j)*ftime(j+1)) < 0.0 )
                I_POS.push_back(j);
        }
    }

    std::sort(I_POS.begin(),I_POS.end());
    I_POS.erase(std::unique(I_POS.begin(),I_POS.end()),I_POS.end());

    Eigen::Map<Eigen::VectorXi> Ipos(I_POS.data(), I_POS.size());

    return Ipos;

}


