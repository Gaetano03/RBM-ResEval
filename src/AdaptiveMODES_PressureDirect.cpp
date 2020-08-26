/*
CODE FOR ERROR EVALUATION OF CONSERVATIVE QUANTITIES USING DIFFERENT RBM TECHNIQUES
!!!!!!(restart files need to be numbered without gaps)!!!!
INPUT ARGUMENTS
Config File RBM
ATTENTION ------> Set DS to the minimum possible in config file
 */

#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Pre-Process.hpp"
#include "System_Calls.hpp"
#include "Post-Process.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"



int main( int argc, char *argv[] ) {

    std::cout << "----------- Direct error on derived pressure starts-----------" << std::endl << std::endl;

    std::cout << "Initializing common variables ... " << std::endl << std::endl;
    prob_settings settings;
    std::string filecfg = argv[1];

    std::vector<Eigen::VectorXd> Err_RBM_Nm_time;
    std::vector<Eigen::VectorXd> ErrP_RBM_Nm_time;
    std::vector<Eigen::VectorXd> EN;

    Read_cfg(filecfg, settings);


    int nC = settings.Cols.size() - 1;
    int Nr;
    Eigen::MatrixXd Coords;
    std::vector<double> t_vec(settings.Ns);
    common_vars(Nr, Coords, t_vec, settings);

    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix(Nr, settings);

    int Nsnap;
    if (settings.Ns % 2 == 0) {
        Nsnap = settings.Ns/2;
        t_vec.pop_back();
    }
    else Nsnap = std::floor(settings.Ns / 2) + 1;

    std::vector<double> t_train(Nsnap);

    t_train[0] = settings.nstart * settings.Dt_cfd;
    for (int i = 1; i < t_train.size(); i++)
        t_train[i] = t_train[i - 1] + settings.Dt_cfd * (double) settings.Ds * 2.0;// + settings.Dt_cfd*(double)settings.Ds/2.0;


    std::cout << "Computing mean/Initial Condition of CFD solution ... " << std::endl;
    //Defining Initial condition
    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC * Nr);

    Eigen::VectorXd svd_cum_sum(settings.Ns);

    if (settings.flag_mean == "IC") {
        std::cout << "Subtracting Initial condition" << std::endl << std::endl;
        Eigen::MatrixXd dummy = sn_set.topRows(nC*Nr);
        Ic = IC(dummy, settings, nC, Nr);
        sn_set.topRows(nC*Nr) = dummy;
    } else {
        std::cout << "Using data without subtracting any reference state" << std::endl << std::endl;
    }

    Eigen::MatrixXd sn_set_p(nC * Nr, Nsnap);
    int count = 0;
    for (int is = 0; is < settings.Ns; is += 2) {
        sn_set_p.col(count) = sn_set.topRows(nC*Nr).col(is);
        count++;
    }

    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(t_vec.size());
    for (int i = 0; i < t_vec.size(); i++)
        norm_sn_set(i) = sn_set.bottomRows(Nr).col(i).norm();

//Defining common scope for POD
    auto methods = settings.flag_method;
    std::vector<std::string>::iterator itPOD;
    itPOD = std::find(methods.begin(), methods.end(), "POD");
    if (itPOD != methods.end()) {

        Eigen::VectorXd lambda(Nsnap);
        Eigen::VectorXd K_pc(Nsnap);
        Eigen::MatrixXd eig_vec(Nsnap, Nsnap);
        std::vector<Eigen::MatrixXd> RecCons(nC);
        int Nm;
        std::vector<int> N_notZero(nC);
        //Check only for POD for now

        for (int ncons = 0; ncons < nC; ncons++) {
            std::cout << "Processing conservative variable " << ncons << std::endl;
            Eigen::MatrixXd Phi = SPOD_basis(sn_set_p.middleRows(ncons * Nr, Nr),
                                    lambda, K_pc, eig_vec,
                                    0,
                                    settings.flag_bc,
                                    settings.flag_filter,
                                    settings.sigma);
            N_notZero[ncons] = Phi.cols();

            std::vector<rbf> surr_coefs = getSurrCoefs(t_train, eig_vec, settings.flag_interp);
            Eigen::MatrixXd coef_t(t_vec.size(), N_notZero[ncons]);

            std::vector<double> tr(1);
            for (int j = 0; j < t_vec.size(); j++) {
                tr[0] = t_vec[j];
                for (int i = 0; i < N_notZero[ncons]; i++)
                    surr_coefs[i].evaluate(tr, coef_t(j, i));
            }

            Nm = std::min(settings.r, N_notZero[ncons]);
            Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm, Nm);
            for (int i = 0; i < Nm; i++)
                Sig(i, i) = std::sqrt(lambda(i));

            RecCons[ncons] = Phi.leftCols(Nm) * Sig * coef_t.leftCols(Nm).transpose();
            if (settings.flag_mean == "IC") {
                for (int it = 0; it < t_vec.size(); it++)
                    RecCons[ncons] += Ic.middleRows(ncons * Nr, Nr);
            }
        }

        std::cout << "Computing pressure from reconstructed vars..." << std::endl;
        Eigen::MatrixXd Rec_pressure(Nr,t_vec.size());
        double u,v,w;

        if ( settings.ndim == 2 ) {
            for (int i = 0; i < Nr; i++) {
                for (int j = 0; j < t_vec.size(); j++) {
                    u = RecCons[1](i,j)/RecCons[0](i,j);
                    v = RecCons[2](i,j)/RecCons[0](i,j);
                    Rec_pressure(i,j) = 0.4*(RecCons[3](i,j) - 0.5*RecCons[0](i,j)*(u*u + v*v));
                }
            }
        } else {
            for (int i = 0; i < Nr; i++) {
                for (int j = 0; j < t_vec.size(); j++) {
                    u = RecCons[1](i,j)/RecCons[0](i,j);
                    v = RecCons[2](i,j)/RecCons[0](i,j);
                    w = RecCons[3](i,j)/RecCons[0](i,j);
                    Rec_pressure(i,j) = 0.4*(RecCons[4](i,j) - 0.5*RecCons[0](i,j)*(u*u + v*v + w*w));
                }
            }
        }

        Eigen::MatrixXd Err_SPOD_map = Eigen::MatrixXd::Zero(Nr, t_vec.size());
        Eigen::VectorXd Err_SPOD_Nm_time = Eigen::VectorXd::Zero(t_vec.size());

        Err_SPOD_map = sn_set.bottomRows(Nr).leftCols(t_vec.size()) - Rec_pressure;

        for (int i = 0; i < t_vec.size(); i++) {
            for (int j = 0; j < Nr; j++)
                Err_SPOD_Nm_time(i) += Err_SPOD_map(j, i) * Err_SPOD_map(j, i);

            Err_SPOD_Nm_time(i) = std::sqrt(Err_SPOD_Nm_time(i)) / norm_sn_set(i);

        }

        std::ofstream errfile;
        std::string file_err_name = "Err_PressurePOD.dat";
        errfile.open(file_err_name);

        for (int nm = 0; nm < t_vec.size(); nm++) {
            errfile << std::setprecision(8) << Err_SPOD_Nm_time(nm) << "\t";
            errfile << std::endl;
        }

        errfile.close();

    }

    //Defining common scope for RDMD
    std::vector<std::string>::iterator itRDMD;
    itRDMD = std::find(methods.begin(), methods.end(), "RDMD");
    if (itRDMD != methods.end()) {

        Eigen::VectorXd lambda =  Eigen::VectorXd::Zero(Nsnap);
        Eigen::VectorXd K_pc(Nsnap);
        Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
        std::vector<Eigen::MatrixXd> RecCons(nC);

        int Nm;

        for (int ncons = 0; ncons < nC; ncons++) {
            std::cout << "Processing conservative variable " << ncons << std::endl;
            Eigen::MatrixXd Phi = RDMD_modes_coefs ( sn_set_p.middleRows(ncons*Nr,Nr),
                                            Coefs,
                                            lambda,
                                            K_pc,
                                            -1,
                                            Nsnap,
                                            settings.En );

            std::vector<rbf> surr_coefs = getSurrCoefs(t_train, Coefs.transpose(), settings.flag_interp);
            Eigen::MatrixXd coef_t(t_vec.size(), Nsnap);

            std::vector<double> tr(1);
            for (int j = 0; j < t_vec.size(); j++) {
                tr[0] = t_vec[j];
                for (int i = 0; i < Nsnap; i++)
                    surr_coefs[i].evaluate(tr, coef_t(j, i));
            }

            int N_max = Phi.cols();
            Nm = std::min(settings.r_RDMD, N_max);

            RecCons[ncons] = Phi.leftCols(Nm) * coef_t.leftCols(Nm).transpose();
            if (settings.flag_mean == "IC") {
                for (int it = 0; it < t_vec.size(); it++)
                    RecCons[ncons] += Ic.middleRows(ncons * Nr, Nr);
            }
        }

        std::cout << "Computing pressure from reconstructed vars..." << std::endl;
        Eigen::MatrixXd Rec_pressure(Nr,t_vec.size());
        double u,v,w;

        if ( settings.ndim == 2 ) {
            for (int i = 0; i < Nr; i++) {
                for (int j = 0; j < t_vec.size(); j++) {
                    u = RecCons[1](i,j)/RecCons[0](i,j);
                    v = RecCons[2](i,j)/RecCons[0](i,j);
                    Rec_pressure(i,j) = 0.4*(RecCons[3](i,j) - 0.5*RecCons[0](i,j)*(u*u + v*v));
                }
            }
        } else {
            for (int i = 0; i < Nr; i++) {
                for (int j = 0; j < t_vec.size(); j++) {
                    u = RecCons[1](i,j)/RecCons[0](i,j);
                    v = RecCons[2](i,j)/RecCons[0](i,j);
                    w = RecCons[3](i,j)/RecCons[0](i,j);
                    Rec_pressure(i,j) = 0.4*(RecCons[4](i,j) - 0.5*RecCons[0](i,j)*(u*u + v*v + w*w));
                }
            }
        }

        Eigen::MatrixXd Err_RDMD_map = Eigen::MatrixXd::Zero(Nr, t_vec.size());
        Eigen::VectorXd Err_RDMD_Nm_time = Eigen::VectorXd::Zero(t_vec.size());

        Err_RDMD_map = sn_set.bottomRows(Nr).leftCols(t_vec.size()) - Rec_pressure;

        for (int i = 0; i < t_vec.size(); i++) {
            for (int j = 0; j < Nr; j++)
                Err_RDMD_Nm_time(i) += Err_RDMD_map(j, i) * Err_RDMD_map(j, i);

            Err_RDMD_Nm_time(i) = std::sqrt(Err_RDMD_Nm_time(i)) / norm_sn_set(i);

        }

        std::ofstream errfile;
        std::string file_err_name = "Err_PressureRDMD.dat";
        errfile.open(file_err_name);

        for (int nm = 0; nm < t_vec.size(); nm++) {
            errfile << std::setprecision(8) << Err_RDMD_Nm_time(nm) << "\t";
            errfile << std::endl;
        }

        errfile.close();

    }

    //Defining common scope for DMD
    std::vector<std::string>::iterator itDMD;
    itDMD = std::find(methods.begin(), methods.end(), "DMD");
    if (itDMD != methods.end()) {

        Eigen::VectorXd lambda_POD;
        Eigen::MatrixXd eig_vec_POD;
        Eigen::VectorXcd lambda_DMD;
        Eigen::MatrixXcd eig_vec_DMD;
        std::vector<Eigen::MatrixXd> RecCons(nC);

        for (int ncons = 0; ncons < nC; ncons++) {

            Eigen::MatrixXcd Phi;
            Eigen::VectorXcd alfa;

            Phi = DMD_basis( sn_set_p.middleRows(ncons*Nr,Nr),
                             lambda_DMD,
                             eig_vec_DMD,
                             lambda_POD,
                             eig_vec_POD,
                             settings.r );

            //         int Nm = Phi.cols();
            //         std::cout << "Number of modes extracted : " << Nm << std::endl;

            Eigen::VectorXcd omega(Phi.cols());
            for ( int i = 0; i < Phi.cols(); i++ )
                omega(i) = std::log(lambda_DMD(i))/(settings.Dt_cfd*settings.Ds*2.0);

            // std::cout << "Calculating coefficients DMD ... " << "\t";
            alfa = Calculate_Coefs_DMD_exact (sn_set_p.middleRows(ncons*Nr,Nr).leftCols(Nsnap-1),
                                              lambda_DMD,
                                              Phi );


            Eigen::MatrixXcd V_and(lambda_DMD.size(), t_vec.size());
            for ( int i = 0; i < lambda_DMD.size(); i++ ) {
                for ( int j = 0; j < t_vec.size(); j++ )
                    V_and(i,j) = std::pow(lambda_DMD(i), (double)j/2.0);
            }

            Eigen::MatrixXcd Psi = Eigen::MatrixXcd::Zero(alfa.size(), t_vec.size());
            for ( int i = 0; i < t_vec.size(); i++ )
                Psi.col(i) = alfa.cwiseProduct(V_and.col(i));

            Eigen::MatrixXcd D_dmd = Phi*Psi;

            RecCons[ncons] = D_dmd.real();
            if (settings.flag_mean == "IC") {
                for (int it = 0; it < t_vec.size(); it++)
                    RecCons[ncons] += Ic.middleRows(ncons * Nr, Nr);
            }
        }

        std::cout << "Computing pressure from reconstructed vars..." << std::endl;
        Eigen::MatrixXd Rec_pressure(Nr,t_vec.size());
        double u,v,w;

        if ( settings.ndim == 2 ) {
            for (int i = 0; i < Nr; i++) {
                for (int j = 0; j < t_vec.size(); j++) {
                    u = RecCons[1](i,j)/RecCons[0](i,j);
                    v = RecCons[2](i,j)/RecCons[0](i,j);
                    Rec_pressure(i,j) = 0.4*(RecCons[3](i,j) - 0.5*RecCons[0](i,j)*(u*u + v*v));
                }
            }
        } else {
            for (int i = 0; i < Nr; i++) {
                for (int j = 0; j < t_vec.size(); j++) {
                    u = RecCons[1](i,j)/RecCons[0](i,j);
                    v = RecCons[2](i,j)/RecCons[0](i,j);
                    w = RecCons[3](i,j)/RecCons[0](i,j);
                    Rec_pressure(i,j) = 0.4*(RecCons[4](i,j) - 0.5*RecCons[0](i,j)*(u*u + v*v + w*w));
                }
            }
        }

        Eigen::MatrixXd Err_DMD_map = Eigen::MatrixXd::Zero(Nr, t_vec.size());
        Eigen::VectorXd Err_DMD_Nm_time = Eigen::VectorXd::Zero(t_vec.size());

        Err_DMD_map = sn_set.bottomRows(Nr).leftCols(t_vec.size()) - Rec_pressure;

        for (int i = 0; i < t_vec.size(); i++) {
            for (int j = 0; j < Nr; j++)
                Err_DMD_Nm_time(i) += Err_DMD_map(j, i) * Err_DMD_map(j, i);

            Err_DMD_Nm_time(i) = std::sqrt(Err_DMD_Nm_time(i)) / norm_sn_set(i);

        }


        std::ofstream errfile;
        std::string file_err_name = "Err_PressureDMD.dat";
        errfile.open(file_err_name);

        for (int nm = 0; nm < t_vec.size(); nm++) {
            errfile << std::setprecision(8) << Err_DMD_Nm_time(nm) << "\t";
            errfile << std::endl;

        }

        errfile.close();

    }


    std::cout << "----------- Direct error on derived pressure starts-----------" << std::endl << std::endl;
    return 0;
}


