//
// Created by haitan on 26/06/2020.
//

#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Pre-Process.hpp"
#include "System_Calls.hpp"
#include "Post-Process.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"


int main( int argc, char *argv[] ) {

    std::cout << "-----------MODES Adaptive Residual Evaluation starts-------------" << std::endl << std::endl;

    std::cout << "Initializing common variables ... " << std::endl << std::endl;

    prob_settings settings;
    std::string filecfg = argv[1];
    Read_cfg(filecfg, settings);

    std::vector<double> t_vec(settings.Ns);
    int Nr;
    Eigen::MatrixXd Coords;

    common_vars(Nr, Coords, t_vec, settings);
// How we upload the snapshot matrix in the most efficient way?
// By now one big igen Matrix where each column has all the vectors of the conservative variables
    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix(Nr, settings);
    std::cout << std::endl;

    std::string binary = "YES";
    Eigen::VectorXd svd_cum_sum(settings.Ns);


//    std::cout << "Computing mean/Initial Condition of CFD solution ... " << std::endl;
//    //Defining Initial condition
//    Eigen::VectorXd mean = sn_set.rowwise().mean();
//    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC * Nr);
//
//    if (settings.flag_mean == "IC") {
//        std::cout << "Subtracting Initial condition" << std::endl << std::endl;
//        Ic = IC(sn_set, settings, nC, Nr);
//    } else {
//        std::cout << "Using data without subtracting any reference state" << std::endl << std::endl;
//    }


    std::cout << "Reading Gradients info ..." << std::endl;
    settings.flag_prob = "GRADIENTS";
    settings.Cols.clear();
//    settings.Cols.push_back(7);
//    settings.Cols.push_back(8);
    settings.Cols.push_back(9);
    settings.Cols.push_back(12);
    Eigen::MatrixXd sn_grad_set = generate_snap_matrix(Nr, settings);

//    std::cout << "Reading velocity info ..." << std::endl;
//    settings.flag_prob = "VELOCITY-2D";
//    settings.Cols.clear();
//    settings.Cols.push_back(3);
//    settings.Cols.push_back(4);
//    settings.Cols.push_back(5);
//    Eigen::MatrixXd sn_vel_set = generate_snap_matrix(Nr, settings);


    std::cout << "Reading matrix for check \n ";
    settings.nstart = 1;
    settings.flag_prob = "SCALAR";
    settings.Cols.clear();
    settings.Cols.push_back(3);
    Eigen::MatrixXd sn_check = generate_snap_matrix(Nr, settings);

    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(sn_check.cols());
    for ( int it = 0; it < sn_check.cols(); it ++)
        norm_sn_set(it) = sn_check.col(it).norm();

    std::cout << "Rearranging matrix of snapshots (with and without gradients)...\n norm of added snapshots" << std::endl;
    //Computing the final matrix with convective contributions
    Eigen::MatrixXd sn_set_improved(Nr, settings.Ns*2);
    double dt_l = 1e-7;//settings.Dt_cfd*settings.Ds/2.0; //As first guess
    std::vector<double> t_improved(2*settings.Ns);
    int count = 0;
    std::cout << "t improved :\n" ;
    for ( int it = 0; it < settings.Ns; it++ ){
        t_improved[count] = t_vec[it];
        count++;
        t_improved[count] = t_vec[it] + dt_l;
        std::cout << t_improved[count-1] << ", " << t_improved[count] << ", ";
        count++;

    }
    std::cout << std::endl;

    int h = 1; //As first guess
    count = 0;
    for ( int it = 0; it < settings.Ns; it++ ){
        sn_set_improved.col(count) = sn_set.col(it);
        count++;
//        sn_set_improved.col(it + settings.Ns) = sn_set.col(it) - dt_l*(sn_vel_set.topRows(Nr).col(it)*sn_grad_set.topRows(Nr).col(it)
//                + sn_vel_set.bottomRows(Nr).col(it)*sn_grad_set.bottomRows(Nr).col(it));
//        sn_set_improved.col(it + settings.Ns) = h*sn_grad_set.topRows(Nr).col(it) + sn_set.col(it);
//        sn_set_improved.col(it + 2*settings.Ns) = h*sn_grad_set.bottomRows(Nr).col(it) + sn_set.col(it);
    Eigen::VectorXd appo = Eigen::VectorXd::Zero(Nr);
        for ( int iPoint = 0; iPoint < Nr; iPoint++ ) {
//            appo(iPoint) = dt_l * (sn_vel_set(iPoint, it) * sn_grad_set(iPoint, it) +
//                                   sn_vel_set(Nr + iPoint, it) * sn_grad_set(Nr + iPoint, it));
            appo(iPoint) = dt_l * (sn_grad_set(iPoint, it) + sn_grad_set(Nr + iPoint, it));
        }
        sn_set_improved.col(count) =  sn_set.col(it) - appo;
        std::cout << sn_set_improved.col(count).norm() << ", " ;
        count++;

    }

    std::cout << std::endl;
    std::cout << "Norms of snapshots in time" << std::endl;
    for ( int it = 0; it < settings.Ns*2; it ++ ){
        std::cout << sn_set_improved.col(it).norm() << ", " ;
    }

    write_Reconstructed_fields ( sn_set_improved.col(1), Coords, settings.out_file, settings.flag_prob, 0);

    std::cout << std::endl;
//    delete sn_grad_set;
//    delete sn_vel_set;

    std::vector<Eigen::VectorXd> Error = {};
    std::vector<double> t(sn_check.cols(),0.0);
    t[0] = settings.nstart*settings.Dt_cfd;
    for ( int it = 1; it < sn_check.cols(); it++ )
        t[it] = t[it-1] + settings.Dt_cfd*settings.Ds;

    {
        std::cout << "Performing POD without gradients" << std::endl;
        Eigen::MatrixXd eig_vec;
        Eigen::VectorXd lambda;
        Eigen::VectorXd K_pc(settings.Ns);

        Eigen::MatrixXd PhiPOD = SPOD_basis(sn_set,
                                            lambda, K_pc, eig_vec,
                                            settings.Nf[0],
                                            settings.flag_bc,
                                            settings.flag_filter,
                                            settings.sigma);

        int N_max = PhiPOD.cols();
        int r = std::min(N_max, settings.r);
        std::cout << "Rank POD : " << r << std::endl;
        Eigen::MatrixXd CoefsPOD = PhiPOD.transpose() * sn_set;
        std::vector<rbf> surr_coefs =  getSurrCoefs (t_vec, CoefsPOD.transpose(), settings.flag_interp);

        Eigen::MatrixXd coef_t(sn_check.cols(), r);

        std::vector<double> tr(1);
        for ( int j = 0; j < sn_check.cols(); j++ ){
            tr[0] = t[j];
            for ( int i = 0; i < r; i++ )
                surr_coefs[i].evaluate(tr, coef_t(j,i));
        }


        Eigen::MatrixXd ErrMap_POD = sn_check - PhiPOD.leftCols(r) * coef_t.transpose();
        Eigen::VectorXd Errtime_POD = Eigen::VectorXd::Zero(sn_check.cols());

        for (int it = 0; it < sn_check.cols(); it++) {
            for (int iPoint = 0; iPoint < Nr; iPoint++) {
                Errtime_POD(it) += ErrMap_POD(iPoint, it) * ErrMap_POD(iPoint, it);
            }

            Errtime_POD(it) = std::sqrt(Errtime_POD(it)) / norm_sn_set(it);
        }

        Error.push_back(Errtime_POD);



//        Eigen::MatrixXd ErrMap_POD = sn_set - PhiPOD.leftCols(settings.r) * CoefsPOD.topRows(settings.r);
//        Eigen::VectorXd Errtime_POD = Eigen::VectorXd::Zero(settings.Ns);
//
//        for (int it = 0; it < settings.Ns; it++) {
//            for (int iPoint = 0; iPoint < Nr; iPoint++) {
//                Errtime_POD(it) += ErrMap_POD(iPoint, it) * ErrMap_POD(iPoint, it);
//            }
//
//            Errtime_POD(it) = std::sqrt(Errtime_POD(it)) / norm_sn_set(it);
//        }
//
//        Error.push_back(Errtime_POD);
    }

    {
        std::cout << "Performing POD with gradients" << std::endl;
        Eigen::MatrixXd eig_vec;
        Eigen::VectorXd lambda;
        Eigen::VectorXd K_pc(settings.Ns*2);

        Eigen::MatrixXd PhiPOD = SPOD_basis(sn_set_improved,
                                            lambda, K_pc, eig_vec,
                                            settings.Nf[0],
                                            settings.flag_bc,
                                            settings.flag_filter,
                                            settings.sigma);

        int N_max = PhiPOD.cols();
        int r = std::min(N_max, settings.r);
        std::cout << "Rank GEPOD : " << r << std::endl;
        Eigen::MatrixXd CoefsPOD = PhiPOD.transpose() * sn_set_improved;
        std::vector<rbf> surr_coefs =  getSurrCoefs (t_improved, CoefsPOD.transpose(), settings.flag_interp);

        Eigen::MatrixXd coef_t(sn_check.cols(), r);

        std::vector<double> tr(1);
        for ( int j = 0; j < sn_check.cols(); j++ ){
            tr[0] = t[j];
            for ( int i = 0; i < r; i++ )
                surr_coefs[i].evaluate(tr, coef_t(j,i));
        }


        Eigen::MatrixXd ErrMap_POD = sn_check - PhiPOD.leftCols(r) * coef_t.transpose();
        Eigen::VectorXd Errtime_POD = Eigen::VectorXd::Zero(sn_check.cols());

        for (int it = 0; it < sn_check.cols(); it++) {
            for (int iPoint = 0; iPoint < Nr; iPoint++) {
                Errtime_POD(it) += ErrMap_POD(iPoint, it) * ErrMap_POD(iPoint, it);
            }

            Errtime_POD(it) = std::sqrt(Errtime_POD(it)) / norm_sn_set(it);
        }

        Error.push_back(Errtime_POD);


//        Eigen::MatrixXd CoefsPOD = PhiPOD.transpose()*sn_set;
//
//        Eigen::MatrixXd ErrMap_POD = sn_set - PhiPOD.leftCols(settings.r)*CoefsPOD.topRows(settings.r);
//        Eigen::VectorXd Errtime_POD = Eigen::VectorXd::Zero(settings.Ns);
//
//        for ( int it = 0; it < settings.Ns; it++ ){
//            for ( int iPoint = 0; iPoint < Nr; iPoint++ ){
//                Errtime_POD(it) += ErrMap_POD(iPoint,it)*ErrMap_POD(iPoint,it);
//            }
//
//            Errtime_POD(it) = std::sqrt(Errtime_POD(it))/norm_sn_set(it);
//        }
//
//        Error.push_back(Errtime_POD);

    }


    std::cout << "Writing error of interpolation and error from projection ... " << std::endl;
    std::ofstream errfile;
    errfile.open("Err_grad_nograd.dat");

    errfile << "T(s), Err_POD, ErrEPOD" << std::endl;

    for ( int it = 0; it < sn_check.cols(); it ++ ){
            errfile << t[it] << ", " << std::setprecision(8) << Error[0](it)
                                    << ", " << std::setprecision(8) << Error[1](it);

        errfile << std::endl;

    }

    errfile.close();


    return 0;
}