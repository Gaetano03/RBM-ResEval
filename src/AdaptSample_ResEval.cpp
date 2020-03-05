/*
CODE FOR RESIDUAL EVALUATION OF DIFFERENT RBM TECHNIQUES USING SU2 CODE
INPUT ARGUMENTS
Config File RBM + Config File SU2
*/

#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Pre-Process.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"


int main( int argc, char *argv[] )
{

    std::cout << "-----------MODES Adaptive Residual Evaluation starts-------------" << std::endl << std::endl;
//
    std::cout << "Initializing common variables ... " << std::endl << std::endl;
    prob_settings settings;
    std::string filecfg = argv[1];
    std::string su2_conf = argv[2];
    std::string decision = argv[3];
    std::string su2dtr_string = "mpirun -np 6 ./SU2_DTR " + su2_conf + " > SU2.log"; // + " > resEval_su2.log";
    int len_s = su2dtr_string.length();
    char su2_sys_call[len_s + 1];
    strcpy(su2_sys_call, su2dtr_string.c_str());

    Read_cfg( filecfg, settings );
//
    std::string root_outputfile;
    root_outputfile.assign ( settings.out_file, 0, settings.out_file.size() - 4);
    std::string root_inputfile;
    root_inputfile.assign ( settings.in_file, 0, settings.in_file.size() - 4);
    std::string input_format;
    input_format.assign ( settings.in_file, settings.in_file.size() - 3, 3);
    std::string rmf_string = "rm -f " + root_outputfile + "_*";
    len_s = rmf_string.length();
    char rmf_sys_call[len_s + 1];
    strcpy(rmf_sys_call, rmf_string.c_str());

    int s_Nf = 5;   //Number of values for the SPOD filter (POD included)
    std::vector<int> Nf(s_Nf);
    Nf[0] = 0;
    Nf[1] = std::ceil(settings.Ns/10.0);
    Nf[2] = std::ceil(settings.Ns/2.0);
    Nf[3] = std::ceil(2.0*settings.Ns/3.0);
    Nf[4] = settings.Ns;

    int nC = settings.Cols.size();
    double alpha = settings.alpha;
    double beta = settings.beta;
    if (settings.ndim == 2) beta = 0.0;

    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(5) << std::to_string(settings.nstart);
    std::string file_1 = root_inputfile + "_" + buffer.str() + "." + input_format;

    // Calculate number of grid points
    int Nr = N_gridpoints ( file_1 );
    std::cout << "Number of grid points : " << Nr << std::endl;

    std::cout << "Reading Coordinates ... \t ";
    Eigen::MatrixXd Coords = read_col( file_1, Nr, settings.Cols_coords );
    std::cout << "Done " << std::endl;

 // How we upload the snapshot matrix in the most efficient way?
 // By now one big igen Matrix where each column has all the vectors of the conservative variables
    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings.Ns, settings.Ds, settings.nstart,
                                                   settings.Cols,
                                                   settings.in_file,
                                                   settings.flag_prob);

    std::cout << "Initializing Vector of times ... " << std::endl;
    std::vector<double> t_vec( settings.Ns );
    t_vec[0] = settings.nstart*settings.Dt_cfd;
    for ( int i = 1; i < settings.Ns; i++ )
        t_vec[i] = t_vec[i-1] + settings.Dt_cfd*settings.Ds;

    // std::vector<double> t_evaluate(2*settings.Ns-1);
    // t_evaluate[0] = settings.nstart*settings.Dt_cfd;
    // for ( int i = 1; i < t_evaluate.size(); i++)
    //     t_evaluate[i] = t_evaluate[i-1] + settings.Dt_cfd*(double)settings.Ds/2.0;

    std::cout << "Computing mean/Initial Condition of CFD solution ... " << std::endl;

    Eigen::VectorXd mean = sn_set.rowwise().mean();
    if ( settings.flag_mean == "YES" ) {
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= mean;
    }

    //Defining Initial condition
    Eigen::VectorXd Ic = IC(settings, nC, Nr);

    std::string binary = "YES";
    Eigen::VectorXd svd_cum_sum(settings.Ns);

    if ( settings.flag_mean == "IC" )
    {
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;
    }


    //Define normalization for conservative variables
    double rho_max, rho_min, rhoU_max, rhoU_min, rhoV_max, rhoV_min, rhoW_max, rhoW_min, rhoE_max, rhoE_min; //add turbulence

    if ( settings.flag_prob == "CONSERVATIVE" ) {
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


    int Nm;
 //Defining common scope for uniform sampling
    {

        //Vector of MatrixXd where to store the evolution in time of conservative variables
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
        Eigen::MatrixXd Phi = Eigen::MatrixXd::Zero(Nr,settings.Ns);
        Eigen::VectorXd lambda = Eigen::VectorXd::Zero(settings.Ns);

        std::vector<smartuq::surrogate::rbf> surr_coefs;
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);
        int N_notZero;
        //Check only for POD for now
        std::cout << "Computing uniform SPOD modes with Nf : " << settings.Nf << "\n";
        Phi = SPOD_basis( sn_set,
                          lambda, K_pc, eig_vec,
                          settings.Nf,
                          settings.flag_bc,
                          settings.flag_filter,
                          settings.sigma);
        N_notZero = Phi.cols();
        if ( settings.r == 0 ) Nm = Nmod(settings.En, K_pc);
        else Nm = std::min(settings.r, N_notZero);
        std::cout << "Number of modes used in reconstruction " << Nm << std::endl;
        surr_coefs = getSurrCoefs(t_vec, eig_vec, settings.flag_interp);

        if ( decision == "-e") {
            for ( int itr = 0; itr < settings.t_res.size(); itr++ ) {
                std::cout << " Computing residuals at time : " << settings.t_res[itr] << std::endl;

                Eigen::MatrixXd coef_t(3, Nm);
                if (settings.t_res[itr] - 2.0 * settings.Dt_res[0] < t_vec[0] ||
                    settings.t_res[itr] > t_vec[t_vec.size() - 1] ) {
                    std::cout
                            << "Define proper Delta_t_res and T_RES vector " << std::endl;
                    exit(EXIT_FAILURE);
                } else {
                    std::vector<double> tr(1);
                    std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[0],
                                                      settings.t_res[itr] - settings.Dt_res[0],
                                                      settings.t_res[itr]};

                    for (int j = 0; j < 3; j++) {
                        tr[0] = t_evaluate[j];
                        for (int i = 0; i < Nm; i++)
                            surr_coefs[i].evaluate(tr, coef_t(j, i));
                    }


                }
                Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm, Nm);
                for (int i = 0; i < Nm; i++)
                    Sig(i, i) = std::sqrt(lambda(i));
                Sn_Cons_time = Phi.leftCols(Nm) * Sig * coef_t.transpose();


                //Introduce an if on the number of conservative variables
                 Sn_Cons_time.middleRows(0, Nr) = Sn_Cons_time.middleRows(0, Nr) * (rho_max - rho_min) + Eigen::MatrixXd::Ones(Nr, settings.Ns)*rho_min;
                 Sn_Cons_time.middleRows(Nr, Nr) = Sn_Cons_time.middleRows(0, Nr) * (rhoU_max - rhoU_min) + Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoU_min;
                 Sn_Cons_time.middleRows(2 * Nr, Nr) = Sn_Cons_time.middleRows(2 * Nr, Nr) * (rhoV_max - rhoV_min) + Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoV_min;
                 Sn_Cons_time.middleRows(3 * Nr, Nr) = Sn_Cons_time.middleRows(3 * Nr, Nr) * (rhoE_max - rhoE_min) + Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoE_min;

                if (settings.flag_mean == "IC") {
                    for (int it = 0; it < 3; it++)
                        Sn_Cons_time.col(it) += Ic;
                } else if (settings.flag_mean == "YES") {
                    for (int it = 0; it < 3; it++)
                        Sn_Cons_time.col(it) += mean;
                }

                std::string mv_string = "mv history_rbm_00002.csv history_pod_uniform_" + std::to_string(itr) + ".csv";
                len_s = mv_string.length();
                char mv_sys_call[len_s + 1];
                strcpy(mv_sys_call, mv_string.c_str());
                std::cout << "Writing time reconstruction " << std::endl;
                // Write_Restart_Cons_Time( Sn_Cons_time, Coords, settings.out_file, t_evaluate.size(), nC, alpha );
                Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, 3, nC, alpha, beta, binary);
                //Executing SU2, removing all useless files, renaming files with residuals
                std::cout << "Calling SU2 for residual evaluation and writing file to history " << std::endl;
                std::system(su2_sys_call);
                std::system(rmf_sys_call);
                std::system(mv_sys_call);
            }
        } else if ( decision == "-r") {
            for ( int nt = 0; nt < settings.t_rec.size(); nt++ ) {

                std::cout << "Reconstructing Momentum at time : " << settings.t_rec[nt] << "\t";
                Eigen::MatrixXd Rec = Reconstruction_S_POD(t_vec,
                                                           K_pc, lambda, eig_vec.transpose(),
                                                           Phi.middleRows(Nr, 2 * Nr), settings.t_rec[nt],
                                                           Nm,
                                                           "VECTOR-2D",
                                                           settings.flag_interp);

                std::cout << "Done" << std::endl;

                Rec.col(0) = Rec.col(0) * (rhoU_max - rhoU_min) + Eigen::MatrixXd::Ones(Nr, 1)*rhoU_min;
                Rec.col(1) = Rec.col(1) * (rhoV_max - rhoV_min) + Eigen::MatrixXd::Ones(Nr, 1)*rhoV_min;

                if (settings.flag_mean == "IC") {
                    Rec.col(0) = Rec.col(0) + Ic.middleRows(Nr, Nr);
                    Rec.col(1) = Rec.col(1) + Ic.middleRows(2*Nr, Nr);
                }
                if (settings.flag_mean == "YES") {
                    Rec.col(0) = Rec.col(0) + mean.middleRows(Nr, Nr);
                    Rec.col(1) = Rec.col(1) + mean.middleRows(2*Nr, Nr);
                }

                std::cout << "Writing reconstructed field ..." << "\t";
                std::string filename = root_outputfile + "_uni.dat";

                write_Reconstructed_fields(Rec, Coords,
                                           filename,
                                           "VECTOR-2D", nt);

                std::cout << "Done" << std::endl << std::endl;
            }

        } else {
            std::cout << "Nothing else to do!" << std::endl;
        }
 //                 Sn_Cons_time.middleRows(ncons*Nr,Nr) = Phi.leftCols(Nm)*Sig*eig_vec.leftCols(Nm).transpose();
    }


//    Defining common scope for adaptive sampling
    {
        int nVar = 12; //that has to contain also first and last snapshot
        Eigen::VectorXi t_pos(nVar);
        std::cout << "Nm =  " << Nm << std::endl;
 //        0, 1, 2, 5, 7, 10, 13, 14, 17, 21, 24, 26, 28, 33, 36, 39, 42, 48, 53, 56, 62, 65, 68, 73, 76, 77,
 //                82, 88, 93, 99;
 //        t_pos << 0, 1, 2, 4, 5, 7, 8, 9, 11, 14, 16, 19, 22, 25, 29, 33, 38, 42, 48, 51, 55, 62, 65, 70, 76, 82, 93, 99;
//        t_pos << 0, 1, 3, 7, 17, 29, 41, 58, 78, 99;
        t_pos << 0, 1, 6, 10, 19, 28, 39, 44, 52, 67, 83,  99;
//        std::vector<double> t_vec( nVar );
 //        for ( int i = 0; i < nVar; i++ )
 //            t_vec[i] = settings.Dt_cfd*settings.Ds*t_pos(i);

        //Vector of MatrixXd where to store the evolution in time of conservative variables
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
        Eigen::MatrixXd Phi = Eigen::MatrixXd::Zero(Nr,nVar);
        Eigen::VectorXd lambda = Eigen::VectorXd::Zero(nVar);

        std::vector<smartuq::surrogate::rbf> surr_coefs ;
        Eigen::VectorXd K_pc(nVar);
        Eigen::MatrixXd eig_vec(nVar, nVar);
        int N_notZero;
        //Check only for POD for now
        std::cout << "Computing adaptive SPOD modes with Nf : " << settings.Nf << "\n";

        Eigen::MatrixXd sub_sn_set = indexing(sn_set, Eigen::ArrayXi::LinSpaced(nC*Nr,0,nC*Nr-1),t_pos);
        Phi = SPOD_basis( sub_sn_set,
                          lambda, K_pc, eig_vec,
                          settings.Nf,
                          settings.flag_bc,
                          settings.flag_filter,
                          settings.sigma);

        Eigen::MatrixXd Coeffs = Phi.transpose()*sn_set;
        surr_coefs = getSurrCoefs(t_vec, Coeffs.transpose(), settings.flag_interp);

 //             Sn_Cons_time.middleRows(ncons*Nr,Nr) = Phi.leftCols(Nm)*Sig*eig_vec.leftCols(Nm).transpose();
        if ( decision == "-e" ) {
            for ( int itr = 0; itr < settings.t_res.size(); itr++ ) {
                std::cout << " Computing residuals at time : " << settings.t_res[itr] << std::endl;

                Eigen::MatrixXd coef_t(3, Nm);
                if (settings.t_res[itr] - 2.0 * settings.Dt_res[0] < t_vec[0] ||
                    settings.t_res[itr] > t_vec[t_vec.size() - 1]) {
                    std::cout
                            << "Define proper Delta_t_res and T_RES vector " << std::endl;
                    exit(EXIT_FAILURE);
                } else {
                    std::vector<double> tr(1);
                    std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[0],
                                                      settings.t_res[itr] - settings.Dt_res[0],
                                                      settings.t_res[itr]};

                    for (int j = 0; j < 3; j++) {
                        tr[0] = t_evaluate[j];
                        for (int i = 0; i < Nm; i++)
                            surr_coefs[i].evaluate(tr, coef_t(j, i));
                    }


                }
                Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm, Nm);
                for (int i = 0; i < Nm; i++)
                    Sig(i, i) = std::sqrt(lambda(i));
                Sn_Cons_time = Phi.leftCols(Nm) * coef_t.transpose();


                //Introduce an if on the number of conservative variables
                Sn_Cons_time.middleRows(0, Nr) = Sn_Cons_time.middleRows(0, Nr) * (rho_max - rho_min) + Eigen::MatrixXd::Ones(Nr, settings.Ns)*rho_min;
                Sn_Cons_time.middleRows(Nr, Nr) = Sn_Cons_time.middleRows(0, Nr) * (rhoU_max - rhoU_min) + Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoU_min;
                Sn_Cons_time.middleRows(2 * Nr, Nr) = Sn_Cons_time.middleRows(2 * Nr, Nr) * (rhoV_max - rhoV_min) + Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoV_min;
                Sn_Cons_time.middleRows(3 * Nr, Nr) = Sn_Cons_time.middleRows(3 * Nr, Nr) * (rhoE_max - rhoE_min) + Eigen::MatrixXd::Ones(Nr, settings.Ns)*rhoE_min;

                if (settings.flag_mean == "IC") {
                    for (int it = 0; it < 3; it++)
                        Sn_Cons_time.col(it) += Ic;
                } else if (settings.flag_mean == "YES") {
                    for (int it = 0; it < 3; it++)
                        Sn_Cons_time.col(it) += mean;
                }

                std::string mv_string = "mv history_rbm_00002.csv history_pod_adaptive_" + std::to_string(itr) + ".csv";
                len_s = mv_string.length();
                char mv_sys_call[len_s + 1];
                strcpy(mv_sys_call, mv_string.c_str());
                std::cout << "Writing time reconstruction " << std::endl;
                // Write_Restart_Cons_Time( Sn_Cons_time, Coords, settings.out_file, t_evaluate.size(), nC, alpha );
                Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, 3, nC, alpha, beta, binary);
                //Executing SU2, removing all useless files, renaming files with residuals
                std::cout << "Calling SU2 for residual evaluation and writing file to history " << std::endl;
                std::system(su2_sys_call);
                std::system(rmf_sys_call);
                std::system(mv_sys_call);
            }
        }
        else if ( decision == "-r" ) {
            for ( int nt = 0; nt < settings.t_rec.size(); nt++ ) {

                std::cout << "Reconstructing Momentum at time : " << settings.t_rec[nt] << "\t";
                Eigen::MatrixXd Rec = Reconstruction_S_POD(t_vec,
                                                           K_pc, lambda, Coeffs,
                                                           Phi.middleRows(Nr, 2 * Nr), settings.t_rec[nt],
                                                           Nm,
                                                           "VECTOR-2D",
                                                           settings.flag_interp);

                std::cout << "Done" << std::endl;

                Rec.col(0) = Rec.col(0) * (rhoU_max - rhoU_min) + Eigen::MatrixXd::Ones(Nr,1) * rhoU_min;
                Rec.col(1) = Rec.col(1) * (rhoV_max - rhoV_min) + Eigen::MatrixXd::Ones(Nr,1) * rhoV_min;

                if (settings.flag_mean == "IC") {
                    Rec.col(0) = Rec.col(0) + Ic.middleRows(Nr, Nr);
                    Rec.col(1) = Rec.col(1) + Ic.middleRows(2*Nr, Nr);
                }
                if (settings.flag_mean == "YES") {
                    Rec.col(0) = Rec.col(0) + mean.middleRows(Nr, Nr);
                    Rec.col(1) = Rec.col(1) + mean.middleRows(2*Nr, Nr);
                }

                std::cout << "Writing reconstructed field ..." << "\t";
                std::string filename = root_outputfile + "_adapt.dat";

                write_Reconstructed_fields(Rec, Coords,
                                           filename,
                                           "VECTOR-2D", nt);

                std::cout << "Done" << std::endl << std::endl;
            }
        }
        else {
            std::cout << "Nothing to do!" << std::endl;
        }
    }

    std::cout << "MODES Adaptive Residual Evaluation ends" << std::endl << std::endl;
    return 0;

}



