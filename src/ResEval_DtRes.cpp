/*
CODE FOR RESIDUAL EVALUATION OF DIFFERENT RBM TECHNIQUES USING SU2 CODE
INPUT ARGUMENTS
Config File RBM + Config File SU2
*/

#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"


int main( int argc, char *argv[] )
{

    std::cout << "-----------MODES Adaptive Residual Evaluation starts-------------" << std::endl << std::endl;

    std::cout << "Initializing common variables ... " << std::endl << std::endl;
    prob_settings settings;
    std::string filecfg = argv[1];
    std::string su2_conf = argv[2];
    std::string su2dtr_string = "mpirun -np 6 ./SU2_DTR " + su2_conf + " > SU2.log"; // + " > resEval_su2.log";
    int len_s = su2dtr_string.length();
    char su2_sys_call[len_s + 1];
    strcpy(su2_sys_call, su2dtr_string.c_str());

    Read_cfg( filecfg, settings );

    std::string root_outputfile;
    root_outputfile.assign ( settings.out_file, 0, settings.out_file.size() - 4);
    std::string root_inputfile;
    root_inputfile.assign ( settings.in_file, 0, settings.in_file.size() - 4);
    std::string input_format;
    input_format.assign ( settings.in_file, settings.in_file.size() - 3, 3);
    std::string rmf_string = "rm -f " + root_outputfile + "_*";
    len_s = rmf_string.length();
    char rmf_sys_call[len_s + 20];
    strcpy(rmf_sys_call, rmf_string.c_str());

    int s_Nf = 1;   //Number of values for the SPOD filter (POD included)
    std::vector<int> Nf(s_Nf);
    Nf[0] = 0;
    //Nf[1] = std::ceil(settings.Ns/10.0);
    //Nf[2] = std::ceil(settings.Ns/2.0);
    //Nf[3] = std::ceil(2.0*settings.Ns/3.0);
//     Nf[1] = settings.Ns;
    int nfj = 0;

    int nC = settings.Cols.size();
    double alpha = settings.alpha;
    double beta = settings.beta;
    std::vector<double> Dt_res = settings.Dt_res;

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
    //Defining Initial condition
    double M = settings.Mach;
    double Re = settings.Re;
    double T = settings.T;
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

    double nuTilde=0.1*mu/rho;				// initial value spalart allmaras model
    double tke = 1.5*n_turb*n_turb*V_magn*V_magn;
    double omega = rho*tke/(std::max(mu*mu_turb2lam_ratio,1.e-25));
    // double rhotke = rho*tke;
    // double rhoomega = rho*omega;
    double rhotke = tke;
    double rhoomega = omega;

    Eigen::VectorXd mean = sn_set.rowwise().mean();
    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC*Nr);

    if ( nC == 4 ) //Laminar 2D Navier-Stokes
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
    }else if (nC== 5) // Turbolent 2D Spalart Allmaras
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(4*Nr, Nr) = nuTilde*Eigen::MatrixXd::Ones(Nr,1);

    }else if ( nC == 6 ) //Turbulent 2D Navier-Stokes
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
    }

    std::string binary = "YES";
    Eigen::VectorXd svd_cum_sum(settings.Ns);

    if ( settings.flag_mean == "IC" ) {
        std::cout << "Subtracting Initial condition" << std::endl << std::endl;
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;
    } else {
        std::cout << "Using data without subtracting any reference state" << std::endl << std::endl;
    }

//Defining common scope for POD-SPOD
    {
        std::cout << "--------------------------------------" << std::endl ;
        std::cout << "-------Performin POD ResEval----------" << std::endl ;
        std::cout << "--------------------------------------" << std::endl ;
        //Vector of MatrixXd where to store the evolution in time of conservative variables
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
        std::vector<Eigen::MatrixXd> Phi(nC);
        std::vector<Eigen::VectorXd> lambda(nC);

        std::vector< std::vector<rbf> > surr_coefs(nC);
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);
        std::vector<int> Nm(nC);
        int N_notZero;
        //Check only for POD for now

        for ( int ncons = 0; ncons < nC; ncons ++ ) {
            std::cout << "Processing conservative variable " << ncons << std::endl;
            Phi[ncons] = SPOD_basis(sn_set.middleRows(ncons * Nr, Nr),
                                    lambda[ncons], K_pc, eig_vec,
                                    Nf[nfj],
                                    settings.flag_bc,
                                    settings.flag_filter,
                                    settings.sigma);
            N_notZero = Phi[ncons].cols();
            if (settings.r == 0) Nm[ncons] = Nmod(settings.En, K_pc);
            else Nm[ncons] = std::min(settings.r, N_notZero);
            std::cout << "Number of modes used in reconstruction " << Nm[ncons] << std::endl;
            surr_coefs[ncons] = getSurrCoefs(t_vec, eig_vec, settings.flag_interp);
        }
        std::cout << std::endl;

//        std::cout << "Computing SPOD " << Nf[nfj] << " reconstruction for each conservative variable ... " << "\n";

        for ( int idtr = 0; idtr < settings.Dt_res.size(); idtr++ ) {
            std::cout << " --------------DT_RES = " << settings.Dt_res[idtr] << "--------------"<< std::endl;
            for (int itr = 0; itr < settings.t_res.size(); itr++) {
                std::cout << " Computing residuals at time t = " << settings.t_res[itr] << std::endl;
                for (int ncons = 0; ncons < nC; ncons++) {

                    Eigen::MatrixXd coef_t(3, Nm[ncons]);
                    if (settings.t_res[itr] - 2.0 * settings.Dt_res[idtr] < t_vec[0] &&
                        settings.t_res[itr] > t_vec[t_vec.size() - 1] && 2 * settings.Dt_res[idtr] > settings.Dt_cfd) {
                        std::cout
                                << "Define proper Delta_t_res and T_RES vector " << std::endl;
                        exit(EXIT_FAILURE);
                    } else {
                        std::vector<double> tr(1);
                        std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[idtr],
                                                          settings.t_res[itr] - settings.Dt_res[idtr],
                                                          settings.t_res[itr]};

                        for (int j = 0; j < 3; j++) {
                            tr[0] = t_evaluate[j];
                            for (int i = 0; i < Nm[ncons]; i++)
                                surr_coefs[ncons][i].evaluate(tr, coef_t(j, i));
                        }


                    }
                    Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm[ncons], Nm[ncons]);
                    for (int i = 0; i < Nm[ncons]; i++)
                        Sig(i, i) = std::sqrt(lambda[ncons](i));
                    Sn_Cons_time.middleRows(ncons * Nr, Nr) = Phi[ncons].leftCols(Nm[ncons]) * Sig * coef_t.transpose();
                }

                if (settings.flag_mean == "IC") {
                    for (int it = 0; it < 3; it++)
                        Sn_Cons_time.col(it) += Ic;
                }

                std::string mv_string;
                if (settings.Ns == settings.r)
                    mv_string = "mv history_rbm_00002.csv history_pod_AllModes_" +
                                std::to_string(settings.Dt_res[idtr]) + "_" + std::to_string(itr) + ".csv";
                else
                    mv_string = "mv history_rbm_00002.csv history_pod_" +
                                std::to_string(Nm[0]) + "_" + std::to_string(settings.Dt_res[idtr]) + "_" + std::to_string(itr) +
                                ".csv";

                len_s = mv_string.length();
                char mv_sys_call[len_s + 10];
                strcpy(mv_sys_call, mv_string.c_str());
                // Write_Restart_Cons_Time( Sn_Cons_time, Coords, settings.out_file, t_evaluate.size(), nC, alpha );
                Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, 3, nC, alpha, beta, binary);
                //Executing SU2, removing all useless files, renaming files with residuals
                std::cout << "Calling SU2 for residual evaluation and writing file to history " << std::endl;
                std::system(su2_sys_call);
                std::system(rmf_sys_call);
                std::system(mv_sys_call);

                std::cout << std::endl;
            }
        }
    }


//Defining scope for DMD ( Rank=-1 preferable, Coeffs = OPT )
    {
        std::cout << "--------------------------------------" << std::endl ;
        std::cout << "-------Performin DMD ResEval----------" << std::endl ;
        std::cout << "--------------------------------------" << std::endl ;
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
        Eigen::VectorXd lambda_POD;
        Eigen::MatrixXd eig_vec_POD;
        std::vector<Eigen::VectorXcd> lambda_DMD(nC);
        Eigen::MatrixXcd eig_vec_DMD;
        std::vector<Eigen::MatrixXcd> Phi(nC);
        std::vector<Eigen::VectorXcd> alfa(nC);
        std::vector<int> Nm(nC);

        for ( int ncons = 0; ncons < nC; ncons++ ) {

            std::cout << "Processing conservative variable " << ncons << std::endl;
//            std::cout << "Extracting basis DMD using rank " << std::min(settings.r, settings.Ns - 1) << "\t";

            if (settings.r == 0) {
                Phi[ncons] = DMD_basis(sn_set.middleRows(ncons * Nr, Nr),
                                       lambda_DMD[ncons],
                                       eig_vec_DMD,
                                       lambda_POD,
                                       eig_vec_POD,
                                       -1);
            } else {
                Phi[ncons] = DMD_basis(sn_set.middleRows(ncons * Nr, Nr),
                                       lambda_DMD[ncons],
                                       eig_vec_DMD,
                                       lambda_POD,
                                       eig_vec_POD,
                                       settings.r);
            }

            //         int Nm = Phi.cols();
            //         std::cout << "Number of modes extracted : " << Nm << std::endl;

            Eigen::VectorXcd omega(Phi[ncons].cols());
            for (int i = 0; i < Phi[ncons].cols(); i++)
                omega(i) = std::log(lambda_DMD[ncons](i)) / (settings.Dt_cfd * settings.Ds);

            // std::cout << "Calculating coefficients DMD ... " << "\t";
            alfa[ncons] = Calculate_Coefs_DMD_exact(sn_set.middleRows(ncons * Nr, Nr).leftCols(settings.Ns - 1),
                                                    lambda_DMD[ncons],
                                                    Phi[ncons]);
            // std::cout << " Done! " << std::endl;

            // std::cout << "Reordering modes DMD ... " << "\t";
            Eigen::VectorXd En = Eigen::VectorXd::Zero(Phi[ncons].cols());
            double T = t_vec[t_vec.size() - 1];

            for (int i = 0; i < Phi[ncons].cols(); i++) {

                double alfa_i = alfa[ncons](i).imag();
                double alfa_r = alfa[ncons](i).real();
                double sigma = omega(i).real();
                En(i) = (alfa_r * alfa_r + alfa_i * alfa_i) * (std::exp(2.0 * sigma * T) - 1.0) / (2.0 * sigma);

            }

            dmd_sort(En, Phi[ncons], lambda_DMD[ncons], alfa[ncons]);
            // std::cout << "Done" << std::endl;

            double sum = 0;
            Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
            for (int i = 0; i < Phi[ncons].cols(); i++) {
                sum += En(i) / En.sum();
                K_pc(i) = sum;
            }
            if (settings.r == 0) {
                Nm[ncons] = Nmod(settings.En, K_pc);
                std::cout << "Number of modes for the desired energetic content : " << Nm[ncons] << std::endl;
            } else {
                Nm[ncons] = std::min(settings.r, settings.Ns - 1);
                std::cout << "Number of modes (fixed) : " << Nm[ncons] << std::endl;
            }

        }

        std::cout << std::endl;

        for ( int idtr = 0; idtr < settings.Dt_res.size(); idtr++) {

            std::cout << " --------------DT_RES = " << settings.Dt_res[idtr] << "--------------"<< std::endl;
            for (int itr = 0; itr < settings.t_res.size(); itr++) {
                std::cout << "Computing residual at time t = " << settings.t_res[itr] << std::endl;
                std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[idtr],
                                                  settings.t_res[itr] - settings.Dt_res[idtr],
                                                  settings.t_res[itr]};

                for (int ncons = 0; ncons < nC; ncons++) {

                    if (((settings.t_res[itr] - 2.0 * settings.Dt_res[idtr]) < t_vec[0]) &&
                        (settings.t_res[itr] > t_vec[t_vec.size() - 1]) && (2 * settings.Dt_res[idtr] > settings.Dt_cfd)) {
                        std::cout << "Define proper Delta_t_res and T_RES vector " << std::endl;
                        exit(EXIT_FAILURE);
                    } else {

                        for (int j = 0; j < 3; j++) {

                            Eigen::MatrixXcd Rec = Reconstruction_DMD(t_evaluate[j],
                                                                      settings.Dt_cfd * settings.Ds,
                                                                      alfa[ncons],
                                                                      Phi[ncons],
                                                                      lambda_DMD[ncons],
                                                                      "SCALAR");

                            Sn_Cons_time.middleRows(ncons * Nr, Nr).col(j) = Rec.real();

                        }
                    }

                }

                if (settings.flag_mean == "IC") {
                    for (int it = 0; it < 3; it++)
                        Sn_Cons_time.col(it) += Ic;
                }

                std::string mv_string;
                if (settings.Ns == settings.r)
                    mv_string = "mv history_rbm_00002.csv history_dmd_AllModes_" + std::to_string(settings.Dt_res[idtr]) + "_" +
                                std::to_string(itr) + ".csv";
                else
                    mv_string = "mv history_rbm_00002.csv history_dmd_" + std::to_string(Nm[0]) + "_" +
                                std::to_string(settings.Dt_res[idtr]) + "_" + std::to_string(itr) + ".csv";

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

                std::cout << std::endl << std::endl;
            }
        }

    }

// // //Defining scope for RDMD
// // //if using the function RDMD_modes_coefs for energybased select
// // //energy level and rank rdmd to zero, for mode based just select
// //rank rdmd to the number of desired modes
    {

        std::cout << "--------------------------------------" << std::endl ;
        std::cout << "-------Performin RDMD ResEval---------" << std::endl ;
        std::cout << "--------------------------------------" << std::endl ;
        //Vector of MatrixXd where to store the evolution in time of conservative variables
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
        std::vector<Eigen::MatrixXd> Phi(nC);
        std::vector< std::vector<rbf> > surr_coefs(nC);
        std::vector<Eigen::VectorXd> lambda(nC);

        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
        int Nm;
        int N_notZero;

        for ( int ncons = 0; ncons < nC; ncons ++ ) {

                std::cout << "Processing conservative variable " << ncons << std::endl;
                Phi[ncons] = RDMD_modes_coefs(sn_set.middleRows(ncons * Nr, Nr),
                                              Coefs,
                                              lambda[ncons],
                                              K_pc,
                                              -1, //Performing singular value hard threshold for DMD reduction at each step
                                              settings.r_RDMD,
                                              settings.En);

                surr_coefs[ncons] = getSurrCoefs(t_vec,
                                                 Coefs.transpose(),
                                                 settings.flag_interp);

                N_notZero = Phi[ncons].cols();
                if (settings.r == 0) Nm = Nmod(settings.En, K_pc);
                else Nm = std::min(settings.r, N_notZero);
                std::cout << "Number of modes used in reconstruction " << Nm << std::endl;

        }


//        for ( int ncons = 0; ncons < nC; ncons++ ) {
//            std::cout << "Processing conservative variable " << ncons << std::endl;
//            Phi[ncons] = RDMD_modes_coefs(sn_set.middleRows(ncons * Nr, Nr),
//                                          Coefs,
//                                          lambda[ncons],
//                                          K_pc,
//                                          -1, //Performing singular value hard threshold for DMD reduction at each step
//                                          settings.r_RDMD,
//                                          settings.En);
//            surr_coefs[ncons] = getSurrCoefs(t_vec,
//                                             Coefs.transpose(),
//                                             settings.flag_interp);
//            N_notZero = Phi[ncons].cols();
//            if (settings.r == 0) Nm = Nmod(settings.En, K_pc);
//            else Nm = std::min(settings.r, N_notZero);
//            std::cout << "Number of modes used in reconstruction " << Nm << std::endl;
//        }

        for ( int idtr = 0; idtr < settings.t_res.size(); idtr++ ) {
            std::cout << " --------------DT_RES = " << settings.Dt_res[idtr] << "--------------"<< std::endl;
            for (int itr = 0; itr < settings.t_res.size(); itr++) {
                std::cout << " Computing residuals at time = " << settings.t_res[itr] << std::endl;

                for (int ncons = 0; ncons < nC; ncons++) {

                    Eigen::MatrixXd coef_t(3, Nm);
                    if (((settings.t_res[itr] - 2.0 * settings.Dt_res[idtr]) < t_vec[0]) &&
                        (settings.t_res[itr] > t_vec[t_vec.size() - 1]) && (2 * settings.Dt_res[idtr] > settings.Dt_cfd)) {
                        std::cout << "Define proper Delta_t_res and T_RES vector " << std::endl;
                        exit(EXIT_FAILURE);
                    } else {
                        std::vector<double> tr(1);
                        std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[idtr],
                                                          settings.t_res[itr] - settings.Dt_res[idtr],
                                                          settings.t_res[itr]};

                        for (int j = 0; j < 3; j++) {
                            tr[0] = t_evaluate[j];
                            for (int i = 0; i < Nm; i++)
                                surr_coefs[ncons][i].evaluate(tr, coef_t(j, i));
                        }
                    }

                    Sn_Cons_time.middleRows(ncons * Nr, Nr) = Phi[ncons].leftCols(Nm) * coef_t.transpose();

                }

                if (settings.flag_mean == "IC") {
                    for (int it = 0; it < 3; it++)
                        Sn_Cons_time.col(it) += Ic;
                }

                std::string mv_string;
                if (settings.Ns == settings.r)
                    mv_string = "mv history_rbm_00002.csv history_rdmd_AllModes_" +
                                std::to_string(settings.Dt_res[idtr]) + "_" + std::to_string(itr) + ".csv";
                else
                    mv_string = "mv history_rbm_00002.csv history_rdmd_" +
                                std::to_string(Nm) + "_" + std::to_string(settings.Dt_res[idtr]) + "_" + std::to_string(itr) +
                                ".csv";

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

                std::cout << std::endl << std::endl;

            }

        }
    }

    std::cout << "MODES Adaptive Residual Evaluation ends" << std::endl << std::endl;
    return 0;

}