
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <utility>

// #include <pagmo/problem.hpp>
// #include <pagmo/types.hpp>
#include "Extract_Basis.hpp"
#include "Pre-Process.hpp"
#include "pagmo.hpp"
#include "Opt_struct.hpp"
#include "Generate_snset.hpp"

int main( int argc, char *argv[] )
{

    std::cout << "-----------MODES Adaptive Sampling starts-------------" << std::endl << std::endl;

    //--------------------------------------------------------------------------------------------//
    //-------------------------------Initializing common Variables--------------------------------//
    //--------------------------------------------------------------------------------------------//

    std::cout << "Initializing common variables ... " << std::endl << std::endl;
    prob_settings settings;
    std::string filecfg = argv[1];
    Read_cfg( filecfg, settings );

    //Reading Number of grid points and Coordinates
    std::cout << "Reading Number of grid points and Coordinates ... " << std::endl;
    std::string root_inputfile;
    root_inputfile.assign ( settings.in_file, 0, settings.in_file.size() - 4);
    std::string input_format;
    input_format.assign ( settings.in_file, settings.in_file.size() - 3, 3);
    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(5) << std::to_string(settings.nstart);
    std::string file_1 = root_inputfile + "_" + buffer.str() + "." + input_format;
    int Np = N_gridpoints ( file_1 );
    int Nr = Np;
    std::cout << "Number of grid points : " << Np << std::endl;
    Eigen::MatrixXd Coords = read_col( file_1, Np, settings.Cols_coords );
    std::cout << "Done " << std::endl;

    //Storing complete snapshot matrix
    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Np, settings);

    std::cout << "Initializing Vector of times ... " << std::endl;
    std::vector<double> t_vec( settings.Ns );
    t_vec[0] = 0.0; //Nstart considered as the initial observation time
    for ( int i = 1; i < settings.Ns; i++ )
        t_vec[i] = t_vec[i-1] + settings.Dt_cfd*settings.Ds;

    std::cout << "Computing mean of CFD solution ... " << std::endl;


    std::cout << "Computing mean/Initial Condition of CFD solution ... " << std::endl;
    //Defining Initial condition
    Eigen::VectorXd mean = sn_set.rowwise().mean();
    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(Np);

    if ( settings.flag_mean == "YES" ) {
        for ( int it = 0; it < settings.Ns; it++ ) sn_set.col(it) -= mean;
    }

    if ( settings.flag_mean == "IC" ) { Ic = IC(sn_set,settings,1,Np); }

    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(settings.Ns);
    for ( int it = 0; it < settings.Ns; it ++)
        norm_sn_set(it) = sn_set.col(it).norm();

    std::cout << "Line 75 " << std::endl;
    Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
    Eigen::VectorXd lambda_POD = Eigen::VectorXd::Zero(settings.Ns);
    Eigen::MatrixXd eig_vec_POD = Eigen::MatrixXd::Zero(settings.Ns,settings.Ns);
    Eigen::VectorXcd lambda_DMD;
    Eigen::MatrixXcd eig_vec_DMD;
    Eigen::VectorXd lambda_RDMD = Eigen::VectorXd::Zero(settings.Ns);
    Eigen::MatrixXd CoefsRDMD = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
//    Eigen::MatrixXcd PhiDMD = DMD_basis(sn_set,
//                                    lambda_DMD,
//                                    eig_vec_DMD,
//                                    lambda_POD,
//                                    eig_vec_POD,
//                                    -1);

    //Try first with SPOD basis (Results shoould be 1,1,1,1,1,1,1,0,..)
    Eigen::MatrixXd PhiPOD = SPOD_basis(sn_set,
                            lambda_POD, K_pc, eig_vec_POD,
                            0,
                            settings.flag_bc,
                            settings.flag_filter,
                            settings.sigma);


    Eigen::MatrixXd PhiRDMD = RDMD_modes_coefs(sn_set,
                                             CoefsRDMD,
                                             lambda_RDMD,
                                             K_pc,
                                             -1, //Performing singular value hard threshold for DMD reduction at each step
                                             settings.r_RDMD,
                                             settings.En);

    Eigen::MatrixXd PHI = PhiPOD.leftCols(settings.r);
    double tol = 1e-10;
    std::cout << "All modes extracted " << std::endl;

    for ( int iAdd = 0; iAdd < settings.r_RDMD; iAdd++ ){
        VectorXd appo = PhiRDMD.col(iAdd) - PHI*(PHI.transpose()*PhiRDMD.col(iAdd));
        std::cout << "Norm of the mode to be added : " << appo.norm() << std::endl;
        if ( appo.norm() < tol ) {
            std::cout << "Mode " << iAdd << " not added" << std::endl;
            continue;
        } else {
            appo = appo/appo.norm();
            PHI.conservativeResize(PHI.rows(),PHI.cols()+1);
            PHI.col(PHI.cols()-1) = appo;
        }
    }

    std::cout << "Number of modes in the final structure : " << PHI.cols() << std::endl;
    std::cout << "Modes should be all orthogonal :\n " << PHI.transpose()*PHI << std::endl;


    Eigen::MatrixXd CoefsCombo = PHI.transpose()*sn_set;
    Eigen::MatrixXd CoefsPOD = PhiPOD.transpose()*sn_set;

    Eigen::MatrixXd ErrMap_Combo = sn_set - PHI*CoefsCombo;
//    Eigen::MatrixXd ErrMap_POD = sn_set - PhiPOD.leftCols(PHI.cols())*CoefsPOD.topRows(PHI.cols());
    Eigen::MatrixXd ErrMap_POD = sn_set - PhiPOD.leftCols(settings.r)*CoefsPOD.topRows(settings.r);
    Eigen::MatrixXd ErrMap_RDMD = sn_set - PhiRDMD.leftCols(settings.r_RDMD)*CoefsRDMD.topRows(settings.r_RDMD);
    Eigen::VectorXd Errtime_Combo = Eigen::VectorXd::Zero(settings.Ns);
    Eigen::VectorXd Errtime_POD = Eigen::VectorXd::Zero(settings.Ns);
    Eigen::VectorXd Errtime_RDMD = Eigen::VectorXd::Zero(settings.Ns);

    for ( int it = 0; it < settings.Ns; it++ ){
        for ( int iPoint = 0; iPoint < Np; iPoint++ ){
            Errtime_Combo(it) += ErrMap_Combo(iPoint,it)*ErrMap_Combo(iPoint,it);
            Errtime_POD(it) += ErrMap_POD(iPoint,it)*ErrMap_POD(iPoint,it);
            Errtime_RDMD(it) += ErrMap_RDMD(iPoint,it)*ErrMap_RDMD(iPoint,it);
        }

        Errtime_Combo(it) = std::sqrt(Errtime_Combo(it))/norm_sn_set(it);
        Errtime_POD(it) = std::sqrt(Errtime_POD(it))/norm_sn_set(it);
        Errtime_RDMD(it) = std::sqrt(Errtime_RDMD(it))/norm_sn_set(it);
    }

    std::cout << "Writing error of interpolation and error from projection ... " << std::endl;

    std::ofstream errfile;
    errfile.open("Err_combo_POD.dat");
    errfile << "T(s), ErrPOD, ErrRDMD, ErrCombo " << std::endl;

    for ( int it = 0; it < settings.Ns; it ++ )
            errfile << t_vec[it] << ", " << std::setprecision(8) << Errtime_POD(it) << ", " << std::setprecision(8) << Errtime_RDMD(it) << ", " << std::setprecision(8) << Errtime_Combo(it) << std::endl;

    errfile.close();


//    std::cout << "Number of DMD modes extracted : " << Phi.cols() << std::endl;

//    Eigen::MatrixXcd PhiTPhi = Phi.transpose()*Phi;
//    Eigen::MatrixXcd Coeffs = PhiTPhi.inverse()*(Phi.transpose()*sn_set);

//    Eigen::MatrixXd Coeffs = Phi.transpose()*sn_set;
//    std::cout << "Performing adaptive sampling with optimization " << std::endl;
//    int nVar_lb = 0;
//    int nVar_ub = 1;
//    int nVar_Opt = Phi.cols();
//    double fVal = 1.0;
//    size_t Npop = 300;
//    int Ngen = 1000 ;
//
//    pagmo::random_device::set_seed(7); //Set seed for reproducible results
//
//    //Variables for output purposes
//    std::vector<double> f_opt_vs_nVar = {};
//    double thr = 1e-5;
//
//    std::string f_name = "Championf_Samples_Rank_" + std::to_string(nVar_Opt) + ".dat";
//    std::ofstream opt_data;
//    opt_data.open(f_name.c_str());
//
//    std::cout << "For nVar = " << nVar_Opt << std::endl;
//    // Define search bounds:
//    std::vector<std::vector<double> > bounds(2, std::vector<double>(nVar_Opt, 0.0));
//    for (int iVar_Opt = 0; iVar_Opt < nVar_Opt; iVar_Opt++) {
//        bounds[0][iVar_Opt] = 0.0;
//        bounds[1][iVar_Opt] = 1.0;
//    }
//
////    pagmo::problem prob{DMD_Best_Modes(bounds, sn_set, Phi, Coeffs, settings)};
//    pagmo::problem prob{ Best_Modes(bounds, sn_set, Phi, Coeffs, settings) };
//    pagmo::gaco uda{1u, 50u, 1.0, 0.0, 0.0, 900u, 24u, 10000000u, 10000000u, 0.0, true, 7};
////    pagmo::sade uda{1u,2u,1u,1e-6,1e-6,true,7};
//    uda.set_verbosity(1u);
//    uda.set_seed(7);
//
//    //I parallalize the population evaluations:
////    uda.set_bfe( pagmo::bfe{} );
//    pagmo::population pop{prob, Npop, 7};
//
//    //Evolving for Ngen iterations
//    for (int i = 0; i < Ngen; i++) {
//
////              pop = algo.evolve(pop);
//        pop=uda.evolve(pop);
//        pagmo::vector_double f = pop.champion_f();
//        pagmo::vector_double T = pop.champion_x();
//
//        opt_data << std::setprecision(8) << f[0];
//
//        int Nmodes = 0;
//        for (int it = 0; it < T.size(); it++) {
//            opt_data << "," << std::setprecision(8) << T[it];
//            Nmodes += static_cast<int> (T[it]);
//        }
//
//        std::cout << "Minimum: " << i << " " << std::setprecision(8) << "f= "
//                  << f[0] << " using a number of modes : " << Nmodes << std::endl;
//        opt_data << std::endl;
//        fVal = f[0];
//
//        if ( f[0] < thr ) break;
//    }
//
//    opt_data.close();
//    pagmo::vector_double T_c = pop.champion_x();
//
//    std::cout << std::endl;
//    std::cout << "--------MODES Adaptive Sampling Ends-----------" << std::endl;
    return 0;

}






















