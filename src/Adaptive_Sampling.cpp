#include <cmath>
#include <initializer_list>
#include <iostream>
#include <utility>

// #include <pagmo/problem.hpp>
// #include <pagmo/types.hpp>
#include "Extract_Basis.hpp"
#include "pagmo.hpp"
#include "Opt_struct.hpp"
#include "Generate_snset.hpp"


std::vector<int> dv_rnd (int min, int max, int N);
double fit (Eigen::MatrixXd _sn_set_, Eigen::VectorXd norm_sn_set, std::vector<int> dv);


int main( int argc, char *argv[] )
{

    std::cout << "-----------MODES Adaptive Sampling starts-------------" << std::endl << std::endl;

    //--------------------------------------------------------------------------------------------//
    //-------------------------------Initializing common Variables--------------------------------//
    //--------------------------------------------------------------------------------------------//
    
    std::cout << "Initializing common variables ... " << std::endl << std::endl;
    prob_settings settings;
    std::string filecfg = argv[1];
    std::string decision = argv[2];
    Read_cfg( filecfg, settings );

    int nVar = settings.Cols.size();

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
    std::cout << "Number of grid points : " << Np << std::endl;
    Eigen::MatrixXd Coords = read_col( file_1, Np, settings.Cols_coords );
    std::cout << "Done " << std::endl;

    //Storing complete snapshot matrix
    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Np, settings.Ns, settings.Ds, settings.nstart,
                                                settings.Cols,
                                                settings.in_file,
                                                settings.flag_prob);

    std::cout << "Initializing Vector of times ... " << std::endl; 
    std::vector<double> t_vec( settings.Ns );
    t_vec[0] = 0.0; //Nstart considered as the initial observation time
    for ( int i = 1; i < settings.Ns; i++ )
        t_vec[i] = t_vec[i-1] + settings.Dt_cfd*settings.Ds;
        
    std::cout << "Computing mean of CFD solution ... " << std::endl;
    Eigen::VectorXd mean = sn_set.rowwise().mean();
    if ( settings.flag_mean == "YES" ) {
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= mean;
    }

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

    double alpha = settings.alpha;
    double beta = settings.beta;
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

    double tke = 1.5*n_turb*n_turb*V_magn*V_magn;
    double omega = rho*tke/(std::max(mu*mu_turb2lam_ratio,1.e-25));
    // double rhotke = rho*tke;
    // double rhoomega = rho*omega;
    double rhotke = tke;
    double rhoomega = omega;

    int nC = settings.Cols.size();
    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC*Np);

    if ( nC == 4 ) //Laminar 2D Navier-Stokes
    {
        Ic.head(Np) = rho*Eigen::MatrixXd::Ones(Np,1);
        Ic.segment(Np, Np) = rhoU*Eigen::MatrixXd::Ones(Np,1);
        Ic.segment(2*Np, Np) = rhoV*Eigen::MatrixXd::Ones(Np,1);
        Ic.segment(3*Np, Np) = rhoE*Eigen::MatrixXd::Ones(Np,1);
    } else if ( nC == 6 ) //Turbulent 2D Navier-Stokes
    {
        Ic.head(Np) = rho*Eigen::MatrixXd::Ones(Np,1);
        Ic.segment(Np, Np) = rhoU*Eigen::MatrixXd::Ones(Np,1);
        Ic.segment(2*Np, Np) = rhoV*Eigen::MatrixXd::Ones(Np,1);
        Ic.segment(3*Np, Np) = rhoE*Eigen::MatrixXd::Ones(Np,1);
        Ic.segment(4*Np, Np) = rhotke*Eigen::MatrixXd::Ones(Np,1);
        Ic.segment(5*Np, Np) = rhoomega*Eigen::MatrixXd::Ones(Np,1);
    } else if ( nC == 7 ) //Turbulent 3D Navier-Stokes
    {
        Ic.head(Np) = rho*Eigen::MatrixXd::Ones(Np,1);
        Ic.segment(Np, Np) = rhoU*Eigen::MatrixXd::Ones(Np,1);
        Ic.segment(2*Np, Np) = rhoV*Eigen::MatrixXd::Ones(Np,1);
        Ic.segment(3*Np, Np) = rhoW*Eigen::MatrixXd::Ones(Np,1); //no sideslip angle
        Ic.segment(4*Np, Np) = rhoE*Eigen::MatrixXd::Ones(Np,1);
        Ic.segment(5*Np, Np) = rhotke*Eigen::MatrixXd::Ones(Np,1);
        Ic.segment(6*Np, Np) = rhoomega*Eigen::MatrixXd::Ones(Np,1);
    }

    if ( settings.flag_mean == "IC" )
    {
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;
    }

    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(settings.Ns);
    for ( int it = 0; it < settings.Ns; it ++)
        norm_sn_set(it) = sn_set.col(it).norm();

    //Define normalization for conservative variables
    double rho_max, rho_min, rhoU_max, rhoU_min, rhoV_max, rhoV_min, rhoW_max, rhoW_min, rhoE_max, rhoE_min; //add turbulence

    if ( settings.flag_prob == "CONSERVATIVE" ) {
        //Introduce an if on the number of conservative variables
        Eigen::ArrayXd temp;

        rho_max = sn_set.middleRows(0, Np).maxCoeff();
        rho_min = sn_set.middleRows(0, Np).minCoeff();
        temp = ( sn_set.middleRows(0, Np).array() - rho_min);
        sn_set.middleRows(0, Np) = temp.matrix()/(rho_max - rho_min);

        rhoU_max = sn_set.middleRows(Np, Np).maxCoeff();
        rhoU_min = sn_set.middleRows(Np, Np).minCoeff();
        temp = ( sn_set.middleRows(0, Np).array() - rhoU_min);
        sn_set.middleRows(Np, Np) = temp.matrix()/(rhoU_max - rhoU_min);

        rhoV_max = sn_set.middleRows(2*Np, Np).maxCoeff();
        rhoV_min = sn_set.middleRows(2*Np, Np).minCoeff();
        temp = ( sn_set.middleRows(2*Np, Np).array() - rhoV_min);
        sn_set.middleRows(2*Np, Np) = temp.matrix()/(rhoV_max - rhoV_min);

        rhoE_max = sn_set.middleRows(3*Np, Np).maxCoeff();
        rhoE_min = sn_set.middleRows(3*Np, Np).minCoeff();
        temp = ( sn_set.middleRows(3*Np, Np).array() - rhoE_min);
        sn_set.middleRows(3*Np, Np) = temp.matrix()/(rhoE_max - rhoE_min);

    }


    //Computing svd of all the snapshots for a first guess of nVar
    Eigen::VectorXd lambda = Eigen::VectorXd::Zero(settings.Ns);
    Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
    Eigen::MatrixXd eig_vec = Eigen::MatrixXd::Zero(settings.Ns,settings.Ns);

    Eigen::MatrixXd Phi = SPOD_basis( sn_set,
                                      lambda, K_pc, eig_vec,
                                      settings.Nf,
                                      settings.flag_bc,
                                      settings.flag_filter,
                                      settings.sigma);
//std::cout << "Vector of energies : " << K_pc << std::endl;
    double thr = 1e-2; //threshold for the projection error
    int Nm = Nmod(1.0 - thr * thr, K_pc);
    std::cout << "Number of samples needed from uniform sampling to have EpsP= " << thr << " : " <<  Nm << std::endl << std::endl;

    //--------------------------------------------------------------------------------------------//
    //-------------------------------Setting Optimization Problem---------------------------------//
    //--------------------------------------------------------------------------------------------//

    if ( decision == "-y") {

        std::cout << "Performing adaptive sampling with optimization " << std::endl;
        int nVar_lb = 2;
        int nVar_ub = settings.Ns;
        int nVar_Opt = Nmod(1.0 - thr * thr, K_pc) - 2;//std::round(settings.Ns/3);
        double fVal = 1.0;
        size_t Npop = 150;
        int Ngen = 1000;
        int Nf = 0; //Filter for POD

        pagmo::random_device::set_seed(7); //Set seed for reproducible results
        std::vector<int> iSamp_pos = {};
        std::string in = "-X-";
        std::string out = "-o-";
        std::string lb = "X-";
        std::string ub = "-X";
        std::vector<std::string> Samp_rep(settings.Ns, out);
        Samp_rep[0] = lb;
        Samp_rep[settings.Ns - 1] = ub;

        //Variables for output purposes
        std::vector<double> f_opt_vs_nVar = {};
        double rtol = .01;
        double abstol = 1e-5;
        double tol = rtol * thr;

        while (std::abs(fVal - thr) > tol && nVar_lb != nVar_ub - 1) {

            std::string f_name = "Championf_Samples_nVar_" + std::to_string(nVar_Opt) + ".dat";
            std::ofstream opt_data;
            opt_data.open(f_name.c_str());

            //Initializing some meaningful individuals for population
            std::vector<double> p1(nVar_Opt);
            p1[0] = settings.Dt_cfd * settings.Ds;
            for (int it = 1; it < nVar_Opt; it++) p1[it] = p1[it - 1] + settings.Dt_cfd * settings.Ds *
                                                                        ((double) settings.Ns - 3.0) /
                                                                        ((double) nVar_Opt - 1.0);

            //Clearing and resetting adaptive sampling vectors for representation
            iSamp_pos.clear();
            for (int it = 1; it < (settings.Ns - 1); it++) Samp_rep[it] = out;

            std::cout << "For nVar = " << nVar_Opt << std::endl;
            // Define search bounds:
            std::vector<std::vector<double> > bounds(2, std::vector<double>(nVar_Opt, 0.0));
            for (int iVar_Opt = 0; iVar_Opt < nVar_Opt; iVar_Opt++) {
                bounds[0][iVar_Opt] = settings.Dt_cfd * settings.Ds;
                bounds[1][iVar_Opt] = settings.Dt_cfd * settings.Ds * ((double) settings.Ns - 2.0);
//            bounds[0][iVar_Opt] = 1.0;
//            bounds[1][iVar_Opt] = (double)settings.Ns-2.0;

            }

            pagmo::problem prob{SPOD_Adapt_Samp(bounds, sn_set, settings, settings.Nf)};
//            pagmo::algorithm algo{pagmo::pso(1u,0.7298,2.05,2.05,0.5,
//                    5u,2u,4u,true,1)};
            pagmo::algorithm algo{pagmo::sade(1u,2u,1u,1e-6,1e-6,
                                             true,7)};
//            pagmo::algorithm algo{pagmo::cmaes( 1u, -1.0, -1.0, -1.0, -1.0, 0.5, 1e-6, 1e-6, true, true, 1)};
            pagmo::population pop{prob, Npop, 7};
            algo.set_verbosity(1);
            //Adding uniform sample
            pop.set_x(0, p1);

            //Evolving for Ngen iterations
            for (int i = 0; i < Ngen; i++) {

                pop = algo.evolve(pop);
                pagmo::vector_double f = pop.champion_f();
                pagmo::vector_double T = pop.champion_x();
                std::cout << "Minimum: " << i << " " << std::setprecision(8) << "f= "
                          << f[0] << std::endl;

                opt_data << std::setprecision(8) << f[0];
                for (int it = 0; it < T.size(); it++) {
                    opt_data << "," << std::setprecision(8) << T[it];
                }

                opt_data << std::endl;
                fVal = f[0];

                if ( f[0] < thr ||
                    std::abs(f[0] - thr) < tol)
                    break;
            }

            opt_data.close();
            pagmo::vector_double T_c = pop.champion_x();
//            fVal = isl.get_population().champion_f()[0];

            //Converting decision vector from double to int to get snapshots column selection
//            iSamp_pos.resize(isl.get_population().champion_x().size() +
//                             2); //+2 takes into account first and last snapshots (always included in the sampling)
            iSamp_pos.resize(T_c.size() +
                             2); //+2 takes into account first and last snapshots (always included in the sampling)
            iSamp_pos[0] = 0;
            iSamp_pos[nVar_Opt + 1] = settings.Ns - 1;

            for (int iVar = 1; iVar < (nVar_Opt + 1); iVar++)
                iSamp_pos[iVar] = static_cast<int>(round(
                        T_c[iVar - 1] / (settings.Dt_cfd * settings.Ds)));

            //Representing Adaptive Sampling on screen
            std::cout << "Adaptive Sampling for nVar = " << nVar_Opt << std::endl;
            std::sort(iSamp_pos.begin(), iSamp_pos.end());
            for (int it = 0; it < iSamp_pos.size(); it++) std::cout << iSamp_pos[it] << '\t';
            iSamp_pos.erase(std::unique(iSamp_pos.begin(), iSamp_pos.end()), iSamp_pos.end());

            int count = 1;
            for (int it = 1; it < settings.Ns - 1; it++) {
                if (it == iSamp_pos[count]) {
                    Samp_rep[it] = in;
                    count++;
                }
            }

            std::cout << std::endl;
            for (std::string appo : Samp_rep) std::cout << appo << std::endl << std::endl;

//        Choosing new nVar_Opt with bisection method
            if (fVal > thr) {
                nVar_lb = nVar_Opt;
                nVar_Opt = nVar_lb + std::round((nVar_ub - nVar_lb) / 2);
            } else {
                nVar_ub = nVar_Opt;
                nVar_Opt = nVar_lb + std::round((nVar_ub - nVar_lb) / 2);
            }

            f_opt_vs_nVar.push_back(fVal);

        }
    }

    //--------------------------------------------------------------------------------------------//
    //-------------------------------Setting Problem without  Opt---------------------------------//
    //--------------------------------------------------------------------------------------------//

    else if ( decision == "-n" ) {

        std::cout << "Performing adaptive sampling with NO optimization " << std::endl;
        double fVal = 1.0;
        int nVar_nOpt = 2;
        std::vector<int> indx = {};
        indx.push_back(0);
        indx.push_back(settings.Ns - 1);
        int max_idx = 0;
        int Np = sn_set.rows();

        while (fVal > thr && nVar_nOpt < settings.Ns) {

            Eigen::VectorXd lambda = Eigen::VectorXd::Zero(nVar_nOpt);
            Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(nVar_nOpt);
            Eigen::MatrixXd eig_vec = Eigen::MatrixXd::Zero(nVar_nOpt, nVar_nOpt);

            std::cout << "New vector of samples : " << std::endl;
            for (int it : indx) std::cout << it << "\t";
            std::cout << std::endl;

            Eigen::Map<Eigen::VectorXi> ci(indx.data(), nVar_nOpt);
            Eigen::MatrixXd sub_sn_set_ = indexing(sn_set, Eigen::ArrayXi::LinSpaced(Np, 0, Np - 1), ci);

            std::cout << "Eigen Sampled vector :\n " << ci.transpose() << std::endl;

            Eigen::MatrixXd Phi_temp = SPOD_basis(sub_sn_set_,
                                                  lambda, K_pc, eig_vec,
                                                  settings.Nf,
                                                  settings.flag_bc,
                                                  settings.flag_filter,
                                                  settings.sigma);

            Eigen::MatrixXd ErrP_SPOD_map = Eigen::MatrixXd::Zero(Np, settings.Ns);
            Eigen::VectorXd ErrP_SPOD_time = Eigen::VectorXd::Zero(settings.Ns);

            Eigen::MatrixXd PhiTPhi = Phi_temp.transpose() * Phi_temp;
            Eigen::MatrixXd dumCoefs = Phi_temp.transpose() * sn_set;
            ErrP_SPOD_map = sn_set - Phi_temp * (PhiTPhi.inverse() * dumCoefs);

            for (int it = 0; it < settings.Ns; it++) {
                int count = 0;
                for (int iP = 0; iP < Np; iP++)
                    ErrP_SPOD_time(it) += ErrP_SPOD_map(iP, it) * ErrP_SPOD_map(iP, it);

                ErrP_SPOD_time(it) = std::sqrt(ErrP_SPOD_time(it)) / norm_sn_set(it);
            }

            //Computing final value of the objective function
            fVal = ErrP_SPOD_time.maxCoeff(&max_idx);

            std::cout << "f = " << fVal << std::endl;

            if (fVal > thr) {
                std::cout << "Adding sample at :" << max_idx << std::endl;
                indx.push_back(max_idx);
            } else {
                std::cout << "Threshold reached\n Number of samples needed : " << indx.size()  << std::endl;
                std::cout << "Sample vector : \n";
                std::sort(indx.begin(), indx.end());
                for (int it : indx) std::cout << it << "\t";
                std::cout << std::endl;
            }

            nVar_nOpt++;
        }
    }

    else if ( decision == "-m") {

        std::cout << "Doing Random search" << std::endl;

        Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(settings.Ns);
        for ( int it = 0; it < settings.Ns; it ++)
            norm_sn_set(it) = sn_set.col(it).norm();

        std::string f_name = "Random_search_" + std::to_string(Nm) + ".dat";
        std::ofstream opt_data;
        opt_data.open(f_name.c_str());
        std::vector<double> f_compare = {};

        for ( int it = 0; it < 2000000; it ++ ) {

            std::vector<int> dv = dv_rnd(1, settings.Ns-2, Nm);
//            std::cout << "Random generated vector at it : " << it << std::endl;
//            for ( int appo : dv ) std::cout << appo << "\t";
//            std::cout << std::endl;
            double fV = fit ( sn_set, norm_sn_set, dv );
            f_compare.push_back(fV);

            std::cout << "Fitness at it " << it << " : " << fV << std::endl;
            std::cout << "Smallest Fitness at it " << it << " : " << *std::min_element(f_compare.begin(), f_compare.end()) << '\n';
            opt_data << std::setprecision(8) << fV;
            for ( int xappo : dv ) opt_data << "," << xappo;

            opt_data << std::endl;
        }

        opt_data.close();
    } else if ( decision == "-s") {

        std::cout << "Performing optimization performance study " << std::endl;
        int nVar_lb = 2;
        int nVar_ub = settings.Ns;
        int nVar_Opt = Nmod(1.0 - thr * thr, K_pc) - 2;//std::round(settings.Ns/3);
        double fVal = 1.0;
        size_t Npop = 150;
        int Ngen = 150;
        int Nf = 0; //Filter for POD

        std::vector<unsigned int> v_seed = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

        std::vector<int> iSamp_pos = {};
        std::string in = "-X-";
        std::string out = "-o-";
        std::string lb = "X-";
        std::string ub = "-X";
        std::vector<std::string> Samp_rep(settings.Ns, out);
        Samp_rep[0] = lb;
        Samp_rep[settings.Ns - 1] = ub;

        //Variables for output purposes
        std::vector<double> f_opt_vs_nVar = {};
        double rtol = .1;
        double abstol = 1e-5;
        double tol = rtol * thr;

        for ( int iseed = 0; iseed < v_seed.size(); iseed++ ) {

            pagmo::random_device::set_seed(v_seed[iseed]);
            std::string f_name = "Championf_Samples_seed_" + std::to_string(v_seed[iseed]) + ".dat";
            std::ofstream opt_data;
            opt_data.open(f_name.c_str());

            //Initializing some meaningful individuals for population
            std::vector<double> p1(nVar_Opt);
            p1[0] = settings.Dt_cfd * settings.Ds;
            for (int it = 1; it < nVar_Opt; it++) p1[it] = p1[it - 1] + settings.Dt_cfd * settings.Ds *
                                                                        ((double) settings.Ns - 3.0) /
                                                                        ((double) nVar_Opt - 1.0);

            //Clearing and resetting adaptive sampling vectors for representation
            iSamp_pos.clear();
            for (int it = 1; it < (settings.Ns - 1); it++) Samp_rep[it] = out;

            std::cout << "For seed " << v_seed[iseed] << std::endl;
            // Define search bounds:
            std::vector<std::vector<double> > bounds(2, std::vector<double>(nVar_Opt, 0.0));
            for (int iVar_Opt = 0; iVar_Opt < nVar_Opt; iVar_Opt++) {
                bounds[0][iVar_Opt] = settings.Dt_cfd * settings.Ds;
                bounds[1][iVar_Opt] = settings.Dt_cfd * settings.Ds * ((double) settings.Ns - 2.0);
//            bounds[0][iVar_Opt] = 1.0;
//            bounds[1][iVar_Opt] = (double)settings.Ns-2.0;

            }

            pagmo::problem prob{SPOD_Adapt_Samp(bounds, sn_set, settings, settings.Nf)};
//            pagmo::algorithm algo{pagmo::pso(1u,0.7298,2.05,2.05,0.5,
//                                             5u,2u,4u,true,1)};
            pagmo::algorithm algo{pagmo::sade(1u,2u,1u,1e-6,1e-6,
                                             true,v_seed[iseed])};
//            pagmo::algorithm algo{pagmo::cmaes( 1u, -1.0, -1.0, -1.0, -1.0, 0.5, 1e-6, 1e-6, true, true, 1)};
            pagmo::population pop{prob, Npop, v_seed[iseed]};
            algo.set_verbosity(1);
            //Adding uniform sample
            pop.set_x(0, p1);

            //Evolving for Ngen iterations
            for (int i = 0; i < Ngen; i++) {

                pop = algo.evolve(pop);
                pagmo::vector_double f = pop.champion_f();
                pagmo::vector_double T = pop.champion_x();
                std::cout << "Minimum: " << i << " " << std::setprecision(8) << "f= "
                          << f[0] << std::endl;

                opt_data << std::setprecision(8) << f[0];
                for (int it = 0; it < T.size(); it++) {
                    opt_data << "," << std::setprecision(8) << T[it];
                }

                opt_data << std::endl;
                fVal = f[0];

                if ( f[0] < thr ||
                     std::abs(f[0] - thr) < tol)
                    break;
            }

            opt_data.close();
            pagmo::vector_double T_c = pop.champion_x();

            //Converting decision vector from double to int to get snapshots column selection
            iSamp_pos.resize(T_c.size() +
                             2); //+2 takes into account first and last snapshots (always included in the sampling)
            iSamp_pos[0] = 0;
            iSamp_pos[nVar_Opt + 1] = settings.Ns - 1;

            for (int iVar = 1; iVar < (nVar_Opt + 1); iVar++)
                iSamp_pos[iVar] = static_cast<int>(round(
                        T_c[iVar - 1] / (settings.Dt_cfd * settings.Ds)));

            //Representing Adaptive Sampling on screen
            std::cout << "Adaptive Sampling for seed = " << v_seed[iseed] << std::endl;
            std::sort(iSamp_pos.begin(), iSamp_pos.end());
            for (int it = 0; it < iSamp_pos.size(); it++) std::cout << iSamp_pos[it] << '\t';
            iSamp_pos.erase(std::unique(iSamp_pos.begin(), iSamp_pos.end()), iSamp_pos.end());

            std::cout << std::endl;
            f_opt_vs_nVar.push_back(fVal);

        }
    } else {
        std::cout << "Nothing to do!" << std::endl;
    }

    std::cout << std::endl;
    std::cout << "--------MODES Adaptive Sampling Ends-----------" << std::endl;
    return 0;

}



//Functions for random search

std::vector<int> dv_rnd (int min, int max, int N) {
    std::vector<int> dv(N);
    std::random_device generator;
    std::uniform_int_distribution<int> distribution(min,max);

    for ( int it = 0; it < N; it ++ ) dv[it] = distribution(generator);

    return dv;

}

double fit (Eigen::MatrixXd _sn_set_, Eigen::VectorXd norm_sn_set, std::vector<int> dv ) {

    int Np = _sn_set_.rows();
    int Ns = _sn_set_.cols();

    //Generating the vector of column indices to read in the full snapshot matrix
    int N_Var = static_cast<int>(dv.size());

    std::vector<int> ci_vec(N_Var+2); //+2 takes into account first and last snapshots (always included in the sampling)
    ci_vec[0] = 0;
    ci_vec[N_Var+1] = Ns-1;

    for ( int iVar = 1; iVar<(N_Var+1); iVar++ )
        ci_vec[iVar] = dv[iVar-1];

//    for (int iappo = 0 ; iappo < N_Var+2; iappo++ )
//        std::cout << "civec[" << iappo << "] = " << ci_vec[iappo] << ", ci_eigen[" << iappo << "] = " << ci(iappo) << std::endl;

    std::sort(ci_vec.begin(),ci_vec.end());
    ci_vec.erase(std::unique(ci_vec.begin(),ci_vec.end()),ci_vec.end());
    N_Var = static_cast<int>(ci_vec.size())-2;
    Eigen::Map<Eigen::VectorXi> ci(ci_vec.data(), N_Var+2);

    Eigen::MatrixXd sn_set_test = indexing(_sn_set_, Eigen::ArrayXi::LinSpaced(Np,0,Np-1),ci);

    //Performing basis extraction
    Eigen::VectorXd lambda = Eigen::VectorXd::Zero(N_Var+2);
    Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(N_Var+2);
    Eigen::MatrixXd eig_vec = Eigen::MatrixXd::Zero(N_Var+2,N_Var+2);

    Eigen::MatrixXd Phi = SPOD_basis( sn_set_test,
                                      lambda, K_pc, eig_vec,
                                      0,
                                      "ZERO",
                                      "BOX",
                                      0.0);
    int Nm = Phi.cols();
    //Computing projection error
    Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm, Nm);

    for ( int iMode = 0; iMode < Nm; iMode++ )
        Sig(iMode,iMode) = std::sqrt(lambda(iMode));

    Eigen::MatrixXd ErrP_SPOD_map = Eigen::MatrixXd::Zero(Np, Ns);
    Eigen::VectorXd ErrP_SPOD_time = Eigen::VectorXd::Zero(Ns);

    Eigen::MatrixXd PhiTPhi = Phi.transpose()*Phi;
    Eigen::MatrixXd dumCoefs = Phi.transpose()*_sn_set_;
    ErrP_SPOD_map = _sn_set_ - Phi*(PhiTPhi.inverse()*dumCoefs);

    for ( int it = 0; it < Ns; it++ ) {
        int count = 0;
        for ( int iP = 0; iP < Np; iP++ )
            ErrP_SPOD_time(it) += ErrP_SPOD_map(iP,it)*ErrP_SPOD_map(iP,it);

        ErrP_SPOD_time(it) = std::sqrt(ErrP_SPOD_time(it))/norm_sn_set(it);
    }

    //Computing final value of the objective function
    double fitness_vector = ErrP_SPOD_time.maxCoeff();
    return fitness_vector;

}




