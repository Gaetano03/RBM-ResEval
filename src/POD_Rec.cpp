/*
Code for adaptive reconstruction based on residual evaluation
Input config file + error file (+ Modes,Coefs and Encontent RDMD if already available)

Output reconstructed field at the desired time instants with the adaptive technique
based on residual evaluation
*/

#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"


int main( int argc, char *argv[] )
{
    std::cout << "Adaptive Reconstruction with ResEval starts " << std::endl;

    prob_settings settings;
    std::string filecfg = argv[1];
    
    //Reading configuration file
    Read_cfg( filecfg, settings );
    double t_0 = settings.nstart*settings.Dt_cfd;
    double alpha = settings.alpha;
    double beta  = settings.beta;

    int s_Nf = 5;
    int Nmethods = s_Nf + 2;
    std::vector<int> Nf(s_Nf);
    Nf[0] = 0;
    Nf[1] = std::ceil(settings.Ns/10.0);
    Nf[2] = std::ceil(settings.Ns/2.0);
    Nf[3] = std::ceil(2.0*settings.Ns/3.0);
    Nf[4] = settings.Ns;
    int Nf_SPOD = 0;

    std::string root_inputfile;
    root_inputfile.assign ( settings.in_file, 0, settings.in_file.size() - 4);
    std::string input_format;
    input_format.assign ( settings.in_file, settings.in_file.size() - 3, 3);

    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(5) << std::to_string(settings.nstart);
    std::string file_1 = root_inputfile + "_" + buffer.str() + "." + input_format;

    int Nr = N_gridpoints ( file_1 );
    std::cout << "Number of grid points : " << Nr << std::endl;

    std::cout << "Reading Coordinates ... \t ";
    Eigen::MatrixXd Coords = read_col( file_1, Nr, settings.Cols_coords );
    std::cout << "Done " << std::endl;

    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings.Ns, settings.Ds, settings.nstart,
                                        settings.Cols,
                                        settings.in_file,
                                        settings.flag_prob);
    // Eigen::MatrixXd sn_set = Eigen::MatrixXd::Zero(settings.ndim*Nr, settings.Ns);
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

    double tke = 1.5*n_turb*n_turb*V_magn*V_magn;
    double omega = rho*tke/(std::max(mu*mu_turb2lam_ratio,1.e-25));
    // double rhotke = rho*tke;
    // double rhoomega = rho*omega;
    double rhotke = tke;
    double rhoomega = omega;

    Eigen::VectorXd mean = sn_set.rowwise().mean();
    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(settings.ndim*Nr);

    if ( settings.ndim == 2 )
    {
        Ic.head(Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
    } else 
    {
        Ic.head(Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoW*Eigen::MatrixXd::Ones(Nr,1);
    }

    if ( settings.flag_mean == "YES" )
    {
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;
    }

    // if ( settings.flag_mean == "YES" )
    // {
    //     for ( int it = 0; it < settings.Ns; it++ )
    //         sn_set.col(it) -= mean;
    // }

//POD reconstruction on each selected time step

    for ( int i = 0; i < settings.t_rec.size(); i++ )
    {

        Eigen::VectorXd Rec_rhoU(Nr);
        Eigen::VectorXd Rec_rhoV(Nr);
        Eigen::VectorXd Rec_rhoW(Nr);

        std::vector<int> pos = {};
        std::cout << " POD reconstruction at time : " << settings.t_rec[i] << std::endl;

        for ( int iDim = 0; iDim < settings.ndim; iDim ++ )
        {


            Eigen::VectorXd lambda(settings.Ns);
            Eigen::VectorXd K_pc(settings.Ns);
            Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);        

            Eigen::MatrixXd Phi = SPOD_basis( sn_set.middleRows(iDim*Nr,Nr),
                                    lambda, K_pc, eig_vec,
                                    Nf_SPOD,
                                    settings.flag_bc, 
                                    settings.flag_filter,  
                                    settings.sigma);

            int Nm;

            if ( settings.r == 0 )
            {
                Nm = Nmod(settings.En, K_pc);
                std::cout << "Number of modes for desired energetic content: " << Nm << std::endl;
            }
            else
            {
                Nm = std::min(settings.r,settings.Ns);
                std::cout << "Number of modes (fixed): " << Nm << std::endl;
            }

            std::vector<double> t_v( settings.Ns );
            t_v[0] = (double)settings.nstart*settings.Dt_cfd;

            for ( int kt = 1; kt < settings.Ns; kt++ )
                t_v[kt] = t_v[kt-1] + settings.Dt_cfd*(double)settings.Ds;

            Eigen::MatrixXd Rec = Reconstruction_S_POD ( t_v,
                                K_pc, lambda, eig_vec.transpose(),
                                Phi, settings.t_rec[i],
                                Nm,
                                "SCALAR",
                                settings.flag_interp ) ;

            if ( iDim == 0 && settings.flag_mean == "YES" ) Rec_rhoU = Rec.col(0) + rhoU*Eigen::MatrixXd::Ones(Nr,1);
            if ( iDim == 1 && settings.flag_mean == "YES" ) Rec_rhoV = Rec.col(0) + rhoV*Eigen::MatrixXd::Ones(Nr,1);
            if ( iDim == 2 && settings.flag_mean == "YES" ) Rec_rhoW = Rec.col(0) + rhoW*Eigen::MatrixXd::Ones(Nr,1);

        }
        
        Eigen::MatrixXd Rec_M(Nr, settings.ndim); 
        
        if ( settings.ndim == 2)
        {
            Rec_M.col(0) = Rec_rhoU;
            Rec_M.col(1) = Rec_rhoV;
        } else
        {
            Rec_M.col(0) = Rec_rhoU;
            Rec_M.col(1) = Rec_rhoV;
            Rec_M.col(2) = Rec_rhoW;
        }
        
        std::cout << "Writing reconstructed field ..." << "\t";
        if (settings.ndim == 2) write_Reconstructed_fields ( Rec_M, Coords, settings.out_file, "VECTOR-2D", i );
        if (settings.ndim == 3) write_Reconstructed_fields ( Rec_M, Coords, settings.out_file, "VECTOR-3D", i );
        std::cout << "Done" << std::endl << std::endl << std::endl;

    }

    std::cout << "POD Reconstruction ends " << std::endl;

    return 0;
}