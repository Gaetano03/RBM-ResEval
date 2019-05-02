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
    std::string su2dtr_string = "mpirun -np 6 ./SU2_DTR " + su2_conf; // + " > resEval_su2.log";
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
    std::vector<double> t_evaluate(settings.Ns*settings.Ds-1);

    t_evaluate[0] = settings.nstart*settings.Dt_cfd;
    for ( int i = 1; i < t_evaluate.size(); i++)
        t_evaluate[i] = t_evaluate[i-1] + settings.Dt_cfd;// + settings.Dt_cfd*(double)settings.Ds/2.0;

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
    double rhoU = rho*V_magn*std::cos(alpha);
    double rhoV = rho*V_magn*std::sin(alpha);
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
    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC*Nr);

    if ( nC == 4 ) //Laminar 2D Navier-Stokes
    {
        Ic.head(Nr) = rho*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(Nr, Nr) = rhoU*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(2*Nr, Nr) = rhoV*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(3*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
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
        Ic.segment(3*Nr, Nr) = Eigen::MatrixXd::Zero(Nr,1); //no sideslip angle
        Ic.segment(4*Nr, Nr) = rhoE*Eigen::MatrixXd::Ones(Nr,1);
        Ic.segment(5*Nr, Nr) = rhotke*Eigen::MatrixXd::Ones(Nr,1); 
        Ic.segment(6*Nr, Nr) = rhoomega*Eigen::MatrixXd::Ones(Nr,1);        
    }

    std::string binary = "NO";
    Eigen::VectorXd svd_cum_sum(settings.Ns);

    if ( settings.flag_mean == "YES" )
    {
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;
    }

//Defining common scope for POD-SPOD
    {

        //Vector of MatrixXd where to store the evolution in time of conservative variables
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, t_evaluate.size());

        Eigen::VectorXd lambda(settings.Ns);
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);
        int Nm;
        int N_notZero;
        //Check only for POD for now
        for ( int nfj = 0; nfj < Nf.size(); nfj++ )
        {
            std::string mv_string = "mv history_rbm.csv history_spod_" + std::to_string(nfj) + ".csv";
            len_s = mv_string.length();
            char mv_sys_call[len_s + 1];
            strcpy(mv_sys_call, mv_string.c_str());
            std::cout << "Computing SPOD " << Nf[nfj] << " reconstruction for each conservative variable ... " << "\n";
            
            for ( int ncons = 0; ncons < nC; ncons ++ )
            {
                std::cout << "Processing conservative variable " << ncons << std::endl;
                Eigen::MatrixXd Phi = SPOD_basis( sn_set.middleRows(ncons*Nr,Nr),
                                        lambda, K_pc, eig_vec,
                                        Nf[nfj],
                                        settings.flag_bc, 
                                        settings.flag_filter,  
                                        settings.sigma);            
                N_notZero = Phi.cols();
                if ( settings.r == 0 ) Nm = Nmod(settings.En, K_pc);
                else Nm = std::min(settings.r, N_notZero);
                std::cout << "Number of modes used in reconstruction " << Nm << std::endl;
                std::vector<rbf> surr_coefs =  getSurrCoefs (t_vec, eig_vec, settings.flag_interp);                
                Eigen::MatrixXd coef_t(t_evaluate.size(), Nm);

                std::vector<double> tr(1);
                for ( int j = 0; j < t_evaluate.size(); j++ )
                {    
                    tr[0] = t_evaluate[j];
                    for ( int i = 0; i < Nm; i++ )
                        surr_coefs[i].evaluate(tr, coef_t(j,i));
                }
                
                Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm, Nm);
                for ( int i = 0; i < Nm; i++ )
                    Sig(i,i) = std::sqrt(lambda(i));
                
                Sn_Cons_time.middleRows(ncons*Nr,Nr) = Phi.leftCols(Nm)*Sig*coef_t.transpose();
                // Sn_Cons_time.middleRows(ncons*Nr,Nr) = Phi.leftCols(Nm)*Sig*eig_vec.leftCols(Nm).transpose();
            }
            //Add mean or Initial condition if it is subtracted
            if ( settings.flag_mean == "YES" )
            {   
                for ( int it = 0; it < t_evaluate.size(); it++ )
                    Sn_Cons_time.col(it) += Ic;
            }
            std::cout << "Writing time reconstruction " << std::endl;
            // Write_Restart_Cons_Time( Sn_Cons_time, Coords, settings.out_file, t_evaluate.size(), nC, alpha );       
            Write_Restart_Cons_Time( Sn_Cons_time, Coords, settings.out_file, t_evaluate.size(), nC, alpha, binary );
            //Executing SU2, removing all useless files, renaming files with residuals
            std::cout << "Calling SU2 for residual evaluation and writing file to history " << std::endl;
            std::system( su2_sys_call );
            std::system( rmf_sys_call );
            std::system( mv_sys_call );

        }
    
    }

//Take care only of this first part for now
//---------------------------------------------------------------------------

// //     for ( int nt = 0; nt < settings.Ns; nt++ )
// //         sn_set.col(nt) += mean;
    
// //     for ( int nt = 0; nt < settings.Ns-1; nt++ )
// //         sn_set_check.col(nt) += mean;

//Defining scope for DMD ( Rank=-1 preferable, Coeffs = OPT )
    {
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, t_evaluate.size());
        Eigen::VectorXd lambda_POD;
        Eigen::MatrixXd eig_vec_POD;
        Eigen::VectorXcd lambda_DMD;
        Eigen::MatrixXcd eig_vec_DMD;
        std::string mv_string = "mv history_rbm.csv history_dmd.csv";
        len_s = mv_string.length();
        char mv_sys_call[len_s + 1];
        strcpy(mv_sys_call, mv_string.c_str());

        for ( int ncons = 0; ncons < nC; ncons ++ )
        {
            std::cout << "Processing conservative variable " << ncons << std::endl;     

            std::cout << "Extracting basis DMD using rank " << std::min(settings.r, settings.Ns-1) << "\t";        
            Eigen::MatrixXcd Phi;
            Eigen::VectorXcd alfa;

            if ( settings.r == 0 )
            {
                Phi = DMD_basis( sn_set.middleRows(ncons*Nr,Nr),
                                lambda_DMD,
                                eig_vec_DMD,
                                lambda_POD,
                                eig_vec_POD,
                                -1 );
            }
            else
            {
                Phi = DMD_basis( sn_set.middleRows(ncons*Nr,Nr),
                                lambda_DMD,
                                eig_vec_DMD,
                                lambda_POD,
                                eig_vec_POD,
                                settings.r );
            }

    //         int Nm = Phi.cols();
    //         std::cout << "Number of modes extracted : " << Nm << std::endl;

            Eigen::VectorXcd omega(Phi.cols());
            for ( int i = 0; i < Phi.cols(); i++ )
                omega(i) = std::log(lambda_DMD(i))/(settings.Dt_cfd*settings.Ds);

            // std::cout << "Calculating coefficients DMD ... " << "\t";            
            alfa = Calculate_Coefs_DMD_exact ( sn_set.middleRows(ncons*Nr,Nr).leftCols(settings.Ns-1),  
                                                                lambda_DMD, 
                                                                Phi );
            // std::cout << " Done! " << std::endl;
            
            // std::cout << "Reordering modes DMD ... " << "\t";
            Eigen::VectorXd En = Eigen::VectorXd::Zero(Phi.cols());
            double T = t_vec[t_vec.size()-1];

            for ( int i = 0 ; i < Phi.cols(); i ++ )
            {

                double alfa_i = alfa(i).imag();
                double alfa_r = alfa(i).real();
                double sigma = omega(i).real();
                En(i) = (alfa_r*alfa_r + alfa_i*alfa_i)*(std::exp(2.0*sigma*T) - 1.0)/(2.0*sigma);

            }

            dmd_sort( En, Phi, lambda_DMD, alfa);
            // std::cout << "Done" << std::endl;

            double sum = 0;
            Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
            for (int i = 0; i < Phi.cols(); i++)
            {
                sum += En(i)/En.sum();
                K_pc(i) = sum;
            }
            int Nm;
            if ( settings.r == 0)
            {
                Nm = Nmod(settings.En, K_pc);
                std::cout << "Number of modes for the desired energetic content : " << Nm << std::endl;   
            }
            else
            {
                Nm = std::min(settings.r,settings.Ns-1);
                std::cout << "Number of modes (fixed) : " << Nm << std::endl;
            }
            

            Eigen::MatrixXcd V_and(lambda_DMD.size(), t_evaluate.size());      
            for ( int i = 0; i < lambda_DMD.size(); i++ )
            {
                for ( int j = 0; j < t_evaluate.size(); j++ )
                    V_and(i,j) = std::pow(lambda_DMD(i), (double)j/(double)settings.Ds);                                                                                         
            }        
            Eigen::MatrixXcd Psi = Eigen::MatrixXcd::Zero(alfa.size(), t_evaluate.size());
            for ( int i = 0; i < t_evaluate.size(); i++ )
                Psi.col(i) = alfa.cwiseProduct(V_and.col(i));


            Eigen::MatrixXcd D_dmd = Phi.leftCols(Nm)*Psi.topRows(Nm);
            Sn_Cons_time.middleRows(ncons*Nr,Nr) = D_dmd.real();
        }
        
        if ( settings.flag_mean == "YES" )
        {   
            for ( int it = 0; it < t_evaluate.size(); it++ )
                Sn_Cons_time.col(it) += Ic;
        }
        std::cout << "Writing time reconstruction " << std::endl;
        // Write_Restart_Cons_Time( Sn_Cons_time, Coords, settings.out_file, t_evaluate.size(), nC, alpha );       
        Write_Restart_Cons_Time( Sn_Cons_time, Coords, settings.out_file, t_evaluate.size(), nC, alpha, binary );
        //Executing SU2, removing all useless files, renaming files with residuals
        std::cout << "Calling SU2 for residual evaluation and writing file to history " << std::endl;
        std::system( su2_sys_call );
        std::system( rmf_sys_call );
        std::system( mv_sys_call );

    }
    

//     for ( int nt = 0; nt < settings.Ns; nt++ )
//         sn_set.col(nt) -= mean;

//     for ( int nt = 0; nt < settings.Ns-1; nt++ )
//         sn_set_check.col(nt) -= mean;

// //Defining scope for RDMD
// //if using the function RDMD_modes_coefs for energybased select 
// //energy level and rank rdmd to zero, for mode based just select
//rank rdmd to the number of desired modes
    {
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, t_evaluate.size());
        std::string mv_string = "mv history_rbm.csv history_rdmd.csv";
        len_s = mv_string.length();
        char mv_sys_call[len_s + 1];
        strcpy(mv_sys_call, mv_string.c_str());

        Eigen::VectorXd lambda = Eigen::VectorXd::Zero(settings.Ns);
        Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
        Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
        Eigen::MatrixXd Phi;

        for ( int ncons = 0; ncons < nC; ncons++ )
        {
            std::cout << "Processing conservative variable " << ncons << std::endl;
            std::cout << "Extracting basis and Coeffs RDMD ... " << "\t";        
        
            Phi = RDMD_modes_coefs ( sn_set.middleRows(ncons*Nr,Nr),
                                    Coefs,
                                    lambda,
                                    K_pc,     
                                    -1, //Performing singular value hard threshold for DMD reduction at each step
                                    settings.r_RDMD,
                                    settings.En );
                
            // else
            // {
            //     std::cout << "Reading basis and extracting Coeffs RDMD ... " << "\t";
            //     std::string file_modes = argv[2];
            //     std::string file_coefs = argv[3];
            //     std::string file_En = argv[4];
            //     Phi = read_modes( file_modes, settings.ndim*Nr, settings.Ns );
            //     Coefs = read_coefs( file_coefs, settings.Ns, settings.Ns );


            //     std::ifstream En_data;
            //     En_data.open( file_En );
            //     if ( !En_data.is_open() )
            //     {
            //         std::cout << "File : " << file_En << " not found" << std::endl;    
            //         exit (EXIT_FAILURE);
            //     }
            //     std::string line_flow_data ;
            //     getline( En_data, line_flow_data );
            //     std::istringstream iss(line_flow_data);
            //     std::string token;

            //     int count = 0;
            //     while( getline( iss, token, ' ') && count < K_pc.size() )
            //     {
            //         K_pc(count) = std::stod(token);
            //         count ++;
            //     } 
            //     En_data.close();
            // }

            int Nm;
            if ( settings.r == 0 )
            {
                Nm = Nmod(settings.En, K_pc);
                std::cout << "number of modes for the desired energetic content " << Nm << std::endl;
            }
            else
            {
                Nm = std::min(settings.r,settings.Ns-1);
                std::cout << "number of modes (fixed) " << Nm << std::endl;
            }

            std::vector<rbf> surr_coefs =  getSurrCoefs (t_vec,
                                                        Coefs.transpose(),
                                                        settings.flag_interp);
            
            Eigen::MatrixXd coef_t(t_evaluate.size(), Nm);

            std::vector<double> tr(1);
            for ( int j = 0; j < t_evaluate.size(); j++ )
            {    
                tr[0] = t_evaluate[j];
                for ( int i = 0; i < Nm; i++ )
                    surr_coefs[i].evaluate(tr, coef_t(j,i));
            }

            // Err_RDMD_map = sn_set_check - Phi.leftCols(Nm)*coef_t.transpose();
            Sn_Cons_time.middleRows(ncons*Nr,Nr) = Phi.leftCols(Nm)*coef_t.transpose();
        }
                    //Add mean or Initial condition if it is subtracted
        if ( settings.flag_mean == "YES" )
        {   
            for ( int it = 0; it < t_evaluate.size(); it++ )
                Sn_Cons_time.col(it) += Ic;
        }
        std::cout << "Writing time reconstruction " << std::endl;
        // Write_Restart_Cons_Time( Sn_Cons_time, Coords, settings.out_file, t_evaluate.size(), nC, alpha );       
        Write_Restart_Cons_Time( Sn_Cons_time, Coords, settings.out_file, t_evaluate.size(), nC, alpha, binary );
        //Executing SU2, removing all useless files, renaming files with residuals
        std::cout << "Calling SU2 for residual evaluation and writing file to history " << std::endl;
        std::system( su2_sys_call );
        std::system( rmf_sys_call );
        std::system( mv_sys_call );

    }

    std::cout << "MODES Adaptive Residual Evaluation ends" << std::endl << std::endl;
    return 0;

}