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
#include "Pre-Process.hpp"

int main( int argc, char *argv[] )
{
    std::cout << "Adaptive Reconstruction with ResEval starts " << std::endl;

    prob_settings settings;
    std::string filecfg = argv[1];
    
    //Reading configuration file
    Read_cfg( filecfg, settings );
    double t_0 = settings.nstart*settings.Dt_cfd;
    double alpha = settings.alpha;
    double beta = settings.beta;
    if (settings.ndim == 2) beta = 0.0;

    int s_Nf = 1;
    int Nmethods = s_Nf + 2;
    std::vector<int> Nf(s_Nf);
    Nf[0] = 0;
//    Nf[1] = std::ceil(settings.Ns/10.0);
//    Nf[2] = std::ceil(settings.Ns/2.0);
//    Nf[3] = std::ceil(2.0*settings.Ns/3.0);
//    Nf[4] = settings.Ns;
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
    //Defining Mean/Initial condition
    Eigen::VectorXd mean = sn_set.rowwise().mean();
    int nC = settings.ndim;
    Eigen::VectorXd Ic = IC( settings, nC, Nr );

    if ( settings.flag_mean == "YES" ) {
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= Ic;
    }

    // if ( settings.flag_mean == "YES" )
    // {
    //     for ( int it = 0; it < settings.Ns; it++ )
    //         sn_set.col(it) -= mean;
    // }

    std::cout << "Reading Residuals ... " << std::endl;

    std::vector<std::string> resfilename = {"history_pod.csv", "history_dmd.csv", "history_rdmd.csv"};
//    std::vector<std::string> resfilename = {"history_spod_0.csv", "history_spod_1.csv", "history_spod_2.csv",
//        "history_spod_3.csv", "history_spod_4.csv", "history_dmd.csv", "history_rdmd.csv"};


    Eigen::MatrixXd Err_RBM = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXd Err_RBM_rhoV = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXd Err_RBM_rhoU = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXd Err_RBM_rhoW = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);

    for ( int i = 0; i < resfilename.size(); i++ )
    {
        std::ifstream file_data;
        file_data.open( resfilename[i] );
            if ( !file_data.is_open() )
        {
            std::cout << "File : " << resfilename[i] << " not found" << std::endl;    
            exit (EXIT_FAILURE);
        }

        std::string line_flow_data ;
        int n_row = 0, count = 0;
        // Reading row of headers
        getline( file_data, line_flow_data );

        while ( getline( file_data, line_flow_data ) && n_row <  settings.Ns*settings.Ds - 1)
        {
            
            std::istringstream iss(line_flow_data);
            std::string token;
            double err;
            count = 0; 
            while( getline( iss, token, ',') )
            {
                err = std::stod(token);

                if ( count == 17 ) Err_RBM_rhoU(n_row, i) = std::pow(10.0, err);
                if ( count == 18 ) Err_RBM_rhoV(n_row, i) = std::pow(10.0, err);
                if ( settings.ndim == 3 && count == 19 ) Err_RBM_rhoW(n_row, i) = std::pow(10.0, err);
                
                count ++;
            } 

            n_row++;
        }

        file_data.close();

    }

//Adaptive reconstruction on each selected time step
    int best_method_idx;

    std::cout << "Initializing Vector of time ... " << std::endl;
    Eigen::VectorXd t_vec( settings.t_res.size());
    for ( int it = 0; it < settings.t_res.size(); it++ ) t_vec(it) = settings.t_res[it];
//    Eigen::VectorXd t_vec( settings.Ns*settings.Ds - 1);
//    t_vec(0) = (double)settings.nstart*settings.Dt_cfd;
//    for ( int i = 1; i < settings.Ns*settings.Ds-1; i++ )
//        t_vec(i) = t_vec(i-1) + settings.Dt_cfd;
    
    double tol = settings.tol;
    int index1, index2;
    Eigen::VectorXd Err_interp(Nmethods);

    for ( int i = 0; i < settings.t_rec.size(); i++ ) {

        Eigen::VectorXd Rec_rhoU(Nr);
        Eigen::VectorXd Rec_rhoV(Nr);
        Eigen::VectorXd Rec_rhoW(Nr);

        std::vector<int> pos = {};
        std::cout << " Adaptive reconstruction at time : " << settings.t_rec[i] << std::endl;

        index1 = 0;
        index2 = 0;
        for ( int nt = 0; nt < t_vec.size()-1; nt ++ ) {
            if ( (settings.t_rec[i] >= t_vec(nt)) && (settings.t_rec[i] <= t_vec(nt+1)) ) {
                index1 = nt;
                index2 = nt+1;
                break;
            }
        }

        if ( index1 == index2 ) {
            std::cout << "Time for reconstruction out of interval!" << std::endl;
            continue;
        }

        for ( int iDim = 0; iDim < settings.ndim; iDim ++ ) {

            if ( iDim == 0 ) Err_RBM = Err_RBM_rhoU;
            if ( iDim == 1 ) Err_RBM = Err_RBM_rhoV;
            if ( iDim == 2 ) Err_RBM = Err_RBM_rhoW;

            int count = 0;
            double Dt = t_vec[index2] - t_vec[index1];
            for ( int k = 0; k < Nmethods; k ++ ) {
                Err_interp(k) = Err_RBM(index1,k) + (Err_RBM(index2,k) - Err_RBM(index1,k))/
                                    Dt*(settings.t_rec[i] - t_vec[index1]);
            }
            std::cout << std::endl;
            double eps = Err_interp.minCoeff( &best_method_idx );
            // std::cout << "Min coeff = " << eps << " in position " << best_method_idx << std::endl;

            //FIX THIS FUNCTION
            std::string method = method_selected ( best_method_idx, Nf_SPOD, Nf );
//            std::cout << "Best method is " << method << " and Nf ( value meaningful only for SPOD ) : " << Nf_SPOD << std::endl;
            std::cout << "Best method is " << method << std::endl;
            std::cout << " Error : " << Err_interp(best_method_idx) << std::endl;

            std::cout << "Computing Reconstruction using selected method " << std::endl;

            if ( method == "SPOD" )
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

                if ( iDim == 0 && settings.flag_mean == "IC" ) Rec_rhoU = Rec.col(0) + Ic.middleRows(0,Nr);
                if ( iDim == 1 && settings.flag_mean == "IC" ) Rec_rhoV = Rec.col(0) + Ic.middleRows(Nr,Nr);
                if ( iDim == 2 && settings.flag_mean == "IC" ) Rec_rhoW = Rec.col(0) + Ic.middleRows(2*Nr,Nr);
                
            }


            if ( method == "DMD" )
            {

                Eigen::VectorXd lambda_POD;
                Eigen::MatrixXd eig_vec_POD;
                Eigen::VectorXcd lambda_DMD;
                Eigen::MatrixXcd eig_vec_DMD;      
                Eigen::MatrixXcd Phi;
                Eigen::VectorXcd alfa;    

                if ( settings.r == 0 )
                {
                    Phi = DMD_basis( sn_set.middleRows(iDim*Nr,Nr),
                                    lambda_DMD,
                                    eig_vec_DMD,
                                    lambda_POD,
                                    eig_vec_POD,
                                    -1 );
                }
                else
                {
                    Phi = DMD_basis( sn_set.middleRows(iDim*Nr,Nr),
                                    lambda_DMD,
                                    eig_vec_DMD,
                                    lambda_POD,
                                    eig_vec_POD,
                                    settings.r );
                }

                alfa = Calculate_Coefs_DMD_exact ( sn_set.middleRows(iDim*Nr,Nr).leftCols(settings.Ns-1),  
                                                                    lambda_DMD,  
                                                                    Phi );                    

                // std::cout << "Reordering modes DMD ... " << "\t";
                Eigen::VectorXd En = Eigen::VectorXd::Zero(Phi.cols());
                double T = t_vec[t_vec.size()];

                Eigen::VectorXcd omega(Phi.cols());
                for ( int idmd = 0; idmd < Phi.cols(); idmd++ )
                    omega(idmd) = std::log(lambda_DMD(idmd))/(settings.Dt_cfd*(double)settings.Ds);


                for ( int idmd = 0 ; idmd < Phi.cols(); idmd ++ )
                {

                    double alfa_i = alfa(idmd).imag();
                    double alfa_r = alfa(idmd).real();
                    double sigma = omega(idmd).real();
                    En(idmd) = (alfa_r*alfa_r + alfa_i*alfa_i)*(std::exp(2.0*sigma*T) - 1.0)/(2.0*sigma);

                }

                dmd_sort( En, Phi, lambda_DMD, alfa);
                // std::cout << "Done" << std::endl;

                double sum = 0;
                Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
                for (int idmd = 0; idmd < Phi.cols(); idmd++)
                {
                    sum += En(idmd)/En.sum();
                    K_pc(idmd) = sum;
                }

                int Nm;

                if ( settings.r == 0)
                {
                    Nm = Nmod(settings.En, K_pc);
                    std::cout << "Number of modes for the desired energetic content : " << Nm << std::endl;
                    
                }
                else
                {
                    Nm = std::min(settings.r, settings.Ns-1);
                    std::cout << "Number of modes (fixed) : " << Nm << std::endl;
                }
            

                Eigen::MatrixXcd Rec = Reconstruction_DMD ( settings.t_rec[i],
                                                        settings.Dt_cfd*settings.Ds,
                                                        alfa.topRows(Nm),
                                                        Phi.leftCols(Nm),
                                                        lambda_DMD.head(Nm),
                                                        "SCALAR" );

                if ( iDim == 0 && settings.flag_mean == "IC" ) Rec_rhoU = Rec.real().col(0) + Ic.middleRows(0,Nr);
                if ( iDim == 1 && settings.flag_mean == "IC" ) Rec_rhoV = Rec.real().col(0) + Ic.middleRows(Nr,Nr);
                if ( iDim == 2 && settings.flag_mean == "IC" ) Rec_rhoW = Rec.real().col(0) + Ic.middleRows(2*Nr,Nr);
            
            }


            if ( method == "RDMD" )
            {
            
                Eigen::VectorXd lambda = Eigen::VectorXd::Zero(settings.Ns);
                Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
                Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
                Eigen::MatrixXd Phi;

            
                
                std::cout << "Extracting basis and Coeffs RDMD ... " << "\t";        
                //You can define rank DMD at each time step from the config file ( use -1 for the adaptive study adviced)
                Phi = RDMD_modes_coefs ( sn_set.middleRows(iDim*Nr,Nr),
                                        Coefs,
                                        lambda,
                                        K_pc,     
                                        -1, //performing DMD with all non-zero eigenvalues
                                        settings.r_RDMD,
                                        settings.En );
                


                int Nm;
                if ( settings.r == 0 )
                {
                    Nm = Nmod(settings.En, K_pc);
                    std::cout << "number of modes for the desired energetic content " << Nm << std::endl;
                }
                else
                {
                    Nm = std::min(settings.r, settings.r_RDMD);
                    std::cout << "number of modes (fixed) " << Nm << std::endl;
                }


                std::vector<double> t_st_vec(settings.Ns);
                t_st_vec[0] = t_0;

                for ( int irdmd = 1; irdmd < settings.Ns; irdmd++ )
                    t_st_vec[irdmd] = t_st_vec[irdmd-1] + settings.Dt_cfd*(double)settings.Ds;


                Eigen::MatrixXd Rec = Reconstruction_RDMD ( settings.t_rec[i],
                                                            t_st_vec,
                                                            Coefs.topRows(Nm),
                                                            Phi.leftCols(Nm),
                                                            "SCALAR",
                                                            settings.flag_interp );

                if ( iDim == 0 && settings.flag_mean == "IC" ) Rec_rhoU = Rec.col(0) + Ic.middleRows(0,Nr);
                if ( iDim == 1 && settings.flag_mean == "IC" ) Rec_rhoV = Rec.col(0) + Ic.middleRows(Nr,Nr);
                if ( iDim == 2 && settings.flag_mean == "IC" ) Rec_rhoW = Rec.col(0) + Ic.middleRows(2*Nr,Nr);

            } 

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

    std::cout << "Adaptive Reconstruction MODES ends " << std::endl;

    return 0;
}
