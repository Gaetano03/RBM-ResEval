/*
Code for adaptive reconstruction based on residual evaluation
Input config file + error file (+ Modes,Coefs and Encontent RDMD if already available)

Output reconstructed field at the desired time instants with the adaptive technique
based on residual evaluation
*/

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
    if ( settings.flag_prob != "CONSERVATIVE"){
        std::cout << "Reconstruction with residual evaluation only implemented for COnservative Variables flag \n "
                     "FLAG_PROB must be CONSERVATIVE\n Exiting ..." << std::endl;
        exit(EXIT_FAILURE);
    }
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
    int nC = settings.Cols.size();
    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC*Nr);

    if ( settings.flag_mean == "YES" ) {
        for ( int it = 0; it < settings.Ns; it++ )
            sn_set.col(it) -= mean;
    }

    if ( settings.flag_mean == "IC" ) {
        Ic = IC(sn_set,settings,nC,Nr);
    }

    std::cout << "Reading Residuals ... " << std::endl;

    std::vector<std::string> resfilename = {"history_pod.csv", "history_dmd.csv", "history_rdmd.csv"};
//    std::vector<std::string> resfilename = {"history_spod_0.csv", "history_spod_1.csv", "history_spod_2.csv",
//        "history_spod_3.csv", "history_spod_4.csv", "history_dmd.csv", "history_rdmd.csv"};

    Eigen::MatrixXd Err_RBM = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXd Err_RBM_rho = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXd Err_RBM_rhoV = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXd Err_RBM_rhoU = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXd Err_RBM_rhoW = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXd Err_RBM_rhoE = Eigen::MatrixXd::Zero(settings.t_res.size(), Nmethods);

    Eigen::MatrixXi Idx_RBM = Eigen::MatrixXi::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXi Idx_RBM_rho = Eigen::MatrixXi::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXi Idx_RBM_rhoV = Eigen::MatrixXi::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXi Idx_RBM_rhoU = Eigen::MatrixXi::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXi Idx_RBM_rhoW = Eigen::MatrixXi::Zero(settings.t_res.size(), Nmethods);
    Eigen::MatrixXi Idx_RBM_rhoE = Eigen::MatrixXi::Zero(settings.t_res.size(), Nmethods);


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

        while ( getline( file_data, line_flow_data ) && n_row <  settings.t_res.size())
        {
            
            std::istringstream iss(line_flow_data);
            std::string token;
            double err;
            int idx;
            count = 0; 
            while( getline( iss, token, ',') )
            {
                err = std::stod(token);
                idx = std::stoi(token);

                if (settings.ndim == 3) {
                    if (count == 1) Err_RBM_rho(n_row, i) = std::pow(10.0, err);
                    if (count == 2) Err_RBM_rhoU(n_row, i) = std::pow(10.0, err);
                    if (count == 3) Err_RBM_rhoV(n_row, i) = std::pow(10.0, err);
                    if (count == 4) Err_RBM_rhoW(n_row, i) = std::pow(10.0, err);
                    if (count == 5) Err_RBM_rhoE(n_row, i) = std::pow(10.0, err);
                    if (count == 6) Idx_RBM_rho(n_row, i) = idx;
                    if (count == 7) Idx_RBM_rhoU(n_row, i) = idx;
                    if (count == 8) Idx_RBM_rhoV(n_row, i) = idx;
                    if (count == 9) Idx_RBM_rhoW(n_row, i) = idx;
                    if (count == 10) Idx_RBM_rhoE(n_row, i) = idx;

                }

                if (settings.ndim == 2) {
                    if (count == 1) Err_RBM_rho(n_row, i) = std::pow(10.0, err);
                    if (count == 2) Err_RBM_rhoU(n_row, i) = std::pow(10.0, err);
                    if (count == 3) Err_RBM_rhoV(n_row, i) = std::pow(10.0, err);
                    if (count == 4) Err_RBM_rhoE(n_row, i) = std::pow(10.0, err);
                    if (count == 5) Idx_RBM_rho(n_row, i) = idx;
                    if (count == 6) Idx_RBM_rhoU(n_row, i) = idx;
                    if (count == 7) Idx_RBM_rhoV(n_row, i) = idx;
                    if (count == 8) Idx_RBM_rhoE(n_row, i) = idx;
                }

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

    std::vector<double> t_pod = {};
    std::vector<double> t_dmd = {};
    std::vector<double> t_rdmd = {};

    for ( int i = 0; i < settings.t_rec.size(); i++ ) {

        Eigen::VectorXd Rec_rho(Nr);
        Eigen::VectorXd Rec_rhoU(Nr);
        Eigen::VectorXd Rec_rhoV(Nr);
        Eigen::VectorXd Rec_rhoW(Nr);
        Eigen::VectorXd Rec_rhoE(Nr);

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

        for ( int iDim = 0; iDim < settings.ndim + 2; iDim ++ ) {

            if (iDim == 0) {
                Err_RBM = Err_RBM_rho;
                Idx_RBM = Idx_RBM_rho;
            }
            if (iDim == 1) {
                Err_RBM = Err_RBM_rhoU;
                Idx_RBM = Idx_RBM_rhoU;
            }
            if (iDim == 2) {
                Err_RBM = Err_RBM_rhoV;
                Idx_RBM = Idx_RBM_rhoV;
            }
            if (iDim == 3 && settings.ndim == 2) {
                Err_RBM = Err_RBM_rhoE;
                Idx_RBM = Idx_RBM_rhoE;
            }
            if (iDim == 3 && settings.ndim == 3) {
                Err_RBM = Err_RBM_rhoW;
                Idx_RBM = Idx_RBM_rhoW;
            }
            if (iDim == 4 ) {
                Err_RBM = Err_RBM_rhoE;
                Idx_RBM = Idx_RBM_rhoE;
            }


            int count = 0;
            double Dt = t_vec[index2] - t_vec[index1];
            for (int k = 0; k < Nmethods; k++) {
                Err_interp(k) = Err_RBM(index1, k) + (Err_RBM(index2, k) - Err_RBM(index1, k)) /
                                                     Dt * (settings.t_rec[i] - t_vec[index1]);
            }
            std::cout << std::endl;
            double eps = Err_interp.minCoeff(&best_method_idx);
            // std::cout << "Min coeff = " << eps << " in position " << best_method_idx << std::endl;

            //FIX THIS FUNCTION
            std::string method = method_selected(best_method_idx, Nf_SPOD, Nf);
//            std::cout << "Best method is " << method << " and Nf ( value meaningful only for SPOD ) : " << Nf_SPOD << std::endl;
            std::cout << "Best method is " << method << " using a number of modes equal to : "
                      << Idx_RBM(index1, best_method_idx) << std::endl;
            std::cout << " Error : " << Err_interp(best_method_idx) << std::endl;

            std::cout << "Computing Reconstruction using selected methods " << std::endl;

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

                int temp2 = Idx_RBM(index1,best_method_idx);
                int temp1 = Phi.cols();

                if ( temp2 > temp1 ) {
                    std::cout << "Something wrong in residual evaluation\n Exiting ... " << std::endl;
                    exit(EXIT_FAILURE);
                }

                int Nm = std::min(temp1, temp2);

//                if ( settings.r == 0 )
//                {
//                    Nm = Nmod(settings.En, K_pc);
//                    std::cout << "Number of modes for desired energetic content: " << Nm << std::endl;
//                }
//                else
//                {
//                    Nm = std::min(settings.r,settings.Ns);
//                    std::cout << "Number of modes (fixed): " << Nm << std::endl;
//                }

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

                if ( iDim == 0 && settings.flag_mean == "IC" ) Rec_rho = Rec.col(0) + Ic.middleRows(0,Nr);
                if ( iDim == 1 && settings.flag_mean == "IC" ) Rec_rhoU = Rec.col(0) + Ic.middleRows(Nr,Nr);
                if ( iDim == 2 && settings.flag_mean == "IC" ) Rec_rhoV = Rec.col(0) + Ic.middleRows(2*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 2 && settings.flag_mean == "IC" ) Rec_rhoE = Rec.col(0) + Ic.middleRows(3*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 3 && settings.flag_mean == "IC" ) Rec_rhoW = Rec.col(0) + Ic.middleRows(3*Nr,Nr);
                if ( iDim == 4 && settings.flag_mean == "IC" ) Rec_rhoE = Rec.col(0) + Ic.middleRows(4*Nr,Nr);

                if ( iDim == 0 && settings.flag_mean == "YES" ) Rec_rho = Rec.col(0) + mean.middleRows(0,Nr);
                if ( iDim == 1 && settings.flag_mean == "YES" ) Rec_rhoU = Rec.col(0) + mean.middleRows(Nr,Nr);
                if ( iDim == 2 && settings.flag_mean == "YES" ) Rec_rhoV = Rec.col(0) + mean.middleRows(2*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 2 && settings.flag_mean == "YES" ) Rec_rhoE = Rec.col(0) + mean.middleRows(3*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 3 && settings.flag_mean == "YES" ) Rec_rhoW = Rec.col(0) + mean.middleRows(3*Nr,Nr);
                if ( iDim == 4 && settings.flag_mean == "YES" ) Rec_rhoE = Rec.col(0) + mean.middleRows(4*Nr,Nr);
            }


            if ( method == "DMD" )
            {

                Eigen::VectorXd lambda_POD;
                Eigen::MatrixXd eig_vec_POD;
                Eigen::VectorXcd lambda_DMD;
                Eigen::MatrixXcd eig_vec_DMD;      
                Eigen::MatrixXcd Phi;
                Eigen::VectorXcd alfa;

                Phi = DMD_basis( sn_set.middleRows(iDim*Nr,Nr),
                                 lambda_DMD,
                                 eig_vec_DMD,
                                 lambda_POD,
                                 eig_vec_POD,
                                 Idx_RBM(index1,best_method_idx));
//                if ( settings.r == 0 )
//                {
//                    Phi = DMD_basis( sn_set.middleRows(iDim*Nr,Nr),
//                                    lambda_DMD,
//                                    eig_vec_DMD,
//                                    lambda_POD,
//                                    eig_vec_POD,
//                                    -1 );
//                }
//                else
//                {
//                    Phi = DMD_basis( sn_set.middleRows(iDim*Nr,Nr),
//                                    lambda_DMD,
//                                    eig_vec_DMD,
//                                    lambda_POD,
//                                    eig_vec_POD,
//                                    settings.r );
//                }

                alfa = Calculate_Coefs_DMD_exact ( sn_set.middleRows(iDim*Nr,Nr).leftCols(settings.Ns-1),  
                                                                    lambda_DMD,  
                                                                    Phi );                    

                // std::cout << "Reordering modes DMD ... " << "\t";
//                Eigen::VectorXd En = Eigen::VectorXd::Zero(Phi.cols());
//                double T = t_vec[t_vec.size()];

//                Eigen::VectorXcd omega(Phi.cols());
//                for ( int idmd = 0; idmd < Phi.cols(); idmd++ )
//                    omega(idmd) = std::log(lambda_DMD(idmd))/(settings.Dt_cfd*(double)settings.Ds);
//
//
//                for ( int idmd = 0 ; idmd < Phi.cols(); idmd ++ )
//                {
//
//                    double alfa_i = alfa(idmd).imag();
//                    double alfa_r = alfa(idmd).real();
//                    double sigma = omega(idmd).real();
//                    En(idmd) = (alfa_r*alfa_r + alfa_i*alfa_i)*(std::exp(2.0*sigma*T) - 1.0)/(2.0*sigma);
//
//                }
//
//                dmd_sort( En, Phi, lambda_DMD, alfa);
//                // std::cout << "Done" << std::endl;
//
//                double sum = 0;
//                Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
//                for (int idmd = 0; idmd < Phi.cols(); idmd++)
//                {
//                    sum += En(idmd)/En.sum();
//                    K_pc(idmd) = sum;
//                }
//
//                int Nm;
//
//                if ( settings.r == 0)
//                {
//                    Nm = Nmod(settings.En, K_pc);
//                    std::cout << "Number of modes for the desired energetic content : " << Nm << std::endl;
//
//                }
//                else
//                {
//                    Nm = std::min(settings.r, settings.Ns-1);
//                    std::cout << "Number of modes (fixed) : " << Nm << std::endl;
//                }
            

//                Eigen::MatrixXcd Rec = Reconstruction_DMD ( settings.t_rec[i],
//                                                        settings.Dt_cfd*settings.Ds,
//                                                        alfa.topRows(Nm),
//                                                        Phi.leftCols(Nm),
//                                                        lambda_DMD.head(Nm),
//                                                        "SCALAR" );

                Eigen::MatrixXcd Rec = Reconstruction_DMD ( settings.t_rec[i],
                                                            settings.Dt_cfd*settings.Ds,
                                                            alfa,
                                                            Phi,
                                                            lambda_DMD,
                                                            "SCALAR" );

                if ( iDim == 0 && settings.flag_mean == "IC" ) Rec_rho = Rec.real().col(0) + Ic.middleRows(0,Nr);
                if ( iDim == 1 && settings.flag_mean == "IC" ) Rec_rhoU = Rec.real().col(0) + Ic.middleRows(Nr,Nr);
                if ( iDim == 2 && settings.flag_mean == "IC" ) Rec_rhoV = Rec.real().col(0) + Ic.middleRows(2*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 2 && settings.flag_mean == "IC" ) Rec_rhoE = Rec.real().col(0) + Ic.middleRows(3*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 3 && settings.flag_mean == "IC" ) Rec_rhoW = Rec.real().col(0) + Ic.middleRows(3*Nr,Nr);
                if ( iDim == 4 && settings.flag_mean == "IC" ) Rec_rhoE = Rec.real().col(0) + Ic.middleRows(4*Nr,Nr);

                if ( iDim == 0 && settings.flag_mean == "YES" ) Rec_rho = Rec.real().col(0) + mean.middleRows(0,Nr);
                if ( iDim == 1 && settings.flag_mean == "YES" ) Rec_rhoU = Rec.real().col(0) + mean.middleRows(Nr,Nr);
                if ( iDim == 2 && settings.flag_mean == "YES" ) Rec_rhoV = Rec.real().col(0) + mean.middleRows(2*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 2 && settings.flag_mean == "YES" ) Rec_rhoE = Rec.real().col(0) + mean.middleRows(3*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 3 && settings.flag_mean == "YES" ) Rec_rhoW = Rec.real().col(0) + mean.middleRows(3*Nr,Nr);
                if ( iDim == 4 && settings.flag_mean == "YES" ) Rec_rhoE = Rec.real().col(0) + mean.middleRows(4*Nr,Nr);

            }


            if ( method == "RDMD" )
            {
            
                Eigen::VectorXd lambda = Eigen::VectorXd::Zero(settings.Ns);
                Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
                Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
                Eigen::MatrixXd Phi;

            
                
                std::cout << "Extracting basis and Coeffs RDMD ... " << "\t";        
                //You can define rank DMD at each time step from the config file ( use -1 for the adaptive study adviced)

                int Nm = Idx_RBM(index1,best_method_idx);
                Phi = RDMD_modes_coefs ( sn_set.middleRows(iDim*Nr,Nr),
                                        Coefs,
                                        lambda,
                                        K_pc,     
                                        -1, //performing DMD with all non-zero eigenvalues
                                        Nm,
                                        settings.En );
                

//                if ( settings.r == 0 )
//                {
//                    Nm = Nmod(settings.En, K_pc);
//                    std::cout << "number of modes for the desired energetic content " << Nm << std::endl;
//                }
//                else
//                {
//                    Nm = std::min(settings.r, settings.r_RDMD);
//                    std::cout << "number of modes (fixed) " << Nm << std::endl;
//                }


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

                if ( iDim == 0 && settings.flag_mean == "IC" ) Rec_rho = Rec.col(0) + Ic.middleRows(0,Nr);
                if ( iDim == 1 && settings.flag_mean == "IC" ) Rec_rhoU = Rec.col(0) + Ic.middleRows(Nr,Nr);
                if ( iDim == 2 && settings.flag_mean == "IC" ) Rec_rhoV = Rec.col(0) + Ic.middleRows(2*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 2 && settings.flag_mean == "IC" ) Rec_rhoE = Rec.col(0) + Ic.middleRows(3*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 3 && settings.flag_mean == "IC" ) Rec_rhoW = Rec.col(0) + Ic.middleRows(3*Nr,Nr);
                if ( iDim == 4 && settings.flag_mean == "IC" ) Rec_rhoE = Rec.col(0) + Ic.middleRows(4*Nr,Nr);

                if ( iDim == 0 && settings.flag_mean == "YES" ) Rec_rho = Rec.col(0) + mean.middleRows(0,Nr);
                if ( iDim == 1 && settings.flag_mean == "YES" ) Rec_rhoU = Rec.col(0) + mean.middleRows(Nr,Nr);
                if ( iDim == 2 && settings.flag_mean == "YES" ) Rec_rhoV = Rec.col(0) + mean.middleRows(2*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 2 && settings.flag_mean == "YES" ) Rec_rhoE = Rec.col(0) + mean.middleRows(3*Nr,Nr);
                if ( iDim == 3 && settings.ndim == 3 && settings.flag_mean == "YES" ) Rec_rhoW = Rec.col(0) + mean.middleRows(3*Nr,Nr);
                if ( iDim == 4 && settings.flag_mean == "YES" ) Rec_rhoE = Rec.col(0) + mean.middleRows(4*Nr,Nr);

            } 

        }
        
        Eigen::MatrixXd Rec_M(Nr, settings.ndim + 2);
        
        if ( settings.ndim == 2)
        {
            Rec_M.col(0) = Rec_rho;
            Rec_M.col(1) = Rec_rhoU;
            Rec_M.col(2) = Rec_rhoV;
            Rec_M.col(3) = Rec_rhoE;
        } else
        {
            Rec_M.col(0) = Rec_rho;
            Rec_M.col(1) = Rec_rhoU;
            Rec_M.col(2) = Rec_rhoV;
            Rec_M.col(3) = Rec_rhoW;
            Rec_M.col(4) = Rec_rhoE;
        }
        
        std::cout << "Writing reconstructed field ..." << "\t";
        write_Reconstructed_fields ( Rec_M, Coords, settings.out_file, "CONSERVATIVE", i );
        std::cout << "Done" << std::endl << std::endl << std::endl;

    }

    std::cout << "Adaptive Reconstruction MODES ends " << std::endl;

    return 0;
}
