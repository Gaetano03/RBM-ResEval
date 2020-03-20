#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"

int main(int argc, char *argv[]) {

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "-----------FLUID-MODES DNS starts-------------" << std::endl << std::endl;

    std::string filecfg = argv[1];
    // std::string mode = argv[2];
    prob_settings settings;

    //Reading configuration file
    Read_cfg( filecfg, settings );
    std::string filename = std::to_string(settings.nstart) + ".q";
    plot3d_info Info = read_plot3d_info (filename);
    std::vector<double> time(settings.Ns);

    std::string file_temp;
    int count = 0;

    std::cout << "Storing Vector of time ... " << "\t";

    for ( int it = settings.nstart; it < (settings.Ns*settings.Ds + settings.nstart); it += settings.Ds )
    {
        file_temp = std::to_string(it) + ".q";
        plot3d_info Info_time = read_plot3d_info (file_temp);
        time[count] = (double)Info_time.T/(settings.Mach*std::sqrt(1.4*settings.P/settings.Rho)); 
        count++;
    }

    double Dt_avg = 0.0;
    double Dt;

    for ( int i = 0; i < time.size()-1; i++ ) Dt_avg += time[i+1] - time[i]; 
    
    Dt_avg = Dt_avg/((double)time.size()-1.0);

    // for ( int i = 0; i < 5; i++ )
    //     std::cout << time[i] << std::endl;

    std::cout << "Complete!" << std::endl;

    int Np = 0;
    //computing number of points of the whole grid
    for ( int iblock = 0; iblock < Info.nblocks; iblock++ )
        Np += Info.ni[iblock]*Info.nj[iblock]*Info.nk[iblock];
    
    Eigen::VectorXd K_pc(settings.Ns);
    
    Config_stream ( settings );

    int Nr = 0;
    // Create matrix of snapshots
    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings.Ns, settings.Ds, settings.nstart,
                                        settings.Cols,
                                        settings.in_file,
                                        settings.flag_prob,
                                        settings.solver);

    //Eigen::VectorXd mean = sn_set.rowwise().mean();

    if ( settings.flag_mean == "YES" )
    {
        Eigen::VectorXd mean = sn_set.rowwise().mean();
        std::cout << "Subtracting mean from snapshots ... " << std::endl << std::endl;
        for ( int i = 0; i < settings.Ns; i++ )
            sn_set.col(i) -= mean;
    }

    if ( settings.flag_mean == "IC" && settings.flag_prob == "VELOCITY-3D" )
    {
        double Ufree = settings.Mach*std::sqrt(1.4*settings.P/settings.Rho);
        Eigen::VectorXd Ic = Eigen::VectorXd::Zero(sn_set.rows());
        Ic.head(sn_set.rows()/3) = Eigen::VectorXd::Ones(sn_set.rows()/3)*Ufree;
        std::cout << "Subtracting IC from snapshots ... " << std::endl << std::endl;
        for ( int i = 0; i < settings.Ns; i++ )
            sn_set.col(i) -= Ic;
    }
    else if ( settings.flag_mean == "IC" && settings.flag_prob != "VELOCITY-3D" )
    {
        std::cout << "Initial condition only implemented for Velocity 3D so far" << std::endl;
        exit( EXIT_FAILURE );
    }


    if ( settings.flag_method[0] == "SPOD")
    {

        Eigen::VectorXd lambda(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);

        std::cout << "Extracting basis SPOD using Nf " << settings.Nf << " ... " << "\t";        

        Eigen::MatrixXd Phi = SPOD_basis( sn_set,
                                lambda, K_pc, eig_vec,
                                settings.Nf,
                                settings.flag_bc, 
                                settings.flag_filter,  
                                settings.sigma);

        std::cout << " Done! " << std::endl << std::endl;

        std::cout << "Number of non-zero modes : " << Phi.cols() << std::endl;
        std::cout << "Eigenvalues : " << lambda << std::endl;
        std::cout << "K_pc : " << K_pc << std::endl;

        int Nrec; 
        if ( settings.r == 0 )
        {    
            Nrec = Nmod( settings.En, K_pc);
            std::cout << "Number of modes for the desired energy content : " << Nrec << std::endl;
        }
        else
        {
            int du = Phi.cols();
            Nrec = std::min(settings.r, du);
            std::cout << " Number of modes : " << Nrec << std::endl;
        }

        if ( settings.flag_wdb_be == "YES" )
        {
            
            std::cout << "Writing modes and Coeffs..." << std::endl;

            if ( settings.flag_prob == "VELOCITY-3D" )
            {
                // Write_Plot3d_Modes( Phi.topRows(Np).leftCols(Nrec), "Modes_U.f", Info );
                // Write_Plot3d_Modes( Phi.middleRows(Np,Np).leftCols(Nrec), "Modes_V.f", Info );
                // Write_Plot3d_Modes( Phi.bottomRows(Np).leftCols(Nrec), "Modes_W.f", Info );
                Eigen::MatrixXd PHI(Np,3*Nrec);
                PHI << Phi.topRows(Np).leftCols(Nrec), 
                        Phi.middleRows(Np,Np).leftCols(Nrec), 
                        Phi.bottomRows(Np).leftCols(Nrec);
                Write_Plot3d_Modes( PHI, "Modes_POD.f", Info );

            }
            
            std::cout << "Writing Coefficients ..." << "\t";
            write_coeffs_sPOD ( eig_vec, time, lambda );
            std::cout << "Complete!" << std::endl;
            std::cout << std::endl;

        }

    } else if ( settings.flag_method[0] == "RDMD")
    {

        Eigen::VectorXd lambda = Eigen::VectorXd::Zero(settings.Ns);
        Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
        
        std::cout << "Extracting basis and Coeffs RDMD ... " << "\t";        
        Eigen::MatrixXd Phi = RDMD_modes_coefs ( sn_set,
                                                Coefs,
                                                lambda,
                                                K_pc,     
                                                settings.r,
                                                settings.r_RDMD,
                                                settings.En );

        int Nrec;

        if ( settings.r_RDMD == 0 )
        {    
            Nrec = Nmod( settings.En, K_pc);
            std::cout << "Number of modes for the desired energy content : " << Nrec << std::endl;
        }
        else
        {
            int du = Phi.cols();
            Nrec = std::min(settings.r_RDMD, du);
            std::cout << " Number of modes : " << Nrec << std::endl;
        }

        std::cout << " Done! " << std::endl << std::endl;

        if ( settings.flag_wdb_be == "YES" )
        {
            std::cout << "Writing modes and coeffs..." << "\t"; 

            if ( settings.flag_prob == "VELOCITY-3D" )
            {
                // Write_Plot3d_Modes( Phi.topRows(Np).leftCols(Nrec), "Modes_U.f", Info );
                // Write_Plot3d_Modes( Phi.middleRows(Np,Np).leftCols(Nrec), "Modes_V.f", Info );
                // Write_Plot3d_Modes( Phi.bottomRows(Np).leftCols(Nrec), "Modes_W.f", Info );
                Eigen::MatrixXd PHI(Np,3*Nrec);
                PHI << Phi.topRows(Np).leftCols(Nrec), 
                        Phi.middleRows(Np,Np).leftCols(Nrec), 
                        Phi.bottomRows(Np).leftCols(Nrec);
                Write_Plot3d_Modes( PHI, "Modes_RDMD.f", Info );

            }

            std::cout << "Writing Coefficients ..." << "\t";
            write_coeffs_sPOD ( Coefs.transpose(), time, lambda );
            std::cout << "Complete!" << std::endl;
            std::cout << std::endl;
        }

    } else if ( settings.flag_method[0] == "DMD")
    {
        double tol = 1e-8;
        double t_0 = time[0];

        Eigen::VectorXd t_vec(time.size()); //Building a vector of uniform times for DMD
        t_vec(0) = t_0;

        for ( int it = 1; it < time.size(); it++ ) t_vec(it) = t_vec(it-1) + Dt_avg; 

        Eigen::VectorXd lambda_POD;
        Eigen::MatrixXd eig_vec_POD;
        Eigen::VectorXcd lambda_DMD;
        Eigen::MatrixXcd eig_vec_DMD;

        std::cout << "Extracting basis ... " << "\t";        
        Eigen::MatrixXcd Phi;
        Eigen::VectorXcd alfa;
        Eigen::MatrixXcd Alfas;
           
        Phi = DMD_basis( sn_set,
                        lambda_DMD,
                        eig_vec_DMD,
                        lambda_POD,
                        eig_vec_POD,
                        settings.r );

        std::cout << " Done! " << std::endl << std::endl;

        int Nm = Phi.cols();
        std::cout << "Number of modes extracted : " << Nm << std::endl;

        Eigen::VectorXcd omega(Nm);
        for ( int i = 0; i < Nm; i++ )
                omega(i) = std::log(lambda_DMD(i))/(Dt_avg);
        
        std::cout << "Calculating coefficients DMD ... " << "\t";

        //Calculating coefficients solving optimization problem
        if ( settings.dmd_coef_flag == "OPT" )
        {
            alfa = Calculate_Coefs_DMD_exact ( sn_set.leftCols(settings.Ns-1),  //matrix of first Ns-1 snaps 
                                                                lambda_DMD,  //slow eigenvalues
                                                                Phi ); //slow exact DMD modes
        }
        //Calculating coefficients with ls
        else if ( settings.dmd_coef_flag == "LS" )
        {
            Eigen::VectorXcd b = Eigen::VectorXcd::Zero(sn_set.rows()); 
            for ( int k = 0; k < sn_set.rows(); k++ )
                b(k).real(sn_set(k,0)); 

            alfa = Phi.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

        }
        //Calculating coefficients with Hybrid method
        else if ( settings.dmd_coef_flag == "HYBRID" )
        {
        
            Alfas = Calculate_Coefs_Matrix_DMD ( sn_set,
                                                Phi,
                                                omega,
                                                t_0,
                                                Dt_avg );


            std::cout << "Writing training points ..." << std::endl;

            std::ofstream train_real;
            train_real.open("train_real.dat");


            for ( int k = 0; k < settings.Ns; k++ )
            {
            
                for( int j = 0; j < Alfas.cols(); j++ ) 
                    train_real << Alfas(k,j).real() << " ";   

            train_real << std::endl;

            }

            train_real.close();


            std::ofstream train_imag;
            train_imag.open("train_imag.dat");


            for ( int k = 0; k < settings.Ns; k++ )
            {
            
                for( int j = 0; j < Alfas.cols(); j++ ) 
                    train_imag << Alfas(k,j).imag() << " ";   

            train_imag << std::endl;

            }

            train_imag.close();

        }
        else
        {
            std::cout << "Method to Calculate DMD coefficients not available! " << std::endl;
            std::cout << "Exiting ... " << std::endl;
            std::exit( EXIT_FAILURE );
        }

        std::cout << " Done! " << std::endl << std::endl;
        

        if ( settings.flag_wdb_be == "YES" )
        {
            
            std::cout << "Writing modes ..." << "\t";

            if ( settings.flag_prob == "VELOCITY-3D" )
            {
                // Write_Plot3d_Modes( Phi.topRows(Np).leftCols(Nm).real(), "Modes_U_r.f", Info );
                // Write_Plot3d_Modes( Phi.middleRows(Np,Np).leftCols(Nm).real(), "Modes_V_r.f", Info );
                // Write_Plot3d_Modes( Phi.bottomRows(Np).leftCols(Nm).real(), "Modes_W_r.f", Info );
                // Write_Plot3d_Modes( Phi.topRows(Np).leftCols(Nm).imag(), "Modes_U_i.f", Info );
                // Write_Plot3d_Modes( Phi.middleRows(Np,Np).leftCols(Nm).imag(), "Modes_V_i.f", Info );
                // Write_Plot3d_Modes( Phi.bottomRows(Np).leftCols(Nm).imag(), "Modes_W_i.f", Info );
                
                Eigen::MatrixXd PHI(Np,3*Nm);
                PHI << Phi.topRows(Np).leftCols(Nm).real(), 
                        Phi.middleRows(Np,Np).leftCols(Nm).real(), 
                        Phi.bottomRows(Np).leftCols(Nm).real();
                Write_Plot3d_Modes( PHI, "Modes_DMD_r.f", Info );
                
                PHI << Phi.topRows(Np).leftCols(Nm).imag(), 
                        Phi.middleRows(Np,Np).leftCols(Nm).imag(), 
                        Phi.bottomRows(Np).leftCols(Nm).imag();
                Write_Plot3d_Modes( PHI, "Modes_DMD_i.f", Info );

            }


            if ( settings.dmd_coef_flag!= "HYBRID" )
            {
                std::cout << "Writing dmd data and time dynamics ..." << "\t";
                write_TimeDynamics_DMD ( omega, alfa, t_vec );
                write_alfa_lam_DMD( alfa, lambda_DMD);
                std::cout << "Complete!" << std::endl;
                std::cout << std::endl;

            }

        }

    } else {

        std::cout << "Only SPOD, DMD and RDMD implemented for CS3D" << std::endl;

    }


    std::cout << std::endl;    
    std::cout << "-----------FLUID MODAL ANALYSIS ends-------------" << std::endl << std::endl;

    return 0;

}
