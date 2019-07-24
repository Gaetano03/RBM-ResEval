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

    Eigen::VectorXd mean = sn_set.rowwise().mean();

    if ( settings.flag_mean == "YES" )
    {
        std::cout << "Subtracting mean from snapshots ... " << std::endl << std::endl;
        for ( int i = 0; i < settings.Ns; i++ )
            sn_set.col(i) -= mean;
    }

    if ( settings.flag_method == "SPOD")
    {

        // std::vector<double> t_vec( settings.Ns );
        // t_vec[0] = 0.0;

        // for ( int i = 1; i < settings.Ns; i++ )
        //     t_vec[i] = t_vec[i-1] + settings.Dt_cfd*settings.Ds;

        // std::cout << std::endl;
        // std::cout << "Initialized vector of times " << std::endl;

        Eigen::VectorXd lambda(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);

        std::cout << "Extracting basis ... " << "\t";        

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
            Nrec = settings.r;
            std::cout << " Number of modes : " << Nrec << std::endl;
        }

        if ( settings.flag_wdb_be == "YES" )
        {
            
            std::cout << "Writing modes and Coeffs..." << std::endl;

            if ( settings.flag_prob == "VELOCITY-3D" )
            {
                Write_Plot3d_Modes( Phi.topRows(Np).leftCols(Nrec), "Modes_U.f", Info );
                Write_Plot3d_Modes( Phi.middleRows(Np,Np).leftCols(Nrec), "Modes_V.f", Info );
                Write_Plot3d_Modes( Phi.bottomRows(Np).leftCols(Nrec), "Modes_W.f", Info );
            }
            
            std::cout << "Writing Coefficients ..." << "\t";
            write_coeffs_sPOD ( eig_vec, time, lambda );
            std::cout << "Complete!" << std::endl;
            std::cout << std::endl;

        }
        

        // if ( settings.flag_rec == "YES" )
        // {

        //     for ( int nt = 0; nt < settings.t_rec.size(); nt++ )
        //     {

        //         std::cout << "Reconstructing field at time : " << settings.t_rec[nt] << "\t";

        //         Eigen::MatrixXd Rec = Reconstruction_S_POD ( t_vec,
        //                             K_pc, lambda, eig_vec.transpose(),
        //                             Phi, settings.t_rec[nt],
        //                             Nrec,
        //                             settings.flag_prob,
        //                             settings.flag_interp ) ;

        //         std::cout << "Done" << std::endl;

        //         if ( settings.flag_mean == "YES" )
        //         {

        //             for ( int i = 0; i < Rec.cols(); i++)
        //                 Rec.col(i) = Rec.col(i) + mean.segment(i*Nr, Nr);

        //         }

        //         std::cout << "Writing reconstructed field ..." << "\t";

        //         write_Reconstructed_fields ( Rec, Coords,
        //                                 settings.out_file,
        //                                 settings.flag_prob, nt );

        //         std::cout << "Done" << std::endl << std::endl;
        //     }

        // }

    } else if ( settings.flag_method == "RDMD")
    {

        Eigen::VectorXd lambda = Eigen::VectorXd::Zero(settings.Ns);
        Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
        Eigen::MatrixXd Phi;

        
        std::cout << "Extracting basis and Coeffs RDMD ... " << "\t";        
    
        Phi = RDMD_modes_coefs ( sn_set,
                                Coefs,
                                lambda,
                                K_pc,     
                                settings.r,
                                settings.r_RDMD,
                                settings.En );
                                
        

        // std::cout << "Check mode orthogonality\n PhiT*Phi :\n " << Phi.transpose()*Phi << std::endl; 

        int Nrec = Phi.cols();

        std::cout << " Done! " << std::endl << std::endl;

        if ( settings.flag_wdb_be == "YES" )
        {
            std::cout << "Writing modes and coeffs..." << "\t"; 

            if ( settings.flag_prob == "VELOCITY-3D" )
            {
                Write_Plot3d_Modes( Phi.topRows(Np).leftCols(Nrec), "Modes_U.f", Info );
                Write_Plot3d_Modes( Phi.middleRows(Np,Np).leftCols(Nrec), "Modes_V.f", Info );
                Write_Plot3d_Modes( Phi.bottomRows(Np).leftCols(Nrec), "Modes_W.f", Info );
            }

            std::cout << "Writing Coefficients ..." << "\t";
            write_coeffs_sPOD ( Coefs.transpose(), time, lambda );
            std::cout << "Complete!" << std::endl;
            std::cout << std::endl;
        }

        // if ( settings.flag_rec == "YES" )
        // {
                               
        //     for ( int nt = 0; nt < settings.t_rec.size(); nt ++)
        //     {
        //         Eigen::MatrixXd Rec;
        //         std::cout << "Reconstructing field at time : " << settings.t_rec[nt] << "\t";


        //         Rec = Reconstruction_RDMD ( settings.t_rec[nt],
        //                                 t_st_vec,
        //                                 Coefs,
        //                                 Phi,
        //                                 settings.flag_prob,
        //                                 settings.flag_interp );


        //         std::cout << "Done" << std::endl;

        //         if ( settings.flag_mean == "YES" )
        //         {

        //             for ( int i = 0; i < Rec.cols(); i++)
        //                 Rec.col(i) = Rec.col(i) + mean.segment(i*Nr, Nr);

        //         }

        //         std::cout << "Writing reconstructed field ..." << "\t";

        //         write_Reconstructed_fields ( Rec, Coords,
        //                                 settings.out_file,
        //                                 settings.flag_prob, nt );

        //         std::cout << "Done" << std::endl << std::endl;

        //     }

    } else {

        std::cout << "Only SPOD and RDMD implemented for CS3D" << std::endl;

    }


    std::cout << std::endl;    
    std::cout << "-----------FLUID MODAL ANALYSIS ends-------------" << std::endl << std::endl;

    return 0;

}
