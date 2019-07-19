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

    if ( settings.flag_method == "SPOD")
    {

        std::vector<double> t_vec( settings.Ns );
        t_vec[0] = 0.0;

        for ( int i = 1; i < settings.Ns; i++ )
            t_vec[i] = t_vec[i-1] + settings.Dt_cfd*settings.Ds;

        std::cout << std::endl;
        std::cout << "Initialized vector of times " << std::endl;

        Eigen::VectorXd lambda(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);

        if ( settings.flag_mean == "YES" )
        {
            std::cout << "Subtracting mean from snapshots ... " << std::endl << std::endl;
            for ( int i = 0; i < settings.Ns; i++ )
                sn_set.col(i) -= mean;
        }


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
            std::cout << "Number of modes for the desired energy content (needs fixes for sPOD) : " << Nrec << std::endl;
        }
        else
        {
            Nrec = settings.r;
            std::cout << " Number of modes (fixed, needs fixes for sPOD) : " << Nrec << std::endl;
        }


        // if ( settings.flag_wdb_be == "YES" )
        // {
            
        //     std::cout << "Writing modes ..." << "\t";
        //     write_modes_sPOD ( Phi.leftCols(Nrec), Coords, settings.flag_prob );
        //     std::cout << "Complete!" << std::endl;

        //     std::cout << "Writing Coefficients ..." << "\t";
        //     write_coeffs_sPOD ( eig_vec.transpose(), t_vec, lambda );
        //     std::cout << "Complete!" << std::endl;
        //     std::cout << std::endl;

        // }

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


    } else {

        std::cout << "Only SPOD implemented so far for DNS" << std::endl;

    }


    std::cout << std::endl;    
    std::cout << "-----------RBM-Clyde end-------------" << std::endl << std::endl;

    return 0;

}
