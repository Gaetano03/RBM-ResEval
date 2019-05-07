#include "write_Outputs.hpp"


void Config_stream ( prob_settings settings )
{

    std::cout << " Number of snapshots : " << settings.Ns << std::endl;
    std::cout << " Delta T between snapshots ( only EQUISPACED for now) : " << settings.Dt_cfd*settings.Ds << std::endl;
    std::cout << " Starting snapshot index : " << settings.nstart << std::endl;
    std::cout << " Data-type to be processed : " << settings.flag_prob << std::endl << std::endl;

    std::cout << "----------- Performing " << settings.flag_method << " --------------" << std::endl;
    std::cout << " Subtracting mean from snapshot = " << settings.flag_mean << std::endl << std::endl;

    if ( settings.flag_method == "SPOD" )
    {
        std::cout << " Filter size : " << settings.Nf << std::endl;
        std::cout << " Filter type : " << settings.flag_filter << std::endl;
        std::cout << " Energy level desired : " << settings.En*100 << "%" << std::endl;
    }   

    if ( settings.flag_method == "DMD" || settings.flag_method == "fbDMD" || settings.flag_method == "HODMD" || settings.flag_method == "mrDMD" || settings.flag_method == "RDMD")
    {
        std::cout << " Rank for the reduced dynamic (if -1 pick all modes, if 0 do SVHT) : " << settings.r << std::endl;
        std::cout << " Method to compute coefficients : " << settings.dmd_coef_flag << std::endl;
    }

    if ( settings.flag_method == "RDMD" )
    {
        std::cout << " Rank for the RDMD (if 0 do the recursion until the desired energetic content ( max = 3*Ns)) : " << settings.r_RDMD << std::endl;
        std::cout << " Energy level desired (to be considered only if rank = 0): " << settings.En*100 << "%" << std::endl;
    }

    if ( settings.flag_method == "mrDMD" )
    {
        std::cout << " Max levels of the multi resolution : " << settings.max_levels << std::endl;
        std::cout << " Number of samples per time bin : " << settings.max_cycles << std::endl;
    }

    if ( settings.flag_method == "HODMD" )
    {
        std::cout << " Number of levels for high order : " << settings.d << std::endl;
    }

    std::cout << std::endl;

}



void write_modes_sPOD ( const Eigen::MatrixXd &Phi_cut, 
                        const Eigen::MatrixXd &Coords, 
                        std::string flag_prob )
{
    
    int Nr = Coords.rows(); 
    std::string filename = "Modes_sPOD.dat";
    std::ofstream flow_data;
    flow_data.open(filename.c_str());

    if ( flag_prob == "SCALAR" )
    {

        // Write row of Headers
        flow_data << "\"PointID\"" << " ";
        flow_data << "\"x\"" << " ";
        flow_data << "\"y\"" << " ";

        if ( Coords.cols() == 3 )
            flow_data << "\"z\"" << " ";

        std::string phi;

        for ( int i = 0; i < Phi_cut.cols(); i++ )
        {

            phi = "\"Phi_" + std::to_string(i+1) + "\""; 
            flow_data << phi << " ";

        }

        flow_data << std::endl;

        //Write fields
        for ( int i = 0; i < Nr; i++ )
        {

            flow_data << i+1 << " ";
            flow_data << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
            flow_data << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
            if ( Coords.cols() == 3 )
                flow_data << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";
                

            for (int j = 0; j < Phi_cut.cols(); j++)
                flow_data << std::setprecision(12) << std::scientific << Phi_cut(i,j) <<  " ";           

        flow_data << std::endl;

        }
    // Close file
    flow_data.close();

    } else if ( flag_prob == "VECTOR-2D" || flag_prob == "VELOCITY-2D" )
    {

        // Write row of Headers
        flow_data << "\"PointID\"" << " ";
        flow_data << "\"x\"" << " ";
        flow_data << "\"y\"" << " ";

        std::string phix;
        std::string phiy;

        for ( int i = 0; i < Phi_cut.cols(); i++ )
        {

            phix = "\"Phi_x_" + std::to_string(i+1) + "\""; 
            flow_data << phix << " ";
            phiy = "\"Phi_y_" + std::to_string(i+1) + "\""; 
            flow_data << phiy << " ";

        }

        flow_data << std::endl;

        //Write fields
        for ( int i = 0; i < Nr; i++ )
        {

            flow_data << i+1 << " ";
            flow_data << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
            flow_data << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";

            for (int j = 0; j < Phi_cut.cols(); j++)
            {

                flow_data << std::setprecision(12) << std::scientific << Phi_cut(i,j) <<  " ";
                flow_data << std::setprecision(12) << std::scientific << Phi_cut(Nr+i,j) <<  " ";            

            }

        flow_data << std::endl;

        }
    // Close file
    flow_data.close();

    } else if ( flag_prob == "VECTOR-3D" || flag_prob == "VELOCITY-3D" )
    {

        // Write row of Headers
        flow_data << "\"PointID\"" << " ";
        flow_data << "\"x\"" << " ";
        flow_data << "\"y\"" << " ";
        flow_data << "\"z\"" << " ";

        std::string phix;
        std::string phiy;
        std::string phiz;

        for ( int i = 0; i < Phi_cut.cols(); i++ )
        {

            phix = "\"Phi_x_" + std::to_string(i+1) + "\""; 
            flow_data << phix << " ";
            phiy = "\"Phi_y_" + std::to_string(i+1) + "\""; 
            flow_data << phiy << " ";
            phiy = "\"Phi_z_" + std::to_string(i+1) + "\""; 
            flow_data << phiz << " ";

        }

        flow_data << std::endl;

        //Write fields
        for ( int i = 0; i < Nr; i++ )
        {

            flow_data << i+1 << " ";
            flow_data << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
            flow_data << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
            flow_data << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";

            for (int j = 0; j < Phi_cut.cols(); j++)
            {

                flow_data << std::setprecision(12) << std::scientific << Phi_cut(i,j) <<  " ";
                flow_data << std::setprecision(12) << std::scientific << Phi_cut(Nr+i,j) <<  " ";
                flow_data << std::setprecision(12) << std::scientific << Phi_cut(2*Nr+i,j) <<  " ";            

            }

        flow_data << std::endl;

        }
    // Close file
    flow_data.close();


    } else 
    {

        std::cout << "Set well problem flag! Exiting ..." << std::endl;
        exit (EXIT_FAILURE);

    }


}


void write_modes_DMD ( const Eigen::MatrixXcd &Phi_cut,
                    const Eigen::MatrixXd &Coords, 
                    std::string flag_prob )
{

        int Nr = Coords.rows(); 
        std::string filenameI = "Modes_DMD_Imag.dat";
        std::string filenameR = "Modes_DMD_Real.dat";
        std::ofstream flow_dataI, flow_dataR;

        flow_dataI.open(filenameI.c_str());
        flow_dataR.open(filenameR.c_str());

        if ( flag_prob == "SCALAR" )
        {

            // Write row of Headers
            flow_dataI << "\"PointID\"" << " ";
            flow_dataI << "\"x\"" << " ";
            flow_dataI << "\"y\"" << " ";

            flow_dataR << "\"PointID\"" << " ";
            flow_dataR << "\"x\"" << " ";
            flow_dataR << "\"y\"" << " ";

            if ( Coords.cols() == 3 )
            {
                flow_dataI << "\"z\"" << " ";
                flow_dataR << "\"z\"" << " ";
            }

            std::string phi;

            for ( int i = 0; i < Phi_cut.cols(); i++ )
            {

                phi = "\"PhiI_" + std::to_string(i+1) + "\""; 
                flow_dataI << phi << " ";
                phi = "\"PhiR_" + std::to_string(i+1) + "\""; 
                flow_dataR << phi << " ";
            }

            flow_dataI << std::endl;
            flow_dataR << std::endl;

            //Write fields
            for ( int i = 0; i < Nr; i++ )
            {

                flow_dataI << i+1 << " ";
                flow_dataI << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                flow_dataI << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
                
                flow_dataR << i+1 << " ";
                flow_dataR << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                flow_dataR << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";

                if ( Coords.cols() == 3 )
                {
                    flow_dataI << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";
                    flow_dataR << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";
                }    

                for (int j = 0; j < Phi_cut.cols(); j++)
                {
                    flow_dataI << std::setprecision(12) << std::scientific << Phi_cut(i,j).imag() <<  " ";
                    flow_dataR << std::setprecision(12) << std::scientific << Phi_cut(i,j).real() <<  " ";           
                }

            flow_dataI << std::endl;
            flow_dataR << std::endl;

            }
        // Close file
        flow_dataI.close();
        flow_dataR.close();

        } else if ( flag_prob == "VECTOR-2D" || flag_prob == "VELOCITY-2D" )
        {

            // Write row of Headers
            flow_dataI << "\"PointID\"" << " ";
            flow_dataI << "\"x\"" << " ";
            flow_dataI << "\"y\"" << " ";

            flow_dataR << "\"PointID\"" << " ";
            flow_dataR << "\"x\"" << " ";
            flow_dataR << "\"y\"" << " ";

            std::string phix;
            std::string phiy;

            for ( int i = 0; i < Phi_cut.cols(); i++ )
            {

                phix = "\"PhiI_x_" + std::to_string(i+1) + "\""; 
                flow_dataI << phix << " ";
                phiy = "\"PhiI_y_" + std::to_string(i+1) + "\""; 
                flow_dataI << phiy << " ";
                phix = "\"PhiR_x_" + std::to_string(i+1) + "\""; 
                flow_dataR << phix << " ";
                phiy = "\"PhiR_y_" + std::to_string(i+1) + "\""; 
                flow_dataR << phiy << " ";

            }

            flow_dataI << std::endl;
            flow_dataR << std::endl;

            //Write fields
            for ( int i = 0; i < Nr; i++ )
            {

                flow_dataI << i+1 << " ";
                flow_dataI << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                flow_dataI << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
                flow_dataR << i+1 << " ";
                flow_dataR << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                flow_dataR << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";

                for (int j = 0; j < Phi_cut.cols(); j++)
                {

                    flow_dataI << std::setprecision(12) << std::scientific << Phi_cut(i,j).imag() <<  " ";
                    flow_dataI << std::setprecision(12) << std::scientific << Phi_cut(Nr+i,j).imag() <<  " ";
                    flow_dataR << std::setprecision(12) << std::scientific << Phi_cut(i,j).real() <<  " ";
                    flow_dataR << std::setprecision(12) << std::scientific << Phi_cut(Nr+i,j).real() <<  " ";            

                }

            flow_dataI << std::endl;
            flow_dataR << std::endl;

            }
        // Close file
        flow_dataI.close();
        flow_dataR.close();

        } else if ( flag_prob == "VECTOR-3D" || flag_prob == "VELOCITY-3D" )
        {

            // Write row of Headers
            flow_dataI << "\"PointID\"" << " ";
            flow_dataI << "\"x\"" << " ";
            flow_dataI << "\"y\"" << " ";
            flow_dataI << "\"z\"" << " ";

            flow_dataR << "\"PointID\"" << " ";
            flow_dataR << "\"x\"" << " ";
            flow_dataR << "\"y\"" << " ";
            flow_dataR << "\"z\"" << " ";

            std::string phix;
            std::string phiy;
            std::string phiz;

            for ( int i = 0; i < Phi_cut.cols(); i++ )
            {

                phix = "\"PhiI_x_" + std::to_string(i+1) + "\""; 
                flow_dataI << phix << " ";
                phiy = "\"PhiI_y_" + std::to_string(i+1) + "\""; 
                flow_dataI << phiy << " ";
                phiy = "\"PhiI_z_" + std::to_string(i+1) + "\""; 
                flow_dataI << phiz << " ";
                phix = "\"PhiR_x_" + std::to_string(i+1) + "\""; 
                flow_dataR << phix << " ";
                phiy = "\"PhiR_y_" + std::to_string(i+1) + "\""; 
                flow_dataR << phiy << " ";
                phiy = "\"PhiR_z_" + std::to_string(i+1) + "\""; 
                flow_dataR << phiz << " ";

            }

            flow_dataI << std::endl;
            flow_dataR << std::endl;

            //Write fields
            for ( int i = 0; i < Nr; i++ )
            {

                flow_dataI << i+1 << " ";
                flow_dataI << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                flow_dataI << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
                flow_dataI << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";

                flow_dataR << i+1 << " ";
                flow_dataR << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
                flow_dataR << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
                flow_dataR << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";

                for (int j = 0; j < Phi_cut.cols(); j++)
                {

                    flow_dataI << std::setprecision(12) << std::scientific << Phi_cut(i,j).imag() <<  " ";
                    flow_dataI << std::setprecision(12) << std::scientific << Phi_cut(Nr+i,j).imag() <<  " ";
                    flow_dataI << std::setprecision(12) << std::scientific << Phi_cut(2*Nr+i,j).imag() <<  " ";

                    flow_dataR << std::setprecision(12) << std::scientific << Phi_cut(i,j).real() <<  " ";
                    flow_dataR << std::setprecision(12) << std::scientific << Phi_cut(Nr+i,j).real() <<  " ";
                    flow_dataR << std::setprecision(12) << std::scientific << Phi_cut(2*Nr+i,j).real() <<  " ";            

                }

            flow_dataI << std::endl;
            flow_dataR << std::endl;

            }
        // Close file
        flow_dataI.close();
        flow_dataR.close();


        } else 
        {

            std::cout << "Set well problem flag! Exiting ..." << std::endl;
            exit (EXIT_FAILURE);

        }

}



void write_coeffs_sPOD ( const Eigen::MatrixXd &Coeffs,
                        const std::vector<double> &t_vec,
                        const Eigen::VectorXd &lam )
{

    int Ns = t_vec.size(); 

    std::string filename = "Coeffs_sPOD.dat";
    std::ofstream flow_data;
    flow_data.open(filename.c_str());

    // Write row of Headers
    flow_data << "\"Time(s)\"" << " ";

    std::string coef;

    for ( int i = 0; i < Coeffs.cols(); i++ )
    {

        coef = "\"Coef_mode_" + std::to_string(i+1) + "\""; 
        flow_data << coef << " ";

    }                     

    flow_data << "\"Lambda\"";
    flow_data << std::endl;

    // Write coefficients
    for ( int i = 0; i < Coeffs.rows(); i++)
    {

        flow_data << std::setprecision(8) << t_vec[i] << " ";

        for ( int j = 0; j < Coeffs.cols(); j++ )
            flow_data << std::setprecision(8) << Coeffs(i,j) << " ";

        flow_data << std::setprecision(8) << lam(i);
        flow_data << std::endl;

    }

    flow_data.close();


}


void write_TimeDynamics_DMD ( const Eigen::VectorXcd omega,
                            const Eigen::VectorXcd alfa,
                            const Eigen::VectorXd t)
{
    
    std::string filename = "TimeDynamics_DMD.dat";
    std::ofstream flow_data;
    flow_data.open(filename);

    int Nt = t.size();
    int Nm = alfa.size();

    for ( int i = 0; i < Nt; i++ )
    {

        flow_data << std::setprecision(12) << std::scientific << t(i) << " ";
        
        for ( int j = 0; j < Nm; j++ )
        {
            std::complex<double> tdyn = alfa(j)*std::exp(omega(j)*t(i));
            flow_data << std::setprecision(12) << std::scientific << tdyn.real() << " ";
        }
        flow_data << std::endl;

    }

}


void write_CoefsDynamics_mrDMD( std::vector<node_mrDMD> &nodes, const int level, const int ns, const int max_levels )
{
    int sum_inf = 0, sum_sup = 0;
    int level_num;

    nodes_mrDMD_sort( nodes );

    for ( int i = 0; i <= max_levels; i++ )
    {
        sum_inf += std::pow(2,i);
        sum_sup = sum_inf + std::pow(2,i+1);

        if ( (nodes.size() == sum_inf) || ((nodes.size() < sum_sup) && (nodes.size() > sum_inf)) )
        {
            level_num = i+1; 
            break;
        }

    }

    if ( level > level_num - 1  )
    {    
        std::cout << " l > l_max\n Not possible to reconstruct the whole dynamics at this level" << std::endl;
        return;
    }

    std::string filename = "CoefsdynLevel_" + std::to_string(level) + ".dat";
    std::ofstream timedyn_data;
    timedyn_data.open(filename);

    int start_index = 0, end_index = 0;

    for ( int i = 0; i < level; i++ )
        start_index += std::pow(2,i);

    end_index = start_index + std::pow(2,level);

    Eigen::VectorXd rank_bins((int)std::pow(2,level));
    for ( int index = start_index, jj = 0; index < end_index; index++, jj++ )
        rank_bins(jj) = nodes[index].Coefs.size();

    int rank_coefs = rank_bins.maxCoeff();

    for ( int i = 0; i < rank_coefs; i++)
    {
        for ( int j = start_index; j < end_index; j++)
        {
            Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(ns, 0, nodes[j].t_end - nodes[j].t_begin);
            //std::cout << "Time bin " << nodes[j].bin_num << std::endl;
            // std::cout << "Coefficients for the time bin : " << nodes[j].Coefs << std::endl;
            // std::cout << "Eigenvalues for the time bin : " << nodes[j].lam << std::endl;
            for ( int k = 0; k < t.size() - 1; k++ )
            {
                if ( i > (nodes[j].Coefs.size() - 1) )
                    timedyn_data << 0.0 << " ";
                else
                {
                    
                    std::complex<double> value = nodes[j].Coefs(i)*std::exp(std::log(nodes[j].lam(i))/nodes[j].dt*t(k)); 
                 //   std::cout << "value " << value << std::endl;
                    timedyn_data << value.real() << " ";
                }
            }
        }
        timedyn_data << std::endl;

    }
        
    timedyn_data.close();

}


void write_Reconstructed_fields ( Eigen::MatrixXd Rec,
                                    const Eigen::MatrixXd &Coords,
                                    std::string outfile,
                                    std::string flag_prob,
                                    const int nt )
{

    int Nr = Coords.rows(); 

    std::string root_outfile;
    root_outfile.assign ( outfile, 0, outfile.size() - 4);
    std::string out_format;
    out_format.assign ( outfile, outfile.size() - 3, 3);

    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(5) << std::to_string(nt);
    std::string filename = root_outfile + "_" + buffer.str() + "." + out_format;
    std::ofstream flow_data;
    flow_data.open(filename);

    if ( flag_prob == "SCALAR" )
    {

        // Write row of Headers
        flow_data << "\"PointID\"" << " ";
        flow_data << "\"x\"" << " ";
        flow_data << "\"y\"" << " ";

        if ( Coords.cols() == 3 )
            flow_data << "z" << " ";


        flow_data << "\"Rec_field\"" << " ";
        flow_data << std::endl;

        //Write fields
        for ( int i = 0; i < Nr; i++ )
        {

            flow_data << i+1 << " ";
            flow_data << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
            flow_data << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
            if ( Coords.cols() == 3 )
                flow_data << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";
                
            flow_data << std::setprecision(12) << std::scientific << Rec(i,0) <<  " ";           

        flow_data << std::endl;

        }
    // Close file
    flow_data.close();

    } else if ( flag_prob == "VECTOR-2D" || flag_prob == "VELOCITY-2D" )
    {

        // Write row of Headers
        flow_data << "\"PointID\"" << " ";
        flow_data << "\"x\"" << " ";
        flow_data << "\"y\"" << " ";

        std::string Recx;
        std::string Recy;

        Recx = "\"rhoU\""; 
        flow_data << Recx << " ";
        Recy = "\"rhoV\""; 
        flow_data << Recy << " ";
        flow_data << std::endl;

        //Write fields
        for ( int i = 0; i < Nr; i++ )
        {
            
            if ( abs(Rec(i,0)) < 1e-12 )    /*It*/
                Rec(i,0) = 0.0;             /*doesn't*/
            if ( abs(Rec(i,1)) < 1e-12 )    /*mean*/
                Rec(i,1) = 0.0;             /*anything*/
            flow_data << i+1 << " ";
            flow_data << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
            flow_data << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
            flow_data << std::setprecision(12) << std::scientific << Rec(i,0) <<  " ";
            flow_data << std::setprecision(12) << std::scientific << Rec(i,1) <<  " ";            

            flow_data << std::endl;

        }
    // Close file
    flow_data.close();

    } else if ( flag_prob == "VECTOR-3D" || flag_prob == "VELOCITY-3D" )
    {

        // Write row of Headers
        flow_data << "\"PointID\"" << " ";
        flow_data << "\"x\"" << " ";
        flow_data << "\"y\"" << " ";
        flow_data << "\"z\"" << " ";

        std::string Recx;
        std::string Recy;
        std::string Recz;

        Recx = "\"rhoU\""; 
        flow_data << Recx << " ";
        Recy = "\"rhoV\""; 
        flow_data << Recy << " ";
        Recz = "\"rhoW\""; 
        flow_data << Recz << " ";
        flow_data << std::endl;

        //Write fields
        for ( int i = 0; i < Nr; i++ )
        {

            flow_data << i+1 << " ";
            flow_data << std::setprecision(12) << std::scientific <<  Coords(i,0)  << " ";
            flow_data << std::setprecision(12) << std::scientific << Coords(i,1)  << " ";
            flow_data << std::setprecision(12) << std::scientific << Coords(i,2)  << " ";
            flow_data << std::setprecision(12) << std::scientific << Rec(i,0) <<  " ";
            flow_data << std::setprecision(12) << std::scientific << Rec(i,1) <<  " ";
            flow_data << std::setprecision(12) << std::scientific << Rec(i,2) <<  " ";            

        flow_data << std::endl;

        }
    // Close file
    flow_data.close();


    } else 
    {

        std::cout << "Set well problem flag! Exiting ..." << std::endl;
        exit (EXIT_FAILURE);

    }


}



void write_modes ( const Eigen::MatrixXd &Phi_cut )
{
    
    int Nr = Phi_cut.rows(); 
    std::string filename = "ModesRDMD.dat";
    std::ofstream flow_data;
    flow_data.open(filename.c_str()); 

    // Write row of Headers
    std::string phi;

    for ( int i = 0; i < Phi_cut.cols(); i++ )
    {
        phi = "\"Phi_" + std::to_string(i+1) + "\""; 
        flow_data << phi << " ";
    }

    flow_data << std::endl;

    //Write fields
    for ( int i = 0; i < Nr; i++ )
    {    
        for (int j = 0; j < Phi_cut.cols(); j++)
            flow_data << std::setprecision(12) << std::scientific << Phi_cut(i,j) <<  " ";           

        flow_data << std::endl;
    }
    // Close file
    flow_data.close();


}


void write_coefs ( const Eigen::MatrixXd &Coefs )
{
    
    int Nt = Coefs.rows(); 
    std::string filename = "CoefsRDMD.dat";
    std::ofstream flow_data;
    flow_data.open(filename.c_str()); 

    // Write row of Headers
    std::string coef;

    for ( int i = 0; i < Coefs.cols(); i++ )
    {
        coef = "\"Coef_t" + std::to_string(i+1) + "\""; 
        flow_data << coef << " ";
    }

    flow_data << std::endl;

    //Write fields
    for ( int i = 0; i < Nt; i++ )
    {    
        for (int j = 0; j < Coefs.cols(); j++)
            flow_data << std::setprecision(12) << std::scientific << Coefs(i,j) <<  " ";           

        flow_data << std::endl;
    }
    // Close file
    flow_data.close();


}



void write_err_j ( const Eigen::MatrixXd data, std::string filename)
{

    std::ofstream datafile;
    datafile.open(filename);
    if ( data.rows() != data.cols() )
    {
        std::cout << "the matrix of data is not well shaped " << std::endl;
    }

    for( int j = 0; j < data.rows(); j++ ) 
    {
        for ( int nm = 0; nm < data.cols(); nm ++ )
            datafile <<  std::setprecision(8) << data(j,nm) << "\t";

        datafile << std::endl;

    }

    datafile.close();
}


void Write_Restart_Cons_Time ( const Eigen::MatrixXd &Rec,
                                    const Eigen::MatrixXd &Coords,
                                    std::string filename,
                                    const int nt,
                                    const int nC,
                                    const double alpha,
                                    const std::string flag )
{

    std::string root_outputfile;
    std::string file_temp;
    root_outputfile.assign ( filename, 0, filename.size() - 4);
    int Nr = Rec.rows()/nC;
    char str_buf[CGNS_STRING_SIZE];

    for ( int it = 0; it < nt; it++ )
    {
        
        std::ofstream restart_file;
        //Numbering the output format for unsteady reconstruction
        std::stringstream buffer;
        buffer << std::setfill('0') << std::setw(5) << std::to_string(it);
        file_temp = root_outputfile + "_" + buffer.str() + ".dat";
        std::cout << "Writing reconstruction file at time instant " << it << std::endl;

        if ( flag == "NO" ) //Write ASCII file
        {        
            restart_file.open(file_temp);
            //Row of Headers
            if ( nC == 4 ) restart_file << "\"PoinID\"\t\"x\"\t\"y\"\t\"Density\"\t\"X-Momentum\"\t\"Y-Momentum\"\t\"Energy\"" << std::endl;
            else if ( nC == 6 ) restart_file << "\"PoinID\"\t\"x\"\t\"y\"\t\"Density\"\t\"X-Momentum\"\t\"Y-Momentum\"\t\"Energy\"\t\"TKE\"\t\"Omega\"" << std::endl;
            else if ( nC == 7 ) restart_file << "\"PoinID\"\t\"x\"\t\"y\"\t\"z\"\t\"Density\"\t\"X-Momentum\"\t\"Y-Momentum\"\t\"Z-Momentum\"\t\"Energy\"\t\"TKE\"\t\"Omega\"" << std::endl;
            else std::cout << "List of headers not compatible with problem definition" << std::endl;

            //Coordinates and fields
            for ( int i = 0; i < Nr; i++ )
            {
                restart_file << i << "\t";
                for ( int iDim = 0; iDim < Coords.cols(); iDim++ )
                    restart_file << std::setprecision(12) << std::scientific << Coords(i,iDim) << "\t";
                for ( int iC = 0; iC < nC; iC++ )
                    restart_file << std::setprecision(12) << std::scientific << Rec(Nr*iC + i,it) << "\t"; 
                restart_file << std::endl;
            }
            //Write Metadata
            double beta = 0.0; //Since none of our simulations have been using sideslip angles so far
            restart_file << "EXT_ITER= " << it + 1 << std::endl;
            restart_file << "AOA= " << std::setprecision(12) << std::scientific << alpha << std::endl;
            restart_file << "SIDESLIP_ANGLE= " << std::setprecision(12) << std::scientific << beta << std::endl;

            restart_file.close();
        } else  //Write bin file
        {

            int nDim = Coords.cols();
            int Np = Nr;
            int var_buf_size = 5;
            int var_buf[5] = {535532, nC+nDim, Np, 1, 8};
            
            //dummy variables
            double *buf = new double[nC+nDim];
            
            //Prepare metadata
            int Restart_ExtIter = it+1;
            double Restart_Metadata[8] = { alpha, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

            //Prepare row of headers
            std::string Variable_Names[nC+nDim];
            if ( nC == 4)
            {
                Variable_Names[0] = "x";
                Variable_Names[1] = "y";
                Variable_Names[2] = "Density";
                Variable_Names[3] = "X-Momentum";
                Variable_Names[4] = "Y-Momentum";
                Variable_Names[5] = "Energy";
            } else if ( nC == 6 )
            {
                Variable_Names[0] = "x";
                Variable_Names[1] = "y";
                Variable_Names[2] = "Density";
                Variable_Names[3] = "X-Momentum";
                Variable_Names[4] = "Y-Momentum";
                Variable_Names[5] = "Energy";
                Variable_Names[6] = "TKE";
                Variable_Names[7] = "Omega";

            } else if ( nC == 7 )
            {
                Variable_Names[0] = "x";
                Variable_Names[1] = "y";
                Variable_Names[2] = "z";
                Variable_Names[3] = "Density";
                Variable_Names[4] = "X-Momentum";
                Variable_Names[5] = "Y-Momentum";
                Variable_Names[6] = "Z-Momentum";
                Variable_Names[7] = "Energy";
                Variable_Names[8] = "TKE";
                Variable_Names[9] = "Omega";   
            } else
            {
                std::cout << "Number of conservative variables not compatible with available cases " << std::endl;
                exit(EXIT_FAILURE);
            }
            
            restart_file.open(file_temp, std::ios::binary);
            if (!restart_file.is_open()) {
                std::cout << "Unable to open bin file " << std::endl;
            }

            restart_file.write((char*) var_buf, var_buf_size*sizeof(int) );
            for ( int iVar = 0; iVar < nC+nDim; iVar++) 
            {
                strncpy(str_buf, Variable_Names[iVar].c_str(), CGNS_STRING_SIZE);
                restart_file.write(str_buf, CGNS_STRING_SIZE*sizeof(char));
            }

            for ( int iPoint = 0; iPoint < Np; iPoint ++)
            { 
                for ( int iVar = 0; iVar < nDim; iVar ++ )
                    buf[iVar] = Coords(iPoint, iVar);
                for ( int iVar = nDim; iVar < (nC + nDim); iVar ++ )
                    buf[iVar] = Rec(iPoint, iVar-nDim);

                restart_file.write((char*) buf, sizeof(double)*(nC+nDim));
            }

            //Writing metadata
            restart_file.write((char*) &Restart_ExtIter, sizeof(int));
            restart_file.write((char*)Restart_Metadata, 8.0*sizeof(double));

            restart_file.close();
            //delete [] buf;

        }
    }
}