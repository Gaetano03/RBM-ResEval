
#include "read_Inputs.hpp"

keywords read_keyword_type( const std::string &key_string )
{
    if( key_string == "ADAPT_POS" )
        return ADAPT_POS;
    if( key_string == "SOLVER" )
        return SOLVER;
    if( key_string == "NS" )
        return NS;
    else if( key_string == "DS" )
        return DS;
    else if( key_string == "EN" )
        return EN;
    else if( key_string == "SIGMA" )
        return SIGMA;
    else if( key_string == "NSTART" )
        return NSTART;
    else if( key_string == "ALPHA" )
        return ALPHA;
    else if( key_string == "BETA" )
        return BETA;
    else if( key_string == "RHO_FREE" )
        return RHO_FREE;
    else if( key_string == "P_FREE" )
        return P_FREE;
    else if( key_string == "MACH" )
        return MACH;
    else if( key_string == "TEMPERATURE" )
        return TEMPERATURE;
    else if( key_string == "REYNOLDS" )
        return REYNOLDS;
    else if( key_string == "VISCOSITY" )
        return VISCOSITY;
    else if( key_string == "NDIM" )
        return NDIM;
    else if( key_string == "DT_CFD" )
        return DT_CFD;
    else if( key_string == "FLAG_DIM" )
        return FLAG_DIM;
    else if( key_string == "FLAG_PROB" )
        return FLAG_PROB;
    else if( key_string == "INPUT_FILE" )
        return INPUT_FILE;
    else if( key_string == "OUTPUT_FILE" )
        return OUTPUT_FILE;
    else if( key_string == "NF" )
        return NF;
    else if( key_string == "COLS_COORDS" )
        return COLS_COORDS;
    else if( key_string == "COLS_FIELDS" )
        return COLS_FIELDS;
    else if( key_string == "FLAG_METHOD" )
        return FLAG_METHOD;
    else if( key_string == "FLAG_MEAN" )
        return FLAG_MEAN;
    else if( key_string == "FLAG_BC" )
        return FLAG_BC;
    else if( key_string == "FLAG_FILTER" )
        return FLAG_FILTER;
    else if( key_string == "FLAG_WDB_BE" )
        return FLAG_WDB_BE;
    else if( key_string == "FLAG_REC" )
        return FLAG_REC;
    else if( key_string == "FLAG_INTERP" )
        return FLAG_INTERP;
    else if( key_string == "T_REC" )
        return T_REC;
    else if( key_string == "T_RES" )
        return T_RES;
    else if( key_string == "DT_RES" )
        return DT_RES;
    else if( key_string == "RANK_RDMD" )
        return RANK_RDMD;
    else if( key_string == "RANK" )
        return RANK;
    else if( key_string == "HO_D" )
        return HO_D;
    else if( key_string == "DMD_COEF_FLAG" )
        return DMD_COEF_FLAG;
    else if( key_string == "MAX_CYCLES" )
        return MAX_CYCLES;
    else if( key_string == "MAX_LEVELS" )
        return MAX_LEVELS;
    else if( key_string == "TOL" )
        return TOL;
    else
    {
        std::cout << key_string << " Not Available ... Something wrong in cfg file" << std::endl;
        exit (EXIT_FAILURE);       
    }

}


void Read_cfg ( const std::string filename, prob_settings &settings )
{

//Initializing values in settings to NULL or default values
    settings.Ns = 0;                         
    settings.Ds = 0;                         
    settings.nstart = 0;                     
    settings.ndim = 0;
    settings.solver = "NONE" ;
    settings.in_file = "NONE";            
    settings.out_file = "NONE";           
    settings.flag_dim = "NONE";           
    settings.flag_prob = "NONE";          
    settings.Cols = {};          
    settings.Cols_coords = {};    
    settings.flag_method = "NONE";        
    settings.flag_wdb_be = "NONE"; 
    settings.P = 101325;
    settings.Rho = 1.228;       
    settings.Dt_cfd = 0.0;                
    settings.alpha = 0.0;                 
    settings.beta = 0.0;                  
    settings.Mach = 0.0;                  
    settings.Re = 0.0;                    
    settings.mu = 0.0;                    
    settings.T = 0.0;                     
    settings.Nf = -1;                    
    settings.En = -1.0;                  
    settings.sigma = -1.0;               
    settings.flag_filter = "NONE";    
    settings.flag_mean = "NO";      
    settings.flag_bc = "NONE";       
    settings.r = -2;                      
    settings.dmd_coef_flag = "NONE"; 
    settings.max_cycles = -1;
    settings.max_levels = -1;
    settings.d = 0;                      
    settings.r_RDMD = -1;                
    settings.flag_rec = "NO";          
    settings.flag_interp = "NONE";       
    settings.t_rec = {};                  
    settings.tol = 0.0;                    


    std::ifstream cFile ( filename );
    if ( cFile.is_open() )
    {

        size_t delimiterPos; 
        std::string name, value;
        std::string line;

        while( getline(cFile, line) )
        {
            line.erase( remove_if(line.begin(), line.end(), isspace),
                                 line.end() );  //include ::  in front of isspace if using namespace std
            if( line[0] == '#' || line.empty() )
                continue;
            delimiterPos = line.find("=");
            name = line.substr(0, delimiterPos);
            value = line.substr(delimiterPos + 1);

            switch (read_keyword_type(name))
            {
                case SOLVER:
                {
                    settings.solver = value;
                    //std::cout << "Problem flag : " << value << std::endl;
                    break;
                }

                case FLAG_PROB:
                {
                    settings.flag_prob = value;
                    //std::cout << "Problem flag : " << value << std::endl;
                    break;
                }

                case FLAG_DIM:
                {
                    settings.flag_dim = value;
                    //std::cout << "Dimension flag : " << value << std::endl;
                    break;
                }

                case NS:
                {
                    settings.Ns = std::stoi(value);
                    //std::cout << "Number of snapshots : " << value << std::endl;
                    break;
                }

                case DS:
                {
                    settings.Ds = std::stod(value);
                    //std::cout << "Delta between selected snapshots (Equispaced) : " << value << std::endl; 
                    break;
                }

                case EN:
                {
                    settings.En = std::stod(value);
                    //std::cout << "Energy content used in the reconstruction : " << value << std::endl;
                    break;
                }

                case SIGMA:
                {
                    settings.sigma = std::stod(value);
                    //std::cout << "Sigma for SPOD gaussian filter : " << value << std::endl;
                    break;
                }

                case TOL:
                {
                    settings.tol = std::stod(value);
                    //std::cout << "Sigma for SPOD gaussian filter : " << value << std::endl;
                    break;
                }                

                case NSTART:
                {
                    settings.nstart = std::stoi(value);
                    //std::cout << "Initial snapshot number : " << value << std::endl;
                    break;
                }

                case DT_RES:
                {
                    std::string str = value;
                    std::stringstream ss(str);

                    double i;

                    while (ss >> i) {

                        settings.Dt_res.push_back(i);

                        if (ss.peek() != ',' || ss.peek() != ' ')
                            ss.ignore();
                    }

                    break;
                }

                case ALPHA:
                {
                    settings.alpha = std::stod(value);
                    //std::cout << "Sigma for SPOD gaussian filter : " << value << std::endl;
                    break;
                }

                case BETA:
                {
                    settings.beta = std::stod(value);
                    //std::cout << "Sigma for SPOD gaussian filter : " << value << std::endl;
                    break;
                }

                case P_FREE:
                {
                    settings.P = std::stod(value);
                    //std::cout << "Sigma for SPOD gaussian filter : " << value << std::endl;
                    break;
                }

                case RHO_FREE:
                {
                    settings.Rho = std::stod(value);
                    //std::cout << "Sigma for SPOD gaussian filter : " << value << std::endl;
                    break;
                }

                case REYNOLDS:
                {
                    settings.Re = std::stod(value);
                    //std::cout << "Sigma for SPOD gaussian filter : " << value << std::endl;
                    break;
                }

                case MACH:
                {
                    settings.Mach = std::stod(value);
                    //std::cout << "Sigma for SPOD gaussian filter : " << value << std::endl;
                    break;
                }

                case TEMPERATURE:
                {
                    settings.T = std::stod(value);
                    //std::cout << "Sigma for SPOD gaussian filter : " << value << std::endl;
                    break;
                }

                case VISCOSITY:
                {
                    settings.mu = std::stod(value);
                    //std::cout << "Sigma for SPOD gaussian filter : " << value << std::endl;
                    break;
                }

                case NDIM:
                {
                    settings.ndim = std::stoi(value);
                    //std::cout << "Initial snapshot number : " << value << std::endl;
                    break;
                }
                
                case DT_CFD:
                {
                    settings.Dt_cfd = std::stod(value);
                    //std::cout << "Dt used in CFD simulation : " << value << std::endl;
                    break;
                }

                case INPUT_FILE:
                {
                    settings.in_file = value;
                    //std::cout << "Input file root name and format : " << value << std::endl;
                    break;
                }

                case OUTPUT_FILE:
                {
                    settings.out_file = value;
                    //std::cout << "Output file root name and format : " << value << std::endl;
                    break;
                }

                case NF:
                {
                    settings.Nf = std::stoi(value);
                    //std::cout << "Filter size for feature extraction : " << value << std::endl;
                    break;
                }

                case FLAG_METHOD:
                {
                    settings.flag_method = value;
                    //std::cout << "Method for feature extraction : " << value << std::endl;
                    break;
                }

                case FLAG_MEAN:
                {
                    settings.flag_mean = value;
                    //std::cout << "Mean subtraction : " << value << std::endl;
                    break;
                }

                case FLAG_BC:
                {
                    settings.flag_bc = value;
                    //std::cout << "Boundary consition for correlation matrix : " << value << std::endl;
                    break;
                }

                case FLAG_FILTER:
                {
                    settings.flag_filter = value;
                    //std::cout << "Filter type for SPOD : " << value << std::endl;
                    break;
                }

                case FLAG_WDB_BE:
                {
                    settings.flag_wdb_be = value;
                    //std::cout << "Write database basis extraction (modes and coefficients) : " << value << std::endl;
                    break;
                }

                case FLAG_REC:
                {
                    settings.flag_rec = value;
                    //std::cout << "Compute reconstructed field : " << value << std::endl;
                    break;
                }

                case FLAG_INTERP:
                {
                    settings.flag_interp = value;
                    //std::cout << "Interpolation technique for rbf : " << value << std::endl;
                    break;
                }

                case RANK:
                {
                    settings.r = std::stoi(value);
                    //std::cout << "DMD rank : " << value << std::endl;
                    break;
                }

                case RANK_RDMD:
                {
                    settings.r_RDMD = std::stoi(value);
                    //std::cout << "Recursive-DMD rank : " << value << std::endl;
                    break;
                }

                case DMD_COEF_FLAG:
                {
                    settings.dmd_coef_flag = value;
                    //std::cout << "DMD coefs method : " << value << std::endl;
                    break;
                }

                case HO_D:
                {
                    settings.d = std::stoi(value);
                    //std::cout << "DMD coefs method : " << value << std::endl;
                    break;
                }

                case MAX_CYCLES:
                {
                    settings.max_cycles = std::stoi(value);
                    //std::cout << "Max_cycles for mrDMD : " << value << std::endl;
                    break;
                }

                case MAX_LEVELS:
                {
                    settings.max_levels = std::stoi(value);
                    //std::cout << "Max_levels for mrDMD : " << value << std::endl;
                    break;
                }

                case T_REC:
                {
                    std::string str = value;
                    std::stringstream ss(str);

                    double i;

                    while ( ss >> i ) {

                        settings.t_rec.push_back(i);

                        if (ss.peek() != ',' || ss.peek() != ' ')
                            ss.ignore();

                    }

                    //std::cout << "Times desired for reconstruction: \t";
                    
                    //for ( i = 0; i < settings.t_rec.size(); i++ )
                        //std::cout << settings.t_rec[i] << "\t";

                    //std::cout << std::endl;

                    break;
                }

                case T_RES:
                {
                    std::string str = value;
                    std::stringstream ss(str);

                    double i;

                    while ( ss >> i ) {

                        settings.t_res.push_back(i);

                        if (ss.peek() != ',' || ss.peek() != ' ')
                            ss.ignore();

                    }

                    break;
                }

                case COLS_COORDS:
                {
                    std::string str = value;
                    std::stringstream ss(str);

                    int i;

                    while ( ss >> i ) {

                        settings.Cols_coords.push_back(i);

                        if (ss.peek() != ',' || ss.peek() != ' ')
                            ss.ignore();

                    }

                    break;
                }

                case ADAPT_POS:
                {
                    std::string str = value;
                    std::stringstream ss(str);

                    int i;

                    while ( ss >> i ) {

                        settings.t_pos.push_back(i);

                        if (ss.peek() != ',' || ss.peek() != ' ')
                            ss.ignore();

                    }

                    break;
                }

                case COLS_FIELDS:
                {
                    std::string str = value;
                    std::stringstream ss(str);

                    int i;

                    while ( ss >> i ) {

                        settings.Cols.push_back(i);

                        if (ss.peek() != ',' || ss.peek() != ' ')
                            ss.ignore();

                    }

                    break;
                }

                default:
                {
                    break;
                }

            }

        }

        cFile.close();

    }
    else
    {
        std::cout << "Unable to open config file! Terminating ... " << '\n';
        exit (EXIT_FAILURE);
    }

}

// Read and change su2 file
void Modify_su2_cfg ( std::string file_in, std::string file_out, double dt_res ) {

    std::ifstream inFile ( file_in );
    std::ofstream outFile ( file_out );

    if ( inFile.is_open() && outFile.is_open() ) {

        size_t delimiterPos;
        std::string name, value;
        std::string line;

        while (getline(inFile, line)) {

            line.erase(remove_if(line.begin(), line.end(), isspace),
                       line.end());  //include ::  in front of isspace if using namespace std
            if (line[0] == '#' || line.empty())
                continue;

            delimiterPos = line.find("=");
            name = line.substr(0, delimiterPos);

            if ( name == "UNST_TIMESTEP" )
                outFile << "UNST_TIMESTEP=" << std::setprecision(16) << dt_res;
            else
                outFile << line;

            outFile << std::endl;

        }

        outFile.close();
        inFile.close();
    } else {
        std::cout << "Problem with SU2 config file " << std::endl;
        exit (EXIT_FAILURE);
    }

}


int N_gridpoints( const std::string file_in) {

    std::string line_flow_data;
    std::ifstream flow_data;
    flow_data.open(file_in.c_str());

    if( !flow_data.is_open() )
    {

        std::cout << " While getting number of grid points, \nFile: " 
            << file_in << " not found " << std::endl;
        exit (EXIT_FAILURE);

    }

    int n_row = 0;


    while(getline( flow_data, line_flow_data )) 
    {
           
        if ( line_flow_data.compare(0,1,"E") == 0 || line_flow_data.compare(0,1,"A") == 0 ) //if reading SU2 native restart file
            break;
            
        n_row++;

   } 

    flow_data.close();

    return (n_row-1);

}




Eigen::MatrixXd read_col( std::string filename, int Nr, std::vector<int> Cols )
{


    Eigen::MatrixXd field (Nr, Cols.size());
    std::ifstream flow_data;
    flow_data.open( filename );    

    if ( !flow_data.is_open() )
    {
        std::cout << "File : " << filename << " not found" << std::endl;    
        exit (EXIT_FAILURE);
    }

    std::string line_flow_data ;

    // Read row of headers
    getline( flow_data, line_flow_data );

    int n_row = 0;

    while ( getline( flow_data, line_flow_data ) )
    {

        Eigen::RowVectorXd point(Cols.size());
        std::istringstream iss(line_flow_data);
        std::string token;
        long double rubbish;
        int count = 0, c = 0; 

        if ( filename.compare( filename.size()-3, 3, "csv") == 0 ) 
        {

            while( getline( iss, token, ',') ) 
            {

                rubbish = std::stold(token);

                //This trick is due to some bad values I found in the restart SU2 file
                if ( (rubbish!=0.0) && ((std::abs(rubbish) < 1.0e-300) || (std::abs(rubbish) > 1.0e300)) )
                {
                    std::cout << " Rubbish value : " << std::setprecision(17) << std::scientific << 
                        rubbish <<  " on the row : "<< n_row << std::endl;
                    rubbish = 0.0;
                }
                
                if ( count == Cols[c])
                {

                    point(c) = rubbish;
                    c++;
                }
                count ++;
            }

        } else if ( filename.compare( filename.size()-3, 3, "dat") == 0 ) 
        {

            while( getline( iss, token, '\t') )
            {

                if ( token.compare(0,1,"E") == 0 || token.compare(0,1,"A") == 0 ) //if reading SU2 restart file
                    break;

                rubbish = std::stold(token);
                

                //This trick is due to some bad values I found in the restart SU2 file
                if ( (rubbish!=0.0) && ((std::abs(rubbish) < 1.0e-300) || (std::abs(rubbish) > 1.0e300)))
                {
                    std::cout << " Rubbish value : " << std::setprecision(17) << std::scientific << 
                        rubbish <<  " on the row : "<< n_row << std::endl;
                    rubbish = 0.0;
                }

                if ( count == Cols[c])
                {

                    point(c) = rubbish;
                    c++;
                }
                count ++;

            } 

        } else
        {
            std::cout << "Bad file format" << std::endl;
            exit (EXIT_FAILURE);
        }

        if ( token.compare(0,1,"E") == 0 || token.compare(0,1,"A") == 0 ) //if reading SU2 restart file
            break;

        field.row(n_row) = point; 
        n_row++;

    }

    flow_data.close();

    return field;
  
}


Eigen::MatrixXd read_modes( std::string filename, int Nr, int r_RDMD )
{
    Eigen::MatrixXd f(Nr, r_RDMD);

    std::ifstream Modes_data;
    Modes_data.open( filename );

        if ( !Modes_data.is_open() )
    {
        std::cout << "File : " << filename << " not found" << std::endl;    
        exit (EXIT_FAILURE);
    }

    std::string line_flow_data ;

    // Read row of headers
    getline( Modes_data, line_flow_data );

    int n_row = 0;

    while ( getline( Modes_data, line_flow_data ) )
    {

        Eigen::RowVectorXd point(r_RDMD);
        std::istringstream iss(line_flow_data);
        std::string token;
        double rubbish;
        int count = 0, c = 0; 

        while( getline( iss, token, ' ') )
        {
            rubbish = std::stod(token);
            if ( count < r_RDMD )
                point(count) = rubbish;
            
            count ++;
        } 

        f.row(n_row) = point; 
        n_row++;

    }

    Modes_data.close();

    return f;

}


Eigen::MatrixXd read_coefs( std::string filename, int Ns, int r_RDMD )
{

    Eigen::MatrixXd f(r_RDMD, Ns);

    std::ifstream Coefs_data;
    Coefs_data.open( filename );

    if ( !Coefs_data.is_open() )
    {
        std::cout << "File : " << filename << " not found" << std::endl;    
        exit (EXIT_FAILURE);
    }

    std::string line_flow_data ;

    // Read row of headers
    getline( Coefs_data, line_flow_data );

    int n_row = 0;

    int count = 0;
    while ( getline( Coefs_data, line_flow_data ) && n_row < r_RDMD )
    {

        Eigen::RowVectorXd point(Ns);
        std::istringstream iss(line_flow_data);
        std::string token;
        double rubbish;
        int count = 0, c = 0; 

        while( getline( iss, token, ' ') && count < Ns )
        {
            rubbish = std::stod(token);
            point(count) = rubbish;
            
            count ++;
        } 

        f.row(n_row) = point; 
        n_row++;

    }

    Coefs_data.close();

    return f;
}



Eigen::MatrixXd read_err_j ( std::string filename, int Ns )
{
        std::ifstream file_data;
        file_data.open( filename );

            if ( !file_data.is_open() )
        {
            std::cout << "File : " << filename << " not found" << std::endl;    
            exit (EXIT_FAILURE);
        }

        std::string line_flow_data ;

        int n_row = 0, count = 0;
        Eigen::MatrixXd Err_map = Eigen::MatrixXd::Zero(Ns, Ns); 

        while ( getline( file_data, line_flow_data ) )
        {

            
            std::istringstream iss(line_flow_data);
            std::string token;
            double err;
            count = 0; 

            while( getline( iss, token, '\t') )
            {
                err = std::stod(token);
                Err_map(n_row, count) = err;

                count ++;
            } 
 
            n_row++;

        }

        file_data.close();

        return Err_map;

}


plot3d_info read_plot3d_info (std::string filename)
{

    plot3d_info Info;
    int nC = 5; // Density, Momentum (3 Components), Total Energy

    // open file in binary mode
    std::ifstream input(filename, std::ios::binary);

    if (input.good())
    {
     // read number of points
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
        input.read((char*) &Info.nblocks, sizeof(int));
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);

        Info.ni.resize(Info.nblocks);
        Info.nj.resize(Info.nblocks);
        Info.nk.resize(Info.nblocks);
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
        
        for ( int i = 0; i < Info.nblocks; i++)
        {
            input.read((char*) &Info.ni[i], sizeof(int));
            input.read((char*) &Info.nj[i], sizeof(int));
            input.read((char*) &Info.nk[i], sizeof(int));
        }

        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
        input.read((char*) &Info.M, sizeof(float));
        input.read((char*) &Info.alpha, sizeof(float));
        input.read((char*) &Info.Re, sizeof(float));
        input.read((char*) &Info.T, sizeof(float));

        input.close();
        // std::cout << "Successfully read plot3d Info" << std::endl;
    }
    else
    {
     std::cout << "Unable to open file .q" << std::endl;
    }

    return Info;
}


std::vector<Eigen::VectorXd> read_plot3d (std::string filename, plot3d_info Info)
{

    int nC = 5; // Density, Momentum (3 Components), Total Energy
    // Read .q file fortran binary

    float dum_f;
    int dum_i;
    // open file in binary mode
    std::ifstream input(filename, std::ios::binary);

    std::vector<Eigen::VectorXd> flow_data(Info.nblocks);

    if (input.good())
    {
     // read number of points
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
        input.read((char*) &dum_i, sizeof(int));
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);

        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
        
        for ( int i = 0; i < Info.nblocks; i++)
        {

            input.read((char*) &dum_i, sizeof(int));
            input.read((char*) &dum_i, sizeof(int));
            input.read((char*) &dum_i, sizeof(int));
        
        }

        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);

        for ( int i = 0; i < Info.nblocks; i++ )
        {
            Eigen::VectorXf data_i(nC*Info.ni[i]*Info.nj[i]*Info.nk[i]);
            input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
            input.read((char*) &dum_f, sizeof(float));
            input.read((char*) &dum_f, sizeof(float));
            input.read((char*) &dum_f, sizeof(float));
            input.read((char*) &dum_f, sizeof(float));
            input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);

            input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
            
            flow_data[i] = Eigen::VectorXd::Zero(nC*Info.ni[i]*Info.nj[i]*Info.nk[i]);
            for ( int j = 0; j < nC*Info.ni[i]*Info.nj[i]*Info.nk[i]; j++ )
                input.read((char*) &data_i(j), sizeof(float));

            flow_data[i] = data_i.cast<double>();
            input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);

        }

        input.close();
        // std::cout << "Successfully read .q file" << std::endl;
        return flow_data;

    }
    else
    {
     std::cout << "Unable to open .q file" << std::endl;
     std::cout << "Terminating ... " << std::endl;
     exit (EXIT_FAILURE);
    }


}

