#ifndef GENERATE_SNSET_HPP
#define GENERATE_SNSET_HPP


#include "read_Inputs.hpp"

Eigen::MatrixXd generate_snap_matrix( const int Nr, const int Ns, const int ds, const int init,
                                        std::vector<int> Cols,
                                        std::string inputfile,
                                        std::string flag_prob = "VELOCITY-2D");



#endif // GENERATE_SNSET_HPP