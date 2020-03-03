	

#include <cmath>
#include <initializer_list>
#include <iostream>
#include <utility>

// #include <pagmo/problem.hpp>
// #include <pagmo/types.hpp>
#include "Extract_Basis.hpp"
#include "pagmo.hpp"
#include "Opt_struct.hpp"

// using namespace Eigen;
// using Eigen::all;


using namespace pagmo;

// Our simple example problem, version 0.
struct HimmelblauFunction {
    // Empty constructor
    // Without an empty constructor the problem is not accepted
    // as a multithreading type
    HimmelblauFunction( ){ }

    //Actual constructor allowing the user to define the boundaries
    HimmelblauFunction( const double x_min, const double x_max, const double y_min,
            const double y_max ) :
        x_min_( x_min ), x_max_( x_max ), y_min_( y_min ), y_max_( y_max )
    { }

    // Mandatory, computes the fitness, i.e. the Himmelblau's function
    std::vector< double > fitness( const std::vector< double > &x ) const {
        std::vector< double > return_value;
//        return_value.push_back( pow( x[0]*x[0] + x[1] - 11.0, 2.0 )
//                + pow( x[0] + x[1]*x[1] - 7.0, 2.0 ) );
std::cout << "Sum variables = " << x[0] + x[1] << std::endl;
        return {pow( x[0]*x[0] + x[1] - 11.0, 2.0 )
                + pow( x[0] + x[1]*x[1] - 7.0, 2.0 ), x[0] + x[1] };
    }

    // Mandatory, returns the box-bounds
    std::pair< std::vector< double >, std::vector< double > > get_bounds( ) const {
        
        std::pair< std::vector< double >, std::vector< double > > box_bounds;
        box_bounds.first.push_back( x_min_ );
        box_bounds.first.push_back( y_min_ );
        box_bounds.second.push_back( x_max_ );
        box_bounds.second.push_back( y_max_ );

        return box_bounds;
    }

    vector_double::size_type get_nec() const {
        return 0;
    }

    vector_double::size_type get_nic() const {
        return 1;
    }

    private:
        double x_min_;
        double y_min_;
        double x_max_;
        double y_max_;


};

int main()
{
   
    //Set seed for reproducible results
    pagmo::random_device::set_seed( 12345 );

    // Define Himmelblau function (range [-5,5]x[-5,5])
    pagmo::problem prob{ HimmelblauFunction( -5, 5, -5, 5) };

    std::cout << "Problem Definition with pagmo \n\n" << prob << std::endl;


    // Solve using DE algorithm
    pagmo::algorithm algo{ pagmo::gaco( ) };

    // Create island with 1000 individuals
    pagmo::island isl = pagmo::island{ algo, prob, 1000 };
    std::cout << "I m here \n\n" << std::endl;
    // Evolve for 1000 generations
    for( int i = 1; i <= 100; i++ ) {
        isl.evolve();
        while (isl.status() != pagmo::evolve_status::idle) {
            isl.wait();
    }

        // printPopulationToFile( isl.get_population( ).get_x( ), "himmelblau_" + std::to_string( i ) , false );
        // printPopulationToFile( isl.get_population( ).get_f( ), "himmelblau_" +  std::to_string( i ) , true );

        // Print current optimum to console
        std::cout << "Minimum: " <<i<<" "<<std::setprecision( 16 ) <<"f= "<< isl.get_population().champion_f()[0] <<", x="<<
                    isl.get_population().champion_x()[0] <<" y="<<isl.get_population().champion_x()[1] <<std::endl;
    }


    return 0;

    // 5 - Output the population
    // std::cout << "The population: \n" << pop;

}














// Example main
// int main()
// {
//     prob_settings settings;
//     std::vector< std::vector< double > > bounds( 2, std::vector< double >( 3, 0.0 ) );

//     // Define search bounds:
//     //I set the bounds of all the variables (...list them here...):
//     bounds[ 0 ][ 0 ] = 0.0; //Lower Bound Var. 1
//     bounds[ 1 ][ 0 ] = 100.0; //Upper Bound Var.1
//     bounds[ 0 ][ 1 ] = 0.0; //Lower Bound Var. 2
//     bounds[ 1 ][ 1 ] = 200.0; //Upper Bound Var. 2
//     bounds[ 0 ][ 2 ] = 0.0; //Lower Bound Var. 3
//     bounds[ 1 ][ 2 ] = 300.0; //Upper Bound Var. 3
//     int Ns = 10;
//     int Np = 10;

//     std::vector<double> v(4, 100.0);
//     double* ptr = &v[0];
//     Eigen::Map<Eigen::VectorXd> my_vect(ptr, 4);

// // https://stackoverflow.com/questions/30181914/copy-data-from-stdvector-to-eigens-matrixxd-in-c

//     std::cout << "Vector mapped my_vec : " << my_vect << std::endl;


//     std::vector<int> dat(4);
//     int i = 0;
//     dat[i] = i + 1; i++;
//     dat[i] = i + 1;
//     dat[i+1] = i + 1; i++;
//     dat[i+1] = i + 1;
//     // typedef Eigen::Matrix<int, -1, -1, Eigen::ColMajor> Cm;
//     Eigen::Map<Eigen::VectorXi> m1(dat.data(), dat.size());
//     std::cout << "Vector std copied in vector c++ : "<< m1 << std::endl;

//     std::cout << "Vector before unique :\n ";
//     for ( int i = 0; i < dat.size(); i++ )
//         std::cout << dat[i] << std::endl;

//     std::sort(dat.begin(),dat.end());
//     dat.erase(std::unique(dat.begin(),dat.end()),dat.end());

//     std::cout << "Vector after unique :\n ";
//     for ( int i = 0; i < dat.size(); i++ )
//         std::cout << dat[i] << std::endl;


//     std::vector<int> ind{1,0};

//     Eigen::MatrixXd sn_set = Eigen::MatrixXd::Zero(Np,Ns);
//     sn_set.col(0) = Eigen::VectorXd::Ones(Np);
//     Eigen::ArrayXi ri = ArrayXi::LinSpaced(Np,0,Np-1);
//     Eigen::ArrayXi ci(3); ci << 1,0,0;

//     Eigen::MatrixXd B = indexing(sn_set, ri, ci);
//     std::cout << "Matrix with selected columns:\n " << B << std::endl;

//     // 2-Problem initialization:
//     pagmo::problem prob{SPOD_Adapt_Samp(bounds, sn_set, settings)};
//     std::cout << "First lower bound: " << prob.get_lb()[0] << std::endl;
//     std::cout << "First upper bound: " << prob.get_ub()[0] << std::endl;
//     std::cout << "Second lower bound: " << prob.get_lb()[1] << std::endl;
//     std::cout << "Second upper bound: " << prob.get_ub()[1] << std::endl;

//     std::cout << "End test " << std::endl;

// }




//Example on website
// using namespace pagmo;


// // Our simple example problem, version 0.
// struct problem_v0 {
//     // Implementation of the objective function.
//     vector_double fitness(const vector_double &deltas) const
//     {
//         return {deltas[0] * deltas[3] * (deltas[0] + deltas[1] + deltas[2]) + deltas[2]};
//     }
//     // Implementation of the box bounds.
//     std::pair<vector_double, vector_double> get_bounds() const
//     {
//         return {{1., 1., 1., 1.}, {5., 5., 5., 5.}};
//     }
// };

// int main()
// {
//     // Construct a pagmo::problem from our example problem.
//     problem p{problem_v0{}};

//     // Compute the value of the objective function
//     // in the point (1, 2, 3, 4).
//     std::cout << "Value of the objfun in (1, 2, 3, 4): " << p.fitness({1, 2, 3, 4})[0] << '\n';

//     // Fetch the lower/upper bounds for the first variable.
//     std::cout << "Lower bounds: [" << p.get_lb()[0] << "]\n";
//     std::cout << "Upper bounds: [" << p.get_ub()[0] << "]\n\n";

//     // Print p to screen.
//     std::cout << p << '\n';
// }







































// // Our simple example problem, version 0.
// struct problem_v0 {

//     problem_v0(){};
//     problem_v0( std::vector< std::vector< double > > &bounds ){

//     }
//     // Implementation of the objective function.
//     vector_double fitness(const vector_double &dv) const
//     {
//         return {dv[0] * dv[3] * (dv[0] + dv[1] + dv[2]) + dv[2]};
//     }
//     // Implementation of the box bounds.
//     std::pair<vector_double, vector_double> get_bounds() const
//     {
//         // return {{1., 1., 1., 1.}, {5., 5., 5., 5.}};
//         return {problemBounds_[0], problemBounds_[1]};
//     }

// private:

//     const std::vector< std::vector< double > > problemBounds_;

// };

// int main()
// {
//     std::vector< std::vector< double > > bounds( 2, std::vector< double >( 4, 0.0 ) );
//     bounds[ 0 ][ 0 ] = 1.0; //Lower Bound Var. 1
//     bounds[ 1 ][ 0 ] = 5.0; //Upper Bound Var.1
//     bounds[ 0 ][ 1 ] = 1.0; //Lower Bound Var. 2
//     bounds[ 1 ][ 1 ] = 5.0; //Upper Bound Var. 2
//     bounds[ 0 ][ 2 ] = 1.0; //Lower Bound Var. 3
//     bounds[ 1 ][ 2 ] = 5.0; //Upper Bound Var. 3
//     bounds[ 0 ][ 3 ] = 1.0; //Lower Bound Var. 4
//     bounds[ 1 ][ 3 ] = 5.0; //Upper Bound Var. 4
//     // Construct a pagmo::problem from our example problem.
    
//     problem_v0 test(bounds);
//     std::cout << "Value of the lowbound for 1st variable: " << test.get_bounds() << '\n';
    
//     // problem p{test};


//     // // Compute the value of the objective function
//     // // in the point (1, 2, 3, 4).
//     // std::cout << "Value of the objfun in (1, 2, 3, 4): " << p.fitness({1, 2, 3, 4})[0] << '\n';

//     // // Fetch the lower/upper bounds for the first variable.
//     // std::cout << "Lower bounds: [" << p.get_lb()[0] << "]\n";
//     // std::cout << "Upper bounds: [" << p.get_ub()[0] << "]\n\n";

//     // // Print p to screen.
//     // std::cout << p << '\n';
// }
