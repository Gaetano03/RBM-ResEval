#ifndef OPT_STRUCT_HPP
#define OPT_STRUCT_HPP

#include <initializer_list>
#include <utility>

#include "read_Inputs.hpp"
#include "Extract_Basis.hpp"
#include "pagmo.hpp"

// Define the problem PaGMO-style
struct SPOD_Adapt_Samp {

    // Empty constructor
    SPOD_Adapt_Samp( ){ }

    SPOD_Adapt_Samp( std::vector< std::vector< double > > &bounds, Eigen::MatrixXd &sn_set, prob_settings &settings,
            int &Nf ) : problemBounds_(bounds), m_sn_set(sn_set), m_settings(settings), m_Nf(Nf) {
        }
    // Fitness:
    std::vector<double> fitness(const std::vector<double> &variables) const;
    // Boundaries of the problem
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const;
    //! Serialization function for Pagmo compatibility
    // template <typename Archive>
    // void serialize(Archive &ar)
    // {
    //     ar(problemBounds_);
    // }

    // //Number of objectives: (that's the max(Eps_P) , projection error)
    // pagmo::vector_double::size_type get_nobj() const {
    //     return 1u; //example in which you have 3 obj
    // }

    // //Overall dimension of the decision vector
    // pagmo::vector_double::size_type get_nx() const {
    //     return problemBounds_[0].size();
    // }

private:

    Eigen::MatrixXd m_sn_set;
    prob_settings m_settings;
    const std::vector< std::vector< double > > problemBounds_;
    int m_Nf;
};


struct SPOD_Adapt_Samp_ {

    // Empty constructor
    SPOD_Adapt_Samp_( ){ }

    SPOD_Adapt_Samp_( std::vector< std::vector< double > > &bounds, Eigen::MatrixXd &sn_set, prob_settings &settings,
                     int &Nf ) : problemBounds_(bounds), m_sn_set(sn_set), m_settings(settings), m_Nf(Nf) {
    }
    // Fitness:
    std::vector<double> fitness(const std::vector<double> &variables) const;
    // Boundaries of the problem
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const;

    //Number of equality and inequality constraints
    pagmo::vector_double::size_type get_nec() const {
        return 0;
    }

    pagmo::vector_double::size_type get_nic() const {
        return 1;
    }

private:

    Eigen::MatrixXd m_sn_set;
    prob_settings m_settings;
    const std::vector< std::vector< double > > problemBounds_;
    int m_Nf;
};



struct SPOD_Adapt_Samp_Int {

    // Empty constructor
    SPOD_Adapt_Samp_Int( ){ }

    SPOD_Adapt_Samp_Int( std::vector< std::vector< double > > &bounds, Eigen::MatrixXd &sn_set, prob_settings &settings,
                      int &Nf ) : problemBounds_(bounds), m_sn_set(sn_set), m_settings(settings), m_Nf(Nf) {
    }

    std::vector<double> fitness(const std::vector<double> &variables) const; // Fitness
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const; // Boundaries of the problem

    //Number of continuous and integer variables
    pagmo::vector_double::size_type get_ncx() const {
        return 0;
    }
    pagmo::vector_double::size_type get_nix() const {
        return problemBounds_[0].size();
    }

private:

    Eigen::MatrixXd m_sn_set;
    prob_settings m_settings;
    const std::vector< std::vector< double > > problemBounds_;
    int m_Nf;
};

struct SPOD_Adapt_Samp_Int_ {

    // Empty constructor
    SPOD_Adapt_Samp_Int_( ){ }

    SPOD_Adapt_Samp_Int_( std::vector< std::vector< double > > &bounds, Eigen::MatrixXd &sn_set, prob_settings &settings,
                         int &Nf ) : problemBounds_(bounds), m_sn_set(sn_set), m_settings(settings), m_Nf(Nf) {
    }

    std::vector<double> fitness(const std::vector<double> &variables) const; // Fitness
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const; // Boundaries of the problem

    //Number of equality and inequality constraints
    pagmo::vector_double::size_type get_nec() const {
        return 0;
    }
    pagmo::vector_double::size_type get_nic() const {
        return 1;
    }
    //Number of continuous and integer variables
    pagmo::vector_double::size_type get_ncx() const {
        return 0;
    }
    pagmo::vector_double::size_type get_nix() const {
        return problemBounds_[0].size();
    }

private:

    Eigen::MatrixXd m_sn_set;
    prob_settings m_settings;
    const std::vector< std::vector< double > > problemBounds_;
    int m_Nf;
};


#endif //OPT_STRUCT_HPP