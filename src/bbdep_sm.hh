/**************************************************************************

   Copyright © 2018 Sergey Poluyan <svpoluyan@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

**************************************************************************/
#ifndef INCLUDED_bbdep_sm_hh
#define INCLUDED_bbdep_sm_hh

#include "bbutils.hh"
#include "dunbrackdata.hh"

#include <core/chemical/AtomType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <numeric/NumericTraits.hh>

namespace pepsgo
{
namespace bbdep
{

struct sm_1d
{
    std::vector<bbdep::Dunbrack_data> lib;
    std::vector<std::vector<int>> libn;
    std::vector<std::pair<double, double>> lib_grid;
    std::vector<std::vector<double>> lib_states;
    std::vector<bbutils::distribution_1d> lib_independent;
    std::vector<std::vector<bbutils::distribution_1d>> lib_cdf_sum_all;
};

struct sm_2d
{
    std::vector<bbdep::Dunbrack_data> lib;
    std::vector<std::vector<int>> libn;
    std::vector<std::pair<double, double>> lib_grid;
    std::vector<std::vector<double>> lib_states_chi1;
    std::vector<bbutils::distribution_1d> lib_independent;
    std::vector<std::vector<bbutils::distribution_1d>> lib_cdf_sum_all;
    std::vector<std::vector<bbutils::distribution_1d>> lib_chi2_depend_chi1;
};

struct sm_3d
{
    std::vector<bbdep::Dunbrack_data> lib;
    std::vector<std::vector<int>> libn;
    std::vector<std::pair<double, double>> lib_grid;
    std::vector<std::vector<double>> lib_states_chi1;
    std::vector<std::vector<std::vector<double>>> lib_states_chi2;
    std::vector<bbutils::distribution_1d> lib_independent;
    std::vector<std::vector<bbutils::distribution_1d>> lib_cdf_sum_all;
    std::vector<std::vector<bbutils::distribution_1d>> lib_chi2_depend_chi1;
    std::vector<std::vector<std::vector<bbutils::distribution_1d>>> lib_chi3_depend_chi12;
};

struct sm_4d
{
    std::vector<bbdep::Dunbrack_data> lib;
    std::vector<std::vector<int>> libn;
    std::vector<std::pair<double, double>> lib_grid;
    std::vector<std::vector<double>> lib_states_chi1;
    std::vector<std::vector<std::vector<double>>> lib_states_chi2;
    std::vector<std::vector<std::vector<std::vector<double>>>> lib_states_chi3;
    std::vector<bbutils::distribution_1d> lib_independent;
    std::vector<std::vector<bbutils::distribution_1d>> lib_cdf_sum_all;
    std::vector<std::vector<bbutils::distribution_1d>> lib_chi2_depend_chi1;
    std::vector<std::vector<std::vector<bbutils::distribution_1d>>> lib_chi3_depend_chi12;
    std::vector<std::vector<std::vector<std::vector<bbutils::distribution_1d>>>> lib_chi4_depend_chi123;
    std::vector<std::vector<std::vector<size_t>>> lib_impossible_conformations;
};

class BBDEP_Dunbrack_sm
{
private:
    std::string path;
    size_t cdf_grid_step;

public:
    std::vector<sm_1d> aa_sm_1d; // 0 ser, 1 val, 2 cys, 3 thr
    std::vector<sm_2d> aa_sm_2d; // 0 trp, 1 his, 2 asn, 3 asp, 4 phe, 5 tyr, 6 ile, 7 leu
    std::vector<sm_3d> aa_sm_3d; // 0 met, 1 glu, 2 gln, 3 pro
    std::vector<sm_4d> aa_sm_4d; // 0 arg, 1 lys

    BBDEP_Dunbrack_sm();
    void set_path(std::string path_to_files);
    void set_step(size_t step);
    void load_data_sm(std::string amino_acids);
    void initialize_all(bool create_cdf_sum, std::string amino_acids, int threads_number);

    void create_cdf_sums(core::chemical::AA amino_acid);

    bbutils::distribution_1d get_chi1_all(std::vector<bbdep::Dunbrack_data> &data) const;
    bbutils::distribution_1d get_chi2_all(std::vector<bbdep::Dunbrack_data> &data) const;
    bbutils::distribution_1d get_chi3_all(std::vector<bbdep::Dunbrack_data> &data) const;
    bbutils::distribution_1d get_chi4_all(std::vector<bbdep::Dunbrack_data> &data) const;
    bbutils::distribution_1d fill_uniformly() const;

    Dunbrack_data get_max(core::chemical::AA amino_acid, double Phi, double Psi);
    size_t get_index_from_phi_psi(const std::vector<std::pair<double, double>> &data, double Phi, double Psi) const;

    void fill_grid_and_states_and_create_cdf_chi1(const std::vector<bbdep::Dunbrack_data> &data,
            std::vector<std::pair<double, double>> &all_chi1_index,
            std::vector<bbutils::distribution_1d> &all_independent,
            std::vector<std::vector<bbutils::distribution_1d>> &all_chi1_value,
            std::vector<std::vector<double>> &states);

    void fill_grid_and_states_and_create_cdf_chi2(const std::vector<bbdep::Dunbrack_data> &data,
            std::vector<std::pair<double, double>> &all_chi1_index,
            std::vector<bbutils::distribution_1d> &all_independent,
            std::vector<std::vector<bbutils::distribution_1d>> &all_chi12_value,
            std::vector<std::vector<bbutils::distribution_1d>> &chi2_depend_chi1,
            std::vector<std::vector<double>> &states);

    void fill_grid_and_states_and_create_cdf_chi3(const std::vector<bbdep::Dunbrack_data> &data,
            std::vector<std::pair<double, double>> &all_chi1_index,
            std::vector<bbutils::distribution_1d> &all_independent,
            std::vector<std::vector<bbutils::distribution_1d>> &all_chi123_value,
            std::vector<std::vector<bbutils::distribution_1d>> &chi2_depend_chi1,
            std::vector<std::vector<std::vector<bbutils::distribution_1d>>> &chi3_depend_chi12,
            std::vector<std::vector<double>> &chi1_states,
            std::vector<std::vector<std::vector<double>>> &chi2_states);

    void fill_grid_and_states_and_create_cdf_chi4(const std::vector<bbdep::Dunbrack_data> &data,
            std::vector<std::pair<double, double>> &all_chi1_index,
            std::vector<bbutils::distribution_1d> &all_independent,
            std::vector<std::vector<bbutils::distribution_1d>> &all_chi1234_value,
            std::vector<std::vector<bbutils::distribution_1d>> &chi2_depend_chi1,
            std::vector<std::vector<std::vector<bbutils::distribution_1d>>> &chi3_depend_chi12,
            std::vector<std::vector<std::vector<std::vector<bbutils::distribution_1d>>>> &chi4_depend_chi123,
            std::vector<std::vector<double>> &chi1_states,
            std::vector<std::vector<std::vector<double>>> &chi2_states,
            std::vector<std::vector<std::vector<std::vector<double>>>> &chi3_states);

    size_t determine_rotamer_state_0_2pi(double degree) const;
    size_t determine_proline_rotamer_state_0_2pi(double degree) const;

    size_t determine_rotamer_state_0_2pi_actual_chi1(size_t index, double degree, core::chemical::AA amino_acid) const;
    size_t determine_rotamer_state_0_2pi_actual_chi2(size_t index, size_t chi1_state, double degree, core::chemical::AA amino_acid) const;
    size_t determine_rotamer_state_0_2pi_actual_chi3(size_t index, size_t chi1_state, size_t chi2_state, double degree, core::chemical::AA amino_acid) const;

    size_t find_index_for_cdf_chi234(core::chemical::AA amino_acid, double Phi, double Psi) const;

    double get_degree_bbdep_from_phi_psi_x01_chi1_dep(size_t index, core::chemical::AA amino_acid, double chi1_positive_degree, double x01) const;
    double get_degree_bbdep_from_phi_psi_x01_chi12_dep(size_t index, core::chemical::AA amino_acid, double chi1_positive_degree, double chi2_positive_degree, double x01) const;
    double get_degree_bbdep_from_phi_psi_x01_chi123_dep(size_t index, core::chemical::AA amino_acid, double chi1_positive_degree, double chi2_positive_degree, double chi3_positive_degree, double x01) const;

    double get_degree_bbdep_from_phi_psi_x01_chi1_dep_actual_states(size_t index, core::chemical::AA amino_acid, double chi1_positive_degree, double x01) const;
    double get_degree_bbdep_from_phi_psi_x01_chi12_dep_actual_states(size_t index, core::chemical::AA amino_acid, double chi1_positive_degree, double chi2_positive_degree, double x01) const;
    double get_degree_bbdep_from_phi_psi_x01_chi123_dep_actual_states(size_t index, core::chemical::AA amino_acid, double chi1_positive_degree, double chi2_positive_degree, double chi3_positive_degree, double x01) const;

    double get_inverse_degree_bbdep_from_phi_psi_x01_chi1_dep(size_t index, core::chemical::AA amino_acid, double chi1_positive_degree, double x01) const;
    double get_inverse_degree_bbdep_from_phi_psi_x01_chi12_dep(size_t index, core::chemical::AA amino_acid, double chi1_positive_degree, double chi2_positive_degree, double x01) const;
    double get_inverse_degree_bbdep_from_phi_psi_x01_chi123_dep(size_t index, core::chemical::AA amino_acid, double chi1_positive_degree, double chi2_positive_degree, double chi3_positive_degree, double x01) const;

    double get_inverse_degree_bbdep_from_phi_psi_x01_chi1_dep_actual_states(size_t index, core::chemical::AA amino_acid, double chi1_positive_degree, double x01) const;
    double get_inverse_degree_bbdep_from_phi_psi_x01_chi12_dep_actual_states(size_t index, core::chemical::AA amino_acid, double chi1_positive_degree, double chi2_positive_degree, double x01) const;
    double get_inverse_degree_bbdep_from_phi_psi_x01_chi123_dep_actual_states(size_t index, core::chemical::AA amino_acid, double chi1_positive_degree, double chi2_positive_degree, double chi3_positive_degree, double x01) const;

    void fill_impossible_conformations(std::vector<bbdep::Dunbrack_data> &data,
                                       std::vector<std::vector<std::vector<size_t>>> &imp_conf);
    bool is_impossible_conformation(std::vector<size_t> conf,
                                    double Phi,
                                    double Psi,
                                    std::vector<std::pair<double, double>> &lys_grid,
                                    std::vector<std::vector<std::vector<size_t>>> &imp_conf);

    bbdep::Dunbrack_data get_first_line(core::chemical::AA amino_acid) const;
    double get_degree_bbind(core::chemical::AA amino_acid, double x01, size_t chinumber) const;
    double get_degree_bbdep_from_phi_psi_x01_chinumber(size_t index, core::chemical::AA amino_acid, double x01, size_t chinumber) const;

};


double pdf_normal_dst(double x, double mu, double sigma);
void plot_chi1_all(pepsgo::bbdep::BBDEP_Dunbrack_sm& bbdep_sm);

}
}

#endif
