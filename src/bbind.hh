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
#ifndef INCLUDED_bbind_hh
#define INCLUDED_bbind_hh

#include <core/chemical/AA.hh>

#include "bbutils.hh"
#include <unordered_map>
#include <any>
#include <boost/multi_array.hpp>

namespace pepsgo
{
namespace bbind
{

// A   G   V   W   H   N   D   F   I   L   C   T   S   M   E   Q   P   Y   R   K
// Ala Gly Val Trp His Asn Asp Phe Ile Leu Cys Thr Ser Met Glu Gln Pro Tyr Arg Lys

// rosetta
// 0 AG
// 1 V
// 2 WHNDFILCTS
// 3 MEQPY
// 4 RK

// bbind top500
// 0 AG ala gly
// 1 SVCT ser val cys thr
// 2 WHNDFYIL trp his asn asp phe tyr ile leu
// 3 MEQP met glu gln pro
// 4 RK arg lys

class BBIND_top
{
public:
    std::string path;
    std::unordered_map<core::chemical::AA, std::any> aa_data;
    std::unordered_map<core::chemical::AA, std::any> aa_dst;
    std::unordered_map<core::chemical::AA, std::vector<std::tuple<double, double, size_t>>> aa_range;

    BBIND_top();
    void set_path(std::string path_to_files);

    void initialize(size_t chi1_step, size_t chi2_step, size_t chi3_step, size_t chi4_step, std::string amino_acids);


    bbutils::distribution_1d make_1d_cdf(core::chemical::AA acid, size_t m);
    bbutils::distribution_2d make_2d_cdf(core::chemical::AA acid, size_t m);
    bbutils::distribution_3d make_3d_cdf(core::chemical::AA acid, size_t m);
    bbutils::distribution_4d make_4d_cdf(core::chemical::AA acid, size_t m);

    void load_data(std::string fname,
                   std::vector<std::vector<double>> &data,
                   std::vector<std::tuple<double, double, size_t>> &range);

    double get_1d(double x, core::chemical::AA acid);

    void load_1d(std::string prefix, std::string acid_name, core::chemical::AA acid);
    void load_2d(std::string prefix, std::string acid_name, core::chemical::AA acid);
    void load_3d(std::string prefix, std::string acid_name, core::chemical::AA acid);
    void load_4d(std::string prefix, std::string acid_name, core::chemical::AA acid);

    void fill_1d(std::string name,
                 std::vector<std::tuple<double, double, size_t>> &range,
                 std::vector<std::vector<double>> &grids,
                 size_t &num_elements,
                 boost::array<int, 1> &grid_sizes,
                 std::vector<double> &f_values);

    void fill_2d(std::string name,
                 std::vector<std::tuple<double, double, size_t>> &range,
                 std::vector<std::vector<double>> &grids,
                 size_t &num_elements,
                 boost::array<int, 2> &grid_sizes,
                 std::vector<double> &f_values);

    void fill_3d(std::string name,
                 std::vector<std::tuple<double, double, size_t>> &range,
                 std::vector<std::vector<double>> &grids,
                 size_t &num_elements,
                 boost::array<int, 3> &grid_sizes,
                 std::vector<double> &f_values);

    void fill_4d(std::string name,
                 std::vector<std::tuple<double, double, size_t>> &range,
                 std::vector<std::vector<double>> &grids,
                 size_t &num_elements,
                 boost::array<int, 4> &grid_sizes,
                 std::vector<double> &f_values);

};


void plot_chi1_all(pepsgo::bbind::BBIND_top& obj);
void plot_chi2(pepsgo::bbind::BBIND_top& obj, core::chemical::AA acid);


}
}

#endif
