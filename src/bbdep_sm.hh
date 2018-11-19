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
    std::string path_to_files;
    size_t cdf_grid_step;

public:
    std::vector<sm_1d> aa_sm_1d; // 0 ser, 1 val, 2 cys, 3 thr
    std::vector<sm_2d> aa_sm_2d; // 0 trp, 1 his, 2 asn, 3 asp, 4 phe, 5 tyr, 6 ile, 7 leu
    std::vector<sm_3d> aa_sm_3d; // 0 met, 1 glu, 2 gln, 3 pro
    std::vector<sm_4d> aa_sm_4d; // 0 arg, 1 lys

    BBDEP_Dunbrack_sm(std::string _path_to_dunbrack_files, size_t _cdf_grid_step);
    void load_data_em();
    void initialize_all(bool create_cdf_sum, std::string amino_acids);
};

}
}

#endif

