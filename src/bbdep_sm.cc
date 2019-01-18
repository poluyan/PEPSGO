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
#include <omp.h>

#include "bbdep_sm.hh"
#include "bbutils.hh"
#include "data_io.hh"

namespace pepsgo
{
namespace bbdep
{

double pdf_normal_dst(double x, double mu, double sigma)
{
    double ps = 2 * sigma * sigma;
    double e = exp(-std::pow(x - mu, 2.0) / ps);
    double C = std::sqrt(2.0 * ps * numeric::NumericTraits<core::Real>::pi());
    return e / C;
}


void plot_chi1_all(pepsgo::bbdep::BBDEP_Dunbrack_sm& bbdep_sm)
{
    auto ser = bbdep_sm.get_chi1_all(bbdep_sm.aa_sm_1d[0].lib);
    auto val = bbdep_sm.get_chi1_all(bbdep_sm.aa_sm_1d[1].lib);
    auto cys = bbdep_sm.get_chi1_all(bbdep_sm.aa_sm_1d[2].lib);
    auto thr = bbdep_sm.get_chi1_all(bbdep_sm.aa_sm_1d[3].lib);

    std::vector<std::vector<double>> to_plot;
    for(size_t i = 0; i != ser.pdf.size(); i++)
    {
        std::vector<double> temp = {ser.grid[i], ser.pdf[i], val.pdf[i], cys.pdf[i], thr.pdf[i]};
        if(temp[0] < 0)
        {
            temp[0] += 360.0;
        }
        to_plot.push_back(temp);
    }
    pepsgo::write_default2d("maps/bbdep/SVCT.dat", to_plot, 4);
}

BBDEP_Dunbrack_sm::BBDEP_Dunbrack_sm() {}


void BBDEP_Dunbrack_sm::set_path(std::string path_to_files)
{
    path = path_to_files;
}
void BBDEP_Dunbrack_sm::set_step(size_t step)
{
    cdf_grid_step = step;
}

void BBDEP_Dunbrack_sm::create_cdf_sums(core::chemical::AA amino_acid)
{
    switch(amino_acid)
    {
        // 0 ser, 1 val, 2 cys, 3 thr
        case core::chemical::aa_ser:
            fill_grid_and_states_and_create_cdf_chi1(aa_sm_1d[0].lib, aa_sm_1d[0].lib_grid,
                    aa_sm_1d[0].lib_independent,
                    aa_sm_1d[0].lib_cdf_sum_all,
                    aa_sm_1d[0].lib_states);
            break;
        case core::chemical::aa_val:
            fill_grid_and_states_and_create_cdf_chi1(aa_sm_1d[1].lib, aa_sm_1d[1].lib_grid,
                    aa_sm_1d[1].lib_independent,
                    aa_sm_1d[1].lib_cdf_sum_all,
                    aa_sm_1d[1].lib_states);
            break;
        case core::chemical::aa_cys:
            fill_grid_and_states_and_create_cdf_chi1(aa_sm_1d[2].lib, aa_sm_1d[2].lib_grid,
                    aa_sm_1d[2].lib_independent,
                    aa_sm_1d[2].lib_cdf_sum_all,
                    aa_sm_1d[2].lib_states);
            break;
        case core::chemical::aa_thr:
            fill_grid_and_states_and_create_cdf_chi1(aa_sm_1d[3].lib, aa_sm_1d[3].lib_grid,
                    aa_sm_1d[3].lib_independent,
                    aa_sm_1d[3].lib_cdf_sum_all,
                    aa_sm_1d[3].lib_states);
            break;


        // 0 trp, 1 his, 2 asn, 3 asp, 4 phe, 5 tyr, 6 ile, 7 leu
        case core::chemical::aa_trp:
            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[0].lib, aa_sm_2d[0].lib_grid,
                    aa_sm_2d[0].lib_independent,
                    aa_sm_2d[0].lib_cdf_sum_all,
                    aa_sm_2d[0].lib_chi2_depend_chi1,
                    aa_sm_2d[0].lib_states_chi1);
            break;
        case core::chemical::aa_his:
            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[1].lib, aa_sm_2d[1].lib_grid,
                    aa_sm_2d[1].lib_independent,
                    aa_sm_2d[1].lib_cdf_sum_all,
                    aa_sm_2d[1].lib_chi2_depend_chi1,
                    aa_sm_2d[1].lib_states_chi1);
            break;
        case core::chemical::aa_asn:
            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[2].lib, aa_sm_2d[2].lib_grid,
                    aa_sm_2d[2].lib_independent,
                    aa_sm_2d[2].lib_cdf_sum_all,
                    aa_sm_2d[2].lib_chi2_depend_chi1,
                    aa_sm_2d[2].lib_states_chi1);
            break;
        case core::chemical::aa_asp:
            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[3].lib, aa_sm_2d[3].lib_grid,
                    aa_sm_2d[3].lib_independent,
                    aa_sm_2d[3].lib_cdf_sum_all,
                    aa_sm_2d[3].lib_chi2_depend_chi1,
                    aa_sm_2d[3].lib_states_chi1);
            break;
        case core::chemical::aa_phe:
            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[4].lib, aa_sm_2d[4].lib_grid,
                    aa_sm_2d[4].lib_independent,
                    aa_sm_2d[4].lib_cdf_sum_all,
                    aa_sm_2d[4].lib_chi2_depend_chi1,
                    aa_sm_2d[4].lib_states_chi1);
            break;
        case core::chemical::aa_tyr:
            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[5].lib, aa_sm_2d[5].lib_grid,
                    aa_sm_2d[5].lib_independent,
                    aa_sm_2d[5].lib_cdf_sum_all,
                    aa_sm_2d[5].lib_chi2_depend_chi1,
                    aa_sm_2d[5].lib_states_chi1);
            break;
        case core::chemical::aa_ile:
            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[6].lib, aa_sm_2d[6].lib_grid,
                    aa_sm_2d[6].lib_independent,
                    aa_sm_2d[6].lib_cdf_sum_all,
                    aa_sm_2d[6].lib_chi2_depend_chi1,
                    aa_sm_2d[6].lib_states_chi1);
            break;
        case core::chemical::aa_leu:
            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[7].lib, aa_sm_2d[7].lib_grid,
                    aa_sm_2d[7].lib_independent,
                    aa_sm_2d[7].lib_cdf_sum_all,
                    aa_sm_2d[7].lib_chi2_depend_chi1,
                    aa_sm_2d[7].lib_states_chi1);
            break;

        // 0 met, 1 glu, 2 gln, 3 pro
        case core::chemical::aa_met:
            fill_grid_and_states_and_create_cdf_chi3(
                aa_sm_3d[0].lib,
                aa_sm_3d[0].lib_grid,
                aa_sm_3d[0].lib_independent,
                aa_sm_3d[0].lib_cdf_sum_all,
                aa_sm_3d[0].lib_chi2_depend_chi1,
                aa_sm_3d[0].lib_chi3_depend_chi12,
                aa_sm_3d[0].lib_states_chi1,
                aa_sm_3d[0].lib_states_chi2);
            break;

        case core::chemical::aa_glu:
            fill_grid_and_states_and_create_cdf_chi3(
                aa_sm_3d[1].lib,
                aa_sm_3d[1].lib_grid,
                aa_sm_3d[1].lib_independent,
                aa_sm_3d[1].lib_cdf_sum_all,
                aa_sm_3d[1].lib_chi2_depend_chi1,
                aa_sm_3d[1].lib_chi3_depend_chi12,
                aa_sm_3d[1].lib_states_chi1,
                aa_sm_3d[1].lib_states_chi2);
            break;

        case core::chemical::aa_gln:
            fill_grid_and_states_and_create_cdf_chi3(
                aa_sm_3d[2].lib,
                aa_sm_3d[2].lib_grid,
                aa_sm_3d[2].lib_independent,
                aa_sm_3d[2].lib_cdf_sum_all,
                aa_sm_3d[2].lib_chi2_depend_chi1,
                aa_sm_3d[2].lib_chi3_depend_chi12,
                aa_sm_3d[2].lib_states_chi1,
                aa_sm_3d[2].lib_states_chi2);
            break;


        case core::chemical::aa_pro:
            fill_grid_and_states_and_create_cdf_chi3(
                aa_sm_3d[3].lib,
                aa_sm_3d[3].lib_grid,
                aa_sm_3d[3].lib_independent,
                aa_sm_3d[3].lib_cdf_sum_all,
                aa_sm_3d[3].lib_chi2_depend_chi1,
                aa_sm_3d[3].lib_chi3_depend_chi12,
                aa_sm_3d[3].lib_states_chi1,
                aa_sm_3d[3].lib_states_chi2);
            break;

        // 0 arg, 1 lys
        case core::chemical::aa_arg:
        {
            fill_grid_and_states_and_create_cdf_chi4(
                aa_sm_4d[0].lib,
                aa_sm_4d[0].lib_grid,
                aa_sm_4d[0].lib_independent,
                aa_sm_4d[0].lib_cdf_sum_all,
                aa_sm_4d[0].lib_chi2_depend_chi1,
                aa_sm_4d[0].lib_chi3_depend_chi12,
                aa_sm_4d[0].lib_chi4_depend_chi123,
                aa_sm_4d[0].lib_states_chi1,
                aa_sm_4d[0].lib_states_chi2,
                aa_sm_4d[0].lib_states_chi3);
            fill_impossible_conformations(aa_sm_4d[0].lib, aa_sm_4d[0].lib_impossible_conformations);
            break;
        }
        case core::chemical::aa_lys:
        {
            fill_grid_and_states_and_create_cdf_chi4(
                aa_sm_4d[1].lib,
                aa_sm_4d[1].lib_grid,
                aa_sm_4d[1].lib_independent,
                aa_sm_4d[1].lib_cdf_sum_all,
                aa_sm_4d[1].lib_chi2_depend_chi1,
                aa_sm_4d[1].lib_chi3_depend_chi12,
                aa_sm_4d[1].lib_chi4_depend_chi123,
                aa_sm_4d[1].lib_states_chi1,
                aa_sm_4d[1].lib_states_chi2,
                aa_sm_4d[1].lib_states_chi3);
            fill_impossible_conformations(aa_sm_4d[1].lib, aa_sm_4d[1].lib_impossible_conformations);
            break;

        }
        default:
            break;
    }
}

// amino acids in 1letter code format
void BBDEP_Dunbrack_sm::initialize_all(bool create_cdf_sum, std::string amino_acids, int threads_number)
{
    amino_acids.erase(std::remove(amino_acids.begin(), amino_acids.end(), 'A'), amino_acids.end());
    amino_acids.erase(std::remove(amino_acids.begin(), amino_acids.end(), 'G'), amino_acids.end());

    load_data_sm(amino_acids);

    if(create_cdf_sum)
    {
        int n = amino_acids.size();
        omp_set_dynamic(0);
        omp_set_num_threads(threads_number);
        #pragma omp parallel for
        for(int i = 0; i < n; ++i)
        {
            create_cdf_sums(core::chemical::aa_from_oneletter_code(amino_acids[i]));
        }      

//        // 0 ser, 1 val, 2 cys, 3 thr
//        if(amino_acids.find("S") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi1(aa_sm_1d[0].lib, aa_sm_1d[0].lib_grid,
//                    aa_sm_1d[0].lib_independent,
//                    aa_sm_1d[0].lib_cdf_sum_all,
//                    aa_sm_1d[0].lib_states);
//        if(amino_acids.find("V") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi1(aa_sm_1d[1].lib, aa_sm_1d[1].lib_grid,
//                    aa_sm_1d[1].lib_independent,
//                    aa_sm_1d[1].lib_cdf_sum_all,
//                    aa_sm_1d[1].lib_states);
//        if(amino_acids.find("C") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi1(aa_sm_1d[2].lib, aa_sm_1d[2].lib_grid,
//                    aa_sm_1d[2].lib_independent,
//                    aa_sm_1d[2].lib_cdf_sum_all,
//                    aa_sm_1d[2].lib_states);
//        if(amino_acids.find("T") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi1(aa_sm_1d[3].lib, aa_sm_1d[3].lib_grid,
//                    aa_sm_1d[3].lib_independent,
//                    aa_sm_1d[3].lib_cdf_sum_all,
//                    aa_sm_1d[3].lib_states);
//
//        // 0 trp, 1 his, 2 asn, 3 asp, 4 phe, 5 tyr, 6 ile, 7 leu
//        if(amino_acids.find("W") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[0].lib, aa_sm_2d[0].lib_grid,
//                    aa_sm_2d[0].lib_independent,
//                    aa_sm_2d[0].lib_cdf_sum_all,
//                    aa_sm_2d[0].lib_chi2_depend_chi1,
//                    aa_sm_2d[0].lib_states_chi1);
//        if(amino_acids.find("H") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[1].lib, aa_sm_2d[1].lib_grid,
//                    aa_sm_2d[1].lib_independent,
//                    aa_sm_2d[1].lib_cdf_sum_all,
//                    aa_sm_2d[1].lib_chi2_depend_chi1,
//                    aa_sm_2d[1].lib_states_chi1);
//        if(amino_acids.find("N") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[2].lib, aa_sm_2d[2].lib_grid,
//                    aa_sm_2d[2].lib_independent,
//                    aa_sm_2d[2].lib_cdf_sum_all,
//                    aa_sm_2d[2].lib_chi2_depend_chi1,
//                    aa_sm_2d[2].lib_states_chi1);
//        if(amino_acids.find("D") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[3].lib, aa_sm_2d[3].lib_grid,
//                    aa_sm_2d[3].lib_independent,
//                    aa_sm_2d[3].lib_cdf_sum_all,
//                    aa_sm_2d[3].lib_chi2_depend_chi1,
//                    aa_sm_2d[3].lib_states_chi1);
//        if(amino_acids.find("F") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[4].lib, aa_sm_2d[4].lib_grid,
//                    aa_sm_2d[4].lib_independent,
//                    aa_sm_2d[4].lib_cdf_sum_all,
//                    aa_sm_2d[4].lib_chi2_depend_chi1,
//                    aa_sm_2d[4].lib_states_chi1);
//        if(amino_acids.find("Y") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[5].lib, aa_sm_2d[5].lib_grid,
//                    aa_sm_2d[5].lib_independent,
//                    aa_sm_2d[5].lib_cdf_sum_all,
//                    aa_sm_2d[5].lib_chi2_depend_chi1,
//                    aa_sm_2d[5].lib_states_chi1);
//        if(amino_acids.find("I") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[6].lib, aa_sm_2d[6].lib_grid,
//                    aa_sm_2d[6].lib_independent,
//                    aa_sm_2d[6].lib_cdf_sum_all,
//                    aa_sm_2d[6].lib_chi2_depend_chi1,
//                    aa_sm_2d[6].lib_states_chi1);
//        if(amino_acids.find("L") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi2(aa_sm_2d[7].lib, aa_sm_2d[7].lib_grid,
//                    aa_sm_2d[7].lib_independent,
//                    aa_sm_2d[7].lib_cdf_sum_all,
//                    aa_sm_2d[7].lib_chi2_depend_chi1,
//                    aa_sm_2d[7].lib_states_chi1);
//
//        // 0 met, 1 glu, 2 gln, 3 pro
//        if(amino_acids.find("M") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi3(
//                aa_sm_3d[0].lib,
//                aa_sm_3d[0].lib_grid,
//                aa_sm_3d[0].lib_independent,
//                aa_sm_3d[0].lib_cdf_sum_all,
//                aa_sm_3d[0].lib_chi2_depend_chi1,
//                aa_sm_3d[0].lib_chi3_depend_chi12,
//                aa_sm_3d[0].lib_states_chi1,
//                aa_sm_3d[0].lib_states_chi2);
//        if(amino_acids.find("E") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi3(
//                aa_sm_3d[1].lib,
//                aa_sm_3d[1].lib_grid,
//                aa_sm_3d[1].lib_independent,
//                aa_sm_3d[1].lib_cdf_sum_all,
//                aa_sm_3d[1].lib_chi2_depend_chi1,
//                aa_sm_3d[1].lib_chi3_depend_chi12,
//                aa_sm_3d[1].lib_states_chi1,
//                aa_sm_3d[1].lib_states_chi2);
//        if(amino_acids.find("Q") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi3(
//                aa_sm_3d[2].lib,
//                aa_sm_3d[2].lib_grid,
//                aa_sm_3d[2].lib_independent,
//                aa_sm_3d[2].lib_cdf_sum_all,
//                aa_sm_3d[2].lib_chi2_depend_chi1,
//                aa_sm_3d[2].lib_chi3_depend_chi12,
//                aa_sm_3d[2].lib_states_chi1,
//                aa_sm_3d[2].lib_states_chi2);
//        if(amino_acids.find("P") != std::string::npos)
//            fill_grid_and_states_and_create_cdf_chi3(
//                aa_sm_3d[3].lib,
//                aa_sm_3d[3].lib_grid,
//                aa_sm_3d[3].lib_independent,
//                aa_sm_3d[3].lib_cdf_sum_all,
//                aa_sm_3d[3].lib_chi2_depend_chi1,
//                aa_sm_3d[3].lib_chi3_depend_chi12,
//                aa_sm_3d[3].lib_states_chi1,
//                aa_sm_3d[3].lib_states_chi2);
//
//        // 0 arg, 1 lys
//        if(amino_acids.find("R") != std::string::npos)
//        {
//            fill_grid_and_states_and_create_cdf_chi4(
//                aa_sm_4d[0].lib,
//                aa_sm_4d[0].lib_grid,
//                aa_sm_4d[0].lib_independent,
//                aa_sm_4d[0].lib_cdf_sum_all,
//                aa_sm_4d[0].lib_chi2_depend_chi1,
//                aa_sm_4d[0].lib_chi3_depend_chi12,
//                aa_sm_4d[0].lib_chi4_depend_chi123,
//                aa_sm_4d[0].lib_states_chi1,
//                aa_sm_4d[0].lib_states_chi2,
//                aa_sm_4d[0].lib_states_chi3);
//            fill_impossible_conformations(aa_sm_4d[0].lib, aa_sm_4d[0].lib_impossible_conformations);
//        }
//        if(amino_acids.find("K") != std::string::npos)
//        {
//            fill_grid_and_states_and_create_cdf_chi4(
//                aa_sm_4d[1].lib,
//                aa_sm_4d[1].lib_grid,
//                aa_sm_4d[1].lib_independent,
//                aa_sm_4d[1].lib_cdf_sum_all,
//                aa_sm_4d[1].lib_chi2_depend_chi1,
//                aa_sm_4d[1].lib_chi3_depend_chi12,
//                aa_sm_4d[1].lib_chi4_depend_chi123,
//                aa_sm_4d[1].lib_states_chi1,
//                aa_sm_4d[1].lib_states_chi2,
//                aa_sm_4d[1].lib_states_chi3);
//            fill_impossible_conformations(aa_sm_4d[1].lib, aa_sm_4d[1].lib_impossible_conformations);
//        }
    }
}


void BBDEP_Dunbrack_sm::load_data_sm(std::string amino_acids)
{
    aa_sm_1d.resize(4); // 0 ser, 1 val, 2 cys, 3 thr
    aa_sm_2d.resize(8); // 0 trp, 1 his, 2 asn, 3 asp, 4 phe, 5 tyr, 6 ile, 7 leu
    aa_sm_3d.resize(4); // 0 met, 1 glu, 2 gln, 3 pro
    aa_sm_4d.resize(2); // 0 arg, 1 lys

    // 0 ser, 1 val, 2 cys, 3 thr
    if(amino_acids.find("S") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "ser.bbdep.rotamers.lib", aa_sm_1d[0].lib, aa_sm_1d[0].libn);
    if(amino_acids.find("V") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "val.bbdep.rotamers.lib", aa_sm_1d[1].lib, aa_sm_1d[1].libn);
    if(amino_acids.find("C") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "cys.bbdep.rotamers.lib", aa_sm_1d[2].lib, aa_sm_1d[2].libn);
    if(amino_acids.find("T") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "thr.bbdep.rotamers.lib", aa_sm_1d[3].lib, aa_sm_1d[3].libn);

    // 0 trp, 1 his, 2 asn, 3 asp, 4 phe, 5 tyr, 6 ile, 7 leu
    if(amino_acids.find("W") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "trp.bbdep.rotamers.lib", aa_sm_2d[0].lib, aa_sm_2d[0].libn);
    if(amino_acids.find("H") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "his.bbdep.rotamers.lib", aa_sm_2d[1].lib, aa_sm_2d[1].libn);
    if(amino_acids.find("N") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "asn.bbdep.rotamers.lib", aa_sm_2d[2].lib, aa_sm_2d[2].libn);
    if(amino_acids.find("D") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "asp.bbdep.rotamers.lib", aa_sm_2d[3].lib, aa_sm_2d[3].libn);
    if(amino_acids.find("F") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "phe.bbdep.rotamers.lib", aa_sm_2d[4].lib, aa_sm_2d[4].libn);
    if(amino_acids.find("Y") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "tyr.bbdep.rotamers.lib", aa_sm_2d[5].lib, aa_sm_2d[5].libn);
    if(amino_acids.find("I") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "ile.bbdep.rotamers.lib", aa_sm_2d[6].lib, aa_sm_2d[6].libn);
    if(amino_acids.find("L") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "leu.bbdep.rotamers.lib", aa_sm_2d[7].lib, aa_sm_2d[7].libn);

    // 0 met, 1 glu, 2 gln, 3 pro
    if(amino_acids.find("M") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "met.bbdep.rotamers.lib", aa_sm_3d[0].lib, aa_sm_3d[0].libn);
    if(amino_acids.find("E") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "glu.bbdep.rotamers.lib", aa_sm_3d[1].lib, aa_sm_3d[1].libn);
    if(amino_acids.find("Q") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "gln.bbdep.rotamers.lib", aa_sm_3d[2].lib, aa_sm_3d[2].libn);
    if(amino_acids.find("P") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "pro.bbdep.rotamers.lib", aa_sm_3d[3].lib, aa_sm_3d[3].libn);

    // 0 arg, 1 lys
    if(amino_acids.find("R") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "arg.bbdep.rotamers.lib", aa_sm_4d[0].lib, aa_sm_4d[0].libn);
    if(amino_acids.find("K") != std::string::npos)
        pepsgo::bbdep::load_data_sm(path + "lys.bbdep.rotamers.lib", aa_sm_4d[1].lib, aa_sm_4d[1].libn);
}

void BBDEP_Dunbrack_sm::fill_grid_and_states_and_create_cdf_chi1(const std::vector<bbdep::Dunbrack_data> &data,
        std::vector<std::pair<double, double>> &all_chi1_index,
        std::vector<bbutils::distribution_1d> &all_independent,
        std::vector<std::vector<bbutils::distribution_1d>> &all_chi1_value,
        std::vector<std::vector<double>> &states)
{
    std::vector<Dunbrack_data> complete_independent(data.begin(), data.end());
    all_independent.push_back(get_chi1_all(complete_independent));

    for(auto i = data.begin(); i != data.end();)
    {
        bbdep::Dunbrack_data search_data;
        search_data.Phi = i->Phi;
        search_data.Psi = i->Psi;

        auto p = std::equal_range(
                     i, data.end(), search_data, [](const bbdep::Dunbrack_data &lhs, const bbdep::Dunbrack_data &rhs) -> bool
        {
            return lhs.Phi < rhs.Phi;
        });

        auto q = std::equal_range(p.first, p.second, search_data,
                                  [](const bbdep::Dunbrack_data &lhs, const bbdep::Dunbrack_data &rhs) -> bool
        {
            return lhs.Psi < rhs.Psi;
        });

        std::vector<bbdep::Dunbrack_data> temp(q.first, q.second);
        all_chi1_index.push_back(std::make_pair(search_data.Phi, search_data.Psi));

        std::vector<bbutils::distribution_1d> temp2;
        temp2.push_back(get_chi1_all(temp));
        all_chi1_value.push_back(temp2);

        //states
        std::vector<double> temp3;
        for(size_t j = 1; j != 4; j++)
        {
            std::vector<bbdep::Dunbrack_data> matches;
            std::copy_if(q.first, q.second, std::back_inserter(matches),
                         [&](const bbdep::Dunbrack_data &value) -> bool {if(j == value.r1) return true; else return false; });

            if(matches.size() != 1)
                std::cout << "wrong fill states 1d" << std::endl;

            temp3.push_back(matches.front().chi1Val < 0 ? matches.front().chi1Val + 360.0 : matches.front().chi1Val);
        }
        states.push_back(temp3);

        i = q.second;
    }
}

void BBDEP_Dunbrack_sm::fill_grid_and_states_and_create_cdf_chi2(const std::vector<bbdep::Dunbrack_data> &data,
        std::vector<std::pair<double, double>> &all_chi1_index,
        std::vector<bbutils::distribution_1d> &all_independent,
        std::vector<std::vector<bbutils::distribution_1d>> &all_chi12_value,
        std::vector<std::vector<bbutils::distribution_1d>> &chi2_depend_chi1,
        std::vector<std::vector<double>> &states)
{
    chi2_depend_chi1.resize(3);

    std::vector<Dunbrack_data> complete_independent(data.begin(), data.end());
    all_independent.push_back(get_chi1_all(complete_independent));
    all_independent.push_back(get_chi2_all(complete_independent));

    for(auto i = data.begin(); i != data.end();)
    {
        bbdep::Dunbrack_data search_data;
        search_data.Phi = i->Phi;
        search_data.Psi = i->Psi;

        auto p = std::equal_range(
                     i, data.end(), search_data, [](const bbdep::Dunbrack_data &lhs, const bbdep::Dunbrack_data &rhs) -> bool
        {
            return lhs.Phi < rhs.Phi;
        });

        auto q = std::equal_range(p.first, p.second, search_data,
                                  [](const bbdep::Dunbrack_data &lhs, const bbdep::Dunbrack_data &rhs) -> bool
        {
            return lhs.Psi < rhs.Psi;
        });

        // range_bounds.push_back(std::make_pair(std::distance(data.begin(),q.first),
        // std::distance(data.begin(),q.second)));

        std::vector<bbdep::Dunbrack_data> temp(q.first, q.second);
        all_chi1_index.push_back(std::make_pair(search_data.Phi, search_data.Psi));

        std::vector<bbutils::distribution_1d> temp2;
        temp2.push_back(get_chi1_all(temp));
        temp2.push_back(get_chi2_all(temp));
        all_chi12_value.push_back(temp2);

        for(size_t j = 1; j != 4; j++)
        {
            std::vector<bbdep::Dunbrack_data> temp3(q.first, q.second);
            temp3.erase(
                std::remove_if(temp3.begin(), temp3.end(), [j](bbdep::Dunbrack_data n)
            {
                return n.r1 != j;
            }),
            temp3.end());

            if(temp3.size())
            {
                chi2_depend_chi1[j - 1].push_back(get_chi2_all(temp3));
            }
            else
            {
                chi2_depend_chi1[j - 1].push_back(fill_uniformly());
            }
        }

        //states
        std::vector<double> temp4;
        for(size_t j = 1; j != 4; j++)
        {
            std::vector<bbdep::Dunbrack_data> matches;
            std::copy_if(q.first, q.second, std::back_inserter(matches), [&](const bbdep::Dunbrack_data &value) -> bool
            {
                if(j == value.r1) return true;
                else return false; });

            double state = 0.0, sum = 0.0;
            for(auto &h : matches)
            {
                state += (h.chi1Val < 0 ? h.chi1Val + 360.0 : h.chi1Val) * h.Probabil;
                sum += h.Probabil;
            }
            state = state / sum;
            temp4.push_back(state);
        }
        states.push_back(temp4);

        i = q.second;
    }
}
void BBDEP_Dunbrack_sm::fill_grid_and_states_and_create_cdf_chi3(const std::vector<bbdep::Dunbrack_data> &data,
        std::vector<std::pair<double, double>> &all_chi1_index,
        std::vector<bbutils::distribution_1d> &all_independent,
        std::vector<std::vector<bbutils::distribution_1d>> &all_chi123_value,
        std::vector<std::vector<bbutils::distribution_1d>> &chi2_depend_chi1,
        std::vector<std::vector<std::vector<bbutils::distribution_1d>>> &chi3_depend_chi12,
        std::vector<std::vector<double>> &chi1_states,
        std::vector<std::vector<std::vector<double>>> &chi2_states)
{
    chi2_depend_chi1.resize(3);

    chi3_depend_chi12.resize(3);
    for(size_t i = 0; i != chi3_depend_chi12.size(); i++)
    {
        chi3_depend_chi12[i].resize(3);
    }

    std::vector<Dunbrack_data> complete_independent(data.begin(), data.end());
    all_independent.push_back(get_chi1_all(complete_independent));
    all_independent.push_back(get_chi2_all(complete_independent));
    all_independent.push_back(get_chi3_all(complete_independent));

    for(auto i = data.begin(); i != data.end();)
    {
        bbdep::Dunbrack_data search_data;
        search_data.Phi = i->Phi;
        search_data.Psi = i->Psi;

        auto p = std::equal_range(
                     i, data.end(), search_data, [](const bbdep::Dunbrack_data &lhs, const bbdep::Dunbrack_data &rhs) -> bool
        {
            return lhs.Phi < rhs.Phi;
        });

        auto q = std::equal_range(p.first, p.second, search_data,
                                  [](const bbdep::Dunbrack_data &lhs, const bbdep::Dunbrack_data &rhs) -> bool
        {
            return lhs.Psi < rhs.Psi;
        });

        std::vector<bbdep::Dunbrack_data> temp(q.first, q.second);
        all_chi1_index.push_back(std::make_pair(search_data.Phi, search_data.Psi));

        std::vector<bbutils::distribution_1d> temp2;
        temp2.push_back(get_chi1_all(temp));
        temp2.push_back(get_chi2_all(temp));
        temp2.push_back(get_chi3_all(temp));
        all_chi123_value.push_back(temp2);

        for(size_t j = 1; j != 4; j++)
        {
            std::vector<bbdep::Dunbrack_data> temp3(q.first, q.second);
            temp3.erase(
                std::remove_if(temp3.begin(), temp3.end(), [j](bbdep::Dunbrack_data n)
            {
                return n.r1 != j;
            }),
            temp3.end());

            chi2_depend_chi1[j - 1].push_back(get_chi2_all(temp3));
        }

        for(size_t j = 1; j != 4; j++)
        {
            for(size_t k = 1; k != 4; k++)
            {
                std::vector<bbdep::Dunbrack_data> temp3(q.first, q.second);
                temp3.erase(std::remove_if(temp3.begin(), temp3.end(),
                                           [j, k](bbdep::Dunbrack_data n)
                {
                    return (n.r1 != j || n.r2 != k);
                }),
                temp3.end());

                if(temp3.size())
                {
                    chi3_depend_chi12[j - 1][k - 1].push_back(get_chi3_all(temp3));
                }
                else
                {
                    chi3_depend_chi12[j - 1][k - 1].push_back(fill_uniformly());
                }
            }
        }

        //states chi1
        std::vector<double> temp4;
        for(size_t j = 1; j != 4; j++)
        {
            std::vector<bbdep::Dunbrack_data> matches;
            std::copy_if(q.first, q.second, std::back_inserter(matches), [&](const bbdep::Dunbrack_data &value) -> bool
            {
                if(j == value.r1) return true;
                else return false; });

            if(!matches.size())
                continue;

            double state = 0.0, sum = 0.0;
            for(auto &h : matches)
            {
                state += (h.chi1Val < 0 ? h.chi1Val + 360.0 : h.chi1Val) * h.Probabil;
                sum += h.Probabil;
            }
            temp4.push_back(state / sum);
        }
        chi1_states.push_back(temp4);

        //states chi2
        std::vector<std::vector<double>> temp5;
        for(size_t j = 1; j != 4; j++)
        {
            std::vector<double> temp6;
            for(size_t k = 1; k != 4; k++)
            {
                std::vector<bbdep::Dunbrack_data> matches;
                std::copy_if(q.first, q.second, std::back_inserter(matches), [&](const bbdep::Dunbrack_data &value) -> bool
                {
                    if(j == value.r1 && k == value.r2) return true;
                    else return false; });

                if(!matches.size())
                    continue;

                double state = 0.0, sum = 0.0;
                for(auto &h : matches)
                {
                    state += (h.chi2Val < 0 ? h.chi2Val + 360.0 : h.chi2Val) * h.Probabil;
                    sum += h.Probabil;
                }
                temp6.push_back(state / sum);
            }
            temp5.push_back(temp6);
        }
        chi2_states.push_back(temp5);

        i = q.second;
    }
}
void BBDEP_Dunbrack_sm::fill_grid_and_states_and_create_cdf_chi4(const std::vector<bbdep::Dunbrack_data> &data,
        std::vector<std::pair<double, double>> &all_chi1_index,
        std::vector<bbutils::distribution_1d> &all_independent,
        std::vector<std::vector<bbutils::distribution_1d>> &all_chi1234_value,
        std::vector<std::vector<bbutils::distribution_1d>> &chi2_depend_chi1,
        std::vector<std::vector<std::vector<bbutils::distribution_1d>>> &chi3_depend_chi12,
        std::vector<std::vector<std::vector<std::vector<bbutils::distribution_1d>>>> &chi4_depend_chi123,
        std::vector<std::vector<double>> &chi1_states,
        std::vector<std::vector<std::vector<double>>> &chi2_states,
        std::vector<std::vector<std::vector<std::vector<double>>>> &chi3_states)
{

    chi2_depend_chi1.resize(3);

    chi3_depend_chi12.resize(3);
    for(size_t i = 0; i != chi3_depend_chi12.size(); i++)
    {
        chi3_depend_chi12[i].resize(3);
    }

    chi4_depend_chi123.resize(3);
    for(size_t i = 0; i != chi4_depend_chi123.size(); i++)
    {
        chi4_depend_chi123[i].resize(3);
        for(size_t j = 0; j != chi4_depend_chi123[i].size(); j++)
        {
            chi4_depend_chi123[i][j].resize(3);
        }
    }

    std::vector<Dunbrack_data> complete_independent(data.begin(), data.end());
    all_independent.push_back(get_chi1_all(complete_independent));
    all_independent.push_back(get_chi2_all(complete_independent));
    all_independent.push_back(get_chi3_all(complete_independent));
    all_independent.push_back(get_chi4_all(complete_independent));

    for(auto i = data.begin(); i != data.end();)
    {
        bbdep::Dunbrack_data search_data;
        search_data.Phi = i->Phi;
        search_data.Psi = i->Psi;

        auto p = std::equal_range(
                     i, data.end(), search_data, [](const bbdep::Dunbrack_data &lhs, const bbdep::Dunbrack_data &rhs) -> bool
        {
            return lhs.Phi < rhs.Phi;
        });

        auto q = std::equal_range(p.first, p.second, search_data,
                                  [](const bbdep::Dunbrack_data &lhs, const bbdep::Dunbrack_data &rhs) -> bool
        {
            return lhs.Psi < rhs.Psi;
        });

        std::vector<bbdep::Dunbrack_data> temp(q.first, q.second);
        all_chi1_index.push_back(std::make_pair(search_data.Phi, search_data.Psi));

        std::vector<bbutils::distribution_1d> temp2;
        temp2.push_back(get_chi1_all(temp));
        temp2.push_back(get_chi2_all(temp));
        temp2.push_back(get_chi3_all(temp));
        temp2.push_back(get_chi4_all(temp));
        all_chi1234_value.push_back(temp2);

        for(size_t j = 1; j != 4; j++)
        {
            std::vector<bbdep::Dunbrack_data> temp3(q.first, q.second);
            temp3.erase(
                std::remove_if(temp3.begin(), temp3.end(), [j](bbdep::Dunbrack_data n)
            {
                return n.r1 != j;
            }),
            temp3.end());

            if(temp3.size())
            {
                chi2_depend_chi1[j - 1].push_back(get_chi2_all(temp3));
            }
            else
            {
                chi2_depend_chi1[j - 1].push_back(fill_uniformly());
            }
        }

        for(size_t j = 1; j != 4; j++)
        {
            for(size_t k = 1; k != 4; k++)
            {
                std::vector<bbdep::Dunbrack_data> temp3(q.first, q.second);
                temp3.erase(std::remove_if(temp3.begin(), temp3.end(),
                                           [j, k](bbdep::Dunbrack_data n)
                {
                    return (n.r1 != j || n.r2 != k);
                }),
                temp3.end());

                if(temp3.size())
                {
                    chi3_depend_chi12[j - 1][k - 1].push_back(get_chi3_all(temp3));
                }
                else
                {
                    chi3_depend_chi12[j - 1][k - 1].push_back(fill_uniformly());
                }
            }
        }

        for(size_t j = 1; j != 4; j++)
        {
            for(size_t k = 1; k != 4; k++)
            {
                for(size_t h = 1; h != 4; h++)
                {
                    std::vector<bbdep::Dunbrack_data> temp3(q.first, q.second);
                    temp3.erase(
                        std::remove_if(temp3.begin(), temp3.end(),
                                       [j, k, h](bbdep::Dunbrack_data n)
                    {
                        return (n.r1 != j || n.r2 != k || n.r3 != h);
                    }),
                    temp3.end());

                    if(temp3.size())
                    {
                        chi4_depend_chi123[j - 1][k - 1][h - 1].push_back(get_chi4_all(temp3));
                    }
                    else
                    {
                        chi4_depend_chi123[j - 1][k - 1][h - 1].push_back(fill_uniformly());
                    }

                    // if(j == 1 && k == 1 && h == 3 && search_data.Phi < 101.0 && search_data.Phi > 99.0 &&
                    // search_data.Psi < -179.0)
                    //{
                    //    for(auto a : temp3)
                    //        std::cout << std::fixed << a << std::endl;
                    //    std::cin.get();
                    //}
                }
            }
        }

        //states chi1
        std::vector<double> temp4;
        for(size_t j = 1; j != 4; j++)
        {
            std::vector<bbdep::Dunbrack_data> matches;
            std::copy_if(q.first, q.second, std::back_inserter(matches), [&](const bbdep::Dunbrack_data &value) -> bool
            {
                if(j == value.r1) return true;
                else return false; });

            if(!matches.size())
                continue;

            double state = 0.0, sum = 0.0;
            for(auto &h : matches)
            {
                state += (h.chi1Val < 0 ? h.chi1Val + 360.0 : h.chi1Val) * h.Probabil;
                sum += h.Probabil;
            }
            temp4.push_back(state / sum);
        }
        chi1_states.push_back(temp4);

        //states chi2
        std::vector<std::vector<double>> temp5;
        for(size_t j = 1; j != 4; j++)
        {
            std::vector<double> temp6;
            for(size_t k = 1; k != 4; k++)
            {
                std::vector<bbdep::Dunbrack_data> matches;
                std::copy_if(q.first, q.second, std::back_inserter(matches), [&](const bbdep::Dunbrack_data &value) -> bool
                {
                    if(j == value.r1 && k == value.r2) return true;
                    else return false; });

                if(!matches.size())
                    continue;

                double state = 0.0, sum = 0.0;
                for(auto &h : matches)
                {
                    state += (h.chi2Val < 0 ? h.chi2Val + 360.0 : h.chi2Val) * h.Probabil;
                    sum += h.Probabil;
                }
                temp6.push_back(state / sum);
            }
            temp5.push_back(temp6);
        }
        chi2_states.push_back(temp5);

        // states chi3
        std::vector<std::vector<std::vector<double>>> temp7;
        for(size_t j = 1; j != 4; j++)
        {
            std::vector<std::vector<double>> temp8;
            for(size_t k = 1; k != 4; k++)
            {
                std::vector<double> temp9;
                for(size_t h = 1; h != 4; h++)
                {
                    std::vector<bbdep::Dunbrack_data> matches;
                    std::copy_if(q.first, q.second, std::back_inserter(matches), [&](const bbdep::Dunbrack_data &value) -> bool
                    {
                        if(j == value.r1 && k == value.r2 && h == value.r3) return true;
                        else return false; });

                    if(!matches.size())
                    {
                        switch(temp9.size())
                        {
                            case 0:
                                temp9.push_back(60.0);
                                break;
                            case 1:
                                temp9.push_back(180.0);
                                break;
                            case 2:
                                temp9.push_back(300.0);
                                break;
                        }
                        continue;
                    }

                    double state = 0.0, sum = 0.0;
                    for(auto &h : matches)
                    {
                        state += (h.chi3Val < 0 ? h.chi3Val + 360.0 : h.chi3Val) * h.Probabil;
                        sum += h.Probabil;
                    }
                    temp9.push_back(state / sum);
                }
                temp8.push_back(temp9);
            }
            temp7.push_back(temp8);
        }
        chi3_states.push_back(temp7);

        i = q.second;
    }
}

void BBDEP_Dunbrack_sm::fill_impossible_conformations(std::vector<bbdep::Dunbrack_data> &data,
        std::vector<std::vector<std::vector<size_t>>> &imp_conf)
{
    std::vector<int> rot_pos = { 1, 2, 3 };
    std::vector<std::vector<size_t>> pp;
    for(size_t i = 1, m = 1; i < 4; i++)
    {
        for(size_t j = 1; j < 4; j++)
        {
            for(size_t k = 1; k < 4; k++)
            {
                for(size_t h = 1; h < 4; h++, m++)
                {
                    // std::cout << i << '\t' << j << '\t' << k << '\t' << h << "\t\t" << m << std::endl;
                    std::vector<size_t> temp;
                    temp.push_back(i);
                    temp.push_back(j);
                    temp.push_back(k);
                    temp.push_back(h);
                    temp.push_back(m);
                    pp.push_back(temp);
                }
            }
        }
    }

    for(auto i = data.begin(); i != data.end();)
    {
        bbdep::Dunbrack_data search_data;
        search_data.Phi = i->Phi;
        search_data.Psi = i->Psi;

        auto p = std::equal_range(
                     i, data.end(), search_data, [](const bbdep::Dunbrack_data &lhs, const bbdep::Dunbrack_data &rhs) -> bool
        {
            return lhs.Phi < rhs.Phi;
        });

        auto q = std::equal_range(p.first, p.second, search_data,
                                  [](const bbdep::Dunbrack_data &lhs, const bbdep::Dunbrack_data &rhs) -> bool
        {
            return lhs.Psi < rhs.Psi;
        });

        std::vector<std::vector<size_t>> temp;
        for(size_t i = 0; i != pp.size(); i++)
        {
            auto a = std::find_if(q.first, q.second, [&](const bbdep::Dunbrack_data &value) -> bool
            {
                if(pp[i][0] == value.r1 && pp[i][1] == value.r2 && pp[i][2] == value.r3 && pp[i][3] == value.r4)
                    return true;
                else
                    return false;
            });
            if(a == q.second)
            {
                // std::cout << pp[i][0] << '\t' << pp[i][1] << '\t' << pp[i][2] << '\t' << pp[i][3] << std::endl;
                std::vector<size_t> rot_string;
                rot_string.push_back(pp[i][0]);
                rot_string.push_back(pp[i][1]);
                rot_string.push_back(pp[i][2]);
                rot_string.push_back(pp[i][3]);
                temp.push_back(rot_string);
            }
        }
        imp_conf.push_back(temp);

        // std::cout << std::endl;

        i = q.second;
    }
    // std::cout << imp_conf.size() << '\t' << imp_conf[0].size() << '\t' << imp_conf[0][1].size() << std::endl;
}


bbutils::distribution_1d BBDEP_Dunbrack_sm::get_chi1_all(std::vector<bbdep::Dunbrack_data> &data) const
{
    std::deque<double> x, pdf;

    size_t step = cdf_grid_step;
    double es = 360;
    double start = -180;

    for(size_t i = 0; i != step; i++)
    {
        x.push_back(start + i * es / step);
        double m = 0;
        for(auto j = data.begin(); j != data.end(); ++j)
            m += j->Probabil * pdf_normal_dst(x.back(), j->chi1Val, j->chi1Sig);
        for(auto j = data.begin(); j != data.end(); ++j)
            m += j->Probabil * pdf_normal_dst(x.back(), j->chi1Val + 360, j->chi1Sig);
        for(auto j = data.begin(); j != data.end(); ++j)
            m += j->Probabil * pdf_normal_dst(x.back(), j->chi1Val - 360, j->chi1Sig);
        // m = data.begin()->Probabil * pdf_normal_dst(x.back(), data.begin()->chi1Val, data.begin()->chi1Sig);
        pdf.push_back(m);
    }

    double S = std::accumulate(pdf.begin(), pdf.end(), 0.0);
    for(auto &a : pdf)
    {
        a /= S;
    }

    std::deque<double> cdf(pdf.size());
    std::partial_sum(pdf.begin(), pdf.end(), cdf.begin(), std::plus<double>());
    cdf.push_front(0);
    x.push_back(180.0);

    bbutils::distribution_1d result;
    result.pdf = pdf;
    result.cdf = cdf;
    result.grid = x;

    return result;
}
bbutils::distribution_1d BBDEP_Dunbrack_sm::get_chi2_all(std::vector<bbdep::Dunbrack_data> &data) const
{
    std::deque<double> x, pdf;

    size_t step = cdf_grid_step;
    double es = 360;
    double start = -180;

    for(size_t i = 0; i != step; i++)
    {
        x.push_back(start + i * es / step);
        double m = 0;
        for(auto j = data.begin(); j != data.end(); ++j)
            m += j->Probabil * pdf_normal_dst(x.back(), j->chi2Val, j->chi2Sig);
        for(auto j = data.begin(); j != data.end(); ++j)
            m += j->Probabil * pdf_normal_dst(x.back(), j->chi2Val + 360, j->chi2Sig);
        for(auto j = data.begin(); j != data.end(); ++j)
            m += j->Probabil * pdf_normal_dst(x.back(), j->chi2Val - 360, j->chi2Sig);
        // m = data.begin()->Probabil * pdf_normal_dst(x.back(), data.begin()->chi2Val, data.begin()->chi2Sig);
        pdf.push_back(m);
    }

    double S = std::accumulate(pdf.begin(), pdf.end(), 0.0);
    for(auto &a : pdf)
    {
        a /= S;
    }

    std::deque<double> cdf(pdf.size());
    std::partial_sum(pdf.begin(), pdf.end(), cdf.begin(), std::plus<double>());
    cdf.push_front(0);
    x.push_back(180.0);

    bbutils::distribution_1d result;
    result.cdf = cdf;
    result.grid = x;

    return result;
}
bbutils::distribution_1d BBDEP_Dunbrack_sm::get_chi3_all(std::vector<bbdep::Dunbrack_data> &data) const
{
    std::deque<double> x, pdf;

    size_t step = cdf_grid_step;
    double es = 360;
    double start = -180;

    for(size_t i = 0; i != step; i++)
    {
        x.push_back(start + i * es / step);
        double m = 0;
        for(auto j = data.begin(); j != data.end(); ++j)
            m += j->Probabil * pdf_normal_dst(x.back(), j->chi3Val, j->chi3Sig);
        for(auto j = data.begin(); j != data.end(); ++j)
            m += j->Probabil * pdf_normal_dst(x.back(), j->chi3Val + 360, j->chi3Sig);
        for(auto j = data.begin(); j != data.end(); ++j)
            m += j->Probabil * pdf_normal_dst(x.back(), j->chi3Val - 360, j->chi3Sig);
        // m = data.begin()->Probabil * pdf_normal_dst(x.back(), data.begin()->chi3Val, data.begin()->chi3Sig);
        pdf.push_back(m);
    }

    double S = std::accumulate(pdf.begin(), pdf.end(), 0.0);
    for(auto &a : pdf)
    {
        a /= S;
    }

    std::deque<double> cdf(pdf.size());
    std::partial_sum(pdf.begin(), pdf.end(), cdf.begin(), std::plus<double>());
    cdf.push_front(0);
    x.push_back(180.0);

    bbutils::distribution_1d result;
    result.cdf = cdf;
    result.grid = x;

    return result;
}
bbutils::distribution_1d BBDEP_Dunbrack_sm::get_chi4_all(std::vector<bbdep::Dunbrack_data> &data) const
{
    std::deque<double> x, pdf;

    size_t step = cdf_grid_step;
    double es = 360;
    double start = -180;

    for(size_t i = 0; i != step; i++)
    {
        x.push_back(start + i * es / step);
        double m = 0;
        for(auto j = data.begin(); j != data.end(); ++j)
            m += j->Probabil * pdf_normal_dst(x.back(), j->chi4Val, j->chi4Sig);
        for(auto j = data.begin(); j != data.end(); ++j)
            m += j->Probabil * pdf_normal_dst(x.back(), j->chi4Val + 360, j->chi4Sig);
        for(auto j = data.begin(); j != data.end(); ++j)
            m += j->Probabil * pdf_normal_dst(x.back(), j->chi4Val - 360, j->chi4Sig);
        // m = data.begin()->Probabil * pdf_normal_dst(x.back(), data.begin()->chi4Val, data.begin()->chi4Sig);
        pdf.push_back(m);
    }

    double S = std::accumulate(pdf.begin(), pdf.end(), 0.0);
    for(auto &a : pdf)
    {
        a /= S;
    }

    std::deque<double> cdf(pdf.size());
    std::partial_sum(pdf.begin(), pdf.end(), cdf.begin(), std::plus<double>());
    cdf.push_front(0);
    x.push_back(180.0);

    bbutils::distribution_1d result;
    result.cdf = cdf;
    result.grid = x;

    return result;
}

bbutils::distribution_1d BBDEP_Dunbrack_sm::fill_uniformly() const
{
    std::deque<double> x, pdf;

    size_t step = cdf_grid_step;
    double es = 360;
    double start = -180;

    for(size_t i = 0; i != step; i++)
    {
        x.push_back(start + i * es / step);
        double m = 1.0;
        // double m = 1.0 * pdf_normal_dst(x.back(), 0.0, 0.000000001);

        pdf.push_back(m);
    }

    double S = std::accumulate(pdf.begin(), pdf.end(), 0.0);
    for(auto &a : pdf)
    {
        a /= S;
    }

    std::deque<double> cdf(pdf.size());
    std::partial_sum(pdf.begin(), pdf.end(), cdf.begin(), std::plus<double>());
    cdf.push_front(0);
    x.push_back(180.0);

    bbutils::distribution_1d result;
    result.cdf = cdf;
    result.grid = x;

    return result;
}

bool BBDEP_Dunbrack_sm::is_impossible_conformation(std::vector<size_t> conf,
        double Phi,
        double Psi,
        std::vector<std::pair<double, double>> &lys_grid,
        std::vector<std::vector<std::vector<size_t>>> &imp_conf)
{
    size_t i = get_index_from_phi_psi(lys_grid, Phi, Psi);

    auto a = std::find_if(imp_conf[i].begin(), imp_conf[i].end(), [&](const std::vector<size_t> &value) -> bool
    {
        if(conf[0] == value[0] && conf[1] == value[1] && conf[2] == value[2] && conf[3] == value[3])
            return true;
        else
            return false;
        return true;
    });
    return (a == imp_conf[i].end()) ? false : true;
}


Dunbrack_data BBDEP_Dunbrack_sm::get_max(core::chemical::AA amino_acid, double Phi, double Psi)
{
    Dunbrack_data result;
    switch(amino_acid)
    {
        // 0 ser, 1 val, 2 cys, 3 thr
        case core::chemical::aa_ser:
            result = get_max_prob_object(aa_sm_1d[0].lib, aa_sm_1d[0].libn, Phi, Psi);
            break;
        case core::chemical::aa_val:
            result = get_max_prob_object(aa_sm_1d[1].lib, aa_sm_1d[1].libn, Phi, Psi);
            break;
        case core::chemical::aa_cys:
            result = get_max_prob_object(aa_sm_1d[2].lib, aa_sm_1d[2].libn, Phi, Psi);
            break;
        case core::chemical::aa_thr:
            result = get_max_prob_object(aa_sm_1d[3].lib, aa_sm_1d[3].libn, Phi, Psi);
            break;

        // 0 trp, 1 his, 2 asn, 3 asp, 4 phe, 5 tyr, 6 ile, 7 leu
        case core::chemical::aa_trp:
            result = get_max_prob_object(aa_sm_2d[0].lib, aa_sm_2d[0].libn, Phi, Psi);
            break;
        case core::chemical::aa_his:
            result = get_max_prob_object(aa_sm_2d[1].lib, aa_sm_2d[1].libn, Phi, Psi);
            break;
        case core::chemical::aa_asn:
            result = get_max_prob_object(aa_sm_2d[2].lib, aa_sm_2d[2].libn, Phi, Psi);
            break;
        case core::chemical::aa_asp:
            result = get_max_prob_object(aa_sm_2d[3].lib, aa_sm_2d[3].libn, Phi, Psi);
            break;
        case core::chemical::aa_phe:
            result = get_max_prob_object(aa_sm_2d[4].lib, aa_sm_2d[4].libn, Phi, Psi);
            break;
        case core::chemical::aa_tyr:
            result = get_max_prob_object(aa_sm_2d[5].lib, aa_sm_2d[5].libn, Phi, Psi);
            break;
        case core::chemical::aa_ile:
            result = get_max_prob_object(aa_sm_2d[6].lib, aa_sm_2d[6].libn, Phi, Psi);
            break;
        case core::chemical::aa_leu:
            result = get_max_prob_object(aa_sm_2d[7].lib, aa_sm_2d[7].libn, Phi, Psi);
            break;

        // 0 met, 1 glu, 2 gln, 3 pro
        case core::chemical::aa_met:
            result = get_max_prob_object(aa_sm_3d[0].lib, aa_sm_3d[0].libn, Phi, Psi);
            break;
        case core::chemical::aa_glu:
            result = get_max_prob_object(aa_sm_3d[1].lib, aa_sm_3d[1].libn, Phi, Psi);
            break;
        case core::chemical::aa_gln:
            result = get_max_prob_object(aa_sm_3d[2].lib, aa_sm_3d[2].libn, Phi, Psi);
            break;
        case core::chemical::aa_pro:
            result = get_max_prob_object(aa_sm_3d[3].lib, aa_sm_3d[3].libn, Phi, Psi);
            break;

        // 0 arg, 1 lys
        case core::chemical::aa_arg:
            result = get_max_prob_object(aa_sm_4d[0].lib, aa_sm_4d[0].libn, Phi, Psi);
            break;
        case core::chemical::aa_lys:
            result = get_max_prob_object(aa_sm_4d[1].lib, aa_sm_4d[1].libn, Phi, Psi);
            break;

        default:
            break;
    }

    return result;
}

size_t BBDEP_Dunbrack_sm::get_index_from_phi_psi(const std::vector<std::pair<double, double>> &data, double Phi, double Psi) const
{
    std::pair<double, double> search_data;
    search_data.first = 10 * std::round(Phi / 10.0);
    search_data.second = 10 * std::round(Psi / 10.0);

    auto p = std::equal_range(data.begin(), data.end(), search_data,
                              [](const std::pair<double, double> &lhs, const std::pair<double, double> &rhs) -> bool
    {
        return lhs.first < rhs.first;
    });

    auto q = std::equal_range(p.first, p.second, search_data,
                              [](const std::pair<double, double> &lhs, const std::pair<double, double> &rhs) -> bool
    {
        return lhs.second < rhs.second;
    });

    return std::distance(data.begin(), q.first);
}

size_t BBDEP_Dunbrack_sm::determine_rotamer_state_0_2pi(double degree) const
{
    std::vector<double> rotameric_states = { 60, 180, 300 };
    auto i = std::min_element(rotameric_states.begin(), rotameric_states.end(),
                              [=](double x, double y)
    {
        return std::abs(x - degree) < std::abs(y - degree);
    });

    return std::distance(rotameric_states.begin(), i);
}

size_t BBDEP_Dunbrack_sm::determine_proline_rotamer_state_0_2pi(double degree) const
{
    std::vector<double> rotameric_states = { 27.3, 334.9 };
    auto i = std::min_element(rotameric_states.begin(), rotameric_states.end(),
                              [=](double x, double y)
    {
        return std::abs(x - degree) < std::abs(y - degree);
    });

    return std::distance(rotameric_states.begin(), i);
}

size_t BBDEP_Dunbrack_sm::determine_rotamer_state_0_2pi_actual_chi1(size_t index, double degree, core::chemical::AA amino_acid) const
{
    std::vector<double> rotameric_states;
    switch(amino_acid)
    {
        // 0 ser, 1 val, 2 cys, 3 thr
        case core::chemical::aa_ser:
            rotameric_states = aa_sm_1d[0].lib_states[index];
            break;
        case core::chemical::aa_val:
            rotameric_states = aa_sm_1d[1].lib_states[index];
            break;
        case core::chemical::aa_cys:
            rotameric_states = aa_sm_1d[2].lib_states[index];
            break;
        case core::chemical::aa_thr:
            rotameric_states = aa_sm_1d[3].lib_states[index];
            break;

        // 0 trp, 1 his, 2 asn, 3 asp, 4 phe, 5 tyr, 6 ile, 7 leu
        case core::chemical::aa_trp:
            rotameric_states = aa_sm_2d[0].lib_states_chi1[index];
            break;
        case core::chemical::aa_his:
            rotameric_states = aa_sm_2d[1].lib_states_chi1[index];
            break;
        case core::chemical::aa_asn:
            rotameric_states = aa_sm_2d[2].lib_states_chi1[index];
            break;
        case core::chemical::aa_asp:
            rotameric_states = aa_sm_2d[3].lib_states_chi1[index];
            break;
        case core::chemical::aa_phe:
            rotameric_states = aa_sm_2d[4].lib_states_chi1[index];
            break;
        case core::chemical::aa_tyr:
            rotameric_states = aa_sm_2d[5].lib_states_chi1[index];
            break;
        case core::chemical::aa_ile:
            rotameric_states = aa_sm_2d[6].lib_states_chi1[index];
            break;
        case core::chemical::aa_leu:
            rotameric_states = aa_sm_2d[7].lib_states_chi1[index];
            break;

        // 0 met, 1 glu, 2 gln, 3 pro
        case core::chemical::aa_met:
            rotameric_states = aa_sm_3d[0].lib_states_chi1[index];
            break;
        case core::chemical::aa_glu:
            rotameric_states = aa_sm_3d[1].lib_states_chi1[index];
            break;
        case core::chemical::aa_gln:
            rotameric_states = aa_sm_3d[2].lib_states_chi1[index];
            break;
        case core::chemical::aa_pro:
            rotameric_states = aa_sm_3d[3].lib_states_chi1[index];
            break;

        // 0 arg, 1 lys
        case core::chemical::aa_arg:
            rotameric_states = aa_sm_4d[0].lib_states_chi1[index];
            break;
        case core::chemical::aa_lys:
            rotameric_states = aa_sm_4d[1].lib_states_chi1[index];
            break;
        default:
            break;
    }
    auto i = std::min_element(rotameric_states.begin(), rotameric_states.end(),
                              [=](double x, double y)
    {
        return std::abs(x - degree) < std::abs(y - degree);
    });
    return std::distance(rotameric_states.begin(), i);
}

size_t BBDEP_Dunbrack_sm::determine_rotamer_state_0_2pi_actual_chi2(size_t index, size_t chi1_state, double degree, core::chemical::AA amino_acid) const
{
    std::vector<double> rotameric_states;
    switch(amino_acid)
    {
        // 0 met, 1 glu, 2 gln, 3 pro
        case core::chemical::aa_met:
            rotameric_states = aa_sm_3d[0].lib_states_chi2[index][chi1_state];
            break;
        case core::chemical::aa_glu:
            rotameric_states = aa_sm_3d[1].lib_states_chi2[index][chi1_state];
            break;
        case core::chemical::aa_gln:
            rotameric_states = aa_sm_3d[2].lib_states_chi2[index][chi1_state];
            break;
        case core::chemical::aa_pro:
            rotameric_states = aa_sm_3d[3].lib_states_chi2[index][chi1_state];
            break;

        // 0 arg, 1 lys
        case core::chemical::aa_arg:
            rotameric_states = aa_sm_4d[0].lib_states_chi2[index][chi1_state];
            break;
        case core::chemical::aa_lys:
            rotameric_states = aa_sm_4d[1].lib_states_chi2[index][chi1_state];
            break;
        default:
            break;
    }
    auto i = std::min_element(rotameric_states.begin(), rotameric_states.end(),
                              [=](double x, double y)
    {
        return std::abs(x - degree) < std::abs(y - degree);
    });
    return std::distance(rotameric_states.begin(), i);
}

size_t BBDEP_Dunbrack_sm::determine_rotamer_state_0_2pi_actual_chi3(size_t index, size_t chi1_state, size_t chi2_state, double degree, core::chemical::AA amino_acid) const
{
    std::vector<double> rotameric_states;
    switch(amino_acid)
    {
        // 0 arg, 1 lys
        case core::chemical::aa_arg:
            rotameric_states = aa_sm_4d[0].lib_states_chi3[index][chi1_state][chi2_state];
            break;
        case core::chemical::aa_lys:
            rotameric_states = aa_sm_4d[1].lib_states_chi3[index][chi1_state][chi2_state];
            break;
        default:
            break;
    }
    auto i = std::min_element(rotameric_states.begin(), rotameric_states.end(),
                              [=](double x, double y)
    {
        return std::abs(x - degree) < std::abs(y - degree);
    });
    return std::distance(rotameric_states.begin(), i);
}

size_t BBDEP_Dunbrack_sm::find_index_for_cdf_chi234(core::chemical::AA amino_acid, double Phi, double Psi) const
{
    size_t index = 0;

    switch(amino_acid)
    {
        // 0 ser, 1 val, 2 cys, 3 thr
        case core::chemical::aa_ser:
            index = get_index_from_phi_psi(aa_sm_1d[0].lib_grid, Phi, Psi);
            break;
        case core::chemical::aa_val:
            index = get_index_from_phi_psi(aa_sm_1d[1].lib_grid, Phi, Psi);
            break;
        case core::chemical::aa_cys:
            index = get_index_from_phi_psi(aa_sm_1d[2].lib_grid, Phi, Psi);
            break;
        case core::chemical::aa_thr:
            index = get_index_from_phi_psi(aa_sm_1d[3].lib_grid, Phi, Psi);
            break;

        // 0 trp, 1 his, 2 asn, 3 asp, 4 phe, 5 tyr, 6 ile, 7 leu
        case core::chemical::aa_trp:
            index = get_index_from_phi_psi(aa_sm_2d[0].lib_grid, Phi, Psi);
            break;
        case core::chemical::aa_his:
            index = get_index_from_phi_psi(aa_sm_2d[1].lib_grid, Phi, Psi);
            break;
        case core::chemical::aa_asn:
            index = get_index_from_phi_psi(aa_sm_2d[2].lib_grid, Phi, Psi);
            break;
        case core::chemical::aa_asp:
            index = get_index_from_phi_psi(aa_sm_2d[3].lib_grid, Phi, Psi);
            break;
        case core::chemical::aa_phe:
            index = get_index_from_phi_psi(aa_sm_2d[4].lib_grid, Phi, Psi);
            break;
        case core::chemical::aa_tyr:
            index = get_index_from_phi_psi(aa_sm_2d[5].lib_grid, Phi, Psi);
            break;
        case core::chemical::aa_ile:
            index = get_index_from_phi_psi(aa_sm_2d[6].lib_grid, Phi, Psi);
            break;
        case core::chemical::aa_leu:
            index = get_index_from_phi_psi(aa_sm_2d[7].lib_grid, Phi, Psi);
            break;

        // 0 met, 1 glu, 2 gln, 3 pro
        case core::chemical::aa_met:
            index = get_index_from_phi_psi(aa_sm_3d[0].lib_grid, Phi, Psi);
            break;
        case core::chemical::aa_glu:
            index = get_index_from_phi_psi(aa_sm_3d[1].lib_grid, Phi, Psi);
            break;
        case core::chemical::aa_gln:
            index = get_index_from_phi_psi(aa_sm_3d[2].lib_grid, Phi, Psi);
            break;
        case core::chemical::aa_pro:
            index = get_index_from_phi_psi(aa_sm_3d[3].lib_grid, Phi, Psi);
            break;

        // 0 arg, 1 lys
        case core::chemical::aa_arg:
            index = get_index_from_phi_psi(aa_sm_4d[0].lib_grid, Phi, Psi);
            break;
        case core::chemical::aa_lys:
            index = get_index_from_phi_psi(aa_sm_4d[1].lib_grid, Phi, Psi);
            break;

        default:
            break;
    }
    return index;
}

double BBDEP_Dunbrack_sm::get_degree_bbdep_from_phi_psi_x01_chi1_dep(size_t index,
        core::chemical::AA amino_acid,
        double chi1_positive_degree,
        double x01) const
{
    double result = 0;

    size_t state = (amino_acid == core::chemical::aa_pro) ? determine_proline_rotamer_state_0_2pi(chi1_positive_degree) : determine_rotamer_state_0_2pi(chi1_positive_degree);

    switch(amino_acid)
    {
        // 0 trp, 1 his, 2 asn, 3 asp, 4 phe, 5 tyr, 6 ile, 7 leu
        case core::chemical::aa_trp:
            result = bbutils::get_1d_from_dst(aa_sm_2d[0].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_his:
            result = bbutils::get_1d_from_dst(aa_sm_2d[1].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_asn:
            result = bbutils::get_1d_from_dst(aa_sm_2d[2].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_asp:
            result = bbutils::get_1d_from_dst(aa_sm_2d[3].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_phe:
            result = bbutils::get_1d_from_dst(aa_sm_2d[4].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_tyr:
            result = bbutils::get_1d_from_dst(aa_sm_2d[5].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_ile:
            result = bbutils::get_1d_from_dst(aa_sm_2d[6].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_leu:
            result = bbutils::get_1d_from_dst(aa_sm_2d[7].lib_chi2_depend_chi1[state][index], x01);
            break;

        // 0 met, 1 glu, 2 gln, 3 pro
        case core::chemical::aa_met:
            result = bbutils::get_1d_from_dst(aa_sm_3d[0].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_glu:
            result = bbutils::get_1d_from_dst(aa_sm_3d[1].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_gln:
            result = bbutils::get_1d_from_dst(aa_sm_3d[2].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_pro:
            result = bbutils::get_1d_from_dst(aa_sm_3d[3].lib_chi2_depend_chi1[state][index], x01);
            break;

        // 0 arg, 1 lys
        case core::chemical::aa_arg:
            result = bbutils::get_1d_from_dst(aa_sm_4d[0].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_lys:
            result = bbutils::get_1d_from_dst(aa_sm_4d[1].lib_chi2_depend_chi1[state][index], x01);
            break;

        default:
            break;
    }
    return result;
}

double BBDEP_Dunbrack_sm::get_degree_bbdep_from_phi_psi_x01_chi12_dep(size_t index,
        core::chemical::AA amino_acid,
        double chi1_positive_degree,
        double chi2_positive_degree,
        double x01) const
{
    double result = 0;

    size_t state1 = (amino_acid == core::chemical::aa_pro) ? determine_proline_rotamer_state_0_2pi(chi1_positive_degree) : determine_rotamer_state_0_2pi(chi1_positive_degree);
    size_t state2 = (amino_acid == core::chemical::aa_pro) ? 0 : determine_rotamer_state_0_2pi(chi2_positive_degree);

    switch(amino_acid)
    {
        case core::chemical::aa_met:
            result = bbutils::get_1d_from_dst(aa_sm_3d[0].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_glu:
            result = bbutils::get_1d_from_dst(aa_sm_3d[1].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_gln:
            result = bbutils::get_1d_from_dst(aa_sm_3d[2].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_pro:
            result = bbutils::get_1d_from_dst(aa_sm_3d[3].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;

        case core::chemical::aa_arg:
            result = bbutils::get_1d_from_dst(aa_sm_4d[0].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_lys:
            result = bbutils::get_1d_from_dst(aa_sm_4d[1].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;

        default:
            break;
    }
    return result;
}

double BBDEP_Dunbrack_sm::get_degree_bbdep_from_phi_psi_x01_chi123_dep(size_t index,
        core::chemical::AA amino_acid,
        double chi1_positive_degree,
        double chi2_positive_degree,
        double chi3_positive_degree,
        double x01) const
{
    double result = 0;

    size_t state1 = determine_rotamer_state_0_2pi(chi1_positive_degree);
    size_t state2 = determine_rotamer_state_0_2pi(chi2_positive_degree);
    size_t state3 = determine_rotamer_state_0_2pi(chi3_positive_degree);

    switch(amino_acid)
    {
        case core::chemical::aa_arg:
            result = bbutils::get_1d_from_dst(aa_sm_4d[0].lib_chi4_depend_chi123[state1][state2][state3][index], x01);
            break;
        case core::chemical::aa_lys:
            result = bbutils::get_1d_from_dst(aa_sm_4d[1].lib_chi4_depend_chi123[state1][state2][state3][index], x01);
            break;
        default:
            break;
    }
    return result;
}

double BBDEP_Dunbrack_sm::get_degree_bbdep_from_phi_psi_x01_chi1_dep_actual_states(size_t index,
        core::chemical::AA amino_acid,
        double chi1_positive_degree,
        double x01) const
{
    double result = 0;

    size_t state = determine_rotamer_state_0_2pi_actual_chi1(index, chi1_positive_degree, amino_acid);

    switch(amino_acid)
    {
        case core::chemical::aa_trp:
            result = bbutils::get_1d_from_dst(aa_sm_2d[0].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_his:
            result = bbutils::get_1d_from_dst(aa_sm_2d[1].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_asn:
            result = bbutils::get_1d_from_dst(aa_sm_2d[2].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_asp:
            result = bbutils::get_1d_from_dst(aa_sm_2d[3].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_phe:
            result = bbutils::get_1d_from_dst(aa_sm_2d[4].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_tyr:
            result = bbutils::get_1d_from_dst(aa_sm_2d[5].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_ile:
            result = bbutils::get_1d_from_dst(aa_sm_2d[6].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_leu:
            result = bbutils::get_1d_from_dst(aa_sm_2d[7].lib_chi2_depend_chi1[state][index], x01);
            break;

        case core::chemical::aa_met:
            result = bbutils::get_1d_from_dst(aa_sm_3d[0].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_glu:
            result = bbutils::get_1d_from_dst(aa_sm_3d[1].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_gln:
            result = bbutils::get_1d_from_dst(aa_sm_3d[2].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_pro:
            result = bbutils::get_1d_from_dst(aa_sm_3d[3].lib_chi2_depend_chi1[state][index], x01);
            break;

        case core::chemical::aa_arg:
            result = bbutils::get_1d_from_dst(aa_sm_4d[0].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_lys:
            result = bbutils::get_1d_from_dst(aa_sm_4d[1].lib_chi2_depend_chi1[state][index], x01);
            break;
        default:
            break;
    }
    return result;
}

double BBDEP_Dunbrack_sm::get_degree_bbdep_from_phi_psi_x01_chi12_dep_actual_states(size_t index,
        core::chemical::AA amino_acid,
        double chi1_positive_degree,
        double chi2_positive_degree,
        double x01) const
{
    double result = 0;

    size_t state1 = determine_rotamer_state_0_2pi_actual_chi1(index, chi1_positive_degree, amino_acid);
    size_t state2 = determine_rotamer_state_0_2pi_actual_chi2(index, state1, chi2_positive_degree, amino_acid);

    switch(amino_acid)
    {
        case core::chemical::aa_met:
            result = bbutils::get_1d_from_dst(aa_sm_3d[0].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_glu:
            result = bbutils::get_1d_from_dst(aa_sm_3d[1].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_gln:
            result = bbutils::get_1d_from_dst(aa_sm_3d[2].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_pro:
            result = bbutils::get_1d_from_dst(aa_sm_3d[3].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;

        case core::chemical::aa_arg:
            result = bbutils::get_1d_from_dst(aa_sm_4d[0].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_lys:
            result = bbutils::get_1d_from_dst(aa_sm_4d[1].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        default:
            break;
    }
    return result;
}

double BBDEP_Dunbrack_sm::get_degree_bbdep_from_phi_psi_x01_chi123_dep_actual_states(size_t index,
        core::chemical::AA amino_acid,
        double chi1_positive_degree,
        double chi2_positive_degree,
        double chi3_positive_degree,
        double x01) const
{
    double result = 0;

    size_t state1 = determine_rotamer_state_0_2pi_actual_chi1(index, chi1_positive_degree, amino_acid);
    size_t state2 = determine_rotamer_state_0_2pi_actual_chi2(index, state1, chi2_positive_degree, amino_acid);
    size_t state3 = determine_rotamer_state_0_2pi_actual_chi3(index, state1, state2, chi3_positive_degree, amino_acid);

    switch(amino_acid)
    {
        case core::chemical::aa_arg:
            result = bbutils::get_1d_from_dst(aa_sm_4d[0].lib_chi4_depend_chi123[state1][state2][state3][index], x01);
            break;
        case core::chemical::aa_lys:
            result = bbutils::get_1d_from_dst(aa_sm_4d[1].lib_chi4_depend_chi123[state1][state2][state3][index], x01);
            break;
        default:
            break;
    }
    return result;
}


double BBDEP_Dunbrack_sm::get_inverse_degree_bbdep_from_phi_psi_x01_chi1_dep(size_t index,
        core::chemical::AA amino_acid,
        double chi1_positive_degree,
        double x01) const
{
    double result = 0;

    size_t state = determine_rotamer_state_0_2pi(chi1_positive_degree);

    switch(amino_acid)
    {
        case core::chemical::aa_trp:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[0].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_his:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[1].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_asn:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[2].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_asp:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[3].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_phe:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[4].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_tyr:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[5].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_ile:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[6].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_leu:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[7].lib_chi2_depend_chi1[state][index], x01);
            break;

        case core::chemical::aa_met:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[0].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_glu:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[1].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_gln:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[2].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_pro:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[3].lib_chi2_depend_chi1[state][index], x01);
            break;

        case core::chemical::aa_arg:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_4d[0].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_lys:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_4d[1].lib_chi2_depend_chi1[state][index], x01);
            break;
        default:
            break;
    }
    return result;
}

double BBDEP_Dunbrack_sm::get_inverse_degree_bbdep_from_phi_psi_x01_chi12_dep(size_t index,
        core::chemical::AA amino_acid,
        double chi1_positive_degree,
        double chi2_positive_degree,
        double x01) const
{
    double result = 0;

    size_t state1 = determine_rotamer_state_0_2pi(chi1_positive_degree);
    size_t state2 = determine_rotamer_state_0_2pi(chi2_positive_degree);

    switch(amino_acid)
    {

        case core::chemical::aa_met:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[0].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_glu:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[1].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_gln:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[2].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_pro:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[3].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;

        case core::chemical::aa_arg:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_4d[0].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_lys:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_4d[1].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        default:
            break;
    }
    return result;
}

double BBDEP_Dunbrack_sm::get_inverse_degree_bbdep_from_phi_psi_x01_chi123_dep(size_t index,
        core::chemical::AA amino_acid,
        double chi1_positive_degree,
        double chi2_positive_degree,
        double chi3_positive_degree,
        double x01) const
{
    double result = 0;

    size_t state1 = determine_rotamer_state_0_2pi(chi1_positive_degree);
    size_t state2 = determine_rotamer_state_0_2pi(chi2_positive_degree);
    size_t state3 = determine_rotamer_state_0_2pi(chi3_positive_degree);

    switch(amino_acid)
    {
        case core::chemical::aa_arg:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_4d[0].lib_chi4_depend_chi123[state1][state2][state3][index], x01);
            break;
        case core::chemical::aa_lys:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_4d[1].lib_chi4_depend_chi123[state1][state2][state3][index], x01);
            break;
        default:
            break;
    }
    return result;
}

double BBDEP_Dunbrack_sm::get_inverse_degree_bbdep_from_phi_psi_x01_chi1_dep_actual_states(size_t index,
        core::chemical::AA amino_acid,
        double chi1_positive_degree,
        double x01) const
{
    double result = 0;

    size_t state = determine_rotamer_state_0_2pi(chi1_positive_degree);

    switch(amino_acid)
    {
        case core::chemical::aa_trp:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[0].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_his:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[1].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_asn:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[2].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_asp:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[3].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_phe:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[4].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_tyr:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[5].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_ile:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[6].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_leu:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_2d[7].lib_chi2_depend_chi1[state][index], x01);
            break;

        case core::chemical::aa_met:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[0].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_glu:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[1].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_gln:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[2].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_pro:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[3].lib_chi2_depend_chi1[state][index], x01);
            break;

        case core::chemical::aa_arg:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_4d[0].lib_chi2_depend_chi1[state][index], x01);
            break;
        case core::chemical::aa_lys:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_4d[1].lib_chi2_depend_chi1[state][index], x01);
            break;
        default:
            break;
    }
    return result;
}

double BBDEP_Dunbrack_sm::get_inverse_degree_bbdep_from_phi_psi_x01_chi12_dep_actual_states(size_t index,
        core::chemical::AA amino_acid,
        double chi1_positive_degree,
        double chi2_positive_degree,
        double x01) const
{
    double result = 0;

    size_t state1 = determine_rotamer_state_0_2pi(chi1_positive_degree);
    size_t state2 = determine_rotamer_state_0_2pi(chi2_positive_degree);

    switch(amino_acid)
    {
        case core::chemical::aa_met:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[0].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_glu:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[1].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_gln:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[2].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_pro:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_3d[3].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;

        case core::chemical::aa_arg:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_4d[0].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        case core::chemical::aa_lys:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_4d[1].lib_chi3_depend_chi12[state1][state2][index], x01);
            break;
        default:
            break;
    }
    return result;
}

double BBDEP_Dunbrack_sm::get_inverse_degree_bbdep_from_phi_psi_x01_chi123_dep_actual_states(size_t index,
        core::chemical::AA amino_acid,
        double chi1_positive_degree,
        double chi2_positive_degree,
        double chi3_positive_degree,
        double x01) const
{
    double result = 0;

    size_t state1 = determine_rotamer_state_0_2pi(chi1_positive_degree);
    size_t state2 = determine_rotamer_state_0_2pi(chi2_positive_degree);
    size_t state3 = determine_rotamer_state_0_2pi(chi3_positive_degree);

    switch(amino_acid)
    {
        case core::chemical::aa_arg:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_4d[0].lib_chi4_depend_chi123[state1][state2][state3][index], x01);
            break;
        case core::chemical::aa_lys:
            result = bbutils::get_inverse_1d_from_dst(aa_sm_4d[1].lib_chi4_depend_chi123[state1][state2][state3][index], x01);
            break;
        default:
            break;
    }
    return result;
}


bbdep::Dunbrack_data BBDEP_Dunbrack_sm::get_first_line(core::chemical::AA amino_acid) const
{
    pepsgo::bbdep::Dunbrack_data result;
    switch(amino_acid)
    {
        case core::chemical::aa_ser:
            result = aa_sm_1d.front().lib.front();
            break;
        case core::chemical::aa_val:
            result = aa_sm_1d[1].lib.front();
            break;
        case core::chemical::aa_cys:
            result = aa_sm_1d[2].lib.front();
            break;
        case core::chemical::aa_thr:
            result = aa_sm_1d.back().lib.front();
            break;

        case core::chemical::aa_trp:
            result = aa_sm_2d.front().lib.front();
            break;
        case core::chemical::aa_his:
            result = aa_sm_2d[1].lib.front();
            break;
        case core::chemical::aa_asn:
            result = aa_sm_2d[2].lib.front();
            break;
        case core::chemical::aa_asp:
            result = aa_sm_2d[3].lib.front();
            break;
        case core::chemical::aa_phe:
            result = aa_sm_2d[4].lib.front();
            break;
        case core::chemical::aa_tyr:
            result = aa_sm_2d[5].lib.front();
            break;
        case core::chemical::aa_ile:
            result = aa_sm_2d[6].lib.front();
            break;
        case core::chemical::aa_leu:
            result = aa_sm_2d.back().lib.front();
            break;

        case core::chemical::aa_met:
            result = aa_sm_3d.front().lib.front();
            break;
        case core::chemical::aa_glu:
            result = aa_sm_3d[1].lib.front();
            break;
        case core::chemical::aa_gln:
            result = aa_sm_3d[2].lib.front();
            break;
        case core::chemical::aa_pro:
            result = aa_sm_3d.back().lib.front();
            break;

        case core::chemical::aa_arg:
            result = aa_sm_4d.front().lib.front();
            break;
        case core::chemical::aa_lys:
            result = aa_sm_4d.back().lib.front();
            break;

        default:
            break;
    }

    return result;
}

double BBDEP_Dunbrack_sm::get_degree_bbind(core::chemical::AA amino_acid, double x01, size_t chinumber) const
{
    double result = 0;

    switch(amino_acid)
    {
        case core::chemical::aa_ser:
            result = bbutils::get_1d_from_dst(aa_sm_1d.front().lib_independent[chinumber], x01);
            break;
        case core::chemical::aa_val:
            result = bbutils::get_1d_from_dst(aa_sm_1d[1].lib_independent[chinumber], x01);
            break;
        case core::chemical::aa_cys:
            result = bbutils::get_1d_from_dst(aa_sm_1d[2].lib_independent[chinumber], x01);
            break;
        case core::chemical::aa_thr:
            result = bbutils::get_1d_from_dst(aa_sm_1d.back().lib_independent[chinumber], x01);
            break;

        case core::chemical::aa_trp:
            result = bbutils::get_1d_from_dst(aa_sm_2d.front().lib_independent[chinumber], x01);
            break;
        case core::chemical::aa_his:
            result = bbutils::get_1d_from_dst(aa_sm_2d[1].lib_independent[chinumber], x01);
            break;
        case core::chemical::aa_asn:
            result = bbutils::get_1d_from_dst(aa_sm_2d[2].lib_independent[chinumber], x01);
            break;
        case core::chemical::aa_asp:
            result = bbutils::get_1d_from_dst(aa_sm_2d[3].lib_independent[chinumber], x01);
            break;
        case core::chemical::aa_phe:
            result = bbutils::get_1d_from_dst(aa_sm_2d[4].lib_independent[chinumber], x01);
            break;
        case core::chemical::aa_tyr:
            result = bbutils::get_1d_from_dst(aa_sm_2d[5].lib_independent[chinumber], x01);
            break;
        case core::chemical::aa_ile:
            result = bbutils::get_1d_from_dst(aa_sm_2d[6].lib_independent[chinumber], x01);
            break;
        case core::chemical::aa_leu:
            result = bbutils::get_1d_from_dst(aa_sm_2d.back().lib_independent[chinumber], x01);
            break;

        case core::chemical::aa_met:
            result = bbutils::get_1d_from_dst(aa_sm_3d.front().lib_independent[chinumber], x01);
            break;
        case core::chemical::aa_glu:
            result = bbutils::get_1d_from_dst(aa_sm_3d[1].lib_independent[chinumber], x01);
            break;
        case core::chemical::aa_gln:
            result = bbutils::get_1d_from_dst(aa_sm_3d[2].lib_independent[chinumber], x01);
            break;
        case core::chemical::aa_pro:
            result = bbutils::get_1d_from_dst(aa_sm_3d.back().lib_independent[chinumber], x01);
            break;

        case core::chemical::aa_arg:
            result = bbutils::get_1d_from_dst(aa_sm_4d.front().lib_independent[chinumber], x01);
            break;
        case core::chemical::aa_lys:
            result = bbutils::get_1d_from_dst(aa_sm_4d.back().lib_independent[chinumber], x01);
            break;

        default:
            break;
    }
    return result;
}

//    std::vector<sm_1d> aa_sm_1d; // 0 ser, 1 val, 2 cys, 3 thr
//    std::vector<sm_2d> aa_sm_2d; // 0 trp, 1 his, 2 asn, 3 asp, 4 phe, 5 tyr, 6 ile, 7 leu
//    std::vector<sm_3d> aa_sm_3d; // 0 met, 1 glu, 2 gln, 3 pro
//    std::vector<sm_4d> aa_sm_4d; // 0 arg, 1 lys

double BBDEP_Dunbrack_sm::get_degree_bbdep_from_phi_psi_x01_chinumber(size_t index,
        core::chemical::AA amino_acid,
        double x01,
        size_t chinumber) const
{
    double result = 0;

    switch(amino_acid)
    {
        case core::chemical::aa_ser:
            result = bbutils::get_1d_from_dst(aa_sm_1d.front().lib_cdf_sum_all[index][chinumber], x01);
            break;
        case core::chemical::aa_val:
            result = bbutils::get_1d_from_dst(aa_sm_1d[1].lib_cdf_sum_all[index][chinumber], x01);
            break;
        case core::chemical::aa_cys:
            result = bbutils::get_1d_from_dst(aa_sm_1d[2].lib_cdf_sum_all[index][chinumber], x01);
            break;
        case core::chemical::aa_thr:
            result = bbutils::get_1d_from_dst(aa_sm_1d.back().lib_cdf_sum_all[index][chinumber], x01);
            break;

        case core::chemical::aa_trp:
            result = bbutils::get_1d_from_dst(aa_sm_2d.front().lib_cdf_sum_all[index][chinumber], x01);
            break;
        case core::chemical::aa_his:
            result = bbutils::get_1d_from_dst(aa_sm_2d[1].lib_cdf_sum_all[index][chinumber], x01);
            break;
        case core::chemical::aa_asn:
            result = bbutils::get_1d_from_dst(aa_sm_2d[2].lib_cdf_sum_all[index][chinumber], x01);
            break;
        case core::chemical::aa_asp:
            result = bbutils::get_1d_from_dst(aa_sm_2d[3].lib_cdf_sum_all[index][chinumber], x01);
            break;
        case core::chemical::aa_phe:
            result = bbutils::get_1d_from_dst(aa_sm_2d[4].lib_cdf_sum_all[index][chinumber], x01);
            break;
        case core::chemical::aa_tyr:
            result = bbutils::get_1d_from_dst(aa_sm_2d[5].lib_cdf_sum_all[index][chinumber], x01);
            break;
        case core::chemical::aa_ile:
            result = bbutils::get_1d_from_dst(aa_sm_2d[6].lib_cdf_sum_all[index][chinumber], x01);
            break;
        case core::chemical::aa_leu:
            result = bbutils::get_1d_from_dst(aa_sm_2d.back().lib_cdf_sum_all[index][chinumber], x01);
            break;

        case core::chemical::aa_met:
            result = bbutils::get_1d_from_dst(aa_sm_3d.front().lib_cdf_sum_all[index][chinumber], x01);
            break;
        case core::chemical::aa_glu:
            result = bbutils::get_1d_from_dst(aa_sm_3d[1].lib_cdf_sum_all[index][chinumber], x01);
            break;
        case core::chemical::aa_gln:
            result = bbutils::get_1d_from_dst(aa_sm_3d[2].lib_cdf_sum_all[index][chinumber], x01);
            break;
        case core::chemical::aa_pro:
            result = bbutils::get_1d_from_dst(aa_sm_3d.back().lib_cdf_sum_all[index][chinumber], x01);
            break;

        case core::chemical::aa_arg:
            result = bbutils::get_1d_from_dst(aa_sm_4d.front().lib_cdf_sum_all[index][chinumber], x01);
            break;
        case core::chemical::aa_lys:
            result = bbutils::get_1d_from_dst(aa_sm_4d.back().lib_cdf_sum_all[index][chinumber], x01);
            break;

        default:
            break;
    }
    return result;
}

}
}
