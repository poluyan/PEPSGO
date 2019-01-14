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
#include "pepsgo.hh"

#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <numeric/NumericTraits.hh>
#include <core/id/TorsionID.hh>
#include <core/conformation/util.hh>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/sequence/util.hh>

#include "data_io.hh"

#include <random>
#include <omp.h>

namespace pepsgo
{

PEPSGO::PEPSGO()
{
    //score_fn = core::scoring::ScoreFunctionFactory::create_score_function("ref2015");
    score_fn = core::scoring::get_score_function();
    std::cout << "score function: " << score_fn->get_name() << std::endl;

    threads_number = 1;
}
void PEPSGO::set_number_of_threads(size_t n)
{
    threads_number = n;
    if(omp_get_num_procs() < threads_number)
    {
        threads_number = omp_get_num_procs();
    }
    std::cout << "number_of_threads: " << threads_number << std::endl;
}
void PEPSGO::extend_peptide()
{
    for(core::Size i = 1; i <= peptide.total_residue(); i++)
    {
        peptide.set_phi(i, -135.0);
        peptide.set_psi(i, 135.0);
        peptide.set_omega(i, 179.8);

        if(peptide.residue(i).type().aa() == core::chemical::aa_pro)
            continue;

        for(core::Size j = 1; j <= peptide.residue(i).nchi(); j++)
        {
            peptide.set_chi(j, i, 90.0);
        }
    }
    peptide.dump_pdb("output_pdb/peptide.pdb");

//    ideal_peptide = peptide;
}
void PEPSGO::set_peptide(std::string _peptide_sequence)
{
    peptide_sequence = _peptide_sequence;
    core::pose::make_pose_from_sequence(peptide, peptide_sequence,
                                        *core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
    extend_peptide();
}
void PEPSGO::set_peptide_from_file()
{
    if(basic::options::option[basic::options::OptionKeys::in::file::fasta].user())
    {
        std::string fasta = basic::options::option[basic::options::OptionKeys::in::file::fasta]()[1];
        std::cout << "reading " << fasta << std::endl;
        peptide_sequence = core::sequence::read_fasta_file_return_str(fasta);
        std::cout << "loaded sequence " << peptide_sequence << std::endl;
    }
    else
    {
        std::cout << "fatal" << std::endl;
    }

    if(peptide_sequence.find("X") != std::string::npos)
    {
        std::cout << "fatal" << std::endl;
    }
    
    core::pose::make_pose_from_sequence(peptide, peptide_sequence,
                                        *core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
    extend_peptide();
}
void PEPSGO::set_bbdep(std::string _bbdep_path)
{
    bbdep_path = _bbdep_path;

    /*std::string lib_path = basic::options::option[basic::options::OptionKeys::in::path::database].value_string() + "rotamer/ExtendedOpt1-5/";
    std::cout << lib_path << std::endl;
    pepsgo::bbdep::BBDEP_Dunbrack_sm bbdep_sm;
    bbdep_sm.set_path(lib_path);
    bbdep_sm.set_step(1000);
    bbdep_sm.initialize_all(true, "SVCT");

    pepsgo::bbdep::plot_chi1_all(bbdep_sm);

    pepsgo::bbind::BBIND_top obj;
    obj.set_path("/ssdwork/ProjectsCPP/mcmc/bbind_top500/");*/

//    obj.initialize(2, 2, 2, 2,"SVCT");
//    obj.initialize(2, 2, 2, 2,"WHNDFYIL");
//    obj.initialize(2, 2, 2, 2,"MEQP");
//    obj.initialize(2, 2, 2, 2,"RK");
}
void PEPSGO::set_bbind(std::string _bbind_path)
{
    bbind_path = _bbind_path;
}
core::Real PEPSGO::objective(const std::vector<double> &x)
{
    for(size_t i = 0, n = opt_vector.size(); i < n; ++i)
    {
        peptide.set_dof(opt_vector[i].dofid, x[i]);
    }
    (*score_fn)(peptide);
    return peptide.energies().total_energy();
}

void PEPSGO::fill_rama2_residue(core::pose::Pose &pep, core::scoring::ScoreFunctionOP &scorefn_rama2b, size_t ind, size_t step)
{
    //core::pose::Pose pep;
    //core::pose::make_pose_from_sequence(pep, subseq,
    //                                    *core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
    core::Size position = 2;

    std::vector<size_t> grid_number = {step, step};

    std::vector<std::vector<double>> grids(grid_number.size());
    std::vector<double> dx(grid_number.size());

    for(size_t i = 0; i != grids.size(); i++)
    {
        std::vector<double> grid(grid_number[i] + 1);
        double startp = -numeric::NumericTraits<core::Real>::pi();
        double endp = numeric::NumericTraits<core::Real>::pi();
        double es = endp - startp;
        for(size_t j = 0; j != grid.size(); j++)
        {
            grid[j] = startp + j*es/double(grid_number[i]);
        }
        grids[i] = grid;
        dx[i] = es/(double(grid_number[i])*2);
    }

    std::vector<std::vector<double>> pdf(grid_number.front(), std::vector<double>(grid_number.back()));
    for(size_t i = 0; i != pdf.size(); i++)
    {
        core::id::TorsionID bb_phi_tid(position, core::id::BB, core::id::phi_torsion);
        pep.set_dof(pep.conformation().dof_id_from_torsion_id(bb_phi_tid), core::Real(grids[0][i] + dx.front()));

        for(size_t j = 0; j != pdf[i].size(); j++)
        {
            core::id::TorsionID bb_psi_tid(position, core::id::BB, core::id::psi_torsion);
            pep.set_dof(pep.conformation().dof_id_from_torsion_id(bb_psi_tid), core::Real(grids[1][j] + dx.back()));

            (*scorefn_rama2b)(pep);
            pdf[i][j] = -pep.energies().total_energy();
        }
    }

    double min_pdf = std::accumulate(pdf.begin(), pdf.end(), pdf[0][0], [](double min, const std::vector<double> &v)
    {
        return std::min(min,
                        *std::max_element(v.begin(),
                                          v.end()));
    });
    double max_pdf = std::accumulate(pdf.begin(), pdf.end(), pdf[0][0], [](double max, const std::vector<double> &v)
    {
        return std::max(max,
                        *std::max_element(v.begin(),
                                          v.end()));
    });
    for(auto &a : pdf)
    {
        for(auto &b : a)
        {
            b = (b - min_pdf) / (max_pdf - min_pdf);
        }
    }

    // only 0 and 1
    for(auto &a : pdf)
    {
        for(auto &b : a)
        {
            if(b > 0)
                b = 2;
        }
    }

    phipsi_rama2_sample[ind] = std::make_shared<trie_based::TrieBased<trie_based::NodeCount<int>,int>>();
    for(size_t i = 0; i != pdf.size(); i++)
    {
        for(size_t j = 0; j != pdf[i].size(); j++)
        {
            if(pdf[i][j] > 1.0)
            {
                std::vector<int> temp = {static_cast<int>(i), static_cast<int>(j)};
                if(!phipsi_rama2_sample[ind]->search(temp))
                {
                    phipsi_rama2_sample[ind]->insert(temp);
                }
            }
        }
    }
}

void PEPSGO::fill_rama2_quantile(size_t step)
{
    // explicit zero in grid points must be avoided
    if(step % 2)
    {
        --step;
    }

    std::set<std::string> phipsi_set;
    for(size_t i = 0, n = peptide.total_residue() - 2; i < n; i++)
    {
        std::string subseq = peptide_sequence.substr(i, 3);
        if(phipsi_set.find(subseq) == phipsi_set.end())
            phipsi_set.insert(subseq);
    }
    std::vector<std::string> phipsi;
    for(const auto &i : phipsi_set)
    {
//        std::cout << i << '\t';
        phipsi.push_back(i);
    }
//    std::cout << std::endl;
    phipsi_rama2_sample.resize(phipsi.size());

    std::vector<core::pose::Pose> peps(phipsi_set.size());
    for(size_t i = 0, n = phipsi.size(); i != n; i++)
    {
        core::pose::make_pose_from_sequence(peps[i], phipsi[i],
                                            *core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
    }

    std::vector<core::scoring::ScoreFunctionOP> scorefns(phipsi_set.size());
    core::scoring::ScoreFunctionOP scorefn_rama2b1(new core::scoring::ScoreFunction);
    scorefn_rama2b1->set_weight(core::scoring::rama2b, 1.0);
    for(auto &i : scorefns)
    {
        i = scorefn_rama2b1;
    }
    int n = phipsi_set.size();
    omp_set_dynamic(0);
    omp_set_num_threads(threads_number);
    #pragma omp parallel for
    for(int i = 0; i < n; ++i)
    {
        fill_rama2_residue(peps[i], scorefns[i], i, step);
    }

    phipsi_rama2_quantile.resize(peptide.total_residue() - 2);
    std::vector<size_t> grid_number = {step, step};
    for(size_t i = 0; i != phipsi_rama2_quantile.size(); i++)
    {
        phipsi_rama2_quantile[i] =
            std::make_shared<empirical_quantile::ImplicitQuantileSorted<int, double>>(
                std::vector<double>(grid_number.size(), -numeric::NumericTraits<core::Real>::pi()),
                std::vector<double>(grid_number.size(), numeric::NumericTraits<core::Real>::pi()),
                grid_number
            );
        auto lower = std::lower_bound(phipsi.begin(), phipsi.end(), peptide_sequence.substr(i, 3));
        phipsi_rama2_quantile[i]->set_sample_shared(phipsi_rama2_sample[std::distance(phipsi.begin(), lower)]);

        // check
//        std::mt19937_64 generator;
//        generator.seed(1);
//        std::uniform_real_distribution<double> ureal01(0.0,1.0);
//
//        std::vector<std::vector<double> > values01;
//        std::vector<std::vector<double> > sampled;
//        std::vector<double> temp1(grid_number.size());
//        std::vector<double> temp2(temp1.size());
//        for(size_t i = 0; i != 1e+6; ++i)
//        {
//            for(size_t j = 0; j != temp1.size(); j++)
//            {
//                temp1[j] = ureal01(generator);
//            }
//            values01.push_back(temp1);
//            sampled.push_back(temp2);
//        }
//        for(size_t j = 0; j != values01.size(); j++)
//            phipsi_rama2_quantile[i]->transform(values01[j], sampled[j]);
//        pepsgo::write_default2d("maps/" + std::to_string(i) + ".dat", sampled, 5);
    }
}


void PEPSGO::fill_opt_vector()
{
    opt_vector.clear();
    // std::string arguments = "phipsi omega allchiexceptpro sctheta12 bbtheta12 mctheta scd12 bbd12 mcd scp12 bbp12 frp12 lrp12 12chi";
    std::string arguments = "phipsi omega allchi";
    pepsgo::insert_to_opt_vector(opt_vector, peptide, arguments, peptide_ranges);
    std::cout << std::get<1>(peptide_ranges.phipsi) << ' ' << std::get<2>(peptide_ranges.phipsi) << '\t';
    std::cout << std::get<1>(peptide_ranges.omega) << ' ' << std::get<2>(peptide_ranges.omega) << '\t';
    std::cout << std::get<1>(peptide_ranges.chi) << ' ' << std::get<2>(peptide_ranges.chi) << std::endl;


}

}
