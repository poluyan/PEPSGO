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
#ifndef INCLUDED_pepsgo_hh
#define INCLUDED_pepsgo_hh

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoreFunction.hh>

#include <string>
#include <vector>

#include "quantile.h"

namespace pepsgo
{

class PEPSGO
{
private:
    std::string bbdep_path;
    std::string bbind_path;

    std::string peptide_sequence;

    core::scoring::ScoreFunctionOP score_fn;

    core::pose::Pose peptide;
    core::pose::Pose ideal_peptide;

    core::pose::Pose abs_min_peptide;
    core::pose::Pose ideal_abs_min_peptide;
    core::pose::Pose all_superposed;
    
    
    /// 
    std::vector<std::shared_ptr<trie_based::TrieBased<trie_based::NodeCount<int>,int>>> phipsi_rama2_sample;
    std::vector<std::shared_ptr<empirical_quantile::ImplicitQuantileSorted<int, double>>> phipsi_rama2_quantile;
    
    ///
    int threads_number;
public:
    PEPSGO();
    void set_number_of_threads(size_t n);
    
    void set_peptide(std::string _peptide_sequence);
    void set_bbdep(std::string _bbdep_path);
    void set_bbind(std::string _bbind_path);
    core::Real objective(const std::vector<double> &x);
    
    // creating rama2 from second to next-to-last residue
    void fill_rama2_residue(core::pose::Pose &pep, core::scoring::ScoreFunctionOP &scorefn_rama2b, size_t ind, size_t step);
    void fill_rama2_sample();
    void fill_rama2_quantile();
};

}
#endif