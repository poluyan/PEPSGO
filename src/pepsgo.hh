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

#include "quantile.hh"
#include "get_dof.hh"
#include "transform.hh"
#include "fragment.hh"

namespace pepsgo
{

class PEPSGO
{
private:
    std::string bbdep_path;
    std::string bbind_path;

    std::string peptide_sequence;
    std::string peptide_amino_acids;

    core::scoring::ScoreFunctionOP score_fn;
    std::vector<core::scoring::ScoreFunctionOP> mt_score_fn;
    
    core::pose::Pose peptide;
    std::vector<core::pose::Pose> mt_peptide;
    
    core::pose::Pose peptide_native;
    core::pose::Pose peptide_native_ideal;
    core::pose::Pose peptide_native_ideal_optimized;
    std::vector<std::uint8_t> native_state;
    
//    core::pose::Pose ideal_peptide;
//    core::pose::Pose abs_min_peptide;
//    core::pose::Pose ideal_abs_min_peptide;
//    core::pose::Pose all_superposed;
    
    ///
    std::vector<pepsgo::opt_element> opt_vector;
    pepsgo::ranges peptide_ranges;

    std::vector<pepsgo::opt_element> opt_vector_native;
    pepsgo::ranges peptide_ranges_native; 
    /// 
    std::vector<std::shared_ptr<trie_based::TrieBased<trie_based::NodeCount<int>,int>>> phipsi_rama2_sample;
    std::vector<std::shared_ptr<empirical_quantile::ImplicitQuantile<int, double>>> phipsi_rama2_quantile;
    
    ///
    int threads_number;
    
    pepsgo::fragment::FragPick frags;
    typedef trie_based::TrieBased<trie_based::NodeCount<std::uint8_t>,std::uint8_t> frag_type;
    std::shared_ptr<frag_type> structures_triebased;
    std::shared_ptr<empirical_quantile::ImplicitQuantile<std::uint8_t, double>> structures_quant;
    
    pepsgo::bbdep::BBDEP_Dunbrack_sm bbdep_sm;
public:
    PEPSGO();
    void set_number_of_threads(size_t n);
    void set_multithread();
    
    void set_peptide(std::string _peptide_sequence);
    void set_peptide_from_file();
    
    void optimize_native();
    void set_native_state();
    bool is_native_in_frag_space();
    
    void create_space_frag(size_t phipsi_step_min, size_t omega_step_min);
    
    void unique_aa();
    
    void set_bbdep(size_t step/*std::string _bbdep_path*/);
    void set_bbind(std::string _bbind_path);
    core::Real objective(const std::vector<double> &x);
    core::Real objective_mt(const std::vector<double> &x, int th_id);
    void write(const std::vector<double> &x);
    
    // creating rama2 from second to next-to-last residue
    void fill_rama2_residue(core::pose::Pose &pep, core::scoring::ScoreFunctionOP &scorefn_rama2b, size_t ind, size_t step);
    void fill_rama2_quantile(size_t step);
    
    void fill_opt_vector();
    
    size_t get_opt_vector_size();
    
    void extend_peptide();
    
    void set_quantile();
};

}
#endif
