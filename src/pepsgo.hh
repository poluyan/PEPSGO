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
#include <filesystem>

#include <mveqf/implicit.h>
#include <get_dof.hh>
#include <transform.hh>
#include <fragment.hh>

namespace pepsgo
{

	class PEPSGO
	{
	private:
		std::string bbdep_path;
		std::string bbind_path;

		std::string peptide_sequence;
		std::string peptide_ss_predicted;
		std::string peptide_amino_acids;

		core::scoring::ScoreFunctionOP score_fn;
		std::vector<core::scoring::ScoreFunctionOP> mt_score_fn;

		core::pose::Pose peptide;
		std::vector<core::pose::Pose> mt_peptide;

		core::pose::Pose peptide_native;
		core::pose::Pose peptide_native_ideal;
		core::pose::Pose peptide_native_ideal_optimized;
		std::vector<std::uint8_t> native_state;

		core::pose::Pose superposed_aa;
		core::pose::Pose superposed_ca;

		///
		std::vector<pepsgo::opt_element> opt_vector;
		pepsgo::ranges peptide_ranges;

		std::vector<pepsgo::opt_element> opt_vector_spec;
		pepsgo::ranges peptide_ranges_spec;

		std::vector<pepsgo::opt_element> opt_vector_native;
		pepsgo::ranges peptide_ranges_native;
		///
		std::vector<std::shared_ptr<mveqf::TrieBased<mveqf::NodeCount<int>,int>>> phipsi_rama2_sample;
		std::vector<std::shared_ptr<mveqf::ImplicitQuantile<int, double>>> phipsi_rama2_quantile;
		std::vector<std::shared_ptr<mveqf::ImplicitQuantile<int, double>>> omega_quantile;

		///
		int threads_number;

		pepsgo::fragment::FragPick frags;
		typedef mveqf::TrieBased<mveqf::NodeCount<std::uint8_t>,std::uint8_t> frag_type;
		std::shared_ptr<frag_type> structures_triebased;
		std::shared_ptr<mveqf::ImplicitQuantile<std::uint8_t, double>> structures_quant;

		pepsgo::bbdep::BBDEP_Dunbrack_sm bbdep_sm;

		size_t task; // optimization task

		/// search space explore
//    size_t extend_space_count;
//    size_t ess_fe_count;
//    std::vector<size_t> th_exss;
//    core::Real best_state_value;
//    std::vector<std::uint8_t> best_state;
	public:
		PEPSGO(int argc, char **argv);
		
		void set_task(size_t task_number);
		
		void set_number_of_threads(size_t n);
		void set_multithread();

		void set_peptide(std::string _peptide_sequence);
		void set_peptide_from_file();

		void optimize_native();
		void set_native_state();
		bool is_native_in_frag_space();

		void create_space_frag(std::pair<std::uint8_t,std::uint8_t> phipsi_minmax, std::pair<std::uint8_t,std::uint8_t> omega_minmax);

		void unique_aa();

		void set_bbdep(size_t step/*std::string _bbdep_path*/);
		void set_bbind(std::string _bbind_path);
		void transform_with_frags(const std::vector<double> &x, std::vector<double> &out) const;
		core::Real objective(const std::vector<double> &x);
		core::Real objective_mo(const std::vector<double> &x, const std::vector<double> &weights, std::vector<double> &res);
		core::Real objective_mt(const std::vector<double> &x, int th_id);
		core::Real objective_mt_mo(const std::vector<double> &x, int th_id, const std::vector<double> &weights, std::vector<double> &res);
		void write(const std::vector<double> &x, std::string fname);
		void write_lbfgs(const std::vector<double> &x, std::string fname);
		core::Real optimize_lbfgs(const std::vector<double> &x);
		core::Real get_CA_rmsd_lbfgs(const std::vector<double> &x);
		core::Real get_AA_rmsd_lbfgs(const std::vector<double> &x);

		core::Real mo_scalarizing(std::vector<double> &criteria_values, const std::vector<double> &weights) const;

		// creating rama2 from second to next-to-last residue
		void fill_rama2_residue(core::pose::Pose &pep, core::scoring::ScoreFunctionOP &scorefn_rama2b, size_t ind, size_t step);
		void fill_rama2_quantile(size_t step);

		void fill_opt_vector();

		size_t get_problem_dimension();

		void extend_peptide();

		void set_quantile();

		void get_native_bb_angles() const;

		void set_task_settings(size_t choise);
		void transform_vec01_to_vecrad(const std::vector<double> &x, std::vector<double> &out) const;
		void transform_simple(const std::vector<double> &x, std::vector<double> &out) const;
		void transform_simple_omega(const std::vector<double> &x, std::vector<double> &out) const;
		void transform_simple_omega_rama2(const std::vector<double> &x, std::vector<double> &out) const;
		void transform_simple_omega_rama2_dun(const std::vector<double> &x, std::vector<double> &out) const;
		void transform_bbfixed_dun(const std::vector<double> &x, std::vector<double> &out) const;
		void transform_bb_chifixed(const std::vector<double> &x, std::vector<double> &out) const;
		void transform_bb_omegachifixed(const std::vector<double> &x, std::vector<double> &out) const;

		core::Real get_CA_rmsd(const std::vector<double> &x);
		core::Real get_AA_rmsd(const std::vector<double> &x);

		void append_to_native_and_superpose_AA(const std::vector<double> &x);
		void append_to_native_and_superpose_CA(const std::vector<double> &x);

		void dumb_superposed();

		void set_omega_quantile(size_t step);

//    void set_extend_search_space(size_t c);
//    void extend_search_space();
//    core::Real objective_mt_extend(const std::vector<double> &x, int th_id);
		bool dir_exists(const std::filesystem::path& p, std::filesystem::file_status s = std::filesystem::file_status{}) const;
	};

}
#endif
