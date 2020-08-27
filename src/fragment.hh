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
#ifndef INCLUDED_fragment_hh
#define INCLUDED_fragment_hh

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <core/pose/Pose.hh>
#include <core/sequence/util.hh>
#include <numeric/conversions.hh>

//#include <core/fragment/SecondaryStructure.hh>
//#include <core/io/external/PsiPredInterface.hh>

#include <string>
#include <vector>

#include <quantile.hh>

namespace pepsgo
{
	namespace fragment
	{

		struct Frag
		{
			char aa;
			core::Real phi;
			core::Real psi;
			core::Real omg;

			Frag();
			Frag(char _aa, core::Real _phi, core::Real _psi, core::Real _omg);
		};

		class FragPick
		{
		private:
			//std::string psipred_fname;
			std::string fragmet_file;
			std::string peptide_seq;
			std::vector<std::uint8_t> native_state;
			std::vector<std::uint8_t> closest;

//    core::fragment::SecondaryStructure ss_profile;
			std::vector<core::Real> confidence;
			std::string ss_predicted;
			std::vector<size_t> step_num_phipsi;
			std::vector<size_t> step_num_omega;

			size_t frag_size;

			std::vector<std::vector<std::vector<Frag>>> all_fragments;

			std::vector<core::Size> phipsi_grid_number;
			std::vector<std::vector<core::Real>> phipsi_grids_radians;
			std::vector<std::vector<core::Real>> phipsi_grids_nodes;
			std::vector<core::Real> phipsi_grids_dx;

			std::vector<core::Size> omega_grid_number;
			std::vector<std::vector<core::Real>> omega_grids_radians;
			std::vector<std::vector<core::Real>> omega_grids_nodes;
			std::vector<core::Real> omega_grids_dx;
			core::Real omega_min_val;
			core::Real omega_max_val;

//    std::vector<std::vector<Frag> > structures;

			typedef mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<std::uint8_t>,std::uint8_t> sample_type;
			std::shared_ptr<sample_type> structures_trie;

			core::pose::Pose peptide;
			core::pose::Pose peptide_centroid;
		public:
			FragPick();
			void set_peptide(const core::pose::Pose &_peptide);

			void set_file();

			void fill_grids();

			size_t get_index_phipsi(core::Real radian, size_t index) const;
			size_t get_index_omega(core::Real radian, size_t index) const;

			void load_frag_file();
			void make_permutations(size_t type);
			void all_possible(size_t residues, size_t n_frags);
			void one_chain(size_t residues, size_t n_frags);
			void multistage(size_t residues, size_t n_frags);
			void one_chain_Von_Neumann(size_t residues, size_t n_frags);
			void coil_chain(size_t residues, size_t n_frags);
//    void write_structures_pdb(size_t delta);
//    void add_phipsi(empirical_quantile::ImplicitQuantile<std::uint8_t,double> &phipsi_quantile);
//
//    void fill_sample(std::shared_ptr<trie_based::TrieBased<trie_based::NodeCount<std::uint8_t>,std::uint8_t>> sample);
//
			void set_storage_shared(std::shared_ptr<sample_type> in_sample);
			std::vector<size_t> get_bounds();
//    void test();
			void set_psipred(std::pair<std::uint8_t,std::uint8_t> phipsi_minmax, std::pair<std::uint8_t,std::uint8_t> omega_minmax);
			void set_native_state(const std::vector<std::uint8_t> &nt);
			std::vector<std::uint8_t> get_closest_to_native();

			std::string get_ss_predicted() const;
		};

	}
}
#endif
