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
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <numeric/NumericTraits.hh>
#include <core/id/TorsionID.hh>
#include <core/conformation/util.hh>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/sequence/util.hh>

#include <core/kinematics/MoveMap.hh>
#include <protocols/minimization_packing/MinMover.hh>

#include <core/scoring/rms_util.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>

#include <pepsgo.hh>
#include <data_io.hh>
#include <transform.hh>
#include <bbtools.hh>

#include <random>
#include <omp.h>

namespace pepsgo
{

	PEPSGO::PEPSGO()
	{
//		score_fn = core::scoring::ScoreFunctionFactory::create_score_function("score12"); // score12 talaris2014 ref2015 talaris2013
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
	void PEPSGO::set_multithread()
	{
		mt_peptide.resize(threads_number);
		for(size_t i = 0; i != mt_peptide.size(); i++)
		{
			mt_peptide[i] = peptide;
		}
		mt_score_fn.resize(threads_number);
		for(size_t i = 0; i != mt_score_fn.size(); i++)
		{
			mt_score_fn[i] = score_fn;
		}
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
		peptide.dump_pdb("output/pdb/extended.pdb");

//    ideal_peptide = peptide;
	}
	void PEPSGO::set_peptide(std::string _peptide_sequence)
	{
		peptide_sequence = _peptide_sequence;
		peptide_ss_predicted = std::string(_peptide_sequence.size(), 'C');
		unique_aa();
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
			peptide_ss_predicted = std::string(peptide_sequence.size(), 'C');
			unique_aa();
			std::cout << "loaded sequence " << peptide_sequence << std::endl;
			std::cout << "unique aa " << peptide_amino_acids << std::endl;
		}
		else
		{
			std::cout << "fatal" << std::endl;
		}

		if(peptide_sequence.find("X") != std::string::npos)
		{
			std::cout << "X in sequence!" << std::endl;
		}

		core::pose::make_pose_from_sequence(peptide, peptide_sequence,
		                                    *core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
		extend_peptide();

		if(basic::options::option[basic::options::OptionKeys::in::file::native].user())
		{
			std::string fname = basic::options::option[basic::options::OptionKeys::in::file::native].value_string();
			std::cout << "native peptide " << fname << std::endl;
			core::import_pose::pose_from_file(peptide_native, fname);
			peptide_native.dump_pdb("output/pdb/peptide_native.pdb");
		}
		else
		{
			std::cout << "native peptide from sequence" << std::endl;
			core::pose::make_pose_from_sequence(peptide_native, peptide_sequence,
			                                    *core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
		}

		if(peptide_sequence != peptide_native.sequence())
		{
			std::cout << "peptide sequence != native sequence" << std::endl;
		}

		pepsgo::bbtools::make_ideal_peptide(peptide_native_ideal, peptide_native);
		peptide_native_ideal_optimized = peptide_native_ideal;

		superposed_aa = peptide_native_ideal;
		superposed_ca = peptide_native_ideal;

		peptide_native_ideal.dump_pdb("output/pdb/peptide_native_ideal.pdb");
	}
	void PEPSGO::set_bbdep(size_t step/*std::string _bbdep_path*/)
	{
		//bbdep_path = _bbdep_path;

		std::string lib_path = basic::options::option[basic::options::OptionKeys::in::path::database].value_string() + "rotamer/ExtendedOpt1-5/";
		std::cout << lib_path << std::endl;

		bbdep_sm.set_path(lib_path);
		//bbdep_sm.set_step(1000);
		bbdep_sm.set_step(step);
		bbdep_sm.initialize_all(true, peptide_amino_acids, threads_number);

//    pepsgo::bbdep::plot_chi1_all(bbdep_sm);

//    pepsgo::bbind::BBIND_top obj;
//    obj.set_path("/ssdwork/ProjectsCPP/mcmc/bbind_top500/");
//    obj.initialize(2, 2, 2, 2,"SVCT");
//    obj.initialize(2, 2, 2, 2,"WHNDFYIL");
//    obj.initialize(2, 2, 2, 2,"MEQP");
//    obj.initialize(2, 2, 2, 2,"RK");
	}
	void PEPSGO::optimize_native()
	{
		core::kinematics::MoveMapOP movemap_minimizer_ = core::kinematics::MoveMapOP(new core::kinematics::MoveMap());
		for(const auto &i : opt_vector_native)
			movemap_minimizer_->set(i.dofid, true);

//    movemap_minimizer_->set_bb(true);
//    movemap_minimizer_->set_chi(true);

		protocols::minimization_packing::MinMover minimizer(movemap_minimizer_, score_fn, "lbfgs", 1e-20, true);
		minimizer.max_iter(5000);
		minimizer.apply(peptide_native_ideal_optimized);
		(*score_fn)(peptide_native_ideal_optimized);

//		peptide_native_ideal_optimized.energies().show_totals(std::cout);
//		peptide_native_ideal_optimized.energies().show(std::cout);
//		peptide_native_ideal_optimized.energies().show_totals(std::cout);
//		std::cout << score_fn->score_by_scoretype(peptide_native_ideal_optimized, core::scoring::fa_elec, false) << std::endl;
//		std::cout << score_fn->score_by_scoretype(peptide_native_ideal_optimized, core::scoring::fa_elec, true) << std::endl;
//		std::cout << peptide_native_ideal_optimized.energies().total_energy() - score_fn->score_by_scoretype(peptide_native_ideal_optimized, core::scoring::fa_elec, true) << std::endl;

		std::cout << "peptide native ideal optimized score " << peptide_native_ideal_optimized.energies().total_energy() << std::endl;
		peptide_native_ideal_optimized.dump_pdb("output/pdb/peptide_native_ideal_optimized.pdb");

		superposed_aa = peptide_native_ideal_optimized;
		superposed_ca = peptide_native_ideal_optimized;
	}
	void PEPSGO::set_native_state()
	{
		native_state.resize(std::get<2>(peptide_ranges_native.omega) + 1);
		for(size_t i = 0; i != std::get<2>(peptide_ranges_native.phipsi) + 1; i++)
		{
			native_state[i] = frags.get_index_phipsi(peptide_native.dof(opt_vector_native[i].dofid), i);
		}
		for(size_t i = std::get<1>(peptide_ranges_native.omega), k = 0; i != std::get<2>(peptide_ranges_native.omega) + 1; i++, k++)
		{
			native_state[i] = frags.get_index_omega(peptide_native.dof(opt_vector_native[i].dofid), k);
		}
		for(size_t i = 0; i != native_state.size(); i++)
		{
			std::cout << int(native_state[i]) << ", ";
		}
		std::cout << std::endl;
	}
	bool PEPSGO::is_native_in_frag_space()
	{
		return structures_triebased->search(native_state);
	}
	void PEPSGO::set_bbind(std::string _bbind_path)
	{
		bbind_path = _bbind_path;
	}
	core::Real PEPSGO::objective(const std::vector<double> &x)
	{
		std::vector<double> vt(x.size());
		transform_vec01_to_vecrad(x, vt);
		for(size_t i = 0, n = opt_vector.size(); i < n; ++i)
		{
			peptide.set_dof(opt_vector[i].dofid, vt[i]);
		}
		(*score_fn)(peptide);
		return peptide.energies().total_energy();
	}

	core::Real PEPSGO::mo_scalarizing(std::vector<double> &criteria_values, const std::vector<double> &weights) const
	{
		/// linear
		core::Real result = 0.0;
		for(size_t i = 0; i != criteria_values.size(); i++)
		{
			result += weights[i]*criteria_values[i];
		}
		return result;
		/// Germeier
//		core::Real make_positive = 1000.0*peptide.size();
//		for(size_t i = 0; i != criteria_values.size(); i++)
//		{
//			criteria_values[i] = weights[i]*(criteria_values[i] + make_positive);
//		}
//		return *std::max_element(criteria_values.begin(), criteria_values.end());
	}

	core::Real PEPSGO::objective_mo(const std::vector<double> &x, const std::vector<double> &weights, std::vector<double> &res)
	{
		std::vector<double> vt(x.size());
		transform_vec01_to_vecrad(x, vt);
		for(size_t i = 0, n = opt_vector.size(); i < n; ++i)
		{
			peptide.set_dof(opt_vector[i].dofid, vt[i]);
		}
		(*score_fn)(peptide);
		core::Real fa_elec = score_fn->score_by_scoretype(peptide, core::scoring::fa_elec, true) +
		                     score_fn->score_by_scoretype(peptide, core::scoring::lk_ball_wtd, true) +
		                     score_fn->score_by_scoretype(peptide, core::scoring::fa_intra_sol, true) +
		                     score_fn->score_by_scoretype(peptide, core::scoring::fa_sol, true)/* +
		                     score_fn->score_by_scoretype(peptide, core::scoring::fa_atr, true) +
		                     score_fn->score_by_scoretype(peptide, core::scoring::fa_rep, true)*/;
//		core::Real fa_nonbound = score_fn->score_by_scoretype(peptide, core::scoring::rama, true);
		std::vector<core::Real> criteria_values = {peptide.energies().total_energy() - fa_elec, fa_elec};
		res[0] = criteria_values[0];
		res[1] = criteria_values[1];
		res[2] = peptide.energies().total_energy();
		return mo_scalarizing(criteria_values, weights);
	}

	core::Real PEPSGO::objective_mt(const std::vector<double> &x, int th_id)
	{
		std::vector<double> vt(x.size());
		transform_vec01_to_vecrad(x, vt);
		for(size_t i = 0, n = opt_vector.size(); i < n; ++i)
		{
			mt_peptide[th_id].set_dof(opt_vector[i].dofid, vt[i]);
		}
		(*mt_score_fn[th_id])(mt_peptide[th_id]);
//    peptide.dump_pdb("output/pdb/peptide.pdb");
		return mt_peptide[th_id].energies().total_energy();
	}

	core::Real PEPSGO::objective_mt_mo(const std::vector<double> &x, int th_id, const std::vector<double> &weights, std::vector<double> &res)
	{
		std::vector<double> vt(x.size());
		transform_vec01_to_vecrad(x, vt);
		for(size_t i = 0, n = opt_vector.size(); i < n; ++i)
		{
			mt_peptide[th_id].set_dof(opt_vector[i].dofid, vt[i]);
		}
		(*mt_score_fn[th_id])(mt_peptide[th_id]);
		core::Real fa_elec = mt_score_fn[th_id]->score_by_scoretype(mt_peptide[th_id], core::scoring::fa_elec, true) +
		                     mt_score_fn[th_id]->score_by_scoretype(mt_peptide[th_id], core::scoring::lk_ball_wtd, true) +
		                     mt_score_fn[th_id]->score_by_scoretype(mt_peptide[th_id], core::scoring::fa_intra_sol, true) +
		                     mt_score_fn[th_id]->score_by_scoretype(mt_peptide[th_id], core::scoring::fa_sol, true)/* +
		                     mt_score_fn[th_id]->score_by_scoretype(mt_peptide[th_id], core::scoring::fa_atr, true) +
		                     mt_score_fn[th_id]->score_by_scoretype(mt_peptide[th_id], core::scoring::fa_rep, true)*/;

//		core::Real fa_nonbound = mt_score_fn[th_id]->score_by_scoretype(mt_peptide[th_id], core::scoring::rama, true);
		std::vector<core::Real> criteria_values = {mt_peptide[th_id].energies().total_energy() - fa_elec, fa_elec};
		res[0] = criteria_values[0];
		res[1] = criteria_values[1];
		res[2] = mt_peptide[th_id].energies().total_energy();
		return mo_scalarizing(criteria_values, weights);
	}

	void PEPSGO::write(const std::vector<double> &x, std::string fname)
	{
		std::vector<double> vt(x.size());
		transform_vec01_to_vecrad(x, vt);
		for(size_t i = 0, n = opt_vector.size(); i < n; ++i)
		{
			peptide.set_dof(opt_vector[i].dofid, vt[i]);
		}
		(*score_fn)(peptide);
		std::cout << "peptide score: " << peptide.energies().total_energy() << std::endl;
		peptide.dump_pdb(fname);
	}

	core::Real PEPSGO::optimize_lbfgs(const std::vector<double> &x)
	{
		std::vector<double> vt(x.size());
		transform_vec01_to_vecrad(x, vt);
		for(size_t i = 0, n = opt_vector.size(); i < n; ++i)
		{
			peptide.set_dof(opt_vector[i].dofid, vt[i]);
		}
		(*score_fn)(peptide);

		core::kinematics::MoveMapOP movemap_minimizer_ = core::kinematics::MoveMapOP(new core::kinematics::MoveMap());
		for(const auto &i : opt_vector_native)
			movemap_minimizer_->set(i.dofid, true);

		protocols::minimization_packing::MinMover minimizer(movemap_minimizer_, score_fn, "lbfgs", 1e-20, true);
		minimizer.max_iter(5000);
		minimizer.apply(peptide);
		(*score_fn)(peptide);
		return peptide.energies().total_energy();
	}

	core::Real PEPSGO::get_CA_rmsd_lbfgs(const std::vector<double> &x)
	{
		std::vector<double> vt(x.size());
		transform_vec01_to_vecrad(x, vt);
		for(size_t i = 0, n = opt_vector.size(); i < n; ++i)
		{
			peptide.set_dof(opt_vector[i].dofid, vt[i]);
		}
		(*score_fn)(peptide);

		core::kinematics::MoveMapOP movemap_minimizer_ = core::kinematics::MoveMapOP(new core::kinematics::MoveMap());
		for(const auto &i : opt_vector_native)
			movemap_minimizer_->set(i.dofid, true);

		protocols::minimization_packing::MinMover minimizer(movemap_minimizer_, score_fn, "lbfgs", 1e-20, true);
		minimizer.max_iter(5000);
		minimizer.apply(peptide);
		(*score_fn)(peptide);
		return core::scoring::CA_rmsd(peptide, peptide_native_ideal_optimized);
	}

	core::Real PEPSGO::get_AA_rmsd_lbfgs(const std::vector<double> &x)
	{
		std::vector<double> vt(x.size());
		transform_vec01_to_vecrad(x, vt);
		for(size_t i = 0, n = opt_vector.size(); i < n; ++i)
		{
			peptide.set_dof(opt_vector[i].dofid, vt[i]);
		}
		(*score_fn)(peptide);

		core::kinematics::MoveMapOP movemap_minimizer_ = core::kinematics::MoveMapOP(new core::kinematics::MoveMap());
		for(const auto &i : opt_vector_native)
			movemap_minimizer_->set(i.dofid, true);

		protocols::minimization_packing::MinMover minimizer(movemap_minimizer_, score_fn, "lbfgs", 1e-20, true);
		minimizer.max_iter(5000);
		minimizer.apply(peptide);
		(*score_fn)(peptide);
		return core::scoring::all_atom_rmsd(peptide, peptide_native_ideal_optimized);
	}

	void PEPSGO::write_lbfgs(const std::vector<double> &x, std::string fname)
	{
		std::vector<double> vt(x.size());
		transform_vec01_to_vecrad(x, vt);
		for(size_t i = 0, n = opt_vector.size(); i < n; ++i)
		{
			peptide.set_dof(opt_vector[i].dofid, vt[i]);
		}
		(*score_fn)(peptide);

		core::kinematics::MoveMapOP movemap_minimizer_ = core::kinematics::MoveMapOP(new core::kinematics::MoveMap());
		for(const auto &i : opt_vector_native)
			movemap_minimizer_->set(i.dofid, true);

		protocols::minimization_packing::MinMover minimizer(movemap_minimizer_, score_fn, "lbfgs", 1e-20, true);
		minimizer.max_iter(5000);
		minimizer.apply(peptide);
		(*score_fn)(peptide);
		std::cout << "peptide score: " << peptide.energies().total_energy() << std::endl;
		peptide.dump_pdb(fname);
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

		phipsi_rama2_sample[ind] = std::make_shared<mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<int>,int>>();
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
			  std::make_shared<mveqf::ImplicitQuantileSorted<int, double>>(
			    std::vector<double>(grid_number.size(), -numeric::NumericTraits<core::Real>::pi()),
			    std::vector<double>(grid_number.size(), numeric::NumericTraits<core::Real>::pi()),
			    grid_number
			  );
			auto lower = std::lower_bound(phipsi.begin(), phipsi.end(), peptide_sequence.substr(i, 3));
			phipsi_rama2_quantile[i]->set_sample_shared_and_fill_count(phipsi_rama2_sample[std::distance(phipsi.begin(), lower)]);

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
		std::string arguments;
		switch(task)
		{
			case 1:
				arguments = "phipsi omega allchi";
				pepsgo::insert_to_opt_vector(opt_vector, peptide, arguments, peptide_ranges);
				std::cout << std::get<1>(peptide_ranges.phipsi) << ' ' << std::get<2>(peptide_ranges.phipsi) << '\t';
				std::cout << std::get<1>(peptide_ranges.omega) << ' ' << std::get<2>(peptide_ranges.omega) << '\t';
				std::cout << std::get<1>(peptide_ranges.chi) << ' ' << std::get<2>(peptide_ranges.chi) << std::endl;
				peptide_ranges.do_chi = true;
				// native
				pepsgo::insert_to_opt_vector(opt_vector_native, peptide_native_ideal_optimized, arguments, peptide_ranges_native);
				set_omega_quantile(72);
				break;
			case 2:
				arguments = "phipsi omega allchi";
				pepsgo::insert_to_opt_vector(opt_vector, peptide, arguments, peptide_ranges);
				std::cout << std::get<1>(peptide_ranges.phipsi) << ' ' << std::get<2>(peptide_ranges.phipsi) << '\t';
				std::cout << std::get<1>(peptide_ranges.omega) << ' ' << std::get<2>(peptide_ranges.omega) << '\t';
				std::cout << std::get<1>(peptide_ranges.chi) << ' ' << std::get<2>(peptide_ranges.chi) << std::endl;
				peptide_ranges.do_chi = true;
				// native
				pepsgo::insert_to_opt_vector(opt_vector_native, peptide_native_ideal_optimized, arguments, peptide_ranges_native);
				fill_rama2_quantile(144);
				set_omega_quantile(72);
				break;
			case 4:
				arguments = "phipsi omega allchi";
				pepsgo::insert_to_opt_vector(opt_vector, peptide, arguments, peptide_ranges);
				std::cout << std::get<1>(peptide_ranges.phipsi) << ' ' << std::get<2>(peptide_ranges.phipsi) << '\t';
				std::cout << std::get<1>(peptide_ranges.omega) << ' ' << std::get<2>(peptide_ranges.omega) << '\t';
				std::cout << std::get<1>(peptide_ranges.chi) << ' ' << std::get<2>(peptide_ranges.chi) << std::endl;
				peptide_ranges.do_chi = true;
				// native
				pepsgo::insert_to_opt_vector(opt_vector_native, peptide_native_ideal_optimized, arguments, peptide_ranges_native);
				fill_rama2_quantile(144);
				set_omega_quantile(72);
				set_bbdep(144);
				break;
			case 5:
				arguments = "phipsi omega allchi";
				pepsgo::insert_to_opt_vector(opt_vector, peptide, arguments, peptide_ranges);
				std::cout << std::get<1>(peptide_ranges.phipsi) << ' ' << std::get<2>(peptide_ranges.phipsi) << '\t';
				std::cout << std::get<1>(peptide_ranges.omega) << ' ' << std::get<2>(peptide_ranges.omega) << '\t';
				std::cout << std::get<1>(peptide_ranges.chi) << ' ' << std::get<2>(peptide_ranges.chi) << std::endl;
				peptide_ranges.do_chi = true;
				// native
				pepsgo::insert_to_opt_vector(opt_vector_native, peptide_native_ideal_optimized, arguments, peptide_ranges_native);

				create_space_frag(std::make_pair(36, 72), std::make_pair(180, 190));
//    obj.create_space_frag(std::make_pair(8, 18), std::make_pair(36, 72));
//    obj.create_space_frag(std::make_pair(8, 12), std::make_pair(36, 72));
//    obj.create_space_frag(std::make_pair(6, 12), std::make_pair(72, 120));
//    obj.create_space_frag(std::make_pair(250, 254), std::make_pair(250, 254));
				fill_rama2_quantile(4);
				set_omega_quantile(72);
				set_bbdep(72);
				break;
			default:
				arguments = "phipsi allchi";
				pepsgo::insert_to_opt_vector(opt_vector, peptide, arguments, peptide_ranges);
				std::cout << std::get<1>(peptide_ranges.phipsi) << ' ' << std::get<2>(peptide_ranges.phipsi) << '\t';
				std::cout << std::get<1>(peptide_ranges.chi) << ' ' << std::get<2>(peptide_ranges.chi) << std::endl;
				peptide_ranges.do_chi = true;
				// native
				pepsgo::insert_to_opt_vector(opt_vector_native, peptide_native_ideal_optimized, arguments, peptide_ranges_native);
				break;
		}
	}

	void PEPSGO::create_space_frag(std::pair<std::uint8_t,std::uint8_t> phipsi_minmax, std::pair<std::uint8_t,std::uint8_t> omega_minmax)
	{
		structures_triebased = std::make_shared<frag_type>();
		frags.set_peptide(peptide);
		frags.set_file();
		frags.set_psipred(phipsi_minmax, omega_minmax);
		frags.fill_grids();
		set_native_state();
		frags.set_storage_shared(structures_triebased);
		frags.load_frag_file();
		frags.set_native_state(native_state);
		frags.make_permutations(1); // 0 - all, 1 - one chain, 2 - one chain Von Neumann, 3 - coil chain
		peptide_ss_predicted = frags.get_ss_predicted();

//		structures_triebased->insert(native_state);

		auto bonds = frags.get_bounds();
		structures_quant = std::make_shared<mveqf::ImplicitQuantile<std::uint8_t, double>>(
		                     std::vector<double>(bonds.size(), -numeric::NumericTraits<core::Real>::pi()),
		                     std::vector<double>(bonds.size(), numeric::NumericTraits<core::Real>::pi()),
		                     bonds);
		structures_quant->set_sample_shared_and_fill_count(structures_triebased);

		auto closest = frags.get_closest_to_native();
		std::cout << "is native in frag space " << is_native_in_frag_space() << std::endl;

//    structures_triebased->insert(closest);

		std::cout << "native vs closest" << std::endl;
		size_t sum = 0;
		for(size_t i = 0; i != closest.size(); i++)
		{
			int t = std::abs(int(native_state[i]) - int(closest[i]));
			sum += std::min(t, std::abs(int(bonds[i]) - t));
			std::cout << bonds[i] << '\t' << int(native_state[i]) << '\t' << int(closest[i]) << '\t'
			          << std::min(t, std::abs(int(bonds[i]) - t)) << std::endl;
		}
		std::cout << "sum = " << sum << std::endl;
	}

	size_t PEPSGO::get_opt_vector_size()
	{
		return opt_vector.size();
	}

	void PEPSGO::unique_aa()
	{
		std::set<char> aa_set;
		for(size_t i = 0; i != peptide_sequence.size(); i++)
		{
			if(aa_set.find(peptide_sequence[i]) == aa_set.end())
				aa_set.insert(peptide_sequence[i]);
		}
		for(const auto &i : aa_set)
		{
			peptide_amino_acids.push_back(i);
		}
	}

	core::Real PEPSGO::get_CA_rmsd(const std::vector<double> &x)
	{
		objective(x);
		return core::scoring::CA_rmsd(peptide, peptide_native_ideal_optimized);
	}
	core::Real PEPSGO::get_AA_rmsd(const std::vector<double> &x)
	{
		objective(x);
		return core::scoring::all_atom_rmsd(peptide, peptide_native_ideal_optimized);
	}

	void PEPSGO::append_to_superposed_aa(const std::vector<double> &x)
	{
		objective(x);
		core::pose::Pose temp_pose = peptide;

		// last bool CA_only
		protocols::simple_moves::SuperimposeMoverOP sm(new protocols::simple_moves::SuperimposeMover(
		      peptide_native_ideal_optimized, 1, peptide_native_ideal_optimized.total_residue(), 1, peptide_native_ideal_optimized.total_residue(), false));
		sm->apply(temp_pose);
		superposed_aa.append_pose_by_jump(temp_pose, superposed_aa.total_residue());
	}

	void PEPSGO::append_to_superposed_ca(const std::vector<double> &x)
	{
		objective(x);
		core::pose::Pose temp_pose = peptide;

		// last bool CA_only
		protocols::simple_moves::SuperimposeMoverOP sm(new protocols::simple_moves::SuperimposeMover(
		      peptide_native_ideal_optimized, 1, peptide_native_ideal_optimized.total_residue(), 1, peptide_native_ideal_optimized.total_residue(), false));
		sm->apply(temp_pose);
		superposed_ca.append_pose_by_jump(temp_pose, superposed_ca.total_residue());
	}

	void PEPSGO::dumb_superposed()
	{
		superposed_aa.dump_pdb("output/pdb/superposed_aa.pdb");
		superposed_ca.dump_pdb("output/pdb/superposed_ca.pdb");
	}

	void PEPSGO::set_task(const size_t choise)
	{
		task = choise > 5 ? 0 : choise;
		std::cout << "Optimization task is " << task << ", ";
		switch(task)
		{
			case 1:
				// only phi/psi/omega/chi angles
				std::cout << "only phi/psi/omega/chi angles" << std::endl;
				break;
			case 2:
				// phi/psi from rama2 (quantile), omega with delta, chi rotamers
				std::cout << "phi/psi from rama2 (quantile), omega with delta, chi rotamers" << std::endl;
				break;
			case 3:
				// phi/psi from rama2 (quantile), omega with delta, chi from TOP (quantile)
				std::cout << "phi/psi from rama2 (quantile), omega with delta, chi from TOP (quantile)" << std::endl;
				break;
			case 4:
				// phi/psi from rama2 (quantile), omega with delta, chi from Dunbrack (quantile)
				std::cout << "phi/psi from rama2 (quantile), omega with delta, chi from Dunbrack (quantile)" << std::endl;
				break;
			case 5:
				// phi/psi from fragments (quantile), omega with delta, chi from Dunbrack (quantile)
				std::cout << "phi/psi/omega from fragments (quantile), chi from Dunbrack (quantile)" << std::endl;
				break;
			default:
				std::cout << "only phi/psi/chi angles as rotamers" << std::endl;
				break;
		}
	}

	void PEPSGO::set_omega_quantile(size_t step)
	{
		omega_quantile.resize(1);
		std::vector<size_t> grid_number = {step};
		omega_quantile.front() = std::make_shared<mveqf::ImplicitQuantileSorted<int, double>>(
		                           std::vector<double>(grid_number.size(), -numeric::NumericTraits<core::Real>::pi()),
		                           std::vector<double>(grid_number.size(), numeric::NumericTraits<core::Real>::pi()),
		                           grid_number
		                         );
		std::vector<std::vector<int>> in_sample = {{0}, {int(step) - 1}};
		omega_quantile.front()->set_sample_and_fill_count(in_sample);
	}

	void PEPSGO::transform_with_frags(const std::vector<double> &x, std::vector<double> &out) const
	{
		std::vector<double> xx = x;
		// transform here
		std::vector<double> frag_vector(std::get<1>(peptide_ranges.chi));
		std::vector<double> x_temp(std::get<1>(peptide_ranges.chi));

		for(size_t i = 0; i != x_temp.size(); i++)
		{
			x_temp[i] = x[i];
		}

		structures_quant->transform(x_temp, frag_vector);

		for(size_t i = 0; i != frag_vector.size(); i++)
		{
			xx[i] = frag_vector[i];
		}

		out = pepsgo::transform::bbdep_experiment_actual_states(
		        xx, opt_vector, peptide_ranges, bbdep_sm, phipsi_rama2_quantile.size());

//    ///  C -> transform
//    // phi/psi
//    std::vector<double> values01(2), sampled(2);
//    for(size_t i = 1, j = 1; i != peptide_ss_predicted.size() - 1; i++, j+=2)
//    {
//        if(peptide_ss_predicted[i] == 'C')
//        {
////            std::cout << "C" << i << '\t' << opt_vector[j].seqpos << '\t' <<
////            opt_vector[j].torsion_name << std::endl;
//            values01[0] = x[j];
//            values01[1] = x[j + 1];
//            phipsi_rama2_quantile[i-1]->transform(values01, sampled);
//            vt[j] = sampled[0];
//            vt[j + 1] = sampled[1];
//        }
//    }
//    // omega
//    std::vector<double> values01_omg(1), sampled_omg(1);
//    for(size_t i = std::get<1>(peptide_ranges.omega), j = 0; i <= std::get<2>(peptide_ranges.omega); i++, j++)
//    {
//        if(peptide_ss_predicted[j] == 'C')
//        {
//            values01_omg[0] = x[i];
//            omega_quantile.front()->transform(values01_omg, sampled_omg);
//            vt[i] = sampled_omg.front();
//        }
//    }
//    ///

//    for(size_t i = 0, n = opt_vector.size(); i < n; ++i)
//    {
//        std::cout << vt[i] << '\t' << i << std::endl;
//    }
	}

	void PEPSGO::transform_vec01_to_vecrad(const std::vector<double> &x, std::vector<double> &out) const
	{
		switch(task)
		{
			case 1:
				transform_simple_omega(x, out);
				break;
			case 2:
				transform_simple_omega_rama2(x, out);
				break;
			case 4:
				transform_simple_omega_rama2_dun(x, out);
				break;
			case 5:
				transform_with_frags(x, out);
				break;
			default:
				transform_simple(x, out);
				break;
		}
	}

	void PEPSGO::transform_simple(const std::vector<double> &x, std::vector<double> &out) const
	{
		for(size_t i = 0; i != x.size(); i++)
		{
			out[i] = -numeric::NumericTraits<core::Real>::pi()*(1.0 - 2.0*x[i]);
		}
	}

	void PEPSGO::transform_simple_omega(const std::vector<double> &x, std::vector<double> &out) const
	{
		for(size_t i = std::get<1>(peptide_ranges.phipsi); i <= std::get<2>(peptide_ranges.phipsi); i++)
		{
			out[i] = -numeric::NumericTraits<core::Real>::pi()*(1.0 - 2.0*x[i]);
		}
		for(size_t i = std::get<1>(peptide_ranges.chi); i <= std::get<2>(peptide_ranges.chi); i++)
		{
			out[i] = -numeric::NumericTraits<core::Real>::pi()*(1.0 - 2.0*x[i]);
		}
		std::vector<double> values01_omg(1), sampled_omg(1);
		for(size_t i = std::get<1>(peptide_ranges.omega); i <= std::get<2>(peptide_ranges.omega); i++)
		{
			values01_omg[0] = x[i];
			omega_quantile.front()->transform(values01_omg, sampled_omg);
			out[i] = sampled_omg.front();
		}
	}

	void PEPSGO::transform_simple_omega_rama2(const std::vector<double> &x, std::vector<double> &out) const
	{
		out[std::get<1>(peptide_ranges.phipsi)] = -numeric::NumericTraits<core::Real>::pi()*(1.0 - 2.0*x[std::get<1>(peptide_ranges.phipsi)]);
		out[std::get<2>(peptide_ranges.phipsi)] = -numeric::NumericTraits<core::Real>::pi()*(1.0 - 2.0*x[std::get<2>(peptide_ranges.phipsi)]);
		std::vector<double> values01(2), sampled(2);
		for(size_t c = std::get<1>(peptide_ranges.phipsi) + 1, pp = 0; pp != phipsi_rama2_quantile.size(); pp++, c += 2)
		{
			values01[0] = x[c];
			values01[1] = x[c + 1];
			phipsi_rama2_quantile[pp]->transform(values01, sampled);
			out[c] = sampled[0];
			out[c + 1] = sampled[1];
		}

		for(size_t i = std::get<1>(peptide_ranges.chi); i <= std::get<2>(peptide_ranges.chi); i++)
		{
			out[i] = -numeric::NumericTraits<core::Real>::pi()*(1.0 - 2.0*x[i]);
		}
		std::vector<double> values01_omg(1), sampled_omg(1);
		for(size_t i = std::get<1>(peptide_ranges.omega); i <= std::get<2>(peptide_ranges.omega); i++)
		{
			values01_omg[0] = x[i];
			omega_quantile.front()->transform(values01_omg, sampled_omg);
			out[i] = sampled_omg.front();
		}
	}

	void PEPSGO::transform_simple_omega_rama2_dun(const std::vector<double> &x, std::vector<double> &out) const
	{
		out[std::get<1>(peptide_ranges.phipsi)] = -numeric::NumericTraits<core::Real>::pi()*(1.0 - 2.0*x[std::get<1>(peptide_ranges.phipsi)]);
		out[std::get<2>(peptide_ranges.phipsi)] = -numeric::NumericTraits<core::Real>::pi()*(1.0 - 2.0*x[std::get<2>(peptide_ranges.phipsi)]);
		std::vector<double> values01(2), sampled(2);
		for(size_t c = std::get<1>(peptide_ranges.phipsi) + 1, pp = 0; pp != phipsi_rama2_quantile.size(); pp++, c += 2)
		{
			values01[0] = x[c];
			values01[1] = x[c + 1];
			phipsi_rama2_quantile[pp]->transform(values01, sampled);
			out[c] = sampled[0];
			out[c + 1] = sampled[1];
		}
		std::vector<double> values01_omg(1), sampled_omg(1);
		for(size_t i = std::get<1>(peptide_ranges.omega); i <= std::get<2>(peptide_ranges.omega); i++)
		{
			values01_omg[0] = x[i];
			omega_quantile.front()->transform(values01_omg, sampled_omg);
			out[i] = sampled_omg.front();
		}
		out = pepsgo::transform::bbdep_experiment_actual_states(x, opt_vector, peptide_ranges, bbdep_sm, phipsi_rama2_quantile.size());
	}

//void PEPSGO::set_extend_search_space(size_t c)
//{
//    extend_space_count = c;
//    ess_fe_count = 0;
//    best_state_value = std::numeric_limits<core::Real>::max();
//    th_exss.resize(threads_number);
//}
//
//void PEPSGO::extend_search_space()
//{
//    auto bonds = frags.get_bounds();
//    auto point = best_state;
//    for(size_t i = 0; i != best_state.size(); i++)
//    {
//        point = best_state;
//        point[i] = point[i] + 1;
//
//        if(point[i] > bonds[i] - 1)
//        {
//            point[i] = 0;
//        }
//        if(!structures_triebased->search(point))
//        {
//            structures_triebased->insert(point);
//        }
//    }
//    for(size_t i = 0; i != best_state.size(); i++)
//    {
//        point = best_state;
//        if(point[i] == 0)
//            point[i] = bonds[i] - 1;
//        else
//            point[i] = point[i] - 1;
//        if(!structures_triebased->search(point))
//        {
//            structures_triebased->insert(point);
//        }
//    }
//}
//
//core::Real PEPSGO::objective_mt_extend(const std::vector<double> &x, int th_id)
//{
//    std::vector<double> vt(x.size());
//    transform_with_frags(x, vt);
//    for(size_t i = 0, n = opt_vector.size(); i < n; ++i)
//    {
//        mt_peptide[th_id].set_dof(opt_vector[i].dofid, vt[i]);
//    }
//    (*mt_score_fn[th_id])(mt_peptide[th_id]);
//
//    std::vector<std::uint8_t> current_state(std::get<2>(peptide_ranges.omega) + 1);
//    for(size_t i = 0; i != std::get<2>(peptide_ranges.phipsi) + 1; i++)
//    {
//        current_state[i] = frags.get_index_phipsi(mt_peptide[th_id].dof(opt_vector_native[i].dofid), i);
//    }
//    for(size_t i = std::get<1>(peptide_ranges.omega), k = 0; i != std::get<2>(peptide_ranges.omega) + 1; i++, k++)
//    {
//        current_state[i] = frags.get_index_omega(mt_peptide[th_id].dof(opt_vector_native[i].dofid), k);
//    }
//
//    core::Real rez = mt_peptide[th_id].energies().total_energy();
//    #pragma omp critical
//    {
//        ess_fe_count++;
//        auto sum = std::accumulate(th_exss.begin(), th_exss.end(), 0);
//        if(rez < best_state_value && ess_fe_count > extend_space_count && sum != threads_number)
//        {
//            th_exss[th_id] = 1;
//            best_state_value = rez;
//            best_state = current_state;
//            extend_search_space();
//            std::cout << "extend it!! " << th_id << std::endl;
//        }
//        if(sum == threads_number)
//        {
//            ess_fe_count = 0;
//            for(auto &i : th_exss)
//                i = 0;
//        }
//    }
//    return rez;
//}

}
