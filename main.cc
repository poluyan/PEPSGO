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
#include <devel/init.hh>

#include <pepsgo.hh>

#include <data_io.hh>
#include "jadecso.h"

struct Example
{
	std::function<core::Real(const std::vector<double>&)> FitnessFunction;
	std::function<core::Real(const std::vector<double>&, int)> FitnessFunctionMultiThread;
};

int main(int argc, char *argv[])
{
	devel::init(argc, argv);
	std::cout << "Start..." << std::endl;

	size_t thread_num = 4;

	pepsgo::PEPSGO obj;
	obj.set_number_of_threads(thread_num);
	obj.set_peptide_from_file();
	obj.set_task(2);
	obj.fill_opt_vector();
	obj.optimize_native();
	obj.set_multithread();

//	std::cout << obj.get_opt_vector_size() << std::endl;
//	std::function<core::Real(const std::vector<double>&)> f = std::bind(&pepsgo::PEPSGO::objective, obj, std::placeholders::_1);
//	std::function<core::Real(const std::vector<double>&, int)> f_mt = std::bind(&pepsgo::PEPSGO::objective_mt, obj, std::placeholders::_1, std::placeholders::_2);
//	std::cout << f(std::vector<core::Real>(obj.get_opt_vector_size(),0.4)) << std::endl;
//	std::cout << f_mt(std::vector<core::Real>(obj.get_opt_vector_size(),0.4),thread_num-1) << std::endl;
//
//	std::cout << obj.get_AA_rmsd(std::vector<core::Real>(obj.get_opt_vector_size(),0.4)) << std::endl;
//	std::cout << obj.get_CA_rmsd(std::vector<core::Real>(obj.get_opt_vector_size(),0.4)) << std::endl;
//
//	Example fObj;
//	fObj.FitnessFunction = std::bind(&pepsgo::PEPSGO::objective, obj, std::placeholders::_1);
//	fObj.FitnessFunctionMultiThread = std::bind(&pepsgo::PEPSGO::objective_mt, obj, std::placeholders::_1, std::placeholders::_2);
//
//	obj.append_to_superposed_aa(std::vector<core::Real>(obj.get_opt_vector_size(),0.4));
//	obj.append_to_superposed_ca(std::vector<core::Real>(obj.get_opt_vector_size(),0.4));
//	obj.dumb_superposed();





std::vector<double> AA_rmsd;
	std::vector<double> CA_rmsd;
	std::vector<double> lAA_rmsd;
	std::vector<double> lCA_rmsd;
	std::vector<double> pAA_rmsd;
	std::vector<double> pCA_rmsd;
	std::vector<double> score_values;
	std::vector<double> lbfgs_score_values;
	std::vector<double> pscore_values;
	std::vector<std::vector<double>> score_front_values;
	size_t BENCH_N = 15;
	size_t fe_max = 2000000;
	for(size_t g = 0; g != BENCH_N; g++)
	{
		jademo::best = std::numeric_limits<double>::max();
		jademo::fe_count = 0;
		jademo::convergence.clear();

		size_t best_index = 0;

		jademo::SubPopulation sub_population(fe_max, 2);
		size_t total_population = 1000, dimension = obj.get_opt_vector_size();
		if(sub_population.Init(total_population, dimension) == jademo::kDone)
		{
			sub_population.set_generator_seed(g + 1);
			sub_population.FitnessFunction = std::bind(&pepsgo::PEPSGO::objective_mt_mo, obj, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
			sub_population.SetAllBoundsVectors(std::vector<double>(dimension, 0.0), std::vector<double>(dimension, 1.0));
			sub_population.SetTargetToMinimum();

			//sub_population.SwitchOffPMCRADE();

			sub_population.SetTotalGenerationsMax(fe_max%total_population ? fe_max/total_population : fe_max/total_population - 1);
			sub_population.SetBestShareP(0.05);
			sub_population.SetAdapitonFrequencyC(0.5);

			sub_population.RunOptimization();

			double temp_best = std::numeric_limits<double>::max();
			for(size_t i = 0; i != sub_population.score_values.size(); i++)
			{
				if(sub_population.score_values[i] < temp_best)
				{
					temp_best = sub_population.score_values[i];
					best_index = i;
				}
			}
			score_values.push_back(sub_population.score_values[best_index]);
			lbfgs_score_values.push_back(obj.optimize_lbfgs(sub_population.f1f2_vec[best_index]));

			obj.write(sub_population.f1f2_vec[best_index], "maps/mo/jade_A12mo/pdb/peptide" + std::to_string(g) + ".pdb");
			obj.write_lbfgs(sub_population.f1f2_vec[best_index], "maps/mo/jade_A12mo/pdb/peptide_lbfgs" + std::to_string(g) + ".pdb");
			obj.append_to_superposed_aa(sub_population.f1f2_vec[best_index]);
			obj.append_to_superposed_ca(sub_population.f1f2_vec[best_index]);
			obj.dumb_superposed();
		}

		AA_rmsd.push_back(obj.get_AA_rmsd(sub_population.f1f2_vec[best_index]));
		CA_rmsd.push_back(obj.get_CA_rmsd(sub_population.f1f2_vec[best_index]));
		lAA_rmsd.push_back(obj.get_AA_rmsd_lbfgs(sub_population.f1f2_vec[best_index]));
		lCA_rmsd.push_back(obj.get_CA_rmsd_lbfgs(sub_population.f1f2_vec[best_index]));
		pepsgo::write_default1d("maps/mo/jade_A12mo/conv" + std::to_string(g) + ".dat", jademo::convergence, 100, 10);

		pepsgo::write_default1d("maps/mo/jade_A12mo/f1_" + std::to_string(g) + ".dat", sub_population.f1_val, 1, 8);
		pepsgo::write_default1d("maps/mo/jade_A12mo/f2_" + std::to_string(g) + ".dat", sub_population.f2_val, 1, 8);

		bool opt_type = 0; // 0 - Pareto, 1 - Slater
		bool type = 0;  // s minimize = 0, l maximize = 1
		std::vector<Datapoint*> dataset; //dataset will contain N K-dim datapoints
	//		size_t idmax_file1; // for sanity check to ensure equal file lengths
	//		std::vector< std::pair<std::string,int> > files;
	//		files.push_back(std::pair<std::string,int>("/work/ProjectsCPP/pareto/f1.dat", type));
	//		files.push_back(std::pair<std::string,int>("/work/ProjectsCPP/pareto/f2.dat", type));


		for(size_t k = 0; k != sub_population.f1_val.size(); k++)
		{
			dataset.push_back(new Datapoint(k));
			dataset[k]->addNumber(sub_population.f1_val[k]);
			dataset[k]->addNumber(sub_population.f2_val[k]);
		}

		std::cout << dataset.size() << std::endl;

		//ParetoAlgo* algo = new BruteforceAlgo(opt_type, type);
		ParetoAlgo* algo = new StablesortAlgo();
		int numPareto = algo->computeFrontier(dataset);

		std::ofstream fOut;
		fOut.open("maps/mo/jade_A12mo/pareto_front_" + std::to_string(g) + ".dat");
		if(!fOut.is_open())
		{
			std::cout << "Error opening file." << std::endl;
			return 1;
		}
		fOut.precision(10);
		for(size_t k = 0; k < dataset.size(); k++)
		{
			if(dataset[k]->getdominationStatus())
			{

	//				if(obj.objective(sub_population.f1f2_vec[k]) < 1)
	//				{
	//					obj.write(sub_population.f1f2_vec[k], "maps/mo/jade_A10mo/pdb/peptide" + std::to_string(k) + ".pdb");
	//					obj.write_lbfgs(sub_population.f1f2_vec[k], "maps/mo/jade_A10mo/pdb/l/peptidel" + std::to_string(k) + ".pdb");
	//				}

				pAA_rmsd.push_back(obj.get_AA_rmsd(sub_population.f1f2_vec[k]));
				pCA_rmsd.push_back(obj.get_CA_rmsd(sub_population.f1f2_vec[k]));
				pscore_values.push_back(sub_population.score_values[k]);

				std::vector<double> tt = dataset[k]->get_vector();
				for(size_t j = 0; j < tt.size(); j++)
				{
					//std::cout << std::scientific << tt[j] << '\t';
					fOut << std::scientific << tt[j] << '\t';
				}
				fOut << std::endl;
			}
		}
		fOut.close();
	}

	pepsgo::write_default1d("maps/mo/jade_A12mo/ca_rmsd.dat", CA_rmsd, 1, 10);
	pepsgo::write_default1d("maps/mo/jade_A12mo/aa_rmsd.dat", AA_rmsd, 1, 10);
	pepsgo::write_default1d("maps/mo/jade_A12mo/score_values.dat", score_values, 1, 10);
	pepsgo::write_default1d("maps/mo/jade_A12mo/lca_rmsd.dat", lCA_rmsd, 1, 10);
	pepsgo::write_default1d("maps/mo/jade_A12mo/laa_rmsd.dat", lAA_rmsd, 1, 10);
	pepsgo::write_default1d("maps/mo/jade_A12mo/lscore_values.dat", lbfgs_score_values, 1, 10);
	pepsgo::write_default1d("maps/mo/jade_A12mo/pca_rmsd.dat", pCA_rmsd, 1, 10);
	pepsgo::write_default1d("maps/mo/jade_A12mo/paa_rmsd.dat", pAA_rmsd, 1, 10);
	pepsgo::write_default1d("maps/mo/jade_A12mo/pscore_values.dat", pscore_values, 1, 10);



}
