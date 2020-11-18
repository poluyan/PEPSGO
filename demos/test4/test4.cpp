#include <devel/init.hh>
#include <pepsgo.hh>

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
	obj.set_task(4);
	obj.fill_opt_vector();
	obj.optimize_native();
	obj.set_multithread();

	std::cout << obj.get_opt_vector_size() << std::endl;
	std::function<core::Real(const std::vector<double>&)> f = std::bind(&pepsgo::PEPSGO::objective, obj, std::placeholders::_1);
	std::function<core::Real(const std::vector<double>&, int)> f_mt = std::bind(&pepsgo::PEPSGO::objective_mt, obj, std::placeholders::_1, std::placeholders::_2);
	std::cout << f(std::vector<core::Real>(obj.get_opt_vector_size(),0.4)) << std::endl;
	std::cout << f_mt(std::vector<core::Real>(obj.get_opt_vector_size(),0.4),thread_num-1) << std::endl;

	std::cout << obj.get_AA_rmsd(std::vector<core::Real>(obj.get_opt_vector_size(),0.4)) << std::endl;
	std::cout << obj.get_CA_rmsd(std::vector<core::Real>(obj.get_opt_vector_size(),0.4)) << std::endl;

	Example fObj;
	fObj.FitnessFunction = std::bind(&pepsgo::PEPSGO::objective, obj, std::placeholders::_1);
	fObj.FitnessFunctionMultiThread = std::bind(&pepsgo::PEPSGO::objective_mt, obj, std::placeholders::_1, std::placeholders::_2);

	obj.append_to_superposed_aa(std::vector<core::Real>(obj.get_opt_vector_size(),0.4));
	obj.append_to_superposed_ca(std::vector<core::Real>(obj.get_opt_vector_size(),0.4));
	obj.dumb_superposed();
}
