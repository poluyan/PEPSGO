#include <pepsgo.hh>

struct Example
{
	std::function<core::Real(const std::vector<double>&)> Function;
	std::function<core::Real(const std::vector<double>&, int)> FunctionMultiThread;
};

int main(int argc, char *argv[])
{
	size_t task_number = 7;     // номер предустановленных параметров
	size_t thread_number = 4;   // количество потоков

	pepsgo::PEPSGO obj(argc, argv);  
	obj.set_task(task_number);       // указание предустановленныйх параметров по номеру

	std::cout << obj.get_problem_dimension() << std::endl; // размер задачи
	std::vector<double> test_vector(obj.get_problem_dimension(), 0.4);
	
	std::function<core::Real(const std::vector<double>&)> f = std::bind(&pepsgo::PEPSGO::objective, obj, std::placeholders::_1);
	std::cout << f(test_vector) << std::endl; // можно double, можно core::Real

	obj.set_number_of_threads(thread_number);  // указание количества потоков для получения многопоточной функции

	std::function<core::Real(const std::vector<double>&, int)> f_multithread = std::bind(&pepsgo::PEPSGO::objective_mt, obj, std::placeholders::_1, std::placeholders::_2);
	std::cout << f_multithread(test_vector, 0) << std::endl;
	std::cout << f_multithread(test_vector, thread_number - 1) << std::endl;

	std::cout << obj.get_AA_rmsd(test_vector) << std::endl; // RMSD отклонение по всем атомам
	std::cout << obj.get_CA_rmsd(test_vector) << std::endl; // RMSD отклонение по атомам главной цепи
	
	obj.write(test_vector, "output/pdb/test.pdb");	         	 // запись структуры в указанный файл
	obj.write_lbfgs(test_vector, "output/pdb/test_lbfgs.pdb"); // локальная оптимизация и запись структуры в указанный файл 

	obj.append_to_native_and_superpose_AA(test_vector);   // добавить к нативной структуре с суперпозицией по всем атомам
	obj.append_to_native_and_superpose_CA(test_vector);   // добавить к нативной структуре с суперпозицией по атомам главной цепи
	obj.dumb_superposed();	                              // записать в output/pdb/superposed.pdb
	
	// пример использования в независимом классе
	Example fObj;
	fObj.Function = std::bind(&pepsgo::PEPSGO::objective, obj, std::placeholders::_1);
	fObj.FunctionMultiThread = std::bind(&pepsgo::PEPSGO::objective_mt, obj, std::placeholders::_1, std::placeholders::_2);
}

/*
#include <pepsgo.hh> /// также требуется динамическая библиотека libpepsgo.so

int main(int argc, char *argv[]) /// в аргументах файл с последовательностью
{
	pepsgo::PEPSGO obj(argc, argv);  /// передача аргументов в объект комплекса 
	obj.set_task(7);       /// выбор предустановленных параметров по номеру [0-8]

	size_t d = obj.get_problem_dimension();  /// размер задачи
	
	/// целевая функция objective, objective_mt, objective_mo, objective_mo_mt
	std::function<core::Real(const std::vector<double>&)> f = 
	std::bind(&pepsgo::PEPSGO::objective, obj, std::placeholders::_1); 
	
	std::vector<double> test_vector(d, 0.4); /// вектор значений
	double score_value = f(test_vector);     /// получение score значения 

	obj.set_number_of_threads(4); /// установка числа потоков, objective_mt(x, th_id)
}*/
