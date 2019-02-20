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

#include "pepsgo.hh"

int main(int argc, char *argv[])
{
    devel::init(argc, argv);
    std::cout << "Start..." << std::endl;

    size_t thread_num = 4;

    pepsgo::PEPSGO obj;
    obj.set_number_of_threads(thread_num);
    obj.set_peptide_from_file();
    obj.fill_opt_vector();
    obj.optimize_native();
//    obj.create_space_frag(std::make_pair(18, 72), std::make_pair(72, 180));
//    obj.create_space_frag(std::make_pair(8, 18), std::make_pair(36, 72));
//    obj.create_space_frag(std::make_pair(8, 12), std::make_pair(36, 72));
    obj.create_space_frag(std::make_pair(8, 10), std::make_pair(54, 72));
    obj.fill_rama2_quantile(4);
    obj.set_omega_quantile(72);
    obj.set_bbdep(360);
    obj.set_multithread();

    std::cout << obj.get_opt_vector_size() << std::endl;
    std::function<core::Real(const std::vector<double>&)> f = std::bind(&pepsgo::PEPSGO::objective, obj, std::placeholders::_1);
    std::function<core::Real(const std::vector<double>&, int)> f_mt = std::bind(&pepsgo::PEPSGO::objective_mt, obj, std::placeholders::_1, std::placeholders::_2);
    std::cout << f(std::vector<core::Real>(obj.get_opt_vector_size(),0.4)) << std::endl;
    std::cout << f_mt(std::vector<core::Real>(obj.get_opt_vector_size(),0.4),thread_num-1) << std::endl;

    std::cout << obj.get_AA_rmsd(std::vector<core::Real>(obj.get_opt_vector_size(),0.4)) << std::endl;
    std::cout << obj.get_CA_rmsd(std::vector<core::Real>(obj.get_opt_vector_size(),0.4)) << std::endl;

    obj.append_to_superposed_aa(std::vector<core::Real>(obj.get_opt_vector_size(),0.4));
    obj.append_to_superposed_ca(std::vector<core::Real>(obj.get_opt_vector_size(),0.4));
    obj.dumb_superposed();
}
