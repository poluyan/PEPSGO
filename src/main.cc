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

#include <basic/options/keys/in.OptionKeys.gen.hh> // Dunbrack lib path
#include <basic/options/option.hh>

//#include "dunbrackdata.hh"
//#include "bbdep_sm.hh"
//#include "data_io.hh"
//#include "bbind.hh"

#include "pepsgo.hh"

int main(int argc, char *argv[])
{
    devel::init(argc, argv);
    std::cout << "Start..." << std::endl;

    pepsgo::PEPSGO obj;
    obj.set_number_of_threads(4);
//    obj.set_peptide("KTWNPATGKWTE");
//    obj.set_peptide("DPCYEVCLQQHGNVKECEEACKHPVE");
    obj.set_peptide_from_file();
    obj.fill_rama2_quantile(10);
    obj.set_bbdep();
    obj.fill_opt_vector();
    obj.create_space_frag();

    std::cout << obj.get_opt_vector_size() << std::endl;
    std::function<core::Real(const std::vector<double>&)> f = std::bind(&pepsgo::PEPSGO::objective, obj, std::placeholders::_1);
    std::cout << f(std::vector<core::Real>(obj.get_opt_vector_size(),0.4)) << std::endl;
}
