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
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include "dunbrackdata.hh"
#include "bbdep_sm.hh"

int main(int argc, char *argv[])
{
    devel::init(argc, argv);
    std::cout << "Start..." << std::endl;

    std::string lib_path = basic::options::option[basic::options::OptionKeys::in::path::database]().front();
    lib_path += "rotamer/ExtendedOpt1-5/";
    std::cout << lib_path << std::endl;
//
//    std::vector<std::vector<int>> num;
//    std::vector<pepsgo::bbdep::Dunbrack_data> lib;
//    pepsgo::bbdep::load_data_sm(lib_path, lib, num);
//    std::cout << lib.size() << std::endl;


    pepsgo::bbdep::BBDEP_Dunbrack_sm bbdep_sm(lib_path, 72);
    bbdep_sm.initialize_all(true, "RK");
}
