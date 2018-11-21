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
#include "data_writing.hh"

void plot_chi1_all(pepsgo::bbdep::BBDEP_Dunbrack_sm& bbdep_sm)
{
    auto ser = bbdep_sm.get_chi1_all(bbdep_sm.aa_sm_1d[0].lib);
    auto val = bbdep_sm.get_chi1_all(bbdep_sm.aa_sm_1d[1].lib);
    auto cys = bbdep_sm.get_chi1_all(bbdep_sm.aa_sm_1d[2].lib);
    auto thr = bbdep_sm.get_chi1_all(bbdep_sm.aa_sm_1d[3].lib);
//    std::cout << rez.pdf.size() << std::endl;
//    std::cout << rez.grid.size() << std::endl;
    std::vector<std::vector<double>> to_plot;
    for(size_t i = 0; i != ser.pdf.size(); i++)
    {
        std::vector<double> temp = {ser.grid[i], ser.pdf[i], val.pdf[i], cys.pdf[i], thr.pdf[i]};
        if(temp[0] < 0)
        {
            temp[0] += 360.0;
        }
        to_plot.push_back(temp);
    }
    write_default2d("maps/bbdep/SVCT.dat", to_plot, 4);
}

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

    pepsgo::bbdep::BBDEP_Dunbrack_sm bbdep_sm(lib_path, 1000);
    bbdep_sm.initialize_all(true, "SVCT");
    
//    plot_chi1_all(bbdep_sm);
}
