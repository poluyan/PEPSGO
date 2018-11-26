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

//#include "dunbrackdata.hh"
//#include "bbdep_sm.hh"
#include "data_io.hh"
#include "bbind.hh"

#include <unordered_map>

int main(int argc, char *argv[])
{
    devel::init(argc, argv);
    std::cout << "Start..." << std::endl;
//
//    std::string lib_path = basic::options::option[basic::options::OptionKeys::in::path::database]().front();
//    lib_path += "rotamer/ExtendedOpt1-5/";
//    std::cout << lib_path << std::endl;
//
//    pepsgo::bbdep::BBDEP_Dunbrack_sm bbdep_sm(lib_path, 1000);
//    bbdep_sm.initialize_all(true, "SVCT");
//
//    plot_chi1_all(bbdep_sm);

    std::unordered_map<core::chemical::AA, std::string> wordMap;

    wordMap.insert({ core::chemical::aa_ala, "ala" });
    wordMap.insert({ core::chemical::aa_ser, "bbb" });


    for(auto element : wordMap)
        std::cout << element.first << " :: " << element.second << std::endl;

    pepsgo::bbind::bbind_top obj("/ssdwork/ProjectsCPP/mcmc/bbind_top500/");
    obj.initialize_all(2,100,100,100,"S");

 /*   auto temp = std::any_cast<std::shared_ptr<linterp::InterpMultilinear<1, double>>>(
                    obj.aa_data[core::chemical::aa_ser]);
                    
    auto dist = std::any_cast<pepsgo::bbutils::distribution_1d>(
                    obj.aa_dst[core::chemical::aa_ser]);
                    
    auto pdf = dist.pdf;
    
    pepsgo::write_default1d("maps/bbind/ser_pdf.dat", pdf, 1, 4);


    for(size_t i = 0; i != 360; i++)
    {
        boost::array<double, 1> args = { double(i) };
//        std::cout << temp->interp(args.begin()) << std::endl;
    }*/

//    std::cout << std::get<2>(obj.aa_range[core::chemical::aa_ser].front()) << std::endl;

//    std::unordered_map<core::chemical::AA, std::any> abc;
//    std::unique_ptr<linterp::InterpMultilinear<1, double>> d;
//    abc.insert( { core::chemical::aa_ala, d });

//	// Insert Few elements in map
//	wordMap.insert( { "First", 1 });
//	wordMap.insert(	{ "Second", 2 });
//	wordMap.insert(	{ "Third", 3 });
//
//	// Overwrite value of an element
//	wordMap["Third"] = 8;
//
//    std::cout << int(core::chemical::aa_ala) << std::endl;
//    std::cout << int(core::chemical::aa_ala) << std::endl;

}
