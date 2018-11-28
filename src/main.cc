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
    obj.set_peptide("KTWNPATGKWTE");

    /*std::string lib_path = basic::options::option[basic::options::OptionKeys::in::path::database].value_string() + "rotamer/ExtendedOpt1-5/";
    std::cout << lib_path << std::endl;
    pepsgo::bbdep::BBDEP_Dunbrack_sm bbdep_sm;
    bbdep_sm.set_path(lib_path);
    bbdep_sm.set_step(1000);
    bbdep_sm.initialize_all(true, "SVCT");

    pepsgo::bbdep::plot_chi1_all(bbdep_sm);

    pepsgo::bbind::BBIND_top obj;
    obj.set_path("/ssdwork/ProjectsCPP/mcmc/bbind_top500/");*/
    
//    obj.initialize(2, 2, 2, 2,"SVCT");
//    obj.initialize(2, 2, 2, 2,"WHNDFYIL");
//    obj.initialize(2, 2, 2, 2,"MEQP");
//    obj.initialize(2, 2, 2, 2,"RK");

    
}
