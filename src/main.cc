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

#include <map>

int main(int argc, char *argv[])
{
    devel::init(argc, argv);
    std::cout << "Start..." << std::endl;

    pepsgo::PEPSGO obj;
    obj.set_number_of_threads(4);
    obj.set_peptide("KTWNPATGKWTE");
    obj.fill_rama2_sample();

    //std::function<core::Real(const std::vector<double>&)> f_display = obj.objective;
}
