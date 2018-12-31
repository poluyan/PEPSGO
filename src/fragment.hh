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
#ifndef INCLUDED_fragment_hh
#define INCLUDED_fragment_hh

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <core/sequence/util.hh>

#include <string>
#include <vector>

namespace pepsgo
{
namespace fragment
{

struct Frag
{
    char aa;
    core::Real phi;
    core::Real psi;
    core::Real omg;
    
    Frag();
    Frag(char _aa, core::Real _phi, core::Real _psi, core::Real _omg);
};

class FragPick
{
private:
    //std::string psipred_fname;
    std::string fragmet_file;
    std::string fasta;
    std::string peptide_seq;
};

}
}
#endif