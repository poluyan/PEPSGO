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

#include <iostream>

namespace pepsgo
{
namespace fragment
{

struct Frag
{
    std::size_t id;
    char aa;
//    core::Real phi;
//    core::Real psi;
//    core::Real omg;
    Frag();
};

}
}
#endif