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
#ifndef INCLUDED_transform_hh
#define INCLUDED_transform_hh

#include <tuple>

namespace pepsgo
{
    
namespace transform
{

struct ranges
{
    std::tuple<bool, size_t, size_t> phipsi;
    std::tuple<bool, size_t, size_t> omega;
    std::tuple<bool, size_t, size_t> chi;

    bool do_phipsi;
    bool do_omega;
    bool do_chi;
};

}
}

#endif