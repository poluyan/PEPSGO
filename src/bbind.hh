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
#include <core/chemical/AA.hh>

#include "bbutils.hh"
#include "linterp.hh"

#include <unordered_map>
#include <any>

namespace pepsgo
{
namespace bbind
{

class bbind_top
{
public:
    std::string path_to_files;
    std::unordered_map<core::chemical::AA, std::any> aa_data;
    std::unordered_map<core::chemical::AA, std::any> aa_dst;
    std::unordered_map<core::chemical::AA, std::tuple<double, double, size_t>> aa_range;

    bbind_top(std::string _path_to_files);

    bbutils::distribution_1d make_1d_cdf(const linterp::InterpMultilinear<1, double> &acid,
                                         const std::vector<std::tuple<double, double, size_t>> &range,
                                         size_t m);

    double get_1d(double x,
                  const linterp::InterpMultilinear<1, double> &f,
                  const std::vector<std::tuple<double, double, size_t>> &range);
};

}
}
