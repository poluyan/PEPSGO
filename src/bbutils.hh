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
#ifndef INCLUDED_bbutils_hh
#define INCLUDED_bbutils_hh

#include <iterator>
#include <deque>
#include <vector>
#include <algorithm>

namespace pepsgo
{
namespace bbutils
{

struct distribution_1d
{
    std::deque<double> pdf;
    std::deque<double> cdf;
    std::deque<double> grid;
};

double get_1d_from_dst(const distribution_1d &dst, double value);
double get_inverse_1d_from_dst(const distribution_1d &dst, double value);

}
}
#endif
