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
#include "bbind.hh"

namespace pepsgo
{
namespace bbind
{


bbind_top::bbind_top(std::string _path_to_files)
    : path_to_files(_path_to_files)
{
}

double bbind_top::get_1d(double x,
                         const linterp::InterpMultilinear<1, double> &f,
                         const std::vector<std::tuple<double, double, size_t>> &range)
{
    if(x < std::get<0>(range[0]) || x > std::get<1>(range[0]))
    {
        return 0;
    }
    boost::array<double, 1> args = { x };
    return f.interp(args.begin());
}

bbutils::distribution_1d bbind_top::make_1d_cdf(const linterp::InterpMultilinear<1, double> &acid,
        const std::vector<std::tuple<double, double, size_t>> &range,
        size_t m)
{
    double start = std::get<0>(range[0]), stop = std::get<1>(range[0]), es = stop - start;
    size_t step = m * std::get<2>(range[0]);

    std::deque<double> pdf, x;
    for(size_t i = 0; i != step; i++)
    {
        x.push_back(start + i * es / step);
        pdf.push_back(get_1d(start + i * es / step, acid, range));
    }
    double S = std::accumulate(pdf.begin(), pdf.end(), 0.0);
    for(auto &a : pdf)
    {
        a /= S;
    }

    std::deque<double> cdf(pdf.size());
    std::partial_sum(pdf.begin(), pdf.end(), cdf.begin(), std::plus<double>());
    cdf.push_front(0);
    x.push_back(std::get<1>(range[0]));

    bbutils::distribution_1d result;
    // result.pdf = pdf;
    result.cdf = cdf;
    result.grid = x;
    return result;
}


}
}
