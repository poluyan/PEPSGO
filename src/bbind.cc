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
#include <fstream>
#include "data_io.hh"

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
    result.pdf = pdf;
    result.cdf = cdf;
    result.grid = x;
    return result;
}

void bbind_top::initialize_all(size_t chi1_step,
                               size_t chi2_step,
                               size_t chi3_step,
                               size_t chi4_step,
                               std::string amino_acids)
{
    std::cout << "\n1d progress";

    if(amino_acids.find("S") != std::string::npos)
    {
        //std::unique_ptr<linterp::InterpMultilinear<1, double>> d;
        std::vector<std::tuple<double, double, size_t>> t;
        linterp::InterpMultilinear<1, double> data;
        std::cout << &data << std::endl;
        load_1d("rota500-", "ser", t, data);
        std::cout << &data << std::endl;




        /*    size_t num_elements;
            std::vector<std::vector<double>> grids;
            boost::array<int, 1> grid_sizes;
            std::vector<double> f_values;
            std::vector<std::tuple<double, double, size_t>> tr;
            fill_1d(path_to_files + "rota500-" + "ser" + ".data", tr, grids, num_elements, grid_sizes, f_values);

            pepsgo::write_default1d("maps/bbind/f_values.dat", f_values, 1, 10);

            auto grid_iter_list = linterp::get_begins_ends(grids.begin(), grids.end());

            //linterp::InterpMultilinear<1, double> dd(grid_iter_list.first.begin(), grid_sizes.begin(), f_values.data(), f_values.data() + num_elements);
            linterp::InterpMultilinear<1, double> dd;
            dd.init(grid_iter_list.first.begin(), grid_sizes.begin(), f_values.data(), f_values.data() + num_elements);

        */
        std::vector<double> tt;
        for(double i = 0; i < 360; i+=0.01)
        {
            boost::array<double, 1> args = { i };
            tt.push_back(data.interp(args.begin()));
        }
        pepsgo::write_default1d("maps/bbind/f_values3.dat", tt, 1, 10);






        /*

        bbutils::distribution_1d dst = make_1d_cdf(d, t, chi1_step);

        aa_dst.insert({ core::chemical::aa_ser, dst });
        //        aa_data.insert({ core::chemical::aa_ser, d });
        aa_range.insert({ core::chemical::aa_ser, t });

        //std::cout << std::get<2>(aa_range[core::chemical::aa_ser].front()) << std::endl;
        std::cout << dst.pdf.size() << std::endl;

        std::cout << ".";*/
    }

}


void bbind_top::load_1d(std::string prefix,
                        std::string acid_name,
                        std::vector<std::tuple<double, double, size_t>> &range,
                        linterp::InterpMultilinear<1, double> &res)
{
    size_t num_elements;
    std::vector<std::vector<double>> grids;
    boost::array<int, 1> grid_sizes;
    std::vector<double> f_values;

    fill_1d(path_to_files + prefix + acid_name + ".data", range, grids, num_elements, grid_sizes, f_values);

    pepsgo::write_default1d("maps/bbind/f_values.dat", f_values, 1, 10);

    auto grid_iter_list = linterp::get_begins_ends(grids.begin(), grids.end());

    res.init(grid_iter_list.first.begin(), grid_sizes.begin(), f_values.data(), f_values.data() + num_elements);

    std::vector<double> tt;
    for(double i = 0; i < 360; i+=0.01)
    {
        boost::array<double, 1> args = { i };
        tt.push_back(res.interp(args.begin()));
    }
    pepsgo::write_default1d("maps/bbind/f_values2.dat", tt, 1, 10);
    std::cout << &res << std::endl;
}

void bbind_top::load_data(std::string fname,
                          std::vector<std::vector<double>> &data,
                          std::vector<std::tuple<double, double, size_t>> &range)
{
    std::ifstream fIn;
    fIn.open(fname.c_str());
    if(!fIn.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }
    size_t dimension = 0;
    std::string temp;
    std::vector<double> push2data;
    bool flag = false;
    std::tuple<double, double, size_t> t;
    while(!fIn.eof())
    {
        if(flag)
        {
            push2data.clear();
            double float_temp = 0.0;
            for(size_t i = 0; i < dimension; i++)
            {
                fIn >> float_temp;
                push2data.push_back(float_temp);
            }
            data.push_back(push2data);
            continue;
        }
        fIn >> temp;
        if("dimensions:" == temp)
        {
            fIn >> dimension;
            dimension++;
        }
        if("x1:" == temp)
        {
            fIn >> std::get<0>(t);
            fIn >> std::get<1>(t);
            fIn >> std::get<2>(t);
            range.push_back(t);
        }
        if("x2:" == temp)
        {
            fIn >> std::get<0>(t);
            fIn >> std::get<1>(t);
            fIn >> std::get<2>(t);
            range.push_back(t);
        }
        if("x3:" == temp)
        {
            fIn >> std::get<0>(t);
            fIn >> std::get<1>(t);
            fIn >> std::get<2>(t);
            range.push_back(t);
        }
        if("x4:" == temp)
        {
            fIn >> std::get<0>(t);
            fIn >> std::get<1>(t);
            fIn >> std::get<2>(t);
            range.push_back(t);
        }
        if("line.)" == temp)
        {
            flag = true;
        }
    }
    fIn.close();
    data.pop_back();
}

void bbind_top::fill_1d(std::string name,
                        std::vector<std::tuple<double, double, size_t>> &range,
                        std::vector<std::vector<double>> &grids,
                        size_t &num_elements,
                        boost::array<int, 1> &grid_sizes,
                        std::vector<double> &f_values)
{

    std::vector<std::vector<double>> data;
    load_data(name, data, range);

    if(range.size() != 1)
    {
        std::cout << "Not a 1d file!" << '\t' << name << std::endl;
        return;
    }

    for(size_t d = 0; d < range.size(); d++)
    {
        std::vector<double> grid1;
        double es = std::get<1>(range[d]) - std::get<0>(range[d]);
        int total = std::get<2>(range[d]);
        int total2 = 2 * total;
        double start_value = data[0][d];
        for(int i = -1; i != total; i++)
        {
            grid1.push_back(es * (1 + 2 * i) / total2);
        }
        grid1.push_back(es * (1 + 2 * total) / total2);
        grids.push_back(grid1);
    }

    for(size_t d = 0; d != grid_sizes.size(); d++)
    {
        grid_sizes[d] = grids[d].size();
    }

    num_elements = 1;
    for(const auto &d : range)
    {
        num_elements *= std::get<2>(d);
    }
    num_elements += range.size() * 2;

    f_values.resize(num_elements);

    f_values.front() = data.back().back();
    for(size_t i = 0; i != data.size(); i++)
    {
        f_values[i + 1] = data[i].back();
    }
    f_values.back() = data.front().back();
}

}
}
