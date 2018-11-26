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


BBIND_top::BBIND_top(std::string _path_to_files)
    : path_to_files(_path_to_files)
{
}

double BBIND_top::get_1d(double x,
                         core::chemical::AA acid)
{
    if(x < std::get<0>(aa_range[acid].front()) || x > std::get<1>(aa_range[acid].front()))
    {
        return 0;
    }
    boost::array<double, 1> args = { x };
    auto f = std::any_cast<std::shared_ptr<linterp::InterpMultilinear<1, double>>>(aa_data[acid]);
    return f->interp(args.begin());
}

bbutils::distribution_1d BBIND_top::make_1d_cdf(core::chemical::AA acid, size_t m)
{
    double start = std::get<0>(aa_range[acid].front()), stop = std::get<1>(aa_range[acid].front()), es = stop - start;
    size_t step = m * std::get<2>(aa_range[acid].front());

    std::deque<double> pdf, x;
    for(size_t i = 0; i != step; i++)
    {
        x.push_back(start + i * es / step);
        pdf.push_back(get_1d(start + i * es / step, acid));
    }
    double S = std::accumulate(pdf.begin(), pdf.end(), 0.0);
    for(auto &a : pdf)
    {
        a /= S;
    }

    std::deque<double> cdf(pdf.size());
    std::partial_sum(pdf.begin(), pdf.end(), cdf.begin(), std::plus<double>());
    cdf.push_front(0);
    x.push_back(std::get<1>(aa_range[acid].front()));

    bbutils::distribution_1d result;
    result.pdf = pdf;
    result.cdf = cdf;
    result.grid = x;
    return result;
}

void BBIND_top::initialize_all(size_t chi1_step,
                               size_t chi2_step,
                               size_t chi3_step,
                               size_t chi4_step,
                               std::string amino_acids)
{
    std::cout << "\n1d progress";

    if(amino_acids.find("S") != std::string::npos)
    {
        std::vector<std::tuple<double, double, size_t>> r;
        linterp::InterpMultilinear<1, double> d;
        load_1d("rota500-", "ser", core::chemical::aa_ser);
        bbutils::distribution_1d dst = make_1d_cdf(core::chemical::aa_ser, chi1_step);
        aa_dst.insert({ core::chemical::aa_ser, dst });
        std::cout << ".";
    }
    if(amino_acids.find("V") != std::string::npos)
    {
        std::vector<std::tuple<double, double, size_t>> r;
        linterp::InterpMultilinear<1, double> d;
        load_1d("rota500-", "val", core::chemical::aa_val);
        bbutils::distribution_1d dst = make_1d_cdf(core::chemical::aa_val, chi1_step);
        aa_dst.insert({ core::chemical::aa_val, dst });
        std::cout << ".";
    }
    if(amino_acids.find("C") != std::string::npos)
    {
        std::vector<std::tuple<double, double, size_t>> r;
        linterp::InterpMultilinear<1, double> d;
        load_1d("rota500-", "cys", core::chemical::aa_cys);
        bbutils::distribution_1d dst = make_1d_cdf(core::chemical::aa_cys, chi1_step);
        aa_dst.insert({ core::chemical::aa_cys, dst });
        std::cout << ".";
    }
    if(amino_acids.find("T") != std::string::npos)
    {
        std::vector<std::tuple<double, double, size_t>> r;
        linterp::InterpMultilinear<1, double> d;
        load_1d("rota500-", "thr", core::chemical::aa_thr);
        bbutils::distribution_1d dst = make_1d_cdf(core::chemical::aa_thr, chi1_step);
        aa_dst.insert({ core::chemical::aa_thr, dst });
        std::cout << ".";
    }

}


void BBIND_top::load_1d(std::string prefix,
                        std::string acid_name, core::chemical::AA acid)
{
    size_t num_elements;
    std::vector<std::vector<double>> grids;
    boost::array<int, 1> grid_sizes;
    std::vector<double> f_values;
    std::vector<std::tuple<double, double, size_t>> range;
    fill_1d(path_to_files + prefix + acid_name + ".data", range, grids, num_elements, grid_sizes, f_values);
    auto grid_iter_list = linterp::get_begins_ends(grids.begin(), grids.end());
    auto res =  std::make_shared<linterp::InterpMultilinear<1, double>>(grid_iter_list.first.begin(), grid_sizes.begin(), f_values.data(), f_values.data() + num_elements);
    aa_data.insert({ acid, std::move(res) });
    aa_range.insert({ acid, range });
}

void BBIND_top::load_data(std::string fname,
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

void BBIND_top::fill_1d(std::string name,
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


void plot_chi1_all(pepsgo::bbind::BBIND_top& obj)
{
    auto s = std::any_cast<std::shared_ptr<linterp::InterpMultilinear<1, double>>>(
                    obj.aa_data[core::chemical::aa_ser]);
    auto v = std::any_cast<std::shared_ptr<linterp::InterpMultilinear<1, double>>>(
                    obj.aa_data[core::chemical::aa_val]);
    auto c = std::any_cast<std::shared_ptr<linterp::InterpMultilinear<1, double>>>(
                    obj.aa_data[core::chemical::aa_cys]);
    auto t = std::any_cast<std::shared_ptr<linterp::InterpMultilinear<1, double>>>(
                    obj.aa_data[core::chemical::aa_thr]);

    std::vector<std::vector<double>> tt;
    for(double i = 0; i < 360; i+=0.1)
    {
        boost::array<double, 1> args = { i };
        tt.push_back(std::vector<double>{i, s->interp(args.begin()),v->interp(args.begin()),c->interp(args.begin()),t->interp(args.begin())});
    }
    pepsgo::write_default2d("maps/bbind/SVCT.dat", tt, 10);
}

}
}
