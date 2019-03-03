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
#include <bbind.hh>
#include <data_io.hh>
#include <linterp.hh>

#include <fstream>

namespace pepsgo
{
namespace bbind
{


BBIND_top::BBIND_top(){}

void BBIND_top::set_path(std::string path_to_files)
{
    path = path_to_files;
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
    //result.pdf = pdf;
    result.cdf = cdf;
    result.grid = x;
    return result;
}

bbutils::distribution_2d BBIND_top::make_2d_cdf(core::chemical::AA acid, size_t m)
{
    std::vector<double> es;
    std::vector<size_t> step;
    for(size_t i = 0; i != aa_range[acid].size(); i++)
    {
        es.push_back(std::get<1>(aa_range[acid][i]) - std::get<0>(aa_range[acid][i]));
        step.push_back(m * std::get<2>(aa_range[acid][i]));
    }

    std::deque<std::deque<double>> pdf;
    std::deque<double> x, y;
    for(size_t i = 0; i != step[0]; i++)
        x.push_back(std::get<0>(aa_range[acid][0]) + i * es[0] / step[0]);
    for(size_t i = 0; i != step[1]; i++)
        y.push_back(std::get<0>(aa_range[acid][1]) + i * es[1] / step[1]);

    x.push_back(std::get<1>(aa_range[acid][0]));
    y.push_back(std::get<1>(aa_range[acid][1]));

    auto f = std::any_cast<std::shared_ptr<linterp::InterpMultilinear<2, double>>>(aa_data[acid]);
    for(size_t i = 0; i != step[0]; i++)
    {
        std::deque<double> temp;
        for(size_t j = 0; j != step[1]; j++)
        {
            boost::array<double, 2> args = { std::get<0>(aa_range[acid][0]) + i * es[0] / step[0],
                                             std::get<0>(aa_range[acid][1]) + j * es[1] / step[1]
                                           };
            double value = f->interp(args.begin());
            temp.push_back(value);
        }
        pdf.push_back(temp);
    }

    // cdf

    std::deque<std::deque<double>> cdf(pdf.size(), std::deque<double>(pdf[0].size()));

    for(size_t i = 0; i != pdf.size(); i++)
    {
        std::partial_sum(pdf[i].begin(), pdf[i].end(), cdf[i].begin(), std::plus<double>());
    }

    std::deque<double> col_pdf;
    for(size_t i = 0; i != cdf.size(); i++)
    {
        col_pdf.push_back(cdf[i].back());
    }

    double S = std::accumulate(col_pdf.begin(), col_pdf.end(), 0.0);
    for(auto &a : col_pdf)
    {
        a /= S;
    }

    std::deque<double> col_cdf(col_pdf.size());
    std::partial_sum(col_pdf.begin(), col_pdf.end(), col_cdf.begin(), std::plus<double>());
    col_cdf.push_front(0);

    std::deque<std::deque<double>> row_dist_bank(pdf.size());
    for(size_t i = 0; i != row_dist_bank.size(); i++)
    {
        row_dist_bank[i] = pdf[i];

        double S = std::accumulate(row_dist_bank[i].begin(), row_dist_bank[i].end(), 0.0);
        for(auto &a : row_dist_bank[i])
        {
            a /= S;
        }

        std::partial_sum(
            row_dist_bank[i].begin(), row_dist_bank[i].end(), row_dist_bank[i].begin(), std::plus<double>());
        row_dist_bank[i].push_front(0);
    }

    bbutils::distribution_2d result;
    // result.pdf = pdf;

    result.row_dist = row_dist_bank;
    result.col_dist = col_cdf;

    result.grid.push_back(x);
    result.grid.push_back(y);
    return result;
}

bbutils::distribution_3d BBIND_top::make_3d_cdf(core::chemical::AA acid, size_t m)
{
    std::vector<double> es;
    std::vector<size_t> step;
    for(size_t i = 0; i != aa_range[acid].size(); i++)
    {
        es.push_back(std::get<1>(aa_range[acid][i]) - std::get<0>(aa_range[acid][i]));
        step.push_back(m * std::get<2>(aa_range[acid][i]));
    }

    std::deque<std::deque<std::deque<double>>> pdf;
    std::deque<double> x, y, z;
    for(size_t i = 0; i != step[0]; i++)
        x.push_back(std::get<0>(aa_range[acid][0]) + i * es[0] / step[0]);
    for(size_t i = 0; i != step[1]; i++)
        y.push_back(std::get<0>(aa_range[acid][1]) + i * es[1] / step[1]);
    for(size_t i = 0; i != step[2]; i++)
        z.push_back(std::get<0>(aa_range[acid][2]) + i * es[2] / step[2]);

    x.push_back(std::get<1>(aa_range[acid][0]));
    y.push_back(std::get<1>(aa_range[acid][1]));
    z.push_back(std::get<1>(aa_range[acid][2]));

    auto f = std::any_cast<std::shared_ptr<linterp::InterpMultilinear<3, double>>>(aa_data[acid]);
    for(size_t i = 0; i != step[0]; i++)
    {
        std::deque<std::deque<double>> temp1;
        for(size_t j = 0; j != step[1]; j++)
        {
            std::deque<double> temp2;
            for(size_t k = 0; k != step[2]; k++)
            {
                boost::array<double, 3> args = { std::get<0>(aa_range[acid][0]) + i * es[0] / step[0],
                                                 std::get<0>(aa_range[acid][1]) + j * es[1] / step[1], std::get<0>(aa_range[acid][2]) + k * es[2] / step[2]
                                               };
                double value = f->interp(args.begin());
                temp2.push_back(value);
            }
            temp1.push_back(temp2);
        }
        pdf.push_back(temp1);
    }

    // std::cout << pdf.size() << '\t' << pdf[0].size() << '\t' << pdf[0][0].size() << std::endl;
    // col cdf

    std::deque<std::deque<double>> first_2d;
    for(size_t i = 0; i != pdf.size(); i++)
    {
        std::deque<double> temp1;
        for(size_t j = 0; j != pdf[i].size(); j++)
        {
            double t = 0.0;
            for(size_t k = 0; k != pdf[i][j].size(); k++)
            {
                t += pdf[i][j][k];
            }
            temp1.push_back(t);
        }
        first_2d.push_back(temp1);
    }

    std::deque<std::deque<double>> first_2d_temp = first_2d;
    for(size_t i = 0; i != pdf.size(); i++)
    {
        std::partial_sum(
            first_2d_temp[i].begin(), first_2d_temp[i].end(), first_2d_temp[i].begin(), std::plus<double>());
    }
    std::deque<double> col_pdf;
    for(size_t i = 0; i != pdf.size(); i++)
    {
        col_pdf.push_back(first_2d_temp[i].back());
    }

    double S = std::accumulate(col_pdf.begin(), col_pdf.end(), 0.0);
    for(auto &a : col_pdf)
    {
        a /= S;
    }

    std::deque<double> col_cdf(col_pdf.size());
    std::partial_sum(col_pdf.begin(), col_pdf.end(), col_cdf.begin(), std::plus<double>());
    col_cdf.push_front(0);

    std::deque<std::deque<double>> row_dist_bank1(pdf.size());
    for(size_t i = 0; i != row_dist_bank1.size(); i++)
    {
        row_dist_bank1[i] = first_2d[i];

        double S = std::accumulate(row_dist_bank1[i].begin(), row_dist_bank1[i].end(), 0.0);
        for(auto &a : row_dist_bank1[i])
        {
            a /= S;
        }

        std::partial_sum(
            row_dist_bank1[i].begin(), row_dist_bank1[i].end(), row_dist_bank1[i].begin(), std::plus<double>());
        row_dist_bank1[i].push_front(0);
    }

    std::deque<std::deque<std::deque<double>>> row_dist_bank2(
        pdf.size(), std::deque<std::deque<double>>(pdf[0].size()));

    for(size_t i = 0; i != row_dist_bank2.size(); i++)
    {
        for(size_t j = 0; j != row_dist_bank2[i].size(); j++)
        {
            row_dist_bank2[i][j] = pdf[i][j];

            double S = std::accumulate(row_dist_bank2[i][j].begin(), row_dist_bank2[i][j].end(), 0.0);
            for(auto &a : row_dist_bank2[i][j])
            {
                a /= S;
            }

            std::partial_sum(row_dist_bank2[i][j].begin(), row_dist_bank2[i][j].end(), row_dist_bank2[i][j].begin(),
                             std::plus<double>());
            row_dist_bank2[i][j].push_front(0);
        }
    }

    bbutils::distribution_3d result;
    // result.pdf = pdf;

    result.row_dist1 = row_dist_bank1;
    result.col_dist = col_cdf;

    result.row_dist2 = row_dist_bank2;

    result.grid.push_back(x);
    result.grid.push_back(y);
    result.grid.push_back(z);

    return result;
}

bbutils::distribution_4d BBIND_top::make_4d_cdf(core::chemical::AA acid, size_t m)
{
    std::vector<double> es;
    std::vector<size_t> step;
    for(size_t i = 0; i != aa_range[acid].size(); i++)
    {
        es.push_back(std::get<1>(aa_range[acid][i]) - std::get<0>(aa_range[acid][i]));
        step.push_back(m * std::get<2>(aa_range[acid][i]));
    }

    std::deque<std::deque<std::deque<std::deque<double>>>> pdf;
    std::deque<double> x, y, z, q;
    for(size_t i = 0; i != step[0]; i++)
        x.push_back(std::get<0>(aa_range[acid][0]) + i * es[0] / step[0]);
    for(size_t i = 0; i != step[1]; i++)
        y.push_back(std::get<0>(aa_range[acid][1]) + i * es[1] / step[1]);
    for(size_t i = 0; i != step[2]; i++)
        z.push_back(std::get<0>(aa_range[acid][2]) + i * es[2] / step[2]);
    for(size_t i = 0; i != step[3]; i++)
        q.push_back(std::get<0>(aa_range[acid][3]) + i * es[3] / step[3]);

    x.push_back(std::get<1>(aa_range[acid][0]));
    y.push_back(std::get<1>(aa_range[acid][1]));
    z.push_back(std::get<1>(aa_range[acid][2]));
    q.push_back(std::get<1>(aa_range[acid][3]));

    auto f = std::any_cast<std::shared_ptr<linterp::InterpMultilinear<4, double>>>(aa_data[acid]);
    for(size_t i = 0; i != step[0]; i++)
    {
        std::deque<std::deque<std::deque<double>>> temp1;
        for(size_t j = 0; j != step[1]; j++)
        {
            std::deque<std::deque<double>> temp2;
            for(size_t k = 0; k != step[2]; k++)
            {
                std::deque<double> temp3;
                for(size_t t = 0; t != step[3]; t++)
                {
                    boost::array<double, 4> args = { std::get<0>(aa_range[acid][0]) + i * es[0] / step[0],
                                                     std::get<0>(aa_range[acid][1]) + j * es[1] / step[1], std::get<0>(aa_range[acid][2]) + k * es[2] / step[2],
                                                     std::get<0>(aa_range[acid][3]) + t * es[3] / step[3]
                                                   };
                    double value = f->interp(args.begin());

                    if(value < std::nextafter(0.0, 1.f))
                    {
                        value = std::nextafter(0.0, 1.f);
                    }
                    temp3.push_back(value);
                }

                /*int allof = 0;
                for( auto a : temp3 )
                {
                    if( a < std::nextafter( 0.0, 1.f ) )
                    {
                        allof++;
                    }
                }
                if( allof == temp3.size() )
                {
                    //std::cout << "THAT IT!\n";
                    temp3.back() = 1.0;
                }*/

                temp2.push_back(temp3);
            }
            temp1.push_back(temp2);
        }
        pdf.push_back(temp1);
    }

    // std::cout << pdf.size() << '\t' << pdf[0].size() << '\t' << pdf[0][0].size() << std::endl;
    // col cdf

    std::deque<std::deque<double>> first_2d;
    for(size_t i = 0; i != pdf.size(); i++)
    {
        std::deque<double> temp1;
        for(size_t j = 0; j != pdf[i].size(); j++)
        {
            temp1.push_back(pdf[i][j][0][0]);
        }
        first_2d.push_back(temp1);
    }

    std::deque<std::deque<double>> first_2d_temp = first_2d;
    for(size_t i = 0; i != pdf.size(); i++)
    {
        std::partial_sum(
            first_2d_temp[i].begin(), first_2d_temp[i].end(), first_2d_temp[i].begin(), std::plus<double>());
    }
    std::deque<double> col_pdf;
    for(size_t i = 0; i != pdf.size(); i++)
    {
        col_pdf.push_back(first_2d_temp[i].back());
    }

    double S = std::accumulate(col_pdf.begin(), col_pdf.end(), 0.0);
    for(auto &a : col_pdf)
    {
        a /= S;
    }

    std::deque<double> col_cdf(col_pdf.size());
    std::partial_sum(col_pdf.begin(), col_pdf.end(), col_cdf.begin(), std::plus<double>());
    col_cdf.push_front(0);

    // std::cout << first_2d.size() << '\t' << first_2d[0].size() << std::endl;

    //

    std::deque<std::deque<double>> row_dist_bank1(pdf.size());
    for(size_t i = 0; i != row_dist_bank1.size(); i++)
    {
        row_dist_bank1[i] = first_2d[i];

        double S = std::accumulate(row_dist_bank1[i].begin(), row_dist_bank1[i].end(), 0.0);
        for(auto &a : row_dist_bank1[i])
        {
            a /= S;
        }

        std::partial_sum(
            row_dist_bank1[i].begin(), row_dist_bank1[i].end(), row_dist_bank1[i].begin(), std::plus<double>());
        row_dist_bank1[i].push_front(0);
    }

    std::deque<std::deque<std::deque<double>>> row_dist_bank2(
        pdf.size(), std::deque<std::deque<double>>(pdf[0].size()));

    for(size_t i = 0; i != row_dist_bank2.size(); i++)
    {
        for(size_t j = 0; j != row_dist_bank2[i].size(); j++)
        {
            row_dist_bank2[i][j] = pdf[i][j][0];

            double S = std::accumulate(row_dist_bank2[i][j].begin(), row_dist_bank2[i][j].end(), 0.0);
            for(auto &a : row_dist_bank2[i][j])
            {
                a /= S;
            }

            std::partial_sum(row_dist_bank2[i][j].begin(), row_dist_bank2[i][j].end(), row_dist_bank2[i][j].begin(),
                             std::plus<double>());
            row_dist_bank2[i][j].push_front(0);
        }
    }

    std::deque<std::deque<std::deque<std::deque<double>>>> row_dist_bank3(
        pdf.size(), std::deque<std::deque<std::deque<double>>>(
            pdf[0].size(), std::deque<std::deque<double>>(pdf[0][0].size())));

    for(size_t i = 0; i != row_dist_bank3.size(); i++)
    {
        for(size_t j = 0; j != row_dist_bank3[i].size(); j++)
        {
            for(size_t k = 0; k != row_dist_bank3[i][j].size(); k++)
            {
                row_dist_bank3[i][j][k] = pdf[i][j][k];

                double S = std::accumulate(row_dist_bank3[i][j][k].begin(), row_dist_bank3[i][j][k].end(), 0.0);
                for(auto &a : row_dist_bank3[i][j][k])
                {
                    a /= S;
                }

                std::partial_sum(row_dist_bank3[i][j][k].begin(), row_dist_bank3[i][j][k].end(),
                                 row_dist_bank3[i][j][k].begin(), std::plus<double>());
                row_dist_bank3[i][j][k].push_front(0);
            }
        }
    }

    bbutils::distribution_4d result;
    // result.pdf = pdf;

    result.row_dist1 = row_dist_bank1;
    result.col_dist = col_cdf;

    result.row_dist2 = row_dist_bank2;
    result.row_dist3 = row_dist_bank3;

    result.grid.push_back(x);
    result.grid.push_back(y);
    result.grid.push_back(z);
    result.grid.push_back(q);

    return result;
}

void BBIND_top::initialize(size_t chi1_step,
                           size_t chi2_step,
                           size_t chi3_step,
                           size_t chi4_step,
                           std::string amino_acids)
{
    std::cout << "\n1d progress";

    if(amino_acids.find("S") != std::string::npos)
    {
        load_1d("rota500-", "ser", core::chemical::aa_ser);
        bbutils::distribution_1d dst = make_1d_cdf(core::chemical::aa_ser, chi1_step);
        aa_dst.insert({ core::chemical::aa_ser, dst });
        std::cout << ".";
    }
    if(amino_acids.find("V") != std::string::npos)
    {
        load_1d("rota500-", "val", core::chemical::aa_val);
        bbutils::distribution_1d dst = make_1d_cdf(core::chemical::aa_val, chi1_step);
        aa_dst.insert({ core::chemical::aa_val, dst });
        std::cout << ".";
    }
    if(amino_acids.find("C") != std::string::npos)
    {
        load_1d("rota500-", "cys", core::chemical::aa_cys);
        bbutils::distribution_1d dst = make_1d_cdf(core::chemical::aa_cys, chi1_step);
        aa_dst.insert({ core::chemical::aa_cys, dst });
        std::cout << ".";
    }
    if(amino_acids.find("T") != std::string::npos)
    {
        load_1d("rota500-", "thr", core::chemical::aa_thr);
        bbutils::distribution_1d dst = make_1d_cdf(core::chemical::aa_thr, chi1_step);
        aa_dst.insert({ core::chemical::aa_thr, dst });
        std::cout << ".";
    }

    std::cout << "\n2d progress";

    if(amino_acids.find("W") != std::string::npos)
    {
        load_2d("rota500-", "trp", core::chemical::aa_trp);
        bbutils::distribution_2d dst = make_2d_cdf(core::chemical::aa_trp, chi2_step);
        aa_dst.insert({ core::chemical::aa_trp, dst });
        std::cout << ".";
    }
    if(amino_acids.find("H") != std::string::npos)
    {
        load_2d("rota500-", "his", core::chemical::aa_his);
        bbutils::distribution_2d dst = make_2d_cdf(core::chemical::aa_his, chi2_step);
        aa_dst.insert({ core::chemical::aa_his, dst });
        std::cout << ".";
    }
    if(amino_acids.find("N") != std::string::npos)
    {
        load_2d("rota500-", "asn", core::chemical::aa_asn);
        bbutils::distribution_2d dst = make_2d_cdf(core::chemical::aa_asn, chi2_step);
        aa_dst.insert({ core::chemical::aa_asn, dst });
        std::cout << ".";
    }
    if(amino_acids.find("D") != std::string::npos)
    {
        load_2d("rota500-", "asp", core::chemical::aa_asp);
        bbutils::distribution_2d dst = make_2d_cdf(core::chemical::aa_asp, chi2_step);
        aa_dst.insert({ core::chemical::aa_asp, dst });
        std::cout << ".";
    }
    if(amino_acids.find("F") != std::string::npos)
    {
        load_2d("rota500-", "phetyr", core::chemical::aa_phe);
        bbutils::distribution_2d dst = make_2d_cdf(core::chemical::aa_phe, chi2_step);
        aa_dst.insert({ core::chemical::aa_phe, dst });
        std::cout << ".";
    }
    if(amino_acids.find("Y") != std::string::npos)
    {
        load_2d("rota500-", "phetyr", core::chemical::aa_tyr);
        bbutils::distribution_2d dst = make_2d_cdf(core::chemical::aa_tyr, chi2_step);
        aa_dst.insert({ core::chemical::aa_tyr, dst });
        std::cout << ".";
    }
    if(amino_acids.find("I") != std::string::npos)
    {
        load_2d("rota500-", "ile", core::chemical::aa_ile);
        bbutils::distribution_2d dst = make_2d_cdf(core::chemical::aa_ile, chi2_step);
        aa_dst.insert({ core::chemical::aa_ile, dst });
        std::cout << ".";
    }
    if(amino_acids.find("L") != std::string::npos)
    {
        load_2d("rota500-", "leu", core::chemical::aa_leu);
        bbutils::distribution_2d dst = make_2d_cdf(core::chemical::aa_leu, chi2_step);
        aa_dst.insert({ core::chemical::aa_leu, dst });
        std::cout << ".";
    }

    std::cout << "\n3d progress";

    if(amino_acids.find("M") != std::string::npos)
    {
        load_3d("rota500-", "met", core::chemical::aa_met);
        bbutils::distribution_3d dst = make_3d_cdf(core::chemical::aa_met, chi3_step);
        aa_dst.insert({ core::chemical::aa_met, dst });
        std::cout << ".";
    }
    if(amino_acids.find("E") != std::string::npos)
    {
        load_3d("rota500-", "glu", core::chemical::aa_glu);
        bbutils::distribution_3d dst = make_3d_cdf(core::chemical::aa_glu, chi3_step);
        aa_dst.insert({ core::chemical::aa_glu, dst });
        std::cout << ".";
    }
    if(amino_acids.find("Q") != std::string::npos)
    {
        load_3d("rota500-", "gln", core::chemical::aa_gln);
        bbutils::distribution_3d dst = make_3d_cdf(core::chemical::aa_gln, chi3_step);
        aa_dst.insert({ core::chemical::aa_gln, dst });
        std::cout << ".";
    }
    if(amino_acids.find("P") != std::string::npos)
    {
        load_3d("rota500-", "pro3d", core::chemical::aa_pro);
        bbutils::distribution_3d dst = make_3d_cdf(core::chemical::aa_pro, chi3_step);
        aa_dst.insert({ core::chemical::aa_pro, dst });
        std::cout << ".";
    }

    std::cout << "\n4d progress";

    if(amino_acids.find("R") != std::string::npos)
    {
        load_4d("rota500-", "arg", core::chemical::aa_arg);
        bbutils::distribution_4d dst = make_4d_cdf(core::chemical::aa_arg, chi4_step);
        aa_dst.insert({ core::chemical::aa_arg, dst });
        std::cout << ".";
    }
    if(amino_acids.find("K") != std::string::npos)
    {
        load_4d("rota500-", "lys", core::chemical::aa_lys);
        bbutils::distribution_4d dst = make_4d_cdf(core::chemical::aa_lys, chi4_step);
        aa_dst.insert({ core::chemical::aa_lys, dst });
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
    fill_1d(path + prefix + acid_name + ".data", range, grids, num_elements, grid_sizes, f_values);
    auto grid_iter_list = linterp::get_begins_ends(grids.begin(), grids.end());
    auto res =  std::make_shared<linterp::InterpMultilinear<1, double>>(grid_iter_list.first.begin(), grid_sizes.begin(), f_values.data(), f_values.data() + num_elements);
    aa_data.insert({ acid, std::move(res) });
    aa_range.insert({ acid, range });
}

void BBIND_top::load_2d(std::string prefix,
                        std::string acid_name, core::chemical::AA acid)
{
    size_t num_elements;
    std::vector<std::vector<double>> grids;
    boost::array<int, 2> grid_sizes;
    std::vector<double> f_values;
    std::vector<std::tuple<double, double, size_t>> range;
    fill_2d(path + prefix + acid_name + ".data", range, grids, num_elements, grid_sizes, f_values);
    auto grid_iter_list = linterp::get_begins_ends(grids.begin(), grids.end());
    auto res =  std::make_shared<linterp::InterpMultilinear<2, double>>(grid_iter_list.first.begin(), grid_sizes.begin(), f_values.data(), f_values.data() + num_elements);
    aa_data.insert({ acid, std::move(res) });
    aa_range.insert({ acid, range });
}

void BBIND_top::load_3d(std::string prefix,
                        std::string acid_name, core::chemical::AA acid)
{
    size_t num_elements;
    std::vector<std::vector<double>> grids;
    boost::array<int, 3> grid_sizes;
    std::vector<double> f_values;
    std::vector<std::tuple<double, double, size_t>> range;
    fill_3d(path + prefix + acid_name + ".data", range, grids, num_elements, grid_sizes, f_values);
    auto grid_iter_list = linterp::get_begins_ends(grids.begin(), grids.end());
    auto res =  std::make_shared<linterp::InterpMultilinear<3, double>>(grid_iter_list.first.begin(), grid_sizes.begin(), f_values.data(), f_values.data() + num_elements);
    aa_data.insert({ acid, std::move(res) });
    aa_range.insert({ acid, range });
}

void BBIND_top::load_4d(std::string prefix,
                        std::string acid_name, core::chemical::AA acid)
{
    size_t num_elements;
    std::vector<std::vector<double>> grids;
    boost::array<int, 4> grid_sizes;
    std::vector<double> f_values;
    std::vector<std::tuple<double, double, size_t>> range;
    fill_4d(path + prefix + acid_name + ".data", range, grids, num_elements, grid_sizes, f_values);
    auto grid_iter_list = linterp::get_begins_ends(grids.begin(), grids.end());
    auto res =  std::make_shared<linterp::InterpMultilinear<4, double>>(grid_iter_list.first.begin(), grid_sizes.begin(), f_values.data(), f_values.data() + num_elements);
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

void BBIND_top::fill_2d(std::string name,
                        std::vector<std::tuple<double, double, size_t>> &range,
                        std::vector<std::vector<double>> &grids,
                        size_t &num_elements,
                        boost::array<int, 2> &grid_sizes,
                        std::vector<double> &f_values)
{
    std::vector<std::vector<double>> data;
    load_data(name, data, range);

    if(range.size() != 2)
    {
        std::cout << "Not a 2d file!" << '\t' << name << std::endl;
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

    std::deque<std::deque<double>> data_2d;
    for(size_t i = 0, a = 0; i != std::get<2>(range[0]); i++)
    {
        std::deque<double> temp;
        for(size_t j = 0; j != std::get<2>(range[1]); j++, a++)
        {
            // j + i*std::get<2>(range[1])
            temp.push_back(data[a].back());
        }
        data_2d.push_back(temp);
    }

    for(size_t i = 0; i != data_2d.size(); i++)
    {
        data_2d[i].push_front(data_2d[i].back());
        data_2d[i].push_back(data_2d[i][1]);
    }
    data_2d.push_front(data_2d.back());
    data_2d.push_back(data_2d[1]);

    for(size_t i = 0; i != grids[0].size(); i++)
    {
        for(size_t j = 0; j != grids[1].size(); j++)
        {
            f_values.push_back(data_2d[i][j]);
        }
    }

    num_elements = f_values.size();
}

void BBIND_top::fill_3d(std::string name,
                        std::vector<std::tuple<double, double, size_t>> &range,
                        std::vector<std::vector<double>> &grids,
                        size_t &num_elements,
                        boost::array<int, 3> &grid_sizes,
                        std::vector<double> &f_values)
{
    std::vector<std::vector<double>> data;
    load_data(name, data, range);

    if(range.size() != 3)
    {
        std::cout << "Not a 3d file!" << '\t' << name << std::endl;
        return;
    }

    //	for (size_t i = 0; i != range.size(); i++)
    //	{
    //		std::cout << std::get < 0 > (range[i]) << '\t' << std::get < 1 > (range[i]) << '\t' << std::get
    //<
    // 2 >
    //(range[i]) << std::endl;
    //	}

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

    std::deque<std::deque<std::deque<double>>> data_3d;
    for(size_t i = 0, a = 0; i != std::get<2>(range[0]); i++)
    {
        std::deque<std::deque<double>> temp1;
        for(size_t j = 0; j != std::get<2>(range[1]); j++)
        {
            std::deque<double> temp2;
            for(size_t k = 0; k != std::get<2>(range[2]); k++, a++)
            {
                // k + j*std::get<2>(range[2]) + i*std::get<2>(range[1])*std::get<2>(range[2])].back()
                temp2.push_back(data[a].back());
            }
            temp1.push_back(temp2);
        }
        data_3d.push_back(temp1);
    }

    for(size_t i = 0; i != data_3d.size(); i++)
    {
        for(size_t j = 0; j != data_3d[i].size(); j++)
        {
            data_3d[i][j].push_front(data_3d[i][j].back());
            data_3d[i][j].push_back(data_3d[i][j][1]);
        }
        data_3d[i].push_front(data_3d[i].back());
        data_3d[i].push_back(data_3d[i][1]);
    }
    data_3d.push_front(data_3d.back());
    data_3d.push_back(data_3d[1]);

    for(size_t i = 0; i != grids[0].size(); i++)
    {
        for(size_t j = 0; j != grids[1].size(); j++)
        {
            for(size_t k = 0; k != grids[2].size(); k++)
            {
                f_values.push_back(data_3d[i][j][k]);
            }
        }
    }

    num_elements = f_values.size();
}

void BBIND_top::fill_4d(std::string name,
                        std::vector<std::tuple<double, double, size_t>> &range,
                        std::vector<std::vector<double>> &grids,
                        size_t &num_elements,
                        boost::array<int, 4> &grid_sizes,
                        std::vector<double> &f_values)
{
    std::vector<std::vector<double>> data;
    load_data(name, data, range);

    if(range.size() != 4)
    {
        std::cout << "Not a 4d file!" << '\t' << name << std::endl;
        return;
    }

    //	for (size_t i = 0; i != range.size(); i++)
    //	{
    //		std::cout << std::get < 0 > (range[i]) << '\t' << std::get < 1 > (range[i]) << '\t' << std::get
    //<
    // 2 >
    //(range[i]) << std::endl;
    //	}

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

    std::deque<std::deque<std::deque<std::deque<double>>>> data_4d;
    for(size_t i = 0, a = 0; i != std::get<2>(range[0]); i++)
    {
        std::deque<std::deque<std::deque<double>>> temp1;
        for(size_t j = 0; j != std::get<2>(range[1]); j++)
        {
            std::deque<std::deque<double>> temp2;
            for(size_t k = 0; k != std::get<2>(range[2]); k++)
            {
                std::deque<double> temp3;
                for(size_t t = 0; t != std::get<2>(range[3]); t++, a++)
                {
                    // t + k*std::get<2>(range[3]) + j*std::get<2>(range[3])*std::get<2>(range[2]) +
                    // i*std::get<2>(range[3])*std::get<2>(range[2])*std::get<2>(range[1])
                    temp3.push_back(
                        data[t + k * std::get<2>(range[3]) + j * std::get<2>(range[3]) * std::get<2>(range[2]) +
                               i * std::get<2>(range[3]) * std::get<2>(range[2]) * std::get<2>(range[1])]
                        .back());
                }
                temp2.push_back(temp3);
            }
            temp1.push_back(temp2);
        }
        data_4d.push_back(temp1);
    }

    for(size_t i = 0; i != data_4d.size(); i++)
    {
        for(size_t j = 0; j != data_4d[i].size(); j++)
        {
            for(size_t k = 0; k != data_4d[i][j].size(); k++)
            {
                data_4d[i][j][k].push_front(data_4d[i][j][k].back());
                data_4d[i][j][k].push_back(data_4d[i][j][k][1]);
            }
            data_4d[i][j].push_front(data_4d[i][j].back());
            data_4d[i][j].push_back(data_4d[i][j][1]);
        }
        data_4d[i].push_front(data_4d[i].back());
        data_4d[i].push_back(data_4d[i][1]);
    }
    data_4d.push_front(data_4d.back());
    data_4d.push_back(data_4d[1]);

    for(size_t i = 0; i != grids[0].size(); i++)
    {
        for(size_t j = 0; j != grids[1].size(); j++)
        {
            for(size_t k = 0; k != grids[2].size(); k++)
            {
                for(size_t t = 0; t != grids[3].size(); t++)
                {
                    f_values.push_back(data_4d[i][j][k][t]);
                }
            }
        }
    }

    num_elements = f_values.size();
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
        tt.push_back(std::vector<double> {i, s->interp(args.begin()),v->interp(args.begin()),c->interp(args.begin()),t->interp(args.begin())});
    }
    pepsgo::write_default2d("maps/bbind/SVCT.dat", tt, 10);
}

void plot_chi2(pepsgo::bbind::BBIND_top& obj, core::chemical::AA acid)
{
    auto f = std::any_cast<std::shared_ptr<linterp::InterpMultilinear<2, double>>>(obj.aa_data[acid]);

    std::vector<double> es;
    std::vector<size_t> step;
    for(size_t i = 0; i != obj.aa_range[acid].size(); i++)
    {
        es.push_back(std::get<1>(obj.aa_range[acid][i]) - std::get<0>(obj.aa_range[acid][i]));
        step.push_back(15*std::get<2>(obj.aa_range[acid][i]));
    }

    std::vector<std::vector<double> > pdf;
    std::vector<std::vector<double> > pdf_list;
    for(size_t i = 0; i != step[0]; i++)
    {
        std::vector<double> temp1;
        for(size_t j = 0; j != step[1]; j++)
        {

            boost::array<double, 2> args = { std::get<0>(obj.aa_range[acid][0]) + i*es[0]/step[0],
                                             std::get<0>(obj.aa_range[acid][1]) + j*es[1]/step[1]
                                           };

            double value = f->interp(args.begin());
            temp1.push_back(value);

            pdf_list.push_back(std::vector<double>{args.front(), args.back(), value});
        }
        pdf.push_back(temp1);
    }
    pepsgo::write_default2d("maps/bbind/2d/" + core::chemical::name_from_aa(acid) + "top500.dat", pdf, 4);//trp his asn asp phe tyr ile leu
    pepsgo::write_default2d("maps/bbind/2d/" + core::chemical::name_from_aa(acid) + "top500_list.dat", pdf_list, 4);//trp his asn asp phe tyr ile leu
}

}
}
