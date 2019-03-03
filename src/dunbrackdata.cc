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
#include <dunbrackdata.hh>

#include <tuple>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>

#include <utility/io/izstream.hh>

namespace pepsgo
{
namespace bbdep
{

std::ostream& operator<<(std::ostream& stream, const Dunbrack_data& obj)
{
    stream << obj.name << ' ' << obj.Phi << ' ' << obj.Psi << ' ' << obj.Count << ' ' << obj.r1 << ' ' << obj.r2 << ' ' << obj.r3 << ' ' << obj.r4;
    stream << ' ' << obj.Probabil << ' ' << obj.chi1Val << ' ' << obj.chi2Val << ' ' << obj.chi3Val << ' ' << obj.chi4Val << ' ';
    stream << obj.chi1Sig << ' ' << obj.chi2Sig << ' ' << obj.chi3Sig << ' ' << obj.chi4Sig;
    return stream;
}

std::string remove_chars(const std::string is, const std::string& chars)
{
    std::string s = is;
    s.erase(std::remove_if(s.begin(), s.end(), [&chars](const char& c)
    {
        return chars.find(c) != std::string::npos;
    }), s.end());
    return s;
}

void load_data_sm(std::string fname, std::vector<Dunbrack_data> &data, std::vector<std::vector<int>> &numbers)
{
    utility::io::izstream fIn(fname.c_str());
//    std::ifstream fIn(fname.c_str());
    if(!fIn.good())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }
    size_t dimension = 16;
    int integer_temp = 0;
    numbers.resize(4);
    std::string temp;
    bool flag = false;
    Dunbrack_data push2data;
    while(!fIn.eof())
    {
        if(flag)
        {
            fIn >> push2data.name;
            fIn >> push2data.Phi;
            fIn >> push2data.Psi;
            fIn >> push2data.Count;
            fIn >> push2data.r1;
            fIn >> push2data.r2;
            fIn >> push2data.r3;
            fIn >> push2data.r4;

            fIn >> temp;
            if(has_only_zeros_or_dots(temp))
            {
                temp = temp + "9";
                push2data.Probabil = std::stod(temp);
            }
            else
            {
                push2data.Probabil = std::stod(temp);
            }


            //fIn >> push2data.Probabil;
            fIn >> push2data.chi1Val;
            fIn >> push2data.chi2Val;
            fIn >> push2data.chi3Val;
            fIn >> push2data.chi4Val;
            fIn >> push2data.chi1Sig;
            fIn >> push2data.chi2Sig;
            fIn >> push2data.chi3Sig;
            fIn >> push2data.chi4Sig;

            data.push_back(push2data);
            continue;
        }
        fIn >> temp;
        if("chi4Sig" == temp)
        {
            fIn >> temp;
            flag = true;

            numbers[1].push_back(numbers[2].size());
            numbers[3].push_back(std::accumulate(numbers[2].begin(), numbers[2].end(), 1, std::multiplies<int>()));
        }
        if("freedom)" == temp)
        {
            fIn >> integer_temp;
            numbers[0].push_back(integer_temp);
        }
        if("angle" == temp)
        {
            while(temp.back() != ']')
            {
                fIn >> temp;
                std::string dst = remove_chars(temp, "[],");
                integer_temp = std::stoi(dst);
                numbers[2].push_back(integer_temp);
            }
        }
    }
    fIn.close();
    data.pop_back();
}

Dunbrack_data get_max_prob_object(std::vector<Dunbrack_data> &data, std::vector<std::vector<int>> &numbers, double Phi, double Psi)
{
    Dunbrack_data search_data;
    search_data.Phi = 10 * std::round(Phi / 10.0);
    search_data.Psi = 10 * std::round(Psi / 10.0);

    auto p = std::equal_range(data.begin(), data.end(), search_data, [](const Dunbrack_data& lhs, const Dunbrack_data& rhs) -> bool
    {	return lhs.Phi < rhs.Phi;});

    size_t index1 = std::distance(data.begin(), std::lower_bound(p.first, p.second, search_data, [](const Dunbrack_data& lhs, const Dunbrack_data& rhs) -> bool
    {	return lhs.Psi < rhs.Psi;}));

    return data[index1];
}

Dunbrack_data get_probabil(std::vector<Dunbrack_data> &data, std::vector<std::vector<int>> &numbers, std::mt19937_64 &generator, double Phi, double Psi)
{
    Dunbrack_data search_data;
    search_data.Phi = 10 * std::round(Phi / 10.0);
    search_data.Psi = 10 * std::round(Psi / 10.0);

    auto p = std::equal_range(data.begin(), data.end(), search_data, [](const Dunbrack_data& lhs, const Dunbrack_data& rhs) -> bool
    {	return lhs.Phi < rhs.Phi;});

    auto q = std::equal_range(p.first, p.second, search_data, [](const Dunbrack_data& lhs, const Dunbrack_data& rhs) -> bool
    {	return lhs.Psi < rhs.Psi;});

    size_t range = std::distance(q.first, q.second) / 4;
    std::uniform_int_distribution<size_t> distribution(0, range);
    return data[std::distance(data.begin(), q.first) + distribution(generator)];
}

bool has_only_zeros_or_dots(const std::string str)
{
    bool flag = true;
    for(std::string::const_iterator it = str.begin(); it != str.end(); ++it)
    {
        if(*it != '0' && *it != '.')
            flag = false;
    }
    return flag;
}

}
}
