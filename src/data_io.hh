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
#ifndef INCLUDED_data_writing_hh
#define INCLUDED_data_writing_hh

#include <fstream>
#include <iostream>

namespace pepsgo
{

template<template<typename, typename...> class Container, class T, typename... Params>
void write_default1d(std::string fname, Container<T, Params...> const& u, size_t step, size_t prec)
{
    std::ofstream fOut;
    fOut.open(fname.c_str());
    if(!fOut.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }
    fOut.precision(prec);
    for(size_t i = 0; i != u.size(); i += step)
    {
        fOut << std::scientific << i << '\t' << u[i] << std::endl;
    }
    fOut.close();
    std::cout << fname << std::endl;
}

template<template<typename, typename...> class Container, class T, typename... Params>
void write_default2d(std::string fname, Container<T, Params...> const& u, size_t prec)
{
    std::ofstream fOut;
    fOut.open(fname.c_str());
    if(!fOut.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }
    fOut.precision(prec);
    for(const auto & i : u)
    {
        for(const auto & j : i)
        {
            fOut << std::scientific << j << '\t';
        }
        fOut << std::endl;
    }
    fOut.close();
    std::cout << fname << std::endl;
}

}

#endif
