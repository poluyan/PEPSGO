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
#include "fragment.hh"

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/picking_old/vall/util.hh>
#include <core/fragment/picking_old/FragmentLibraryManager.hh>

#include <cmath>

namespace pepsgo
{
namespace fragment
{

Frag::Frag(): aa(' ')
{

}
Frag::Frag(char _aa, core::Real _phi, core::Real _psi, core::Real _omg)
{
    aa = _aa;
    phi = _phi;
    psi = _psi;
    omg = _omg;
}


template <typename T>
bool increase(const std::vector<std::vector<T>> &v, std::vector<size_t> &it)
{
    for(size_t i = 0, size = it.size(); i != size; i++)
    {
        const size_t index = size - 1 - i;
        ++it[index];
        if(it[index] == v[index].size())
        {
            it[index] = 0;
        }
        else
        {
            return true;
        }
    }
    return false;
}

template <typename T>
std::vector<T> get_line(const std::vector<std::vector<T>> &v, std::vector<size_t> &it)
{
    std::vector<T> rez(v.size());
    for(size_t i = 0, size = v.size(); i != size; i++)
    {
        rez[i] = v[i][it[i]];
    }
    return rez;
}

template <typename T>
std::vector<std::vector<T>> iterate(const std::vector<std::vector<T>> &v)
{
    std::vector<size_t> it(v.size(), 0);
    std::vector<std::vector<T>> values;
    do
    {
        values.push_back(get_line(v, it));
    }
    while(increase(v, it));
    return values;
}

template <typename T>
size_t how_much_will_be_iterate(const std::vector<std::vector<T>> &v)
{
    size_t count = 0;
    std::vector<size_t> it(v.size(), 0);
    std::vector<std::vector<T>> values;
    do
    {
        get_line(v, it);
        count++;
    }
    while(increase(v, it));
    return count;
}

std::vector<size_t> get_permut_sizes(size_t total_resiudes, size_t residues, size_t frag_size)
{
    std::vector<std::vector<size_t>> matrix(residues, std::vector<size_t>(total_resiudes, 0));
    for(size_t i = 0; i != matrix.size(); i++)
    {
        for(size_t j = i; j != frag_size + i; j++)
        {
            matrix[i][j] = 1;
        }
    }

    for(size_t i = 0; i != matrix.size(); i++)
    {
        for(size_t j = 0; j != matrix[i].size(); j++)
        {
            std::cout << matrix[i][j] << ' ';
        }
        std::cout << std::endl;
    }

    std::vector<size_t> sums(total_resiudes);
    for(size_t i = 0; i != sums.size(); i++)
    {
        sums[i] = 0;
        for(size_t j = 0; j != matrix.size(); j++)
        {
            sums[i] += matrix[j][i];
        }
    }

    for(size_t i = 0; i != sums.size(); i++)
    {
        std::cout << sums[i] << ' ';
    }
    std::cout << std::endl;

    return sums;
}

std::vector<size_t> get_permut_sizes2(size_t total_resiudes, size_t residues, size_t frag_size)
{
    size_t matrix_size = total_resiudes%frag_size ? total_resiudes/frag_size + 1 : total_resiudes/frag_size;
    std::vector<std::vector<size_t>> matrix(matrix_size, std::vector<size_t>(total_resiudes, 0));
    size_t delta = frag_size;
    for(size_t i = 0, t = delta, k = 0; i != matrix.size(); i++, t+=delta)
    {
        for(size_t j = k; (j != t) && (j != matrix[i].size()); j++)
        {
            matrix[i][j] = 1;
            k = j + 1;
        }
    }
    for(size_t i = 0; i != matrix.size(); i++)
    {
        for(size_t j = 0; j != matrix[i].size(); j++)
        {
            std::cout << matrix[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    std::vector<size_t> sums(total_resiudes);
    for(size_t i = 0; i != sums.size(); i++)
    {
        sums[i] = 0;
        for(size_t j = 0; j != matrix.size(); j++)
        {
            sums[i] += matrix[j][i];
        }
    }
    for(size_t i = 0; i != sums.size(); i++)
    {
        std::cout << sums[i] << ' ';
    }
    std::cout << std::endl;

    return sums;
}

std::vector<std::vector<Frag>> get_frag_mix(std::vector<std::vector<Frag>> res, size_t total_resiudes, size_t residues, size_t frag_size)
{
    std::vector<std::vector<Frag>> matrix(residues, std::vector<Frag>(total_resiudes));
    for(size_t i = 0; i != matrix.size(); i++)
    {
        for(size_t j = i, t = 0; j != frag_size + i; j++, t++)
        {
            matrix[i][j] = res[i][t];
        }
    }
    std::vector<std::vector<Frag>> sums;
    for(size_t i = 0; i != matrix.front().size(); i++)
    {
        std::vector<Frag> temp;
        for(size_t j = 0; j != matrix.size(); j++)
        {
            if(matrix[j][i].aa != ' ')
                temp.push_back(matrix[j][i]);
        }
        sums.push_back(temp);
    }
    return sums;
}

std::vector<std::vector<Frag>> get_frag_mix2(std::vector<std::vector<Frag>> res, size_t total_resiudes, size_t residues, size_t frag_size)
{
    size_t matrix_size = total_resiudes%frag_size ? total_resiudes/frag_size + 1 : total_resiudes/frag_size;
    std::vector<std::vector<Frag>> matrix(matrix_size, std::vector<Frag>(total_resiudes));
    size_t delta = frag_size;

    for(size_t i = 0, t = delta, k = 0; i != matrix.size(); i++, t+=delta)
    {
        for(size_t j = k, m = 0; (j != t) && (j != matrix[i].size()); j++, m++)
        {
            matrix[i][j] = res[i][m];
            k = j + 1;
        }
    }
    std::vector<std::vector<Frag>> sums;
    for(size_t i = 0; i != matrix.front().size(); i++)
    {
        std::vector<Frag> temp;
        for(size_t j = 0; j != matrix.size(); j++)
        {
            if(matrix[j][i].aa != ' ')
                temp.push_back(matrix[j][i]);
        }
        sums.push_back(temp);
    }
    return sums;
}


FragPick::FragPick()
{

}

void FragPick::set_peptide(const core::pose::Pose &_peptide)
{
    peptide = _peptide;
    peptide_seq = peptide.sequence();
}

void FragPick::set_file()
{
    frag_size = basic::options::option[basic::options::OptionKeys::frags::frag_sizes]()[1];
    std::cout << "frag_sizes " << frag_size << std::endl;
    std::cout << "fragmet_file " << fragmet_file << std::endl;
}

void FragPick::fill_grids(size_t phipsi_step, size_t omega_step)
{
    phipsi_grid_number = std::vector<core::Size>(2*peptide_seq.size() - 2, phipsi_step);
    omega_grid_number = std::vector<core::Size>(peptide_seq.size() - 1, omega_step);

    phipsi_grids_dx.resize(phipsi_grid_number.size());
    phipsi_grids_radians.resize(phipsi_grid_number.size());
    phipsi_grids_nodes.resize(phipsi_grid_number.size());
    for(size_t i = 0; i != phipsi_grids_radians.size(); i++)
    {
        std::vector<core::Real> grid(phipsi_grid_number[i] + 1);
        std::vector<core::Real> grid_nodes(phipsi_grid_number[i]);
        core::Real startp = -numeric::NumericTraits<core::Real>::pi();
        core::Real endp = numeric::NumericTraits<core::Real>::pi();
        core::Real es = endp - startp;
        for(size_t j = 0; j != grid.size(); j++)
        {
            grid[j] = startp + j*es/core::Real(phipsi_grid_number[i]);
        }

        phipsi_grids_radians[i] = grid;
        phipsi_grids_dx[i] = es/(core::Real(phipsi_grid_number[i])*2);

        for(size_t j = 0; j != grid_nodes.size(); j++)
        {
            grid_nodes[j] = grid[j] + phipsi_grids_dx[i];
        }
        phipsi_grids_nodes[i] = grid_nodes;
    }

    omega_grids_dx.resize(omega_grid_number.size());
    omega_grids_radians.resize(omega_grid_number.size());
    omega_grids_nodes.resize(omega_grid_number.size());
    for(size_t i = 0; i != omega_grids_radians.size(); i++)
    {
        std::vector<core::Real> grid(omega_grid_number[i] + 1);
        std::vector<core::Real> grid_nodes(omega_grid_number[i]);
        core::Real startp = -numeric::NumericTraits<core::Real>::pi();
        core::Real endp = numeric::NumericTraits<core::Real>::pi();
        core::Real es = endp - startp;
        for(size_t j = 0; j != grid.size(); j++)
        {
            grid[j] = startp + j*es/core::Real(omega_grid_number[i]);
        }
        omega_grids_radians[i] = grid;
        omega_grids_dx[i] = es/(core::Real(omega_grid_number[i])*2);
        for(size_t j = 0; j != grid_nodes.size(); j++)
        {
            grid_nodes[j] = grid[j] + omega_grids_dx[i];
        }
        omega_grids_nodes[i] = grid_nodes;
    }
}

size_t FragPick::get_index_phipsi(core::Real radian, size_t index) const
{
    auto pos = std::lower_bound(phipsi_grids_radians[index].begin(), phipsi_grids_radians[index].end(), radian);
    size_t dist = std::distance(phipsi_grids_radians[index].begin(), pos);
    dist = dist > phipsi_grids_radians[index].size() - 1 ? dist - 1 : dist;
    return dist > 0 ? dist - 1 : 0;
}

size_t FragPick::get_index_omega(core::Real radian, size_t index) const
{
    auto pos = std::lower_bound(omega_grids_radians[index].begin(), omega_grids_radians[index].end(), radian);
    size_t dist = std::distance(omega_grids_radians[index].begin(), pos);
    dist = dist > omega_grids_radians[index].size() - 1 ? dist - 1 : dist;
    return dist > 0 ? dist - 1 : 0;
}


void FragPick::set_storage_shared(std::shared_ptr<sample_type> in_sample)
{
    structures_trie = std::move(in_sample);
}

void FragPick::load_frag_file()
{
    std::cout << "loading fragment file...";

    core::fragment::ConstantLengthFragSet fragsetNmer(frag_size);//(new core::fragment::ConstantLengthFragSet( 3 ));
    fragsetNmer.read_fragment_file(fragmet_file);

    for(size_t i = 1, iend = fragsetNmer.nr_frames(); i <= iend; ++i)
    {
        std::vector<std::vector<Frag>> pos;
        core::fragment::FrameList frames;
        fragsetNmer.frames(i, frames);
        for(size_t j = 1, jend = frames.front()->nr_frags(); j <= jend; ++j)
        {
            std::vector<Frag> one_frag;

            //core::fragment::FragData fd = frames.front()->fragment(j);
            for(size_t k = 1, kend = frames.front()->fragment(j).size(); k <= kend; ++k)
            {
                //frames.front()->fragment(j)->show(std::cout);
                core::fragment::BBTorsionSRFD const & frag_i =
                    dynamic_cast< core::fragment::BBTorsionSRFD const &>(*(frames.front()->fragment(j).get_residue(k)));

                one_frag.push_back(Frag(frag_i.sequence(), frag_i.torsion(1), frag_i.torsion(2), frag_i.torsion(3)));
            }
            pos.push_back(one_frag);
        }
        all_fragments.push_back(pos);
    }
    std::cout << "ok" << std::endl;
    std::cout << "residues " << all_fragments.size() << std::endl;
    std::cout << "n_frags " << all_fragments.front().size() << std::endl;
}

void FragPick::make_permutations()
{
    bool full = false;

    size_t residues = all_fragments.size(); // not total_residues!
    size_t n_frags = all_fragments.front().size(); // n_frags ^ residues

    if(full)
    {
        std::vector<size_t> sizes = get_permut_sizes(peptide.total_residue(), residues, frag_size);

        std::vector<std::vector<int>> for_selected(sizes.size());
        for(size_t i = 0, k = 0; i != for_selected.size(); i++, k+=frag_size)
        {
            for_selected[i].resize(sizes[i]);
            for(size_t j = 0; j != sizes[i]; j++)
            {
                for_selected[i][j] = j;
            }
        }
        std::cout << "it will be " << std::fixed <<
                  n_frags << " ^ " << residues << " * " << how_much_will_be_iterate(for_selected) << " = "
                  << size_t(std::pow(n_frags, residues)*how_much_will_be_iterate(for_selected)) << std::endl;
        std::cin.get();

        std::vector<std::vector<int>> permut_for_selected = iterate(for_selected);

        std::vector<std::vector<int>> variable_values(residues, std::vector<int>(n_frags));
        for(size_t i = 0; i != variable_values.size(); i++)
        {
            for(size_t j = 0; j != n_frags; j++)
            {
                variable_values[i][j] = j;
                std::cout << variable_values[i][j] << ' ';
            }
            std::cout << std::endl;
        }
        std::vector<std::vector<int>> permut = iterate(variable_values);    // n_frags^residues
//        structures.clear();
        for(size_t i = 0; i != permut.size(); i++)
        {
            std::vector<std::vector<Frag>> to_distr;
            for(size_t j = 0; j != permut[i].size(); j++)
            {
                to_distr.push_back(all_fragments[j][permut[i][j]]);
            }
//                std::cout << "to_distr " << to_distr.size() << std::endl;
            auto distr = get_frag_mix(to_distr, peptide.total_residue(), residues, frag_size);

            for(size_t j = 0; j != permut_for_selected.size(); j++)
            {
                std::vector<Frag> temp;
                for(size_t k = 0; k != permut_for_selected[j].size(); k++)
                {
                    temp.push_back(distr[k][permut_for_selected[j][k]]);
                }
//                structures.push_back(temp);


                std::vector<double> phipsi_values_radians(phipsi_grid_number.size());
                std::vector<double> omega_values_radians(omega_grid_number.size());
                std::vector<std::uint8_t> to_trie(phipsi_values_radians.size() + omega_values_radians.size());

                phipsi_values_radians.front() = numeric::conversions::radians(temp.front().psi);

                phipsi_values_radians.front() = numeric::conversions::radians(temp.front().psi);
                omega_values_radians.front() = numeric::conversions::radians(temp.front().omg);
                for(size_t j = 1, k = 1; j != temp.size() - 1; j++, k+=2)
                {
                    phipsi_values_radians[k] = numeric::conversions::radians(temp[j].phi);
                    phipsi_values_radians[k + 1] = numeric::conversions::radians(temp[j].psi);
                    omega_values_radians[j] = numeric::conversions::radians(temp[j].omg);
                }
                phipsi_values_radians.back() = numeric::conversions::radians(temp.back().phi);

                for(size_t j = 0; j != phipsi_values_radians.size(); j++)
                {
                    to_trie[j] = get_index_phipsi(phipsi_values_radians[j], j);
                }
                for(size_t j = phipsi_values_radians.size(), k = 0; j != to_trie.size(); j++, k++)
                {
                    to_trie[j] = get_index_omega(omega_values_radians[k], k);
                }
                if(!structures_trie->search(to_trie))
                    structures_trie->insert(to_trie);
            }
        }
        std::cout << permut.size() << std::endl;
//        std::cout << structures.size() << std::endl;
    }
    else
    {

        std::vector<size_t> sizes = get_permut_sizes2(peptide.total_residue(), residues, frag_size);

        size_t matrix_size = peptide.total_residue()%frag_size ? peptide.total_residue()/frag_size + 1 : peptide.total_residue()/frag_size;

        std::vector<std::vector<int>> for_selected(sizes.size());
        for(size_t i = 0, k = 0; i != for_selected.size(); i++, k+=frag_size)
        {
            for_selected[i].resize(sizes[i]);
            for(size_t j = 0; j != sizes[i]; j++)
            {
                for_selected[i][j] = j;
            }
        }
        std::cout << "it will be " << std::fixed <<
                  n_frags << " ^ " << matrix_size << " * " << how_much_will_be_iterate(for_selected) << " = "
                  << size_t(std::pow(n_frags, matrix_size)*how_much_will_be_iterate(for_selected)) << std::endl;
        std::cin.get();

        std::vector<std::vector<int>> permut_for_selected = iterate(for_selected);

        std::vector<std::vector<int>> variable_values(matrix_size, std::vector<int>(n_frags));
        for(size_t i = 0; i != variable_values.size(); i++)
        {
            for(size_t j = 0; j != n_frags; j++)
            {
                variable_values[i][j] = j;
                //std::cout << variable_values[i][j] << ' ';
            }
            //std::cout << std::endl;
        }
        std::vector<std::vector<int>> permut = iterate(variable_values);    // n_frags^residues
        std::cout << permut.size() << std::endl;
        std::cout << all_fragments.size() << std::endl;
//        structures.clear();
        for(size_t i = 0; i != permut.size(); i++)
        {
            std::vector<std::vector<Frag>> to_distr;
            size_t delta = frag_size;
            for(size_t j = 0, t = 0; j != permut[i].size(); j++, t+=delta)
            {
                to_distr.push_back(all_fragments[t][permut[i][j]]);
                if(t + delta >= all_fragments.size())
                    delta = 1;
            }
            //std::cout << "to_distr " << to_distr.size() << std::endl;
            auto distr = get_frag_mix2(to_distr, peptide.total_residue(), residues, frag_size);

            for(size_t j = 0; j != permut_for_selected.size(); j++)
            {
                std::vector<Frag> temp;
                for(size_t k = 0; k != permut_for_selected[j].size(); k++)
                {
                    temp.push_back(distr[k][permut_for_selected[j][k]]);
                }
//                structures.push_back(temp);

                std::vector<double> phipsi_values_radians(phipsi_grid_number.size());
                std::vector<double> omega_values_radians(omega_grid_number.size());
                std::vector<std::uint8_t> to_trie(phipsi_values_radians.size() + omega_values_radians.size());

                phipsi_values_radians.front() = numeric::conversions::radians(temp.front().psi);

                phipsi_values_radians.front() = numeric::conversions::radians(temp.front().psi);
                omega_values_radians.front() = numeric::conversions::radians(temp.front().omg);
                for(size_t j = 1, k = 1; j != temp.size() - 1; j++, k+=2)
                {
                    phipsi_values_radians[k] = numeric::conversions::radians(temp[j].phi);
                    phipsi_values_radians[k + 1] = numeric::conversions::radians(temp[j].psi);
                    omega_values_radians[j] = numeric::conversions::radians(temp[j].omg);
                }
                phipsi_values_radians.back() = numeric::conversions::radians(temp.back().phi);

                for(size_t j = 0; j != phipsi_values_radians.size(); j++)
                {
                    to_trie[j] = get_index_phipsi(phipsi_values_radians[j], j);
                }
                for(size_t j = phipsi_values_radians.size(), k = 0; j != to_trie.size(); j++, k++)
                {
                    to_trie[j] = get_index_omega(omega_values_radians[k], k);
                }
                if(!structures_trie->search(to_trie))
                    structures_trie->insert(to_trie);


            }
        }
        std::cout << permut.size() << std::endl;
//        std::cout << structures.size() << std::endl;
    }
}

}
}
