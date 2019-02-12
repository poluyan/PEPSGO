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

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/kinematics/MoveMap.hh>
#include <protocols/minimization_packing/MinMover.hh>

#include <cmath>
#include <fstream>
#include "bbtools.hh"

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

size_t column_count(const std::vector<std::vector<size_t>> &matrix, size_t c, size_t filln)
{
    size_t count = 0;
    for(size_t i = 0; i != matrix.size(); i++)
    {
        if(matrix[i][c] != filln)
            count += matrix[i][c];
        else
            count++;
    }
    return count;
}

std::vector<size_t> get_permut_sizes3(size_t total_resiudes, size_t residues, size_t frag_size, std::string ss_pred, size_t max_in_column, size_t &matrix_size, std::vector<std::vector<size_t>> &full_matrix)
{
    size_t filln = 1e6;
    std::vector<std::vector<size_t>> matrix(residues, std::vector<size_t>(total_resiudes, 0));
    for(size_t i = 0; i != matrix.size(); i++)
    {
        for(size_t j = i; j != frag_size + i; j++)
        {
            matrix[i][j] = 1;
        }
    }
    size_t delta = frag_size, k = 0;
    for(size_t i = 0, t = delta; i < matrix.size(); i+=delta, t+=delta)
    {
        for(size_t j = k; (j != t) && (j != matrix[i].size()); j++)
        {
            matrix[i][j] = filln;
            k = j + 1;
        }
    }
    for(size_t i = k; i < matrix.back().size(); i++)
    {
        matrix.back()[i] = filln;
    }
    for(size_t i = 0; i != matrix.size(); i++)
    {
        for(size_t j = 0; j != matrix[i].size(); j++)
        {
            if(matrix[i][j] == filln)
            {
                matrix[i][j] = 1;
                continue;
            }
            if(matrix[i][j] == 0)
                continue;
            if(ss_pred[j] == 'C' && column_count(matrix, j, filln) <= max_in_column)
                //if(j > 1 && j < matrix[i].size() - 2 && ss_pred[j] == 'C' && column_count(matrix, j, filln) <= max_in_column)
                matrix[i][j] = 1;
            else
                matrix[i][j] = 0;
        }
    }
    full_matrix = matrix;

    auto empty_rows = std::remove_if(matrix.begin(), matrix.end(), [](const std::vector<size_t>& row)
    {
        return std::accumulate(row.begin(), row.end(), 0) == 0;
    });
    matrix.erase(empty_rows, matrix.end());
    for(size_t i = 0; i != matrix.size(); i++)
    {
        for(size_t j = 0; j != matrix[i].size(); j++)
        {
            std::cout << matrix[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

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

    matrix_size = matrix.size();
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

std::vector<std::vector<Frag>> get_frag_mix3(std::vector<std::vector<Frag>> res,
                            size_t total_resiudes, size_t residues, size_t frag_size,
                            const std::vector<std::vector<size_t>> &imatrix)
{
//    std::cout << residues << "x" << total_resiudes << std::endl;
//    std::cout << res.size() << "x" << res.front().size() << std::endl;
//    std::cin.get();
    std::vector<std::vector<Frag>> matrix;
    for(size_t i = 0, l = 0; i != imatrix.size(); i++)
    {
        std::vector<Frag> t(total_resiudes);
        bool is_pushed = false;
        for(size_t j = i, k = 0; j != frag_size + i; j++, k++)
        {
            if(imatrix[i][j] == 1)
            {
                t[j] = res[l][k];
                is_pushed = true;
            }
        }
        if(is_pushed)
        {
            matrix.push_back(t);
            l++;
        }
    }
//    std::cout << residues << "x" << total_resiudes << std::endl;
//    for(size_t i = 0; i != matrix.size(); i++)
//    {
//        for(size_t j = 0; j != matrix[i].size(); j++)
//        {
//            std::cout << matrix[i][j].aa << ' ';
//        }
//        std::cout << std::endl;
//    }

//    std::cout << "asv" << std::endl;
//    std::cin.get();

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
    peptide_centroid = _peptide;
    bbtools::to_centroid(peptide_centroid);
    peptide_seq = peptide.sequence();
}

void FragPick::set_file()
{
    if(basic::options::option[basic::options::OptionKeys::in::file::frag_files].user())
    {
        fragmet_file = basic::options::option[basic::options::OptionKeys::in::file::frag_files]()[1];
    }
    else
    {
        std::cout << "fatal" << std::endl;
    }

    frag_size = basic::options::option[basic::options::OptionKeys::frags::frag_sizes]()[1];
    std::cout << "frag_sizes " << frag_size << std::endl;
    std::cout << "fragmet_file " << fragmet_file << std::endl;
}

void FragPick::fill_grids()
{
//    phipsi_grid_number = std::vector<core::Size>(2*peptide_seq.size() - 2, phipsi_step);
//    omega_grid_number = std::vector<core::Size>(peptide_seq.size() - 1, omega_step);
    phipsi_grid_number.resize(2*peptide_seq.size() - 2);
    phipsi_grid_number.front() = step_num_phipsi.front();
    for(size_t i = 1, k = 1; k != step_num_phipsi.size() - 1; i+=2, k++)
    {
        phipsi_grid_number[i] = step_num_phipsi[k];
        phipsi_grid_number[i + 1] = step_num_phipsi[k];
    }
    phipsi_grid_number.back() = step_num_phipsi.back();
    omega_grid_number.resize(peptide_seq.size() - 1);
    for(size_t i = 0; i != omega_grid_number.size(); i++)
    {
        omega_grid_number[i] = step_num_omega[i];
    }

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

void FragPick::all_possible(size_t residues, size_t n_frags)
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
//    std::vector<std::vector<int>> permut = iterate(variable_values);    // n_frags^residues
    //        structures.clear();
    std::vector<size_t> it(variable_values.size(), 0);
    do
    {
        std::vector<int> permut = get_line(variable_values, it);
        std::vector<std::vector<Frag>> to_distr;
        for(size_t j = 0; j != permut.size(); j++)
        {
            to_distr.push_back(all_fragments[j][permut[j]]);
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

            std::vector<core::Real> phipsi_values_radians(phipsi_grid_number.size());
            std::vector<core::Real> omega_values_radians(omega_grid_number.size());
            std::vector<std::uint8_t> to_trie(phipsi_values_radians.size() + omega_values_radians.size());

            phipsi_values_radians.front() = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp.front().psi));
            omega_values_radians.front() = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp.front().omg));
            for(size_t j = 1, k = 1; j != temp.size() - 1; j++, k+=2)
            {
                phipsi_values_radians[k] = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp[j].phi));
                phipsi_values_radians[k + 1] = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp[j].psi));
                omega_values_radians[j] = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp[j].omg));
            }
            phipsi_values_radians.back() = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp.back().phi));

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
    while(increase(variable_values, it));

//    std::cout << permut.size() << std::endl;
    //        std::cout << structures.size() << std::endl;
}

void FragPick::one_chain(size_t residues, size_t n_frags)
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
//    std::vector<std::vector<int>> permut = iterate(variable_values);    // n_frags^residues
//    std::cout << permut.size() << std::endl;
    std::cout << all_fragments.size() << std::endl;
    //        structures.clear();
    std::ofstream fOut("sample.dat");
    std::set<std::vector<std::uint8_t>> sample;
    size_t min_sum = 256*native_state.size();
    auto bonds = get_bounds();
    size_t count = 0;
    std::vector<size_t> it(variable_values.size(), 0);
    do
    {
        count++;
        if(!(count%1000000))
            std::cout << "count " << count << '\t' << sample.size() << '\t' << min_sum << std::endl;
        std::vector<int> permut = get_line(variable_values, it);
        std::vector<std::vector<Frag>> to_distr;
        size_t delta = frag_size;
        for(size_t j = 0, t = 0; j != permut.size(); j++, t+=delta)
        {
            to_distr.push_back(all_fragments[t][permut[j]]);
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

            std::vector<core::Real> phipsi_values_radians(phipsi_grid_number.size());
            std::vector<core::Real> omega_values_radians(omega_grid_number.size());
            std::vector<std::uint8_t> to_trie(phipsi_values_radians.size() + omega_values_radians.size());

            phipsi_values_radians.front() = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp.front().psi));
            omega_values_radians.front() = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp.front().omg));
            for(size_t j = 1, k = 1; j != temp.size() - 1; j++, k+=2)
            {
                phipsi_values_radians[k] = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp[j].phi));
                phipsi_values_radians[k + 1] = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp[j].psi));
                omega_values_radians[j] = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp[j].omg));
            }
            phipsi_values_radians.back() = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp.back().phi));

            for(size_t j = 0; j != phipsi_values_radians.size(); j++)
            {
                to_trie[j] = get_index_phipsi(phipsi_values_radians[j], j);
            }
            for(size_t j = phipsi_values_radians.size(), k = 0; j != to_trie.size(); j++, k++)
            {
                to_trie[j] = get_index_omega(omega_values_radians[k], k);
            }

            ///

            /*peptide_centroid.set_psi(1, temp.front().psi);
            peptide_centroid.set_omega(1, temp.front().omg);

            for(size_t j = 1; j != temp.size() - 1; j++)
            {
                peptide_centroid.set_phi(j + 1, temp[j].phi);
                peptide_centroid.set_psi(j + 1, temp[j].psi);
                peptide_centroid.set_omega(j + 1, temp[j].omg);
            }
            peptide_centroid.set_phi(peptide_centroid.total_residue(), temp.back().phi);

            //            core::scoring::ScoreFunctionOP score_cent = core::scoring::ScoreFunctionFactory::create_score_function("cen_std");
            //            core::kinematics::MoveMapOP movemap_minimizer_ = core::kinematics::MoveMapOP(new core::kinematics::MoveMap());
            //            movemap_minimizer_->set_bb(true);
            //            protocols::minimization_packing::MinMover minimizer(movemap_minimizer_, score_cent, "lbfgs", 1e-20, true);
            //            minimizer.max_iter(5000);
            //            minimizer.apply(peptide_centroid);

            //peptide_centroid.dump_pdb("1.pdb");
            //std::cin.get();
            //            numeric::Real current_distance(peptide_centroid.residue(1).xyz("CA").distance( peptide_centroid.residue(peptide_centroid.total_residue()).xyz("CA") ));
            //            if(current_distance > 10)
            //                continue;
            //            std::cout << "hit it" << std::endl;
            //            peptide_centroid.dump_pdb("1.pdb");
            //            std::cin.get();

            std::vector<core::Real> phipsi_values_radians(phipsi_grid_number.size());
            std::vector<core::Real> omega_values_radians(omega_grid_number.size());
            std::vector<std::uint8_t> to_trie(phipsi_values_radians.size() + omega_values_radians.size());

            phipsi_values_radians.front() = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(peptide_centroid.psi(1)));
            omega_values_radians.front() = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(peptide_centroid.omega(1)));
            for(size_t j = 1, k = 1, seqpos = 2; j != temp.size() - 1; j++, k+=2, seqpos++)
            {
                phipsi_values_radians[k] = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(peptide_centroid.phi(seqpos)));
                phipsi_values_radians[k + 1] = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(peptide_centroid.psi(seqpos)));
                omega_values_radians[j] = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(peptide_centroid.omega(seqpos)));
            }
            phipsi_values_radians.back() = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(peptide_centroid.phi(peptide_centroid.total_residue())));

            for(size_t j = 0; j != phipsi_values_radians.size(); j++)
            {
                to_trie[j] = get_index_phipsi(phipsi_values_radians[j], j);
            }
            for(size_t j = phipsi_values_radians.size(), k = 0; j != to_trie.size(); j++, k++)
            {
                to_trie[j] = get_index_omega(omega_values_radians[k], k);
            }*/

            ///


//            if(!structures_trie->search(to_trie))
            if(sample.find(to_trie) == sample.end())
            {
//                structures_trie->insert(to_trie);
                sample.insert(to_trie);

                //
                size_t sum = 0;
                std::vector<std::uint8_t> dist(native_state.size());
                for(size_t k = 0; k != dist.size(); k++)
                {
                    int t = std::abs(static_cast<int>(native_state[k]) - static_cast<int>(to_trie[k]));
                    sum += std::min(t, std::abs(int(bonds[k]) - t));
                }
                if(min_sum > sum)
                {
                    min_sum = sum;
                    closest = to_trie;
                    std::cout << min_sum << std::endl;
                }

//                for(const auto &k : to_trie)
//                {
//                    fOut << int(k) << '\t';
//                }
//                fOut << std::endl;
            }
//            if(structures_trie->get_total_count() > 1e6)
//                break;
        }
//        if(structures_trie->get_total_count() > 1e6)
//            break;
    }
    while(increase(variable_values, it));
    fOut.close();
//    std::cout << permut.size() << std::endl;
//    std::cout << structures.size() << std::endl;
}

void FragPick::one_chain_Von_Neumann(size_t residues, size_t n_frags)
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
//    std::vector<std::vector<int>> permut = iterate(variable_values);    // n_frags^residues
//    std::cout << permut.size() << std::endl;
    std::cout << all_fragments.size() << std::endl;
    //        structures.clear();
    std::ofstream fOut("sample.dat");
    std::set<std::vector<std::uint8_t>> sample;
    size_t min_sum = 256*native_state.size();
    auto bonds = get_bounds();
    size_t count = 0;
    std::vector<size_t> it(variable_values.size(), 0);
    do
    {
        count++;
        if(!(count%1000000))
            std::cout << "count " << count << '\t' << sample.size() << '\t' << min_sum << std::endl;
        std::vector<int> permut = get_line(variable_values, it);
        std::vector<std::vector<Frag>> to_distr;
        size_t delta = frag_size;
        for(size_t j = 0, t = 0; j != permut.size(); j++, t+=delta)
        {
            to_distr.push_back(all_fragments[t][permut[j]]);
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

            std::vector<core::Real> phipsi_values_radians(phipsi_grid_number.size());
            std::vector<core::Real> omega_values_radians(omega_grid_number.size());
            std::vector<std::uint8_t> to_trie(phipsi_values_radians.size() + omega_values_radians.size());

            phipsi_values_radians.front() = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp.front().psi));
            omega_values_radians.front() = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp.front().omg));
            for(size_t j = 1, k = 1; j != temp.size() - 1; j++, k+=2)
            {
                phipsi_values_radians[k] = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp[j].phi));
                phipsi_values_radians[k + 1] = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp[j].psi));
                omega_values_radians[j] = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp[j].omg));
            }
            phipsi_values_radians.back() = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp.back().phi));

            for(size_t j = 0; j != phipsi_values_radians.size(); j++)
            {
                to_trie[j] = get_index_phipsi(phipsi_values_radians[j], j);
            }
            for(size_t j = phipsi_values_radians.size(), k = 0; j != to_trie.size(); j++, k++)
            {
                to_trie[j] = get_index_omega(omega_values_radians[k], k);
            }

            auto init_point = to_trie;
            auto point = to_trie;
            for(size_t i = 0; i != init_point.size(); i++)
            {
                point = init_point;
                point[i] = point[i] + 1;

                if(point[i] > static_cast<std::uint8_t>(bonds[i] - 1))
                {
                    point[i] = 0;
                }
                if(!structures_trie->search(point))
//            if(sample.find(point) == sample.end())
                {
                    structures_trie->insert(point);
//                sample.insert(point);

                    //
                    size_t sum = 0;
                    std::vector<std::uint8_t> dist(native_state.size());
                    for(size_t k = 0; k != dist.size(); k++)
                    {
                        int t = std::abs(static_cast<int>(native_state[k]) - static_cast<int>(point[k]));
                        sum += std::min(t, std::abs(int(bonds[k]) - t));
                    }
                    if(min_sum > sum)
                    {
                        min_sum = sum;
                        closest = point;
                        std::cout << min_sum << std::endl;
                    }

//                for(const auto &k : point)
//                {
//                    fOut << int(k) << '\t';
//                }
//                fOut << std::endl;
                }
            }
            for(size_t i = 0; i != point.size(); i++)
            {
                point = init_point;
                if(point[i] != 0)
                    point[i] = point[i] - 1;
                else
                    point[i] = point[i] - 1;


                if(!structures_trie->search(point))
//            if(sample.find(point) == sample.end())
                {
                    structures_trie->insert(point);
//                sample.insert(point);

                    //
                    size_t sum = 0;
                    std::vector<std::uint8_t> dist(native_state.size());
                    for(size_t k = 0; k != dist.size(); k++)
                    {
                        int t = std::abs(static_cast<int>(native_state[k]) - static_cast<int>(point[k]));
                        sum += std::min(t, std::abs(int(bonds[k]) - t));
                    }
                    if(min_sum > sum)
                    {
                        min_sum = sum;
                        closest = point;
                        std::cout << min_sum << std::endl;
                    }

//                for(const auto &k : point)
//                {
//                    fOut << int(k) << '\t';
//                }
//                fOut << std::endl;
                }
            }

            if(!structures_trie->search(to_trie))
//            if(sample.find(to_trie) == sample.end())
            {
                structures_trie->insert(to_trie);
//                sample.insert(to_trie);

                //
                size_t sum = 0;
                std::vector<std::uint8_t> dist(native_state.size());
                for(size_t k = 0; k != dist.size(); k++)
                {
                    int t = std::abs(static_cast<int>(native_state[k]) - static_cast<int>(to_trie[k]));
                    sum += std::min(t, std::abs(int(bonds[k]) - t));
                }
                if(min_sum > sum)
                {
                    min_sum = sum;
                    closest = to_trie;
                    std::cout << min_sum << std::endl;
                }

//                for(const auto &k : to_trie)
//                {
//                    fOut << int(k) << '\t';
//                }
//                fOut << std::endl;
            }
        }
    }
    while(increase(variable_values, it));
    fOut.close();
//    std::cout << permut.size() << std::endl;
//    std::cout << structures.size() << std::endl;
}


void FragPick::coil_chain(size_t residues, size_t n_frags)
{
    size_t col_max = 2;
    size_t matrix_size = 0;
    std::vector<std::vector<size_t>> full_matrix;
    std::vector<size_t> sizes = get_permut_sizes3(peptide.total_residue(), residues, frag_size, ss_predicted, col_max, matrix_size, full_matrix);

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
//    for(size_t i = 0; i != variable_values.size(); i++)
//    {
//        for(size_t j = 0; j != n_frags; j++)
//        {
//            variable_values[i][j] = j;
//            std::cout << variable_values[i][j] << ' ';
//        }
//        std::cout << std::endl;
//    }
//    std::cin.get();
    size_t count = 0;
//    std::vector<std::vector<int>> permut = iterate(variable_values);    // n_frags^residues
    //        structures.clear();
    std::set<std::vector<std::uint8_t>> sample;
    size_t min_sum = 256*native_state.size();
    auto bonds = get_bounds();
    std::vector<size_t> it(variable_values.size(), 0);
    do
    {
        std::vector<int> permut = get_line(variable_values, it);
        std::vector<std::vector<Frag>> to_distr;
        for(size_t j = 0; j != permut.size(); j++)
        {
            to_distr.push_back(all_fragments[j][permut[j]]);
        }
//        std::cout << "to_distr " << to_distr.size() << std::endl;
//        std::cin.get();
        auto distr = get_frag_mix3(to_distr, peptide.total_residue(), residues, frag_size, full_matrix);

        for(size_t j = 0; j != permut_for_selected.size(); j++)
        {
            std::vector<Frag> temp;
            for(size_t k = 0; k != permut_for_selected[j].size(); k++)
            {
                temp.push_back(distr[k][permut_for_selected[j][k]]);
            }
            //                structures.push_back(temp);

            std::vector<core::Real> phipsi_values_radians(phipsi_grid_number.size());
            std::vector<core::Real> omega_values_radians(omega_grid_number.size());
            std::vector<std::uint8_t> to_trie(phipsi_values_radians.size() + omega_values_radians.size());

            phipsi_values_radians.front() = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp.front().psi));
            omega_values_radians.front() = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp.front().omg));
            for(size_t j = 1, k = 1; j != temp.size() - 1; j++, k+=2)
            {
                phipsi_values_radians[k] = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp[j].phi));
                phipsi_values_radians[k + 1] = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp[j].psi));
                omega_values_radians[j] = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp[j].omg));
            }
            phipsi_values_radians.back() = bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(temp.back().phi));

            for(size_t j = 0; j != phipsi_values_radians.size(); j++)
            {
                to_trie[j] = get_index_phipsi(phipsi_values_radians[j], j);
            }
            for(size_t j = phipsi_values_radians.size(), k = 0; j != to_trie.size(); j++, k++)
            {
                to_trie[j] = get_index_omega(omega_values_radians[k], k);
            }
            count++;
            if(!(count%1000000))
                std::cout << "count " << count << '\t' << sample.size() << '\t' << min_sum << std::endl;
//            if(!structures_trie->search(to_trie))
            if(sample.find(to_trie) == sample.end())
            {
//                structures_trie->insert(to_trie);
                sample.insert(to_trie);

                size_t sum = 0;
                std::vector<std::uint8_t> dist(native_state.size());
                for(size_t k = 0; k != dist.size(); k++)
                {
                    int t = std::abs(static_cast<int>(native_state[k]) - static_cast<int>(to_trie[k]));
                    sum += std::min(t, std::abs(int(bonds[k]) - t));
                }
                if(min_sum > sum)
                {
                    min_sum = sum;
                    closest = to_trie;
                    std::cout << min_sum << std::endl;
                }
            }
        }
    }
    while(increase(variable_values, it));
//    std::cout << permut.size() << std::endl;
    std::cout << "count " << count << std::endl;
    //        std::cout << structures.size() << std::endl;
}

void FragPick::make_permutations(size_t type)
{
    size_t residues = all_fragments.size(); // not total_residues!
    size_t n_frags = all_fragments.front().size(); // n_frags ^ residues

    switch(type)
    {
        case 0:
            all_possible(residues, n_frags);
            break;
        case 1:
            one_chain(residues, n_frags);
            break;
        case 2:
            one_chain_Von_Neumann(residues, n_frags);
            break;
        case 3:
            coil_chain(residues, n_frags);
            break;
        default:
            break;
    }
}

std::vector<size_t> FragPick::get_bounds()
{
    std::cout << "get_total_count " <<  structures_trie->get_total_count() << std::endl;
    std::vector<size_t> phipsiomg_grid_number;

    phipsiomg_grid_number.insert(phipsiomg_grid_number.end(), phipsi_grid_number.begin(), phipsi_grid_number.end());
    phipsiomg_grid_number.insert(phipsiomg_grid_number.end(), omega_grid_number.begin(), omega_grid_number.end());

    std::cout << phipsiomg_grid_number.size() << std::endl;
    return phipsiomg_grid_number;
}

void FragPick::set_psipred(std::pair<std::uint8_t,std::uint8_t> phipsi_minmax, std::pair<std::uint8_t,std::uint8_t> omega_minmax)
{
    std::string horiz_fname;
    if(basic::options::option[basic::options::OptionKeys::in::file::psipred_ss2].user())
    {
        horiz_fname = basic::options::option[basic::options::OptionKeys::in::file::psipred_ss2].value_string();
        std::cout << "psipred.horiz " << horiz_fname << std::endl;
    }

//    confidence.resize(ss_profile.total_residue());
//    for(core::Size i = 1; i <= ss_profile.total_residue(); i++)
//    {
//        confidence[i - 1] = ss_profile.helix_fraction(i);
//        if(ss_profile.sheet_fraction(i) > confidence[i - 1])
//            confidence[i - 1] = ss_profile.sheet_fraction(i);
//        if(ss_profile.loop_fraction(i) > confidence[i - 1])
//            confidence[i - 1] = ss_profile.loop_fraction(i);
//    }

    ss_predicted.clear();
    confidence.clear();
    std::ifstream fPsiHoriz;
    fPsiHoriz.open(horiz_fname.c_str());
    if(!fPsiHoriz.is_open())
    {
        std::cout << "not open" << std::endl;
    }
    std::string line;
    while(getline(fPsiHoriz, line))
    {
        if(line.find("#") != std::string::npos || line.empty())
            continue;
        if(line.find("Conf") != std::string::npos)
        {
            for(size_t i = 0; i != line.size(); i++)
            {
                if(std::isdigit(line[i]) == 0)
                    continue;
                confidence.push_back(std::stoi(line.substr(i, 1)));
            }
        }
        if(line.find("Pred") != std::string::npos)
        {
            for(size_t i = 0; i != line.size(); i++)
            {
                if(line[i] == 'C' || line[i] == 'H' || line[i] == 'E')
                {
                    ss_predicted.push_back(line[i]);
                }
            }
        }
    }
    fPsiHoriz.close();

    for(const auto &i : ss_predicted)
    {
        std::cout << i;
    }
    std::cout << std::endl;
    for(const auto &i : peptide_seq)
    {
        std::cout << i;
    }
    std::cout << std::endl;
    for(const auto &i : confidence)
    {
        std::cout << i;
    }
    std::cout << std::endl;
    if(confidence.size() != peptide.total_residue())
    {
        std::cout << "confidence size != peptide total residue" << std::endl;
    }

    double phipsi = phipsi_minmax.second - phipsi_minmax.first;
    double omega = omega_minmax.second - omega_minmax.first;
    step_num_phipsi.resize(confidence.size());
    step_num_omega.resize(confidence.size());
    for(size_t i = 0; i != confidence.size(); i++)
    {
        step_num_phipsi[i] = phipsi_minmax.first + phipsi*(confidence[i] + 1)/10.0;
        if(ss_predicted[i] == 'C')
        {
//            step_num_phipsi[i] /= 2;
//            if(i > 1 && i < confidence.size() - 2)
//            {
//                step_num_phipsi[i] /= 2;
//                step_num_phipsi[i] = 1;
//            }
//            step_num_phipsi[i] = 36;
        }
        else
        {
            /*if(i > 1 && i < confidence.size() - 2)
            {
                if(ss_predicted[i + 1] == 'C')
                {
                    step_num_phipsi[i] /= 2;
                }
                if(ss_predicted[i - 1] == 'C')
                {
                    step_num_phipsi[i] /= 2;
                }
            }*/
        }
        step_num_omega[i] = omega_minmax.first + omega*(confidence[i] + 1)/10.0;
        std::cout << step_num_phipsi[i] << '\t' << step_num_omega[i] << std::endl;
    }
}

void FragPick::set_native_state(const std::vector<std::uint8_t> &nt)
{
    native_state = nt;
}
std::vector<std::uint8_t> FragPick::get_closest_to_native()
{
    return closest;
}

std::string FragPick::get_ss_predicted() const
{
    return ss_predicted;
}

}
}
