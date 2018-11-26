
//
// Copyright (c) 2012 Ronaldo Carpio
//                                     
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee,
// provided that the above copyright notice appear in all copies and   
// that both that copyright notice and this permission notice appear
// in supporting documentation.  The authors make no representations
// about the suitability of this software for any purpose.          
// It is provided "as is" without express or implied warranty.
//               

/*
This is a C++ header-only library for N-dimensional linear interpolation on a rectangular grid. Implements two methods:
* Multilinear: Interpolate using the N-dimensional hypercube containing the point. Interpolation step is O(2^N) 
* Simplicial: Interpolate using the N-dimensional simplex containing the point. Interpolation step is O(N log N), but less accurate.
Requires boost/multi_array library.

For a description of the algorithms, see:
* Weiser & Zarantonello (1988), "A Note on Piecewise Linear and Multilinear Table Interpolation in Many Dimensions", _Mathematics of Computation_ 50 (181), p. 189-196
* Davies (1996), "Multidimensional Triangulation and Interpolation for Reinforcement Learning", _Proceedings of Neural Information Processing Systems 1996_
*/

// Orgiginal file from linterp library: github.com/rncarpio/linterp/blob/master/src/linterp.h 
// Modified by Sergey Poluyan
// There are no significant changes. Everything except linear interpolation was removed. Some code simplification was done.

// Usage: 
//    std::vector<double> f_values = {0.5, 1.0, 0.75};
//    std::vector<std::vector<double>> grids = { {1, 2, 3} };
//    boost::array<int, 1> grid_sizes = { 3 };
//    auto grid_iter_list = linterp::get_begins_ends( grids.begin(), grids.end() );
//    auto obj = linterp::InterpMultilinear<1, double>(grid_iter_list.first.begin(), grid_sizes.begin(), f_values.data(), f_values.data() + f_values.size() );
//    boost::array<double, 1> args = { 1.25 };
//    std::cout << obj.interp( args.begin() ) << std::endl;

#ifndef linterp_hh
#define linterp_hh

#include <vector>
#include <array>

#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace linterp
{

template <int N, class T>
class NDInterpolator
{
public:
    typedef boost::numeric::ublas::array_adaptor<T> grid_type;
    typedef boost::const_multi_array_ref<T, N> array_type;
    typedef std::unique_ptr<array_type> array_type_ptr;

    array_type_ptr m_pF;

    std::vector<grid_type> m_grid_list;

    // constructors assume that [f_begin, f_end) is a contiguous array in C-order
    // non ref-counted constructor.
    template <class IterT1, class IterT2, class IterT3>
    NDInterpolator(IterT1 grids_begin, IterT2 grids_len_begin, IterT3 f_begin, IterT3 f_end)
    {
        init(grids_begin, grids_len_begin, f_begin, f_end);
    }

    template <class IterT1, class IterT2, class IterT3>
    void init(IterT1 grids_begin, IterT2 grids_len_begin, IterT3 f_begin, IterT3 f_end)
    {
        set_grids(grids_begin, grids_len_begin);
        set_f_array(f_begin, f_end);
    }

    template <class IterT1, class IterT2>
    void set_grids(IterT1 grids_begin, IterT2 grids_len_begin)
    {
        m_grid_list.clear();
        for(int i=0; i<N; i++)
        {
            int gridLength = grids_len_begin[i];
            T const *grid_ptr = &(*grids_begin[i]);
            m_grid_list.push_back(grid_type(gridLength, (T*) grid_ptr));
        }
    }

    // assumes that [f_begin, f_end) is a contiguous array in C-order
    template <class IterT>
    void set_f_array(IterT f_begin, IterT f_end)
    {
        unsigned int nGridPoints = 1;
        std::array<int,N> sizes;
        for(unsigned int i=0; i<m_grid_list.size(); i++)
        {
            sizes[i] = m_grid_list[i].size();
            nGridPoints *= sizes[i];
        }
        int f_len = f_end - f_begin;
        if(f_len != static_cast<int>(nGridPoints))
        {
            throw std::invalid_argument("f has wrong size");
        }
        m_pF.reset(new array_type(f_begin, sizes));

    }

    // -1 is before the first grid point
    // N-1 (where grid.size() == N) is after the last grid point
    int find_cell(int dim, T x) const
    {
        grid_type const &grid(m_grid_list[dim]);
        if(x < *(grid.begin())) return -1;
        else if(x >= *(grid.end()-1)) return grid.size()-1;
        else
        {
            auto i_upper = std::upper_bound(grid.begin(), grid.end(), x);
            return i_upper - grid.begin() - 1;
        }
    }

    // return the value of f at the given cell and vertex
    T get_f_val(std::array<int,N> const &cell_index, std::array<int,N> const &v_index) const
    {
        std::array<int,N> f_index;

        for(int i=0; i<N; i++)
        {
            if(cell_index[i] < 0)
            {
                f_index[i] = 0;
            }
            else if(cell_index[i] >= static_cast<int>(m_grid_list[i].size()-1))
            {
                f_index[i] = m_grid_list[i].size()-1;
            }
            else
            {
                f_index[i] = cell_index[i] + v_index[i];
            }
        }
        return (*m_pF)(f_index);
    }

    T get_f_val(std::array<int,N> const &cell_index, int v) const
    {
        std::array<int,N> v_index;
        for(int dim=0; dim<N; dim++)
        {
            v_index[dim] = (v >> (N-dim-1)) & 1;						// test if the i-th bit is set
        }
        return get_f_val(cell_index, v_index);
    }
};

template <int N, class T>
class InterpMultilinear : public NDInterpolator<N, T>
{
public:
    typedef NDInterpolator<N, T> super;

    template <class IterT1, class IterT2, class IterT3>
    InterpMultilinear(IterT1 grids_begin, IterT2 grids_len_begin, IterT3 f_begin, IterT3 f_end)
        : super(grids_begin, grids_len_begin, f_begin, f_end)
    {}

    template <class IterT1, class IterT2>
    static T linterp_nd_unitcube(IterT1 f_begin, IterT1 f_end, IterT2 xi_begin, IterT2 xi_end)
    {
        int n = xi_end - xi_begin;
        int f_len = f_end - f_begin;
        assert(1 << n == f_len);
        T sub_lower, sub_upper;
        if(n == 1)
        {
            sub_lower = f_begin[0];
            sub_upper = f_begin[1];
        }
        else
        {
            sub_lower = linterp_nd_unitcube(f_begin, f_begin + (f_len/2), xi_begin + 1, xi_end);
            sub_upper = linterp_nd_unitcube(f_begin + (f_len/2), f_end, xi_begin + 1, xi_end);
        }
        T result = sub_lower + (*xi_begin)*(sub_upper - sub_lower);
        return result;
    }

    template <class IterT>
    T interp(IterT x_begin) const
    {
        std::array<T,1> result;
        std::array< std::array<T,1>, N > coord_iter;
        for(int i=0; i<N; i++)
        {
            coord_iter[i][0] = x_begin[i];
        }
        interp_vec(1, coord_iter.begin(), coord_iter.end(), result.begin());
        return result[0];
    }

    template <class IterT1, class IterT2>
    void interp_vec(int n, IterT1 coord_iter_begin, IterT1 coord_iter_end, IterT2 i_result) const
    {
        assert(N == coord_iter_end - coord_iter_begin);
        std::array<int,N> index;
        int c;
        T y;
        std::vector<T> f(1 << N);
        std::array<T,N> x;

        for(int i=0; i<n; i++)  								// loop over each point
        {
            for(int dim=0; dim<N; dim++)  						// loop over each dimension
            {
                auto const &grid(super::m_grid_list[dim]);
                c = this->find_cell(dim, coord_iter_begin[dim][i]);
                if(c == -1)  					// before first grid point
                {
                    y = 1.0;
                }
                else if(c == static_cast<int>(grid.size()-1))  	// after last grid point
                {
                    y = 0.0;
                }
                else
                {
                    y = (coord_iter_begin[dim][i] - grid[c]) / (grid[c + 1] - grid[c]);
                    if(y < 0.0) y=0.0;
                    else if(y > 1.0) y=1.0;
                }
                index[dim] = c;
                x[dim] = y;
            }
            // copy f values at vertices
            for(int v=0; v < (1 << N); v++)  					// loop over each vertex of hypercube
            {
                f[v] = this->get_f_val(index, v);
            }
            *i_result++ = linterp_nd_unitcube(f.begin(), f.end(), x.begin(), x.end());
        }
    }
};

template <class IterT>
std::pair<std::vector<typename IterT::value_type::const_iterator>, std::vector<typename IterT::value_type::const_iterator> > get_begins_ends(IterT iters_begin, IterT iters_end)
{
    typedef typename IterT::value_type T;
    typedef std::vector<typename T::const_iterator> VecT;
    int N = iters_end - iters_begin;
    std::pair<VecT, VecT> result;
    result.first.resize(N);
    result.second.resize(N);
    for(int i=0; i<N; i++)
    {
        result.first[i] = iters_begin[i].begin();
        result.second[i] = iters_begin[i].end();
    }
    return result;
}

}


#endif //_linterp_h