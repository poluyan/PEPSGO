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

#ifndef QUANTILE_HH
#define QUANTILE_HH

#include <trie_based.hh>

#include <numeric>

namespace empirical_quantile
{

template <typename T, typename U>
class Quantile
{
protected:
    std::vector<U> lb;
    std::vector<U> ub;
    std::vector<U> dx;
    std::vector<size_t> grid_number;
    std::vector<std::vector<U>> grids;
public:
    Quantile();
    Quantile(std::vector<U> in_lb, std::vector<U> in_ub, std::vector<size_t> in_gridn);
    void set_grid_and_gridn(std::vector<U> in_lb, std::vector<U> in_ub, std::vector<size_t> in_gridn);
};

template <typename T, typename U>
Quantile<T, U>::Quantile()
{
}

template <typename T, typename U>
Quantile<T, U>::Quantile(std::vector<U> in_lb,
                         std::vector<U> in_ub,
                         std::vector<size_t> in_gridn)
{
    lb = in_lb;
    ub = in_ub;
    grid_number = in_gridn;

    dx.resize(grid_number.size());
    grids.resize(grid_number.size());
    for(size_t i = 0; i != grids.size(); i++)
    {
        std::vector<U> grid(grid_number[i] + 1);
        U startp = lb[i];
        U endp = ub[i];
        U es = endp - startp;
        for(size_t j = 0; j != grid.size(); j++)
        {
            grid[j] = startp + j*es/U(grid_number[i]);
        }
        grids[i] = grid;
        dx[i] = es/(U(grid_number[i])*2);
    }
}

template <typename T, typename U>
void Quantile<T, U>::set_grid_and_gridn(std::vector<U> in_lb, std::vector<U> in_ub, std::vector<size_t> in_gridn)
{
    Quantile<T, U>(in_lb, in_ub, in_gridn);
}


template <typename T, typename U>
class ExplicitQuantile : public Quantile<T, U>
{
protected:
    using Quantile<T, U>::grids;
    using Quantile<T, U>::dx;

    typedef std::vector<std::vector<U>> sample_type;
    std::shared_ptr<sample_type> sample;

    size_t count_less(const std::vector<U> &layer, U target) const;
    std::pair<size_t, U> quantile_transform(const std::vector<U> &layer, size_t ind, U val01) const;
public:
    ExplicitQuantile();
    ExplicitQuantile(std::vector<U> in_lb, std::vector<U> in_ub, std::vector<size_t> in_gridn);
    using Quantile<T, U>::set_grid_and_gridn;
    void set_sample(const std::vector<std::vector<T>> &in_sample);
    void set_sample_shared(std::shared_ptr<sample_type> in_sample);
    void transform(const std::vector<U>& in01, std::vector<U>& out) const;
};

template <typename T, typename U>
ExplicitQuantile<T, U>::ExplicitQuantile(): Quantile<T, U>()
{
}

template <typename T, typename U>
ExplicitQuantile<T, U>::ExplicitQuantile(std::vector<U> in_lb,
        std::vector<U> in_ub,
        std::vector<size_t> in_gridn): Quantile<T, U>(in_lb, in_ub, in_gridn)
{
}

template <typename T, typename U>
void ExplicitQuantile<T, U>::set_sample(const std::vector<std::vector<T>> &in_sample)
{
    sample = std::make_shared< std::vector<std::vector<U>> >();
    for(size_t i = 0; i != in_sample.size(); ++i)
    {
        std::vector<U> temp;
        for(size_t j = 0; j != in_sample[i].size(); ++j)
        {
            temp.push_back(grids[j][in_sample[i][j]] + dx[j]);
        }
        sample->push_back(temp);
    }
}
template <typename T, typename U>
void ExplicitQuantile<T, U>::set_sample_shared(std::shared_ptr<sample_type> in_sample)
{
    sample = std::move(in_sample);
}

template <typename T, typename U>
void ExplicitQuantile<T, U>::transform(const std::vector<U>& in01, std::vector<U>& out) const
{
    std::vector<size_t> m(grids.size());
    for(size_t i = 0, g = in01.size(); i != g; i++)
    {
        std::vector<U> row(sample->size());
        size_t index = 0;
        for(size_t j = 0, n = sample->size(); j != n; j++)
        {
            bool flag = true;
            for(size_t k = 0; k != i; k++)
            {
                if(!((*sample)[j][k] > grids[k][m[k]] && (*sample)[j][k] < grids[k][m[k] + 1]))
                {
                    flag = false;
                    break;
                }
            }
            if(flag)
            {
                row[index] = (*sample)[j][i];
                ++index;
            }
        }
        row.resize(index);

        auto rez = quantile_transform(row, i, in01[i]);
        out[i] = rez.second;
        m[i] = rez.first;
    }
}

template <typename T, typename U>
size_t ExplicitQuantile<T, U>::count_less(const std::vector<U> &layer, U target) const
{
    return std::count_if(layer.begin(), layer.end(), [&](const U &v)
    {
        return v < target;
    });
}

template <typename T, typename U>
std::pair<size_t, U> ExplicitQuantile<T, U>::quantile_transform(const std::vector<U> &layer, size_t ind, U val01) const
{
    size_t count = grids[ind].size() - 1, step, c1 = 0, c2 = 0, m = 0;
    U f1 = 0.0, f2 = 0.0, n = layer.size();
    auto first = grids[ind].begin();
    auto it = grids[ind].begin();
    while(count > 0)
    {
        it = first;
        step = count / 2;
        std::advance(it, step);
        m = std::distance(grids[ind].begin(), it);

        c1 = count_less(layer, grids[ind][m]);
        f1 = c1/n;

        if(f1 < val01)
        {
            c2 = count_less(layer, grids[ind][m + 1]);
            f2 = c2/n;

            if(val01 < f2)
                break;

            first = ++it;
            count -= step + 1;
        }
        else
            count = step;
    }

    if(count == 0)
    {
        c2 = count_less(layer, grids[ind][m + 1]);
        f2 = c2/n;
    }

    if(c1 == c2)
    {
        if(c1 == 0)
        {
            auto min_val = *std::min_element(layer.begin(), layer.end()) - 2.0*dx[ind];
            auto lb_min = std::lower_bound(grids[ind].begin(), grids[ind].end(), min_val);
            size_t min_ind = std::distance(grids[ind].begin(), lb_min);
            return std::make_pair(min_ind, grids[ind][min_ind] + 2.0*val01*dx[ind]);
        }
        if(c1 == layer.size())
        {
            auto max_val = *std::max_element(layer.begin(), layer.end()) - 2.0*dx[ind];
            auto lb_max = std::lower_bound(grids[ind].begin(), grids[ind].end(), max_val);
            size_t max_ind = std::distance(grids[ind].begin(), lb_max);
            return std::make_pair(max_ind, grids[ind][max_ind] + 2.0*val01*dx[ind]);
        }

        U target = grids[ind][m];
        U diff = std::numeric_limits<U>::max();
        size_t index = 0;
        for(size_t i = 1; i != layer.size(); ++i)
        {
            U curr = std::abs(layer[i] - target);
            if(diff > curr)
            {
                diff = curr;
                index = i;
            }
        }

        auto lb = std::lower_bound(grids[ind].begin(), grids[ind].end(), layer[index] - 2.0*dx[ind]);
        size_t lb_ind = std::distance(grids[ind].begin(), lb);
        return std::make_pair(lb_ind, grids[ind][lb_ind] + 2.0*val01*dx[ind]);
    }
    return std::make_pair(m, grids[ind][m] + (val01 - f1) * (grids[ind][m + 1] - grids[ind][m]) / (f2 - f1));
}



template <typename T, typename U>
class ImplicitQuantile : public Quantile<T, U>
{
protected:
    typedef trie_based::TrieBased<trie_based::NodeCount<T>,T> sample_type;
    std::shared_ptr<sample_type> sample;

    using Quantile<T, U>::grids;
    using Quantile<T, U>::dx;

    size_t count_less(trie_based::NodeCount<T> *layer, size_t m) const;
    std::pair<size_t, U> quantile_transform(trie_based::NodeCount<T> *layer, size_t ind, U val01) const;
public:
    ImplicitQuantile();
    ImplicitQuantile(std::vector<U> in_lb, std::vector<U> in_ub, std::vector<size_t> in_gridn);
    void set_sample_and_fill_count(const std::vector<std::vector<T>> &in_sample);
    void set_sample_shared_and_fill_count(std::shared_ptr<sample_type> in_sample);
    void set_sample_shared(std::shared_ptr<sample_type> in_sample);
    void transform(const std::vector<U>& in01, std::vector<U>& out) const;
    size_t get_node_count() const;
    size_t get_link_count() const;
};

template <typename T, typename U>
ImplicitQuantile<T, U>::ImplicitQuantile()
{
}

template <typename T, typename U>
ImplicitQuantile<T, U>::ImplicitQuantile(std::vector<U> in_lb,
        std::vector<U> in_ub,
        std::vector<size_t> in_gridn) : Quantile<T, U>(in_lb, in_ub, in_gridn)
{

}

template <typename T, typename U>
size_t ImplicitQuantile<T, U>::get_node_count() const
{
    return sample->get_node_count();
}

template <typename T, typename U>
size_t ImplicitQuantile<T, U>::get_link_count() const
{
    return sample->get_link_count();
}

template <typename T, typename U>
void ImplicitQuantile<T, U>::set_sample_and_fill_count(const std::vector<std::vector<T>> &in_sample)
{
    sample = std::make_shared<sample_type>();
    sample->set_dimension(grids.size());
    for(const auto & i : in_sample)
        sample->insert(i);
    sample->fill_tree_count();
}

template <typename T, typename U>
void ImplicitQuantile<T, U>::set_sample_shared_and_fill_count(std::shared_ptr<sample_type> in_sample)
{
    sample = std::move(in_sample);
    sample->fill_tree_count();
}

template <typename T, typename U>
void ImplicitQuantile<T, U>::set_sample_shared(std::shared_ptr<sample_type> in_sample)
{
    sample = std::move(in_sample);
}

template <typename T, typename U>
size_t ImplicitQuantile<T, U>::count_less(trie_based::NodeCount<T> *layer, size_t r) const
{
    size_t c = 0;
    for(size_t i = 0; i != layer->children.size(); ++i)
    {
        if(static_cast<size_t>(layer->children[i]->index) < r)
        {
            c += layer->children[i]->count;
        }
    }
    return c;
}
template <typename T, typename U>
void ImplicitQuantile<T, U>::transform(const std::vector<U>& in01, std::vector<U>& out) const
{
    auto p = sample->root.get();
    for(size_t i = 0; i != in01.size(); i++)
    {
        auto rez = quantile_transform(p, i, in01[i]);
        out[i] = rez.second;
        p = p->children[rez.first].get();
    }
}
template <typename T, typename U>
std::pair<size_t, U> ImplicitQuantile<T, U>::quantile_transform(trie_based::NodeCount<T> *layer, size_t ind, U val01) const
{
    size_t m = 0, count = grids[ind].size() - 1, step, c1 = 0, c2 = 0;
    U f1 = 0.0, f2 = 0.0, sample_size_u = static_cast<U>(layer->count);
    auto first = grids[ind].begin();
    auto it = grids[ind].begin();

    while(count > 0)
    {
        it = first;
        step = count / 2;
        std::advance(it, step);
        m = std::distance(grids[ind].begin(), it);

        c1 = count_less(layer, m);
        f1 = c1/sample_size_u;

        if(f1 < val01)
        {
            c2 = count_less(layer, m + 1);
            f2 = c2/sample_size_u;

            if(val01 < f2)
                break;

            first = ++it;
            count -= step + 1;
        }
        else
            count = step;
    }
    if(count == 0)
    {
        c2 = count_less(layer, m + 1);
        f2 = c2/sample_size_u;
    }
    if(c1 == c2)
    {
        if(c1 == 0)
        {
            auto min_val_it = std::min_element(layer->children.begin(), layer->children.end(),
                                               [](const std::shared_ptr<trie_based::NodeCount<T>> &l,
                                                  const std::shared_ptr<trie_based::NodeCount<T>> &r)
            {
                return l->index < r->index;
            });
            size_t min_ind = std::distance(layer->children.begin(), min_val_it);
            return std::make_pair(min_ind, grids[ind][layer->children[min_ind]->index] + 2.0*val01*dx[ind]);
        }
        if(c1 == layer->count)
        {
            auto max_val_it = std::max_element(layer->children.begin(), layer->children.end(),
                                               [](const std::shared_ptr<trie_based::NodeCount<T>> &l,
                                                  const std::shared_ptr<trie_based::NodeCount<T>> &r)
            {
                return l->index < r->index;
            });
            size_t max_ind = std::distance(layer->children.begin(), max_val_it);
            return std::make_pair(max_ind, grids[ind][layer->children[max_ind]->index] + 2.0*val01*dx[ind]);
        }
        int diff = std::numeric_limits<int>::max();
        size_t index = 0;
        int min_ind = static_cast<int>(layer->children[index]->index);
        for(size_t i = 1; i != layer->children.size(); ++i)
        {
            int t = static_cast<int>(layer->children[i]->index);
            int curr = std::abs(t - static_cast<int>(m));
            if(diff > curr)
            {
                diff = curr;
                index = i;
                min_ind = t;
            }
            else if(diff == curr)
            {
                if(min_ind > t)
                {
                    min_ind = t;
                    index = i;
                }
            }
        }
        return std::make_pair(index, grids[ind][layer->children[index]->index] + 2.0*val01*dx[ind]);
    }
    size_t index = 0;
    T target = m;
    for(size_t j = 1; j < layer->children.size(); j++)
    {
        if(layer->children[j]->index == target)
        {
            index = j;
            break;
        }
    }
    return std::make_pair(index, grids[ind][m] + (val01 - f1) * (grids[ind][m + 1] - grids[ind][m]) / (f2 - f1));
}

template <typename T, typename U>
class ImplicitQuantileSorted : public ImplicitQuantile<T, U>
{
protected:
    using ImplicitQuantile<T, U>::grids;
    using ImplicitQuantile<T, U>::sample;
    using ImplicitQuantile<T, U>::dx;

    using sample_type = typename ImplicitQuantile<T, U>::sample_type;
    void sort_layer(trie_based::NodeCount<T> *p);

    size_t count_less_binary(trie_based::NodeCount<T> *layer, T target) const;
    std::pair<size_t, U> quantile_transform(trie_based::NodeCount<T> *layer, const std::vector<size_t> &psum, size_t ind, U val01) const;
public:
    ImplicitQuantileSorted();
    ImplicitQuantileSorted(std::vector<U> in_lb, std::vector<U> in_ub, std::vector<size_t> in_gridn);

    void set_sample_and_fill_count(const std::vector<std::vector<T>> &in_sample);
    void set_sample_shared_and_fill_count(std::shared_ptr<sample_type> in_sample);

    void sort();
    void transform(const std::vector<U>& in01, std::vector<U>& out) const;

};

template <typename T, typename U>
ImplicitQuantileSorted<T, U>::ImplicitQuantileSorted()
{
}

template <typename T, typename U>
ImplicitQuantileSorted<T, U>::ImplicitQuantileSorted(std::vector<U> in_lb,
        std::vector<U> in_ub,
        std::vector<size_t> in_gridn) : ImplicitQuantile<T, U>(in_lb, in_ub, in_gridn)
{
}

template <typename T, typename U>
void ImplicitQuantileSorted<T, U>::set_sample_and_fill_count(const std::vector<std::vector<T>> &in_sample)
{
    sample = std::make_shared<sample_type>();
    sample->set_dimension(grids.size());
    for(const auto & i : in_sample)
        sample->insert(i);
    sample->fill_tree_count();
    sort();
}


template <typename T, typename U>
void ImplicitQuantileSorted<T, U>::set_sample_shared_and_fill_count(std::shared_ptr<sample_type> in_sample)
{
    sample = std::move(in_sample);
    sample->fill_tree_count();
    sort();
}

template <typename T, typename U>
void ImplicitQuantileSorted<T, U>::sort()
{
    sort_layer(sample->root.get());
    std::sort(sample->last_layer.begin(), sample->last_layer.end(),
              [](const std::shared_ptr<trie_based::NodeCount<T>> &l, const std::shared_ptr<trie_based::NodeCount<T>> &r)
    {
        return l->index < r->index;
    });
}


template <typename T, typename U>
void ImplicitQuantileSorted<T,U>::sort_layer(trie_based::NodeCount<T> *p)
{
    bool must_sort = !std::is_sorted(p->children.begin(), p->children.end(),
                                     [](const std::shared_ptr<trie_based::NodeCount<T>> &l,
                                        const std::shared_ptr<trie_based::NodeCount<T>> &r)
    {
        return l->index < r->index;
    });
    if(must_sort)
    {
        std::sort(p->children.begin(), p->children.end(),
                  [](const std::shared_ptr<trie_based::NodeCount<T>> &l, const std::shared_ptr<trie_based::NodeCount<T>> &r)
        {
            return l->index < r->index;
        });

    }
    if(p->children != sample->last_layer) // bad comparison here
    {
        for(auto &i : p->children)
        {
            sort_layer(i.get());
        }
    }
}


template <typename T, typename U>
void ImplicitQuantileSorted<T, U>::transform(const std::vector<U>& in01, std::vector<U>& out) const
{
    auto *p = sample->root.get();
    for(size_t i = 0; i != in01.size(); ++i)
    {
        std::vector<size_t> psum(p->children.size() + 1, 0);
        for(size_t j = 1, k = 0; j != p->children.size(); ++j)
        {
            k += p->children[j-1]->count;
            psum[j] = k;
        }
        psum[p->children.size()] = p->count;

        auto rez = quantile_transform(p, psum, i, in01[i]);
        out[i] = rez.second;
        p = p->children[rez.first].get();
    }
}

template <typename T, typename U>
size_t ImplicitQuantileSorted<T, U>::count_less_binary(trie_based::NodeCount<T> *layer, T target) const
{
    auto lb = std::lower_bound(layer->children.begin(), layer->children.end(), target,
                               [](const std::shared_ptr<trie_based::NodeCount<T>> &l,
                                  const T &r)
    {
        return l->index < r;
    });
    size_t pos = std::distance(layer->children.begin(), lb);
    if(lb == layer->children.end())
        pos = layer->children.size(); // to psum! which is layer->children.size() + 1
    return pos; // to psum!
}

template <typename T, typename U>
std::pair<size_t, U> ImplicitQuantileSorted<T, U>::quantile_transform(trie_based::NodeCount<T> *layer, const std::vector<size_t> &psum, size_t ind, U val01) const
{
    size_t m = 0, count = grids[ind].size() - 1, step, c1 = 0, c2 = 0;
    U f1 = 0.0, f2 = 0.0, sample_size_u = static_cast<U>(layer->count);
    auto first = grids[ind].begin();
    auto it = grids[ind].begin();

    while(count > 0)
    {
        it = first;
        step = count / 2;
        std::advance(it, step);
        m = std::distance(grids[ind].begin(), it);

        c1 = psum[count_less_binary(layer, m)];
        f1 = c1/sample_size_u;

        if(f1 < val01)
        {
            c2 = psum[count_less_binary(layer, m + 1)];
            f2 = c2/sample_size_u;

            if(val01 < f2)
                break;

            first = ++it;
            count -= step + 1;
        }
        else
        {
            count = step;
        }
    }

    if(count == 0)
    {
        c2 = psum[count_less_binary(layer, m + 1)];
        f2 = c2/sample_size_u;
    }

    if(c1 == c2)
    {
        if(c1 == 0)
        {
            return std::make_pair(0, grids[ind][layer->children.front()->index] + 2.0*val01*dx[ind]);
        }
        if(c1 == layer->count)
        {
            return std::make_pair(layer->children.size() - 1, grids[ind][layer->children.back()->index] + 2.0*val01*dx[ind]);
        }

        T target = m;
        auto pos = std::lower_bound(layer->children.begin(), layer->children.end(), target,
                                    [](const std::shared_ptr<trie_based::NodeCount<T>> &l,
                                       const T &r)
        {
            return l->index < r;
        });
        size_t index = std::distance(layer->children.begin(), pos);

        if(index > 0)
        {
            int curr1 = std::abs(static_cast<int>(layer->children[index]->index) - static_cast<int>(m));
            int curr2 = std::abs(static_cast<int>(layer->children[index - 1]->index) - static_cast<int>(m));

            if(curr1 < curr2)
            {
                return std::make_pair(index, grids[ind][layer->children[index]->index] + 2.0*val01*dx[ind]);
            }
            else if(curr1 == curr2)
            {
                if(layer->children[index - 1]->index < layer->children[index]->index)
                    return std::make_pair(index - 1, grids[ind][layer->children[index - 1]->index] + 2.0*val01*dx[ind]);
                else
                    return std::make_pair(index, grids[ind][layer->children[index]->index] + 2.0*val01*dx[ind]);
            }
            else
            {
                return std::make_pair(index - 1, grids[ind][layer->children[index - 1]->index] + 2.0*val01*dx[ind]);
            }

        }
        return std::make_pair(index, grids[ind][layer->children[index]->index] + 2.0*val01*dx[ind]);
    }
    T target = m;
    auto pos = std::lower_bound(layer->children.begin(), layer->children.end(), target,
                                [](const std::shared_ptr<trie_based::NodeCount<T>> &l,
                                   const T &r)
    {
        return l->index < r;
    });
    T index = std::distance(layer->children.begin(), pos);
    return std::make_pair(index, grids[ind][m] + (val01 - f1) * (grids[ind][m + 1] - grids[ind][m]) / (f2 - f1));
}

}

#endif
