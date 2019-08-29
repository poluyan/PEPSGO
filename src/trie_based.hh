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

#ifndef TRIE_BASED_HH
#define TRIE_BASED_HH

#include <vector>
#include <algorithm>
#include <memory>

namespace trie_based
{

template <template <typename> class T, typename I>
struct TrieNode
{
    I index;
    std::vector<std::shared_ptr<T<I>>> children;
    TrieNode() : index(0) { }
    TrieNode(I ind) : index(ind) { }
};

template <typename I>
struct Node: public TrieNode<Node, I>
{
    Node() : TrieNode<Node, I>() {}
    Node(I ind) : TrieNode<Node, I>(ind) {}
};

template <typename I>
struct NodeCount: public TrieNode<NodeCount, I>
{
    size_t count;
    NodeCount() : TrieNode<NodeCount, I>(), count(0) {}
    NodeCount(I ind) : TrieNode<NodeCount, I>(ind), count(0) {}
};

template <typename T, typename I>
class TrieBased
{
protected:
    size_t dimension;
public:
    std::shared_ptr<T> root;
    std::vector<std::shared_ptr<T>> last_layer;
    TrieBased();
    TrieBased(size_t dim);
    ~TrieBased();
    void set_dimension(size_t dim);
    size_t get_dimension() const;
    void insert(const std::vector<I> &key);
    void insert(const std::vector<I> &key, size_t number);
    bool search(const std::vector<I> &key) const;
    void fill_tree_count();
    bool empty() const;
    void remove_tree();
    size_t get_total_count() const;
    std::vector<I> get_and_remove_last();
    size_t get_link_count() const;
    size_t get_node_count() const;
protected:
    void get_link_number(T *p, size_t &count) const;
    void get_node_number(T *p, size_t &count) const;
    void fill_tree_count(T *p);
    void get_number(T *p, size_t &count) const;
    void is_all_empty(T *p) const;
    void delete_last(int dim);
};
template <typename T, typename I>
TrieBased<T,I>::TrieBased()
{
    root = std::make_shared<T>();
}
template <typename T, typename I>
TrieBased<T,I>::TrieBased(size_t dim) : dimension(dim)
{
    root = std::make_shared<T>();
}
template <typename T, typename I>
TrieBased<T,I>::~TrieBased() {}
template <typename T, typename I>
void TrieBased<T,I>::set_dimension(size_t dim)
{
    dimension = dim;
}
template <typename T, typename I>
size_t TrieBased<T,I>::get_dimension() const
{
    return dimension;
}
template <typename T, typename I>
bool TrieBased<T,I>::empty() const
{
    return root->children.empty();
}
template <typename T, typename I>
size_t TrieBased<T,I>::get_total_count() const
{
    size_t count = 0;
    for(auto &i : last_layer)
    {
        count += i.use_count();
    }
    return count - last_layer.size();
}
template <typename T, typename I>
void TrieBased<T,I>::get_node_number(T *p, size_t &count) const
{
    for(const auto &i : p->children)
    {
        ++count;
        get_node_number(i.get(), count);
    }
}
template <typename T, typename I>
size_t TrieBased<T,I>::get_node_count() const
{
    size_t count = 0;
    get_node_number(root.get(), count);
    for(const auto &i : last_layer)
        count -= i.use_count();
    return count + 2*last_layer.size() + 1;
}
template <typename T, typename I>
void TrieBased<T,I>::get_link_number(T *p, size_t &count) const
{
    for(const auto &i : p->children)
    {
        ++count;
        get_node_number(i.get(), count);
    }
}
template <typename T, typename I>
size_t TrieBased<T,I>::get_link_count() const
{
    size_t count = 0;
    get_node_number(root.get(), count);
//    for(const auto &i : last_layer)
//        count -= i.use_count();
    return count;
}
template <typename T, typename I>
void TrieBased<T,I>::insert(const std::vector<I> &key)
{
    auto p = root.get();
    for(size_t i = 0; i != key.size() - 1; i++)
    {
        auto value = key[i];
        auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<T> &obj)
        {
            return obj->index == value;
        });
        if(it == p->children.end())
        {
            p->children.emplace_back(std::make_shared<T>(value));
            p->children.shrink_to_fit();
            p = p->children.back().get();
        }
        else
        {
            p = p->children[std::distance(p->children.begin(), it)].get();
        }
    }
    auto value = key.back();
    auto it = std::find_if(last_layer.begin(), last_layer.end(), [&value](const std::shared_ptr<T> &obj)
    {
        return obj->index == value;
    });
    size_t dist = 0;
    if(it == last_layer.end())
    {
        last_layer.emplace_back(std::make_shared<T>(value));
        last_layer.shrink_to_fit();
        dist = last_layer.size() - 1;
    }
    else
    {
        dist = std::distance(last_layer.begin(), it);
    }

    it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<T> &obj)
    {
        return obj->index == value;
    });
    if(it == p->children.end())
    {
        std::shared_ptr<T> ptr(last_layer[dist]);
        p->children.emplace_back(ptr);
        p->children.shrink_to_fit();
    }
}
template <typename T, typename I>
void TrieBased<T,I>::insert(const std::vector<I> &key, size_t count)
{
    auto p = root.get();
    for(size_t i = 0; i != key.size() - 1; i++)
    {
        p->count += count;
        auto value = key[i];
        auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<T> &obj)
        {
            return obj->index == value;
        });
        if(it == p->children.end())
        {
            p->children.emplace_back(std::make_shared<T>(value));
            p->children.shrink_to_fit();
            p = p->children.back().get();
        }
        else
        {
            p = p->children[std::distance(p->children.begin(), it)].get();
        }
    }
    auto value = key.back();
    auto it = std::find_if(last_layer.begin(), last_layer.end(), [&value](const std::shared_ptr<T> &obj)
    {
        return obj->index == value;
    });
    size_t dist = 0;
    if(it == last_layer.end())
    {
        last_layer.emplace_back(std::make_shared<T>(value));
//        last_layer.back()->count += count;
        last_layer.back()->count = 1;
        last_layer.shrink_to_fit();
        dist = last_layer.size() - 1;
    }
    else
    {
        dist = std::distance(last_layer.begin(), it);
        last_layer[dist]->count = 1;
//        last_layer[dist]->count += count;
    }
    p->count += count;
    it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<T> &obj)
    {
        return obj->index == value;
    });
    if(it == p->children.end())
    {
        std::shared_ptr<T> ptr(last_layer[dist]);
        p->children.emplace_back(ptr);
        p->children.shrink_to_fit();
    }
}
template <typename T, typename I>
bool TrieBased<T,I>::search(const std::vector<I> &key) const
{
    auto p = root.get();
    for(size_t i = 0; i != key.size(); i++)
    {
        auto value = key[i];
        auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<T> &obj)
        {
            return obj->index == value;
        });
        if(it == p->children.end())
        {
            return false;
        }
        else
        {
            p = p->children[std::distance(p->children.begin(), it)].get();
        }
    }
    return true;
}
template <typename T, typename I>
std::vector<I> TrieBased<T,I>::get_and_remove_last()
{
    std::vector<I> sample;

    auto p = root.get();
    if(p->children.empty())
        return sample;

    for(size_t i = 0; i != dimension; ++i)
    {
        p->count--;
        sample.push_back(p->children.back()->index);
        p = p->children.back().get();
    }

    size_t dim = sample.size() - 1;

    delete_last(dim);

    last_layer.erase(
        std::remove_if(last_layer.begin(), last_layer.end(),
                       [&sample](const std::shared_ptr<T> &obj)
    {
        if(obj->index != sample.back())
            return false;
        if(obj.use_count() == 1)
            return true;
        else
            return false;
    }),
    last_layer.end());

    return sample;
}
template <typename T, typename I>
void TrieBased<T,I>::remove_tree()
{
    auto p = root.get();
    while(!p->children.empty())
    {
        auto t = get_and_remove_last();
    }
}
template <typename T, typename I>
void TrieBased<T,I>::delete_last(int dim)
{
    if(dim < 0)
        return;

    auto p = root.get();
    if(p->children.empty())
        return;

    for(int i = 0; i != dim; ++i)
    {
        p = p->children.back().get();
    }
    p->children.pop_back();

    if(p->children.empty())
    {
        dim = dim - 1;
        delete_last(dim);
    }
}
template <typename T, typename I>
void TrieBased<T,I>::fill_tree_count()
{
    fill_tree_count(root.get());
    size_t count = 0;
    for(auto &i : root->children)
    {
        count += i->count;
    }
    root->count = count;
}
template <typename T, typename I>
void TrieBased<T,I>::fill_tree_count(T *p)
{
    for(auto &i : p->children)
    {
        size_t count = 0;
        get_number(i.get(), count);
        i->count = count > 0 ? count : 1;
        fill_tree_count(i.get());
    }
}
template <typename T, typename I>
void TrieBased<T,I>::get_number(T *p, size_t &count) const
{
    for(auto &i : p->children)
    {
        if(i->children.empty())
            ++count;
        get_number(i.get(), count);
    }
}

}

#endif
