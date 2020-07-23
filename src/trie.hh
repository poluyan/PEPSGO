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

#ifndef TRIE_HH
#define TRIE_HH

#include <memory>
#include <vector>
#include <algorithm>

namespace mveqf
{
	namespace trie_based
	{
		template <typename TNode, typename TIndex>
		class Trie
		{
		protected:
			size_t dimension;
		public:
			std::shared_ptr<TNode> root;
			Trie();
			Trie(size_t dim);
			~Trie();
			void set_dimension(size_t dim);
			size_t get_dimension() const;
			bool empty() const;
			void insert(const std::vector<TIndex> &key, size_t number);
			bool search(const std::vector<TIndex> &key) const;
		};
		template <typename TNode, typename TIndex>
		Trie<TNode,TIndex>::Trie()
		{
			root = std::make_shared<TNode>();
		}
		template <typename TNode, typename TIndex>
		Trie<TNode,TIndex>::Trie(size_t dim) : dimension(dim)
		{
			root = std::make_shared<TNode>();
		}
		template <typename TNode, typename TIndex>
		Trie<TNode,TIndex>::~Trie() {}
		template <typename TNode, typename TIndex>
		void Trie<TNode,TIndex>::set_dimension(size_t dim)
		{
			dimension = dim;
		}
		template <typename TNode, typename TIndex>
		size_t Trie<TNode,TIndex>::get_dimension() const
		{
			return dimension;
		}
		template <typename TNode, typename TIndex>
		bool Trie<TNode,TIndex>::empty() const
		{
			return root->children.empty();
		}
		template <typename TNode, typename TIndex>
		void Trie<TNode,TIndex>::insert(const std::vector<TIndex> &key, size_t count)
		{
			auto p = root.get();
			for(const auto &i : key)
			{
				p->count += count;
				auto it = std::find_if(p->children.begin(), p->children.end(), [&i](const std::shared_ptr<TNode> &obj)
				{
					return obj->index == i;
				});
				if(it == p->children.end())
				{
					p->children.emplace_back(std::make_shared<TNode>(i));
					p->children.shrink_to_fit();
					p = p->children.back().get();
				}
				else
				{
					p = p->children[std::distance(p->children.begin(), it)].get();
				}
			}
			p->count += count;
		}
		template <typename TNode, typename TIndex>
		bool Trie<TNode,TIndex>::search(const std::vector<TIndex> &key) const
		{
			auto p = root.get();
			for(const auto &i : key)
			{
				auto it = std::find_if(p->children.begin(), p->children.end(), [&i](const std::shared_ptr<TNode> &obj)
				{
					return obj->index == i;
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
	}
}

#endif
