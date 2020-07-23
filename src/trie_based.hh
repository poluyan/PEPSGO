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
#include <map>
#include <unordered_map>
#include <algorithm>
#include <memory>

namespace mveqf
{
	namespace trie_based
	{
		template <typename TNode, typename TIndex>
		class TrieBased
		{
		protected:
			size_t dimension;
		public:
			std::shared_ptr<TNode> root;
			std::vector<std::shared_ptr<TNode>> last_layer;
			TrieBased();
			TrieBased(size_t dim);
			~TrieBased();
			void set_dimension(size_t dim);
			size_t get_dimension() const;
			void insert(const std::vector<TIndex> &key);
			void insert(const std::vector<TIndex> &key, size_t number);
			bool search(const std::vector<TIndex> &key) const;
			void fill_tree_count();
			bool empty() const;
			void remove_tree();
			size_t get_total_count() const;
			std::vector<TIndex> get_and_remove_last();
			size_t get_link_count() const;
			size_t get_node_count() const;
			std::map<size_t,std::vector<TIndex>> get_layer_count() const;
		protected:
			void layer_count(size_t current_layer, TNode *p, std::map<size_t,std::vector<TIndex>> &layers) const;
			void get_link_number(TNode *p, size_t &count) const;
			void get_node_number(TNode *p, size_t &count) const;
			void fill_tree_count(TNode *p);
			void get_number(TNode *p, size_t &count) const;
			void is_all_empty(TNode *p) const;
			void delete_last(int dim);
		};
		template <typename TNode, typename TIndex>
		TrieBased<TNode,TIndex>::TrieBased()
		{
			root = std::make_shared<TNode>();
		}
		template <typename TNode, typename TIndex>
		TrieBased<TNode,TIndex>::TrieBased(size_t dim) : dimension(dim)
		{
			root = std::make_shared<TNode>();
		}
		template <typename TNode, typename TIndex>
		TrieBased<TNode,TIndex>::~TrieBased() {}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::set_dimension(size_t dim)
		{
			dimension = dim;
		}
		template <typename TNode, typename TIndex>
		size_t TrieBased<TNode,TIndex>::get_dimension() const
		{
			return dimension;
		}
		template <typename TNode, typename TIndex>
		bool TrieBased<TNode,TIndex>::empty() const
		{
			return root->children.empty();
		}
		template <typename TNode, typename TIndex>
		size_t TrieBased<TNode,TIndex>::get_total_count() const
		{
			size_t count = 0;
			for(auto &i : last_layer)
			{
				count += i.use_count();
			}
			return count - last_layer.size();
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::get_node_number(TNode *p, size_t &count) const
		{
			for(const auto &i : p->children)
			{
				++count;
				get_node_number(i.get(), count);
			}
		}
		template <typename TNode, typename TIndex>
		size_t TrieBased<TNode,TIndex>::get_node_count() const
		{
			size_t count = 0;
			get_node_number(root.get(), count);
			for(const auto &i : last_layer)
				count -= i.use_count();
			return count + 2*last_layer.size() + 1;
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::get_link_number(TNode *p, size_t &count) const
		{
			for(const auto &i : p->children)
			{
				++count;
				get_node_number(i.get(), count);
			}
		}
		template <typename TNode, typename TIndex>
		size_t TrieBased<TNode,TIndex>::get_link_count() const
		{
			size_t count = 0;
			get_node_number(root.get(), count);
//    for(const auto &i : last_layer)
//        count -= i.use_count();
			return count;
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::insert(const std::vector<TIndex> &key)
		{
			auto p = root.get();
			for(size_t i = 0; i != key.size() - 1; i++)
			{
				auto value = key[i];
				auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<TNode> &obj)
				{
					return obj->index == value;
				});
				if(it == p->children.end())
				{
					p->children.emplace_back(std::make_shared<TNode>(value));
					p->children.shrink_to_fit();
					p = p->children.back().get();
				}
				else
				{
					p = p->children[std::distance(p->children.begin(), it)].get();
				}
			}
			auto value = key.back();
			auto it = std::find_if(last_layer.begin(), last_layer.end(), [&value](const std::shared_ptr<TNode> &obj)
			{
				return obj->index == value;
			});
			size_t dist = 0;
			if(it == last_layer.end())
			{
				last_layer.emplace_back(std::make_shared<TNode>(value));
				last_layer.shrink_to_fit();
				dist = last_layer.size() - 1;
			}
			else
			{
				dist = std::distance(last_layer.begin(), it);
			}

			it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<TNode> &obj)
			{
				return obj->index == value;
			});
			if(it == p->children.end())
			{
				std::shared_ptr<TNode> ptr(last_layer[dist]);
				p->children.emplace_back(ptr);
				p->children.shrink_to_fit();
			}
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::insert(const std::vector<TIndex> &key, size_t count)
		{
			auto p = root.get();
			for(size_t i = 0; i != key.size() - 1; i++)
			{
				p->count += count;
				auto value = key[i];
				auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<TNode> &obj)
				{
					return obj->index == value;
				});
				if(it == p->children.end())
				{
					p->children.emplace_back(std::make_shared<TNode>(value));
					p->children.shrink_to_fit();
					p = p->children.back().get();
				}
				else
				{
					p = p->children[std::distance(p->children.begin(), it)].get();
				}
			}
			auto value = key.back();
			auto it = std::find_if(last_layer.begin(), last_layer.end(), [&value](const std::shared_ptr<TNode> &obj)
			{
				return obj->index == value;
			});
			size_t dist = 0;
			if(it == last_layer.end())
			{
				last_layer.emplace_back(std::make_shared<TNode>(value));
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
			it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<TNode> &obj)
			{
				return obj->index == value;
			});
			if(it == p->children.end())
			{
				std::shared_ptr<TNode> ptr(last_layer[dist]);
				p->children.emplace_back(ptr);
				p->children.shrink_to_fit();
			}
		}
		template <typename TNode, typename TIndex>
		bool TrieBased<TNode,TIndex>::search(const std::vector<TIndex> &key) const
		{
			auto p = root.get();
			for(size_t i = 0; i != key.size(); i++)
			{
				auto value = key[i];
				auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<TNode> &obj)
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
		template <typename TNode, typename TIndex>
		std::vector<TIndex> TrieBased<TNode,TIndex>::get_and_remove_last()
		{
			std::vector<TIndex> sample;

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
			                 [&sample](const std::shared_ptr<TNode> &obj)
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
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::remove_tree()
		{
			auto p = root.get();
			while(!p->children.empty())
			{
				auto t = get_and_remove_last();
			}
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::delete_last(int dim)
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
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::fill_tree_count()
		{
			fill_tree_count(root.get());
			size_t count = 0;
			for(auto &i : root->children)
			{
				count += i->count;
			}
			root->count = count;
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::fill_tree_count(TNode *p)
		{
			for(auto &i : p->children)
			{
				size_t count = 0;
				get_number(i.get(), count);
				i->count = count > 0 ? count : 1;
				fill_tree_count(i.get());
			}
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::get_number(TNode *p, size_t &count) const
		{
			for(auto &i : p->children)
			{
				if(i->children.empty())
					++count;
				get_number(i.get(), count);
			}
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::layer_count(size_t current_layer, TNode *p, std::map<size_t, std::vector<TIndex>> &layers) const
		{
			auto it = layers.find(current_layer);
			if(it != layers.end())
			{
				it->second.push_back(p->children.size());
			}
			else
			{
				layers.insert(std::pair<size_t, std::vector<TIndex>>(current_layer, std::vector<TIndex>(1, p->children.size())));
			}
			for(const auto &i : p->children)
			{
				layer_count(current_layer + 1, i.get(), layers);
			}
		}
		template <typename TNode, typename TIndex>
		std::map<size_t, std::vector<TIndex>> TrieBased<TNode,TIndex>::get_layer_count() const
		{
			std::map<size_t, std::vector<TIndex>> layers;
			size_t cur_layer = 0;
			layer_count(cur_layer, root.get(), layers);
			layers.erase(std::prev(layers.end()));
			return layers;
		}

		template <typename I>
		class TrieLayer
		{
		public:
			size_t dimension;
			std::vector< std::multimap<I,std::vector<I>> > layers;
			typename std::multimap<I,std::vector<I>>::iterator root;
			bool sorted;
		public:
			TrieLayer();
			TrieLayer(size_t dim);
			void set_dimension(size_t dim);
			size_t get_dimension() const;
			bool empty() const;
			void insert(const std::vector<I> &key);
			bool search(const std::vector<I> &key) const;
			void sort();
//    void print() const;
		};

		template <typename I>
		TrieLayer<I>::TrieLayer() {}

		template <typename I>
		TrieLayer<I>::TrieLayer(size_t dim): dimension(dim)
		{
			layers.resize(dimension);
			layers.front().insert(std::pair<I, std::vector<I>>(0, std::vector<I>()));
			root = layers.front().find(0);
		}
		template <typename I>
		void TrieLayer<I>::set_dimension(size_t dim)
		{
			dimension = dim;
			layers.resize(dimension);
			layers.front().insert(std::pair<I, std::vector<I>>(0, std::vector<I>()));
			root = layers.front().find(0);
		}
		template <typename I>
		size_t TrieLayer<I>::get_dimension() const
		{
			return dimension;
		}
		template <typename I>
		bool TrieLayer<I>::empty() const
		{
			return root->second.empty();
		}
		template <typename I>
		void TrieLayer<I>::insert(const std::vector<I> &key)
		{
			sorted = false;

			auto pos = std::find(root->second.begin(), root->second.end(), key.front());
			if(pos == root->second.end())
				root->second.push_back(key.front());

			for(size_t i = 0; i != dimension - 1; i++)
			{
				auto it = layers[i+1].find(key[i]);
				if(it == layers[i+1].end())
				{
					layers[i+1].insert(std::pair<I, std::vector<I>>(key[i],std::vector<I>(1, key[i+1])));
				}
				else
				{
					auto pos = std::find(it->second.begin(), it->second.end(), key[i + 1]);
					if(pos == it->second.end())
						it->second.push_back(key[i + 1]);
				}
			}
		}
		template <typename I>
		void TrieLayer<I>::sort()
		{
			for(size_t i = 0; i != layers.size(); i++)
			{
				for(auto it = layers[i].begin(); it != layers[i].end(); ++it)
				{
					std::sort(it->second.begin(), it->second.end());
				}
			}
			sorted = true;
		}
//template <typename I>
//void TrieLayer<I>::print() const
//{
//    for(size_t i = 0; i != layers.size(); i++)
//    {
//        for(auto it = layers[i].begin(); it != layers[i].end(); ++it)
//        {
//            std::cout << int(it->first) << " => ";
//            for(size_t j = 0; j != it->second.size(); j++)
//            {
//                std::cout << int(it->second[j]) << ' ';
//            }
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//    }
//}
		template <typename I>
		bool TrieLayer<I>::search(const std::vector<I> &key) const
		{
			if(sorted)
			{
				auto pos = std::lower_bound(root->second.begin(), root->second.end(), key.front());
				if(pos == root->second.end())
					return false;
				if(*pos != key.front())
					return false;

				for(size_t i = 0; i != dimension - 1; i++)
				{
					auto it = layers[i + 1].find(key[i]);
					if(it == layers[i + 1].end())
					{
						return false;
					}
					else
					{
						auto pos = std::lower_bound(it->second.begin(), it->second.end(), key[i + 1]);
						if(pos == it->second.end())
							return false;
						if(*pos != key[i + 1])
							return false;
					}
				}
				return true;
			}
			else
			{
				auto pos = std::find(root->second.begin(), root->second.end(), key.front());
				if(pos == root->second.end())
					return false;
				for(size_t i = 0; i != dimension - 1; i++)
				{
					auto it = layers[i + 1].find(key[i]);
					if(it == layers[i + 1].end())
					{
						return false;
					}
					else
					{
						auto pos = std::find(it->second.begin(), it->second.end(), key[i + 1]);
						if(pos == it->second.end())
							return false;
					}
				}
				return true;
			}
		}

		struct invect
		{
			char vname;
			int index;
			size_t count;
			invect(char _vname, int _index, int _count)
			{
				vname = _vname;
				index = _index;
				count = _count;
			}
		};

		template <typename I>
		class Graph
		{
		public:
			std::unordered_map< char, std::vector<invect> > layers;
			typename std::unordered_map<char, std::vector<invect>>::iterator root;
		public:
			Graph();
			void print() const;
		};

		template <typename I>
		Graph<I>::Graph()
		{
			std::vector<invect> y, a,b,c,d,e, f,g,h,j,k,l, z;
			y.push_back(invect('a',0,4));
			y.push_back(invect('b',1,1));
			y.push_back(invect('e',2,4));
			y.push_back(invect('d',3,2));
			y.push_back(invect('c',4,3));
			layers.insert(std::make_pair('y',y));

			a.push_back(invect('g',0,1));
			a.push_back(invect('h',2,1));
			a.push_back(invect('f',3,2));
			layers.insert(std::make_pair('a',a));

			b.push_back(invect('h',0,1));
			layers.insert(std::make_pair('b',b));

			c.push_back(invect('h',0,1));
			c.push_back(invect('h',3,1));
			c.push_back(invect('h',4,1));
			layers.insert(std::make_pair('c',c));

			d.push_back(invect('j',0,1));
			d.push_back(invect('h',3,1));
			layers.insert(std::make_pair('d',d));

			e.push_back(invect('k',0,2));
			e.push_back(invect('l',1,1));
			e.push_back(invect('l',2,1));
			layers.insert(std::make_pair('e',e));


			f.push_back(invect('z',2,1));
			f.push_back(invect('z',3,1));
			layers.insert(std::make_pair('f',f));

			g.push_back(invect('z',1,1));
			layers.insert(std::make_pair('g',g));

			h.push_back(invect('z',0,1));
			layers.insert(std::make_pair('h',h));

			j.push_back(invect('z',2,1));
			layers.insert(std::make_pair('j',j));

			k.push_back(invect('z',0,1));
			k.push_back(invect('z',4,1));
			layers.insert(std::make_pair('k',k));

			l.push_back(invect('z',4,1));
			layers.insert(std::make_pair('l',l));

			layers.insert(std::make_pair('z',z));

			root = layers.find('y');

			/*std::vector<invect> y, a,b,c,d,e, f,g,h,j,k,l, z;
			y.push_back(invect('a',0,1));
			y.push_back(invect('a',1,1));
			y.push_back(invect('a',2,1));
			y.push_back(invect('a',3,1));
			y.push_back(invect('a',4,1));
			layers.insert(std::make_pair('y',y));

			a.push_back(invect('b',0,1));
			a.push_back(invect('b',1,1));
			a.push_back(invect('b',2,1));
			a.push_back(invect('b',3,1));
			a.push_back(invect('b',4,1));
			layers.insert(std::make_pair('a',a));

			b.push_back(invect('z',0,1));
			b.push_back(invect('z',1,1));
			b.push_back(invect('z',2,1));
			b.push_back(invect('z',3,1));
			b.push_back(invect('z',4,1));
			layers.insert(std::make_pair('b',b));

			layers.insert(std::make_pair('z',z));

			root = layers.find('y');*/
		}


//template <typename I>
//void Graph<I>::print() const
//{
//    for(const auto & i : layers)
//    {
//
//        std::cout << i.first << " => ";
//        for(auto it = i.second.begin(); it != i.second.end(); ++it)
//        {
//            std::cout << it->vname << ' ' << it->index << ' ' << it->count << '\t';
//        }
//        std::cout << std::endl;
//    }
//}

	}

}
#endif
