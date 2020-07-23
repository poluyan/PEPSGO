/**************************************************************************

   Copyright © 2020 Sergey Poluyan <svpoluyan@gmail.com>

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
#ifndef IMPLICIT_HH
#define IMPLICIT_HH

#include <trie_node.hh>
#include <quantile.hh>

namespace mveqf
{
	template <typename TIndex, typename TFloat>
	class ImplicitQuantile : public Quantile<TIndex, TFloat>
	{
	protected:
		typedef trie_based::TrieBased<trie_based::NodeCount<TIndex>,TIndex> sample_type;
		std::shared_ptr<sample_type> sample;

		//using Quantile<TIndex, TFloat>::grids;
		using Quantile<TIndex, TFloat>::grid_number;
		using Quantile<TIndex, TFloat>::dx;
		using Quantile<TIndex, TFloat>::lb;
		using Quantile<TIndex, TFloat>::ub;

		using Quantile<TIndex, TFloat>::get_grid_value;

		std::pair<size_t, size_t> count_less(trie_based::NodeCount<TIndex> *layer, const size_t &r) const;
		std::pair<size_t, TFloat> quantile_transform(trie_based::NodeCount<TIndex> *layer, size_t ind, TFloat val01) const;
	public:
		ImplicitQuantile() = default;
		ImplicitQuantile(std::vector<TFloat> in_lb, std::vector<TFloat> in_ub, std::vector<size_t> in_gridn);
		ImplicitQuantile(const ImplicitQuantile&) = delete;
		ImplicitQuantile& operator=(const ImplicitQuantile&) = delete;
		void set_sample(const std::vector<std::vector<TIndex>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights) override;
		void set_sample_and_fill_count(const std::vector<std::vector<TIndex>> &in_sample);
		void set_sample_shared_and_fill_count(std::shared_ptr<sample_type> in_sample);
		void set_sample_shared(std::shared_ptr<sample_type> in_sample);
		void transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const override;
		void transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const override;
		size_t get_node_count() const;
		size_t get_link_count() const;
		using Quantile<TIndex, TFloat>::get_the_closest_grid_node_to_the_value;
		using Quantile<TIndex, TFloat>::get_real_node_values;
		~ImplicitQuantile();
	};

	template <typename TIndex, typename TFloat>
	ImplicitQuantile<TIndex, TFloat>::ImplicitQuantile(std::vector<TFloat> in_lb,
	    std::vector<TFloat> in_ub,
	    std::vector<size_t> in_gridn) : Quantile<TIndex, TFloat>(in_lb, in_ub, in_gridn)
	{}

	template <typename TIndex, typename TFloat>
	ImplicitQuantile<TIndex, TFloat>::~ImplicitQuantile()
	{}

	template <typename TIndex, typename TFloat>
	size_t ImplicitQuantile<TIndex, TFloat>::get_node_count() const
	{
		return sample->get_node_count();
	}

	template <typename TIndex, typename TFloat>
	size_t ImplicitQuantile<TIndex, TFloat>::get_link_count() const
	{
		return sample->get_link_count();
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantile<TIndex, TFloat>::set_sample_and_fill_count(const std::vector<std::vector<TIndex>> &in_sample)
	{
		sample = std::make_shared<sample_type>();
		sample->set_dimension(grid_number.size());
		for(const auto & i : in_sample)
			sample->insert(i);
		sample->fill_tree_count();
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TIndex>> &in_sample)
	{
		set_sample_and_fill_count(in_sample);
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample)
	{
		sample = std::make_shared<sample_type>();
		sample->set_dimension(grid_number.size());
		for(size_t i = 0; i != in_sample.size(); ++i)
		{
			std::vector<TIndex> temp(in_sample[i].size());
			for(size_t j = 0; j != in_sample[i].size(); ++j)
			{
				temp[j] = get_the_closest_grid_node_to_the_value(lb[j], ub[j], grid_number[j], in_sample[i][j]);
			}
//			for(size_t j = 0; j != temp.size(); j++)
//			{
//				std::cout << temp[j] << ' ';
//			}
//			std::cout << std::endl;
			sample->insert(temp);
		}
		sample->fill_tree_count();
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights)
	{
		set_sample(in_sample);
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantile<TIndex, TFloat>::set_sample_shared_and_fill_count(std::shared_ptr<sample_type> in_sample)
	{
		sample = std::move(in_sample);
		sample->fill_tree_count();
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantile<TIndex, TFloat>::set_sample_shared(std::shared_ptr<sample_type> in_sample)
	{
		sample = std::move(in_sample);
	}

	template <typename TIndex, typename TFloat>
	std::pair<size_t, size_t> ImplicitQuantile<TIndex, TFloat>::count_less(trie_based::NodeCount<TIndex> *layer, const size_t &r) const
	{
		std::pair<size_t, size_t> res;
		for(const auto &i : layer->children)
		{
			size_t j = static_cast<size_t>(i->index);
			if(j < r + 1)
			{
				res.second += i->count;
				if(j < r)
				{
					res.first += i->count;
				}
			}
		}
		return res;
	}
	template <typename TIndex, typename TFloat>
	void ImplicitQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const
	{
		auto p = sample->root.get();
		for(size_t i = 0, k; i != in01.size(); ++i)
		{
			std::tie(k, out[i]) = quantile_transform(p, i, in01[i]);
			p = p->children[k].get();
		}
	}
	template <typename TIndex, typename TFloat>
	void ImplicitQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const
	{
		auto p = sample->root.get();
		for(size_t i = 0; i != in01.size(); ++i)
		{
			//std::tie(k, out[i]) = quantile_transform(p, i, in01[i]);
			auto [k, result] = quantile_transform(p, i, in01[i]);
			out[i] = p->children[k]->index;
			p = p->children[k].get();
		}
	}
	template <typename TIndex, typename TFloat>
	std::pair<size_t, TFloat> ImplicitQuantile<TIndex, TFloat>::quantile_transform(trie_based::NodeCount<TIndex> *layer, size_t ind, TFloat val01) const
	{
		size_t m = 0, count = grid_number[ind], step, a = 0, b = 0;
		TFloat x = 0.0, y = 0.0, p = static_cast<TFloat>(layer->count);
		//auto first = grids[ind].begin();
		//auto it = grids[ind].begin();
		size_t it = 0, first = 0;

		while(count > 0)
		{
			it = first;
			step = count / 2;
			it += step;
			m = it;
			//std::advance(it, step);
			//m = std::distance(grids[ind].begin(), it);

			std::tie(a, b) = count_less(layer, m);
			x = a/p;

			if(x < val01)
			{
				y = b/p;
				if(val01 < y)
					break;

				first = ++it;
				count -= step + 1;
			}
			else
				count = step;
		}
		if(count == 0)
		{
			y = b/p;
		}
		if(a == b)
		{
			if(a == 0)
			{
				auto min_val_it = std::min_element(layer->children.begin(), layer->children.end(),
				                                   [](const std::shared_ptr<trie_based::NodeCount<TIndex>> &l,
				                                      const std::shared_ptr<trie_based::NodeCount<TIndex>> &r)
				{
					return l->index < r->index;
				});
				size_t min_ind = std::distance(layer->children.begin(), min_val_it);
				//return std::make_pair(min_ind, grids[ind][layer->children[min_ind]->index] + 2.0*val01*dx[ind]);
				return std::make_pair(min_ind, get_grid_value(ind, layer->children[min_ind]->index) + 2.0*val01*dx[ind]);
			}
			if(a == layer->count)
			{
				auto max_val_it = std::max_element(layer->children.begin(), layer->children.end(),
				                                   [](const std::shared_ptr<trie_based::NodeCount<TIndex>> &l,
				                                      const std::shared_ptr<trie_based::NodeCount<TIndex>> &r)
				{
					return l->index < r->index;
				});
				size_t max_ind = std::distance(layer->children.begin(), max_val_it);
				//return std::make_pair(max_ind, grids[ind][layer->children[max_ind]->index] + 2.0*val01*dx[ind]);
				return std::make_pair(max_ind, get_grid_value(ind, layer->children[max_ind]->index) + 2.0*val01*dx[ind]);
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
			//return std::make_pair(index, grids[ind][layer->children[index]->index] + 2.0*val01*dx[ind]);
			return std::make_pair(index, get_grid_value(ind, layer->children[index]->index) + 2.0*val01*dx[ind]);
		}
		size_t index = 0;
		TIndex target = m;
		for(size_t j = 1; j < layer->children.size(); j++)
		{
			if(layer->children[j]->index == target)
			{
				index = j;
				break;
			}
		}
		//return std::make_pair(index, grids[ind][m] + (val01 - x) * (grids[ind][m + 1] - grids[ind][m]) / (y - x));
		return std::make_pair(index, get_grid_value(ind, m) + (val01 - x) * (get_grid_value(ind, m + 1) - get_grid_value(ind, m)) / (y - x));
	}

	template <typename TIndex, typename TFloat>
	class ImplicitQuantileSorted : public ImplicitQuantile<TIndex, TFloat>
	{
	protected:
//		using ImplicitQuantile<TIndex, TFloat>::grids;
		using ImplicitQuantile<TIndex, TFloat>::grid_number;
		using ImplicitQuantile<TIndex, TFloat>::sample;
		using ImplicitQuantile<TIndex, TFloat>::dx;
		using ImplicitQuantile<TIndex, TFloat>::lb;
		using ImplicitQuantile<TIndex, TFloat>::ub;

		using ImplicitQuantile<TIndex, TFloat>::get_grid_value;

		using ImplicitQuantile<TIndex, TFloat>::get_the_closest_grid_node_to_the_value;

		using sample_type = typename ImplicitQuantile<TIndex, TFloat>::sample_type;
		void sort_layer(trie_based::NodeCount<TIndex> *p);

		size_t count_less_binary(trie_based::NodeCount<TIndex> *layer, TIndex target) const;
		std::pair<size_t, TFloat> quantile_transform(trie_based::NodeCount<TIndex> *layer, const std::vector<size_t> &psum, size_t ind, TFloat val01) const;
	public:
		ImplicitQuantileSorted() = default;
		ImplicitQuantileSorted(std::vector<TFloat> in_lb, std::vector<TFloat> in_ub, std::vector<size_t> in_gridn);
		ImplicitQuantileSorted(const ImplicitQuantileSorted&) = delete;
		ImplicitQuantileSorted& operator=(const ImplicitQuantileSorted&) = delete;
		void set_sample(const std::vector<std::vector<TIndex>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights) override;
		void set_sample_and_fill_count(const std::vector<std::vector<TIndex>> &in_sample);
		void set_sample_shared_and_fill_count(std::shared_ptr<sample_type> in_sample);
		void sort();
		void transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const override;
		void transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const override;
	};

	template <typename TIndex, typename TFloat>
	ImplicitQuantileSorted<TIndex, TFloat>::ImplicitQuantileSorted(std::vector<TFloat> in_lb,
	    std::vector<TFloat> in_ub,
	    std::vector<size_t> in_gridn) : ImplicitQuantile<TIndex, TFloat>(in_lb, in_ub, in_gridn)
	{
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantileSorted<TIndex, TFloat>::set_sample_and_fill_count(const std::vector<std::vector<TIndex>> &in_sample)
	{
		sample = std::make_shared<sample_type>();
		sample->set_dimension(grid_number.size());
		for(const auto & i : in_sample)
			sample->insert(i);
		sample->fill_tree_count();
		sort();
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantileSorted<TIndex, TFloat>::set_sample(const std::vector<std::vector<TIndex>> &in_sample)
	{
		set_sample_and_fill_count(in_sample);
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantileSorted<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample)
	{
		sample = std::make_shared<sample_type>();
		sample->set_dimension(grid_number.size());
		for(size_t i = 0; i != in_sample.size(); ++i)
		{
			std::vector<TIndex> temp(in_sample[i].size());
			for(size_t j = 0; j != in_sample[i].size(); ++j)
			{
				temp[j] = get_the_closest_grid_node_to_the_value(lb[j], ub[j], grid_number[j], in_sample[i][j]);
			}
			sample->insert(temp);
		}
		sample->fill_tree_count();
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantileSorted<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights)
	{
		set_sample(in_sample);
	}


	template <typename TIndex, typename TFloat>
	void ImplicitQuantileSorted<TIndex, TFloat>::set_sample_shared_and_fill_count(std::shared_ptr<sample_type> in_sample)
	{
		sample = std::move(in_sample);
		sample->fill_tree_count();
		sort();
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantileSorted<TIndex, TFloat>::sort()
	{
		sort_layer(sample->root.get());
		std::sort(sample->last_layer.begin(), sample->last_layer.end(),
		          [](const std::shared_ptr<trie_based::NodeCount<TIndex>> &l, const std::shared_ptr<trie_based::NodeCount<TIndex>> &r)
		{
			return l->index < r->index;
		});
	}


	template <typename TIndex, typename TFloat>
	void ImplicitQuantileSorted<TIndex,TFloat>::sort_layer(trie_based::NodeCount<TIndex> *p)
	{
		bool must_sort = !std::is_sorted(p->children.begin(), p->children.end(),
		                                 [](const std::shared_ptr<trie_based::NodeCount<TIndex>> &l,
		                                    const std::shared_ptr<trie_based::NodeCount<TIndex>> &r)
		{
			return l->index < r->index;
		});
		if(must_sort)
		{
			std::sort(p->children.begin(), p->children.end(),
			          [](const std::shared_ptr<trie_based::NodeCount<TIndex>> &l, const std::shared_ptr<trie_based::NodeCount<TIndex>> &r)
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


	template <typename TIndex, typename TFloat>
	void ImplicitQuantileSorted<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const
	{
		auto *p = sample->root.get();
		for(size_t i = 0, k; i != in01.size(); ++i)
		{
			std::vector<size_t> psum(p->children.size() + 1, 0);
			for(size_t j = 1, m = 0; j != p->children.size(); ++j)
			{
				m += p->children[j-1]->count;
				psum[j] = m;
			}
			psum[p->children.size()] = p->count;

			//auto [k, res] = quantile_transform(p, psum, i, in01[i]);
			//out[i] = res;
			std::tie(k, out[i]) = quantile_transform(p, psum, i, in01[i]);
			p = p->children[k].get();
		}
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantileSorted<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const
	{
		auto *p = sample->root.get();
		for(size_t i = 0; i != in01.size(); ++i)
		{
			std::vector<size_t> psum(p->children.size() + 1, 0);
			for(size_t j = 1, m = 0; j != p->children.size(); ++j)
			{
				m += p->children[j-1]->count;
				psum[j] = m;
			}
			psum[p->children.size()] = p->count;

			auto [k, res] = quantile_transform(p, psum, i, in01[i]);
			out[i] = p->children[k]->index;
//			std::tie(k, out[i]) = quantile_transform(p, psum, i, in01[i]);
			p = p->children[k].get();
		}
	}

	template <typename TIndex, typename TFloat>
	size_t ImplicitQuantileSorted<TIndex, TFloat>::count_less_binary(trie_based::NodeCount<TIndex> *layer, TIndex target) const
	{
		auto lb = std::lower_bound(layer->children.begin(), layer->children.end(), target,
		                           [](const std::shared_ptr<trie_based::NodeCount<TIndex>> &l,
		                              const TIndex &r)
		{
			return l->index < r;
		});
		size_t pos = std::distance(layer->children.begin(), lb);
		if(lb == layer->children.end())
			pos = layer->children.size(); // to psum! which is layer->children.size() + 1
		return pos; // to psum!
	}

	template <typename TIndex, typename TFloat>
	std::pair<size_t, TFloat> ImplicitQuantileSorted<TIndex, TFloat>::quantile_transform(trie_based::NodeCount<TIndex> *layer, const std::vector<size_t> &psum, size_t ind, TFloat val01) const
	{
		size_t m = 0, count = grid_number[ind], step, c1 = 0, c2 = 0;
		TFloat f1 = 0.0, f2 = 0.0, sample_size_u = static_cast<TFloat>(layer->count);
		//auto first = grids[ind].begin();
		//auto it = grids[ind].begin();
		size_t it = 0, first = 0;

		while(count > 0)
		{
			it = first;
			step = count / 2;
			it += step;
			m = it;

			//std::advance(it, step);
			//m = std::distance(grids[ind].begin(), it);

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
				return std::make_pair(0, get_grid_value(ind, layer->children.front()->index) + 2.0*val01*dx[ind]);
			}
			if(c1 == layer->count)
			{
				return std::make_pair(layer->children.size() - 1, get_grid_value(ind, layer->children.back()->index) + 2.0*val01*dx[ind]);
			}

			TIndex target = m;
			auto pos = std::lower_bound(layer->children.begin(), layer->children.end(), target,
			                            [](const std::shared_ptr<trie_based::NodeCount<TIndex>> &l,
			                               const TIndex &r)
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
					return std::make_pair(index, get_grid_value(ind, layer->children[index]->index) + 2.0*val01*dx[ind]);
				}
				else if(curr1 == curr2)
				{
					if(layer->children[index - 1]->index < layer->children[index]->index)
						return std::make_pair(index - 1, get_grid_value(ind, layer->children[index - 1]->index) + 2.0*val01*dx[ind]);
					else
						return std::make_pair(index, get_grid_value(ind, layer->children[index]->index) + 2.0*val01*dx[ind]);
				}
				else
				{
					return std::make_pair(index - 1, get_grid_value(ind, layer->children[index - 1]->index) + 2.0*val01*dx[ind]);
				}

			}
			return std::make_pair(index, get_grid_value(ind, layer->children[index]->index) + 2.0*val01*dx[ind]);
		}
		TIndex target = m;
		auto pos = std::lower_bound(layer->children.begin(), layer->children.end(), target,
		                            [](const std::shared_ptr<trie_based::NodeCount<TIndex>> &l,
		                               const TIndex &r)
		{
			return l->index < r;
		});
		TIndex index = std::distance(layer->children.begin(), pos);
		return std::make_pair(index, get_grid_value(ind, m) + (val01 - f1) * (get_grid_value(ind, m + 1) - get_grid_value(ind, m)) / (f2 - f1));
	}


	template <typename TIndex, typename TFloat>
	class ImplicitQuantileSortedInterp : public ImplicitQuantile<TIndex, TFloat>
	{
	protected:
//		using ImplicitQuantile<TIndex, TFloat>::grids;
		using ImplicitQuantile<TIndex, TFloat>::grid_number;
		using ImplicitQuantile<TIndex, TFloat>::sample;
		using ImplicitQuantile<TIndex, TFloat>::dx;

		using ImplicitQuantile<TIndex, TFloat>::get_grid_value;

		using sample_type = typename ImplicitQuantile<TIndex, TFloat>::sample_type;

		std::pair<size_t, size_t> count_less(trie_based::NodeCount<TIndex> *layer, const size_t &r) const;
		std::pair<size_t, TFloat> quantile_transform(trie_based::NodeCount<TIndex> *layer, size_t ind, TFloat val01) const;
	public:
		ImplicitQuantileSortedInterp() = default;
		ImplicitQuantileSortedInterp(std::vector<TFloat> in_lb, std::vector<TFloat> in_ub, std::vector<size_t> in_gridn);
		ImplicitQuantileSortedInterp(const ImplicitQuantileSortedInterp&) = delete;
		ImplicitQuantileSortedInterp& operator=(const ImplicitQuantileSortedInterp&) = delete;
		void transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const override;
		void transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const override;
	};

	template <typename TIndex, typename TFloat>
	ImplicitQuantileSortedInterp<TIndex, TFloat>::ImplicitQuantileSortedInterp(std::vector<TFloat> in_lb,
	    std::vector<TFloat> in_ub,
	    std::vector<size_t> in_gridn) : ImplicitQuantile<TIndex, TFloat>(in_lb, in_ub, in_gridn) { }

	template <typename TIndex, typename TFloat>
	std::pair<size_t, size_t> ImplicitQuantileSortedInterp<TIndex, TFloat>::count_less(trie_based::NodeCount<TIndex> *layer, const size_t& r) const
	{
		std::pair<size_t, size_t> res;
		for(size_t i = 0; i != layer->children.size(); ++i)
		{
			size_t j = static_cast<size_t>(layer->children[i]->index);
			if(j < r + 1)
			{
				res.second += layer->children[i]->count;
				if(j < r)
				{
					res.first += layer->children[i]->count;
				}
			}
		}
		return res;
	}
	template <typename TIndex, typename TFloat>
	void ImplicitQuantileSortedInterp<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const
	{
		auto p = sample->root.get();
		for(size_t i = 0, k; i != in01.size(); ++i)
		{
			std::tie(k, out[i]) = quantile_transform(p, i, in01[i]);
			p = p->children[k].get();
		}
	}
	template <typename TIndex, typename TFloat>
	void ImplicitQuantileSortedInterp<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const
	{
		auto p = sample->root.get();
		for(size_t i = 0; i != in01.size(); ++i)
		{
			//std::tie(k, out[i]) = quantile_transform(p, i, in01[i]);
			auto [k, result] = quantile_transform(p, i, in01[i]);
			out[i] = p->children[k]->index;
			p = p->children[k].get();
		}
	}
	template <typename TIndex, typename TFloat>
	std::pair<size_t, TFloat> ImplicitQuantileSortedInterp<TIndex, TFloat>::quantile_transform(trie_based::NodeCount<TIndex> *layer, size_t ind, TFloat val01) const
	{
		size_t m = 0, count = grid_number[ind], step, a = 0, b = 0;
		TFloat f1 = 0.0, f2 = 0.0, sample_size_u = static_cast<TFloat>(layer->count);
		//auto first = grids[ind].begin();
		//auto it = grids[ind].begin();
		size_t it = 0, first = 0;

		while(count > 0)
		{
			it = first;
			step = count / 2;
			it += step;
			m = it;

			//std::advance(it, step);
			//m = std::distance(grids[ind].begin(), it);

			std::tie(a, b) = count_less(layer, m);
			f1 = a/sample_size_u;

			if(f1 < val01)
			{
				f2 = b/sample_size_u;

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
			f2 = b/sample_size_u;
		}
		if(a == b)
		{
			if(a == 0)
			{
				auto min_val_it = std::min_element(layer->children.begin(), layer->children.end(),
				                                   [](const std::shared_ptr<trie_based::NodeCount<TIndex>> &l,
				                                      const std::shared_ptr<trie_based::NodeCount<TIndex>> &r)
				{
					return l->index < r->index;
				});
				size_t min_ind = std::distance(layer->children.begin(), min_val_it);
				return std::make_pair(min_ind, get_grid_value(ind, layer->children[min_ind]->index) + 2.0*val01*dx[ind]);
			}
			if(a == layer->count)
			{
				auto max_val_it = std::max_element(layer->children.begin(), layer->children.end(),
				                                   [](const std::shared_ptr<trie_based::NodeCount<TIndex>> &l,
				                                      const std::shared_ptr<trie_based::NodeCount<TIndex>> &r)
				{
					return l->index < r->index;
				});
				size_t max_ind = std::distance(layer->children.begin(), max_val_it);
				return std::make_pair(max_ind, get_grid_value(ind, layer->children[max_ind]->index) + 2.0*val01*dx[ind]);
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
			return std::make_pair(index, get_grid_value(ind, layer->children[index]->index) + 2.0*val01*dx[ind]);
		}
		size_t index = 0;
		TIndex target = m;
		for(size_t j = 1; j < layer->children.size(); j++)
		{
			if(layer->children[j]->index == target)
			{
				index = j;
				break;
			}
		}
		return std::make_pair(index, get_grid_value(ind, m) + (val01 - f1) * (get_grid_value(ind, m + 1) - get_grid_value(ind, m)) / (f2 - f1));
	}




	template <typename TIndex, typename TFloat>
	class ImplicitGraphQuantile : public Quantile<TIndex, TFloat>
	{
	protected:
		//using Quantile<TIndex, TFloat>::grids;
		using Quantile<TIndex, TFloat>::grid_number;
		using Quantile<TIndex, TFloat>::dx;
		using Quantile<TIndex, TFloat>::lb;
		using Quantile<TIndex, TFloat>::ub;

		using Quantile<TIndex, TFloat>::get_grid_value;

		typedef trie_based::TrieLayer<TIndex> sample_type;
		std::shared_ptr<sample_type> sample;

		size_t count_less_binary(const std::vector<TIndex> &layer, TIndex target) const;
		std::pair<size_t, TFloat> quantile_transform(const std::vector<TIndex> &layer, const std::vector<size_t> &psum, size_t ind, TFloat val01) const;
	public:
		ImplicitGraphQuantile() = default;
		ImplicitGraphQuantile(std::vector<TFloat> in_lb, std::vector<TFloat> in_ub, std::vector<size_t> in_gridn);
		using Quantile<TIndex, TFloat>::set_grid_and_gridn;
		void set_sample(const std::vector<std::vector<TIndex>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights) override;
		void set_sample(std::shared_ptr<sample_type> in_sample);
		void transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const override;
		void transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const override;
	};

	template <typename TIndex, typename TFloat>
	ImplicitGraphQuantile<TIndex, TFloat>::ImplicitGraphQuantile(std::vector<TFloat> in_lb,
	    std::vector<TFloat> in_ub,
	    std::vector<size_t> in_gridn) : Quantile<TIndex, TFloat>(in_lb, in_ub, in_gridn)
	{

	}
	template <typename TIndex, typename TFloat>
	void ImplicitGraphQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TIndex>> &in_sample)
	{
		sample = std::make_shared<sample_type>();
		sample->set_dimension(grid_number.size());
		for(const auto & i : in_sample)
			sample->insert(i);
		sample->sort();
	}
	template <typename TIndex, typename TFloat>
	void ImplicitGraphQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample)
	{
		sample = std::make_shared<sample_type>();
		sample->set_dimension(grid_number.size());
		for(size_t i = 0; i != in_sample.size(); ++i)
		{
			std::vector<TIndex> temp(in_sample[i].size());
			for(size_t j = 0; j != in_sample[i].size(); ++j)
			{
				temp[j] = get_the_closest_grid_node_to_the_value(lb[j], ub[j], grid_number[j], in_sample[i][j]);
			}
			sample->insert(temp);
		}
		sample->sort();
	}

	template <typename TIndex, typename TFloat>
	void ImplicitGraphQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights)
	{
		set_sample(in_sample);
	}

	template <typename TIndex, typename TFloat>
	void ImplicitGraphQuantile<TIndex, TFloat>::set_sample(std::shared_ptr<sample_type> in_sample)
	{
		sample = std::move(in_sample);
		sample->sort();
	}
	template <typename TIndex, typename TFloat>
	void ImplicitGraphQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const
	{
//    // first step
//    std::vector<size_t> psum(sample->root->second.size() + 1, 0);
//    for(size_t j = 1, k = 0; j != sample->root->second.size(); ++j)
//    {
//        k += 1;
//        psum[j] = k;
//    }
//    psum[sample->root->second.size()] = sample->root->second.size();
//
//    auto rez = quantile_transform(sample->root->second, psum, 0, in01[0]);
//    out[0] = rez.second;
//    //p = p->children[rez.first].get();


		auto p = sample->root;
		for(size_t i = 0; i != in01.size(); ++i)
		{
			std::vector<size_t> psum(p->second.size() + 1, 0);
			for(size_t j = 1, k = 0; j != p->second.size(); ++j)
			{
				k += 1;
				psum[j] = k;
			}
			psum[p->second.size()] = p->second.size();

			auto rez = quantile_transform(p->second, psum, i, in01[i]);
			out[i] = rez.second;
			if(i < in01.size() - 1)
				p = sample->layers[i + 1].find(p->second[rez.first]);
		}
	}
	template <typename TIndex, typename TFloat>
	void ImplicitGraphQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const
	{
		auto p = sample->root;
		for(size_t i = 0; i != in01.size(); ++i)
		{
			std::vector<size_t> psum(p->second.size() + 1, 0);
			for(size_t j = 1, k = 0; j != p->second.size(); ++j)
			{
				k += 1;
				psum[j] = k;
			}
			psum[p->second.size()] = p->second.size();

			auto rez = quantile_transform(p->second, psum, i, in01[i]);
			out[i] = rez.first;
			if(i < in01.size() - 1)
				p = sample->layers[i + 1].find(p->second[rez.first]);
		}
	}
	template <typename TIndex, typename TFloat>
	std::pair<size_t, TFloat> ImplicitGraphQuantile<TIndex, TFloat>::quantile_transform(const std::vector<TIndex> &layer, const std::vector<size_t> &psum, size_t ind, TFloat val01) const
	{
		size_t m = 0, count = grid_number[ind], step, c1 = 0, c2 = 0;
		TFloat f1 = 0.0, f2 = 0.0, sample_size_u = static_cast<TFloat>(layer.size());
		//auto first = grids[ind].begin();
		//auto it = grids[ind].begin();
		size_t it = 0, first = 0;

		while(count > 0)
		{
			it = first;
			step = count / 2;
			it += step;
			m = it;
			//std::advance(it, step);
			//m = std::distance(grids[ind].begin(), it);

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
				return std::make_pair(0, get_grid_value(ind, layer.front()) + 2.0*val01*dx[ind]);
			}
			if(c1 == layer.size())
			{
				return std::make_pair(layer.size() - 1, get_grid_value(ind, layer.back()) + 2.0*val01*dx[ind]);
			}

			TIndex target = m;
			auto pos = std::lower_bound(layer.begin(), layer.end(), target);
			size_t index = std::distance(layer.begin(), pos);

			if(index > 0)
			{
				int curr1 = std::abs(static_cast<int>(layer[index]) - static_cast<int>(m));
				int curr2 = std::abs(static_cast<int>(layer[index - 1]) - static_cast<int>(m));

				if(curr1 < curr2)
				{
					return std::make_pair(index, get_grid_value(ind, layer[index]) + 2.0*val01*dx[ind]);
				}
				else if(curr1 == curr2)
				{
					if(layer[index - 1] < layer[index])
						return std::make_pair(index - 1, get_grid_value(ind, layer[index - 1]) + 2.0*val01*dx[ind]);
					else
						return std::make_pair(index, get_grid_value(ind, layer[index]) + 2.0*val01*dx[ind]);
				}
				else
				{
					return std::make_pair(index - 1, get_grid_value(ind, layer[index - 1]) + 2.0*val01*dx[ind]);
				}

			}
			return std::make_pair(index, get_grid_value(ind, layer[index]) + 2.0*val01*dx[ind]);
		}
		TIndex target = m;
		auto pos = std::lower_bound(layer.begin(), layer.end(), target);
		TIndex index = std::distance(layer.begin(), pos);
		return std::make_pair(index, get_grid_value(ind, m) + (val01 - f1) * (get_grid_value(ind, m + 1) - get_grid_value(ind, m)) / (f2 - f1));
	}
	template <typename TIndex, typename TFloat>
	size_t ImplicitGraphQuantile<TIndex, TFloat>::count_less_binary(const std::vector<TIndex> &layer, TIndex target) const
	{
		auto lb = std::lower_bound(layer.begin(), layer.end(), target);
		size_t pos = std::distance(layer.begin(), lb);
		if(lb == layer.end())
			pos = layer.size(); // to psum! which is layer->children.size() + 1
		return pos; // to psum!
	}


// fsa

	template <typename TIndex, typename TFloat>
	class GraphQuantile : public Quantile<TIndex, TFloat>
	{
	protected:
//		using Quantile<TIndex, TFloat>::grids;
		using Quantile<TIndex, TFloat>::grid_number;
		using Quantile<TIndex, TFloat>::dx;

		using Quantile<TIndex, TFloat>::get_grid_value;

		typedef trie_based::Graph<TIndex> sample_type;
		std::shared_ptr<sample_type> sample;

		size_t count_less_binary(const std::vector<trie_based::invect> &layer, TIndex target) const;
		std::pair<size_t, TFloat> quantile_transform(const std::vector<trie_based::invect> &layer, const std::vector<size_t> &psum, size_t ind, TFloat val01) const;
	public:
		GraphQuantile() = default;
		GraphQuantile(std::vector<TFloat> in_lb, std::vector<TFloat> in_ub, std::vector<size_t> in_gridn);
		using Quantile<TIndex, TFloat>::set_grid_and_gridn;
		void set_sample(const std::vector<std::vector<TIndex>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights) override;
		void transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const override;
		void transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const override;
	};

	template <typename TIndex, typename TFloat>
	GraphQuantile<TIndex, TFloat>::GraphQuantile(std::vector<TFloat> in_lb,
	    std::vector<TFloat> in_ub,
	    std::vector<size_t> in_gridn) : Quantile<TIndex, TFloat>(in_lb, in_ub, in_gridn)
	{
		sample = std::make_shared<sample_type>();
	}

	template <typename TIndex, typename TFloat>
	void GraphQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const
	{
		auto p = sample->root;
		for(size_t i = 0, k; i != in01.size(); ++i)
		{
			std::vector<size_t> psum(p->second.size() + 1, 0);
			for(size_t j = 1, m = 0; j != p->second.size(); ++j)
			{
				m += p->second[j-1].count;
				psum[j] = m;
				psum.back() += p->second[j-1].count;
			}
			psum.back() += p->second.back().count;
			std::tie(k, out[i]) = quantile_transform(p->second, psum, i, in01[i]);
			p = sample->layers.find(p->second[k].vname);
		}
	}
	template <typename TIndex, typename TFloat>
	void GraphQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TIndex>> &in_sample)
	{

	}
	template <typename TIndex, typename TFloat>
	void GraphQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample)
	{

	}
	template <typename TIndex, typename TFloat>
	void GraphQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights)
	{

	}
	template <typename TIndex, typename TFloat>
	void GraphQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const
	{
		auto p = sample->root;
		for(size_t i = 0; i != in01.size(); ++i)
		{
			std::vector<size_t> psum(p->second.size() + 1, 0);
			for(size_t j = 1, m = 0; j != p->second.size(); ++j)
			{
				m += p->second[j-1].count;
				psum[j] = m;
				psum.back() += p->second[j-1].count;
			}
			psum.back() += p->second.back().count;
			//std::tie(k, out[i]) = quantile_transform(p->second, psum, i, in01[i]);
			auto [k, result] = quantile_transform(p->second, psum, i, in01[i]);
			p = sample->layers.find(p->second[k].vname);
		}
	}
	template <typename TIndex, typename TFloat>
	std::pair<size_t, TFloat> GraphQuantile<TIndex, TFloat>::quantile_transform(const std::vector<trie_based::invect> &layer, const std::vector<size_t> &psum, size_t ind, TFloat val01) const
	{
		size_t m = 0, count = grid_number[ind], step, c1 = 0, c2 = 0;
		TFloat f1 = 0.0, f2 = 0.0, sample_size_u = static_cast<TFloat>(psum.back());
		//auto first = grids[ind].begin();
		//auto it = grids[ind].begin();
		size_t it = 0, first = 0;

		while(count > 0)
		{
			it = first;
			step = count / 2;
			it += step;
			m = it;
//			std::advance(it, step);
//			m = std::distance(grids[ind].begin(), it);

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
				return std::make_pair(0, get_grid_value(ind, layer.front().index) + 2.0*val01*dx[ind]);
			}
			if(c1 == psum.back())
			{
				return std::make_pair(layer.size() - 1, get_grid_value(ind, layer.back().index) + 2.0*val01*dx[ind]);
			}

			TIndex target = m;
			auto pos = std::lower_bound(layer.begin(), layer.end(), target,
			                            [](const trie_based::invect &l,
			                               const TIndex &r)
			{
				return l.index < r;
			});
			size_t index = std::distance(layer.begin(), pos);

			if(index > 0)
			{
				int curr1 = std::abs(static_cast<int>(layer[index].index) - static_cast<int>(m));
				int curr2 = std::abs(static_cast<int>(layer[index - 1].index) - static_cast<int>(m));

				if(curr1 < curr2)
				{
					return std::make_pair(index, get_grid_value(ind, layer[index].index) + 2.0*val01*dx[ind]);
				}
				else if(curr1 == curr2)
				{
					if(layer[index - 1].index < layer[index].index)
						return std::make_pair(index - 1, get_grid_value(ind, layer[index - 1].index) + 2.0*val01*dx[ind]);
					else
						return std::make_pair(index, get_grid_value(ind, layer[index].index) + 2.0*val01*dx[ind]);
				}
				else
				{
					return std::make_pair(index - 1, get_grid_value(ind, layer[index - 1].index) + 2.0*val01*dx[ind]);
				}

			}
			return std::make_pair(index, get_grid_value(ind, layer[index].index) + 2.0*val01*dx[ind]);
		}
		TIndex target = m;
		auto pos = std::lower_bound(layer.begin(), layer.end(), target,
		                            [](const trie_based::invect &l,
		                               const TIndex &r)
		{
			return l.index < r;
		});
		TIndex index = std::distance(layer.begin(), pos);
		return std::make_pair(index, get_grid_value(ind, m) + (val01 - f1) * (get_grid_value(ind, m + 1) - get_grid_value(ind, m)) / (f2 - f1));
	}
	template <typename TIndex, typename TFloat>
	size_t GraphQuantile<TIndex, TFloat>::count_less_binary(const std::vector<trie_based::invect> &layer, TIndex target) const
	{
		auto lb = std::lower_bound(layer.begin(), layer.end(), target,
		                           [](const trie_based::invect &l,
		                              const TIndex &r)
		{
			return l.index < r;
		});
		size_t pos = std::distance(layer.begin(), lb);
		if(lb == layer.end())
			pos = layer.size(); // to psum! which is layer->children.size() + 1
		return pos; // to psum!
	}

	template <typename TIndex, typename TFloat>
	class ImplicitTrieQuantile : public ImplicitQuantile<TIndex, TFloat>
	{
	protected:
//		using ImplicitQuantile<TIndex, TFloat>::grids;
//		using ImplicitQuantile<TIndex, TFloat>::dx;
		using ImplicitQuantile<TIndex, TFloat>::grid_number;
		using ImplicitQuantile<TIndex, TFloat>::lb;
		using ImplicitQuantile<TIndex, TFloat>::ub;

		typedef trie_based::Trie<trie_based::NodeCount<TIndex>,TIndex> trie_type;
		std::shared_ptr<trie_type> sample;

//		using ImplicitQuantile<TIndex, TFloat>::count_less;
		using ImplicitQuantile<TIndex, TFloat>::quantile_transform;
//		const std::vector<size_t> weights;
	public:
		ImplicitTrieQuantile() = default;
		ImplicitTrieQuantile(std::vector<TFloat> in_lb, std::vector<TFloat> in_ub, std::vector<size_t> in_gridn);
		void set_sample_shared(std::shared_ptr<trie_type> in_sample);
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights) override;
		void transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const override;
		void transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const override;
		using ImplicitQuantile<TIndex, TFloat>::get_the_closest_grid_node_to_the_value;
	};

	template <typename TIndex, typename TFloat>
	ImplicitTrieQuantile<TIndex, TFloat>::ImplicitTrieQuantile(std::vector<TFloat> in_lb,
	    std::vector<TFloat> in_ub,
	    std::vector<size_t> in_gridn) : ImplicitQuantile<TIndex, TFloat>(in_lb, in_ub, in_gridn)
	{
	}
	template <typename TIndex, typename TFloat>
	void ImplicitTrieQuantile<TIndex, TFloat>::set_sample_shared(std::shared_ptr<trie_type> in_sample)
	{
		sample = std::move(in_sample);
	}

	template <typename TIndex, typename TFloat>
	void ImplicitTrieQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const
	{
		auto p = sample->root.get();
		for(size_t i = 0, k; i != in01.size(); i++)
		{
			std::tie(k, out[i]) = quantile_transform(p, i, in01[i]);
			p = p->children[k].get();
		}
	}

	template <typename TIndex, typename TFloat>
	void ImplicitTrieQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const
	{
		auto p = sample->root.get();
		for(size_t i = 0; i != in01.size(); i++)
		{
			//std::tie(k, out[i]) = quantile_transform(p, i, in01[i]);
			auto [k, result] = quantile_transform(p, i, in01[i]);
			out[i] = p->children[k]->index;
			p = p->children[k].get();
		}
	}

	template <typename TIndex, typename TFloat>
	void ImplicitTrieQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample)
	{
		sample = std::make_shared<trie_type>();
		sample->set_dimension(grid_number.size());
		for(size_t i = 0; i != in_sample.size(); ++i)
		{
			std::vector<TIndex> temp(in_sample[i].size());
			for(size_t j = 0; j != in_sample[i].size(); ++j)
			{
				temp[j] = get_the_closest_grid_node_to_the_value(lb[j], ub[j], grid_number[j], in_sample[i][j]);
			}
//			for(size_t j = 0; j != temp.size(); j++)
//			{
//				std::cout << temp[j] << ' ';
//			}
//			std::cout << std::endl;
			sample->insert(temp, 1);
		}
		//sample->fill_tree_count();
	}

	template <typename TIndex, typename TFloat>
	void ImplicitTrieQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights)
	{
		sample = std::make_shared<trie_type>();
		sample->set_dimension(grid_number.size());
		for(size_t i = 0; i != in_sample.size(); ++i)
		{
			std::vector<TIndex> temp(in_sample[i].size());
			for(size_t j = 0; j != in_sample[i].size(); ++j)
			{
				temp[j] = get_the_closest_grid_node_to_the_value(lb[j], ub[j], grid_number[j], in_sample[i][j]);
			}
//			for(size_t j = 0; j != temp.size(); j++)
//			{
//				std::cout << temp[j] << ' ';
//			}
//			std::cout << std::endl;
			sample->insert(temp, weights[i]);
		}
		//sample->fill_tree_count();
	}
}

#endif
