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


#include <trie.hh>
#include <trie_based.hh>
#include <trie_node.hh>

#include <numeric>

namespace mveqf
{
//	template <typename T>
//	bool increase(const std::vector<std::vector<T>>& v, std::vector<std::size_t>& it)
//	{
//		for(std::size_t i = 0, size = it.size(); i != size; ++i)
//		{
//			const std::size_t index = size - 1 - i;
//			++it[index];
//			if(it[index] == v[index].size())
//				it[index] = 0;
//			else
//				return true;
//		}
//		return false;
//	}
//	template <typename T>
//	std::vector<T> do_job(const std::vector<std::vector<T>>& v, std::vector<std::size_t>& it)
//	{
//		std::vector<T> rez(v.size());
//		for(std::size_t i = 0, size = v.size(); i != size; ++i)
//			rez[i] = v[i][it[i]];
//		return rez;
//	}

	template <typename T>
	bool recur_one(size_t n,
	               size_t places_left,
	               bool hit_max,
	               std::vector<T> &sequnce,
	               std::vector<T> &result,
	               size_t& c,
	               const size_t count)
	{
		if(c == count)
			return false;
		if(places_left == 0)
		{
			++c;
			if(c == count)
			{
				result = sequnce;
				return true;
			}
			return false;
		}
		if(places_left == 1 && !hit_max)
		{
			sequnce.push_back(n);
			recur_one(n, 0, true, sequnce, result, c, count);
			sequnce.pop_back();
			return false;
		}
		for(size_t i = 0; i < n + 1; i++)
		{
			sequnce.push_back(i);
			recur_one(n, places_left - 1, i < n ? hit_max : true, sequnce, result, c, count);
			sequnce.pop_back();
		}
		return false;
	}

	template <typename T>
	std::vector<T> return_permutation(size_t n, size_t k, size_t count)
	{
		std::vector<T> result;
		size_t c = 0;
		for(size_t i = 0; i != n; i++)
		{
			std::vector<T> t;
			if(recur_one(i, k, false, t, result, c, count))
				break;
		}
		return result;
	}

	template <typename TIndex, typename TFloat>
	class Quantile
	{
	protected:
		std::vector<TFloat> lb;
		std::vector<TFloat> ub;
		std::vector<TFloat> dx;
		std::vector<size_t> grid_number;
		std::vector<TFloat> grid_ranges;
//		std::vector<std::vector<TFloat>> grids;

		size_t get_lb(TFloat lb, TFloat ub, size_t gridn, const TFloat &value) const;
		TFloat get_min_delta_from_grid_node(TFloat lb, TFloat ub, size_t gridn, TFloat value);
		size_t get_optimal_gridn_linear(const TFloat lb, const TFloat ub, const TFloat value, const TFloat delta);
	public:
		explicit Quantile();
		explicit Quantile(std::vector<TFloat> in_lb, std::vector<TFloat> in_ub, std::vector<size_t> in_gridn);
		void set_grid_and_gridn(std::vector<TFloat> in_lb, std::vector<TFloat> in_ub, std::vector<size_t> in_gridn);
		void set_grid_from_sample(std::vector<TFloat> in_lb, std::vector<TFloat> in_ub, const std::vector<std::vector<TFloat>> &in_sample);
		virtual void transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const = 0;
		virtual void transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const = 0;
		virtual void set_sample(const std::vector<std::vector<TIndex>> &in_sample) = 0;
		virtual void set_sample(const std::vector<std::vector<TFloat>> &in_sample) = 0;
		virtual void set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights) = 0;
		inline TFloat get_grid_value(size_t current_dimension, size_t index) const;
		std::vector<size_t> get_grid_number() const;
		std::vector<std::vector<TFloat>> get_real_node_values(const std::vector<std::vector<TFloat>> &in_sample) const;
		size_t get_the_closest_grid_node_to_the_value(TFloat lb, TFloat ub, size_t gridn, TFloat value) const;
		virtual ~Quantile() = 0;
	};

	template <typename TIndex, typename TFloat>
	Quantile<TIndex, TFloat>::Quantile()
	{
	}

	template <typename TIndex, typename TFloat>
	Quantile<TIndex, TFloat>::~Quantile()
	{
	}

	template <typename TIndex, typename TFloat>
	Quantile<TIndex, TFloat>::Quantile(std::vector<TFloat> in_lb,
	                                   std::vector<TFloat> in_ub,
	                                   std::vector<size_t> in_gridn)
	{
		set_grid_and_gridn(in_lb, in_ub, in_gridn);
	}

	template <typename TIndex, typename TFloat>
	std::vector<std::vector<TFloat>> Quantile<TIndex, TFloat>::get_real_node_values(const std::vector<std::vector<TFloat>> &in_sample) const
	{
		std::vector<std::vector<TFloat>> values = in_sample;
		for(size_t i = 0; i != in_sample.size(); i++)
		{
			for(size_t j = 0; j != in_sample[i].size(); j++)
			{
				values[i][j] = get_grid_value(j, get_the_closest_grid_node_to_the_value(lb[j], ub[j], grid_number[j], in_sample[i][j])) + dx[j];
			}
		}
		return values;
	}

	template <typename TIndex, typename TFloat>
	void Quantile<TIndex, TFloat>::set_grid_and_gridn(std::vector<TFloat> in_lb, std::vector<TFloat> in_ub, std::vector<size_t> in_gridn)
	{
		lb = in_lb;
		ub = in_ub;
		grid_number = in_gridn;

		dx.resize(grid_number.size());

//		grids.resize(grid_number.size());
//		for(size_t i = 0; i != grids.size(); i++)
//		{
//			std::vector<TFloat> grid(grid_number[i] + 1);
//			TFloat startp = lb[i];
//			TFloat endp = ub[i];
//			TFloat es = endp - startp;
//			for(size_t j = 0; j != grid.size(); j++)
//			{
//				grid[j] = startp + j*es/TFloat(grid_number[i]);
//			}
//			grids[i] = grid;
//			dx[i] = es/(TFloat(grid_number[i])*2);
//		}

		grid_ranges.resize(grid_number.size());
		for(size_t i = 0; i != grid_number.size(); i++)
		{
			grid_ranges[i] = ub[i] - lb[i];
			dx[i] = grid_ranges[i]/(TFloat(grid_number[i])*2);
		}
	}


	template <typename TIndex, typename TFloat>
	void Quantile<TIndex, TFloat>::set_grid_from_sample(const std::vector<TFloat> in_lb, const std::vector<TFloat> in_ub, const std::vector<std::vector<TFloat>> &in_sample)
	{
		lb = in_lb;
		ub = in_ub;
		grid_number.resize(lb.size(), 1);

//		TFloat min_range = std::numeric_limits<TFloat>::max();
//		for(size_t i = 0; i != lb.size(); i++)
//		{
//			TFloat range = ub[i] - lb[i];
//			if(min_range < range)
//				min_range = range;
//		}
//		TFloat init_delta = min_range/10;
//		for(size_t i = 0; i != in_sample.size(); i++)
//		{
//			for(size_t j = 0; j != in_sample[i].size(); j++)
//			{
//				size_t current_grid = get_optimal_gridn_linear(lb[j], ub[j], in_sample[i][j], init_delta);
//				if(current_grid > grid_number[j])
//					grid_number[j] = current_grid;
//			}
//		}

		//std::cout << "init grid" << std::endl;
		//for(size_t i = 0; i != grid_number.size(); i++)
		//{
		//	std::cout << grid_number[i] << ' ';
		//}
		//std::cout << std::endl;


//		std::vector<size_t> init_grid_number(grid_number.size());
//		for(size_t i = 0; i != init_grid_number.size(); i++)
//			init_grid_number[i] = 2000;//std::numeric_limits<TIndex>::max() - 1;
//		std::vector<std::vector<int>> variable_values(init_grid_number.size());
//		for(size_t i = 0; i != variable_values.size(); i++)
//		{
//			variable_values[i].resize(init_grid_number[i]);
//			for(size_t j = 0; j != variable_values[i].size(); j++)
//			{
//				variable_values[i][j] = j;
//			}
//		}
//		std::vector<std::size_t> it(variable_values.size(), 0);
//		std::vector<std::vector<int>> temp;
//
//		do
//		{
//			auto val = do_job(variable_values, it);
//			auto current_grid = grid_number;
//			for(size_t i = 0; i != current_grid.size(); i++)
//			{
//				current_grid[i] += val[i];
//				std::cout << current_grid[i] << ' ';
//			}
//			std::cout << std::endl;
////			std::cin.get();
//
//			bool unique = true;
//			auto sample = std::make_shared<trie_based::TrieBased<trie_based::NodeCount<int>,int>>();
//			sample->set_dimension(current_grid.size());
//			for(size_t i = 0; i != in_sample.size(); ++i)
//			{
//				std::vector<TIndex> temp(in_sample[i].size());
//				for(size_t j = 0; j != in_sample[i].size(); ++j)
//				{
//					temp[j] = get_the_closest_grid_node_to_the_value(lb[j], ub[j], current_grid[j], in_sample[i][j]);
//				}
//				if(!sample->search(temp))
//					sample->insert(temp);
//				else
//				{
//					unique = false;
//					break;
//				}
//			}
//			if(unique)
//			{
//				grid_number = current_grid;
//
//				std::cout << "new grid" << std::endl;
//				for(size_t i = 0; i != grid_number.size(); i++)
//				{
//					std::cout << grid_number[i] << ' ';
//				}
//				std::cout << std::endl;
//
//				break;
//			}
//		}
//		while(increase(variable_values, it));

		std::vector<size_t> init_grid_number(grid_number.size());
		for(size_t i = 0; i != init_grid_number.size(); i++)
			init_grid_number[i] = std::numeric_limits<TIndex>::max() > 256 ? 1024 : 256;
		size_t search_range = std::numeric_limits<TIndex>::max() > 256 ? 1024 : 256;
		bool all_unique = false;
		for(size_t count = 1; !all_unique; count++)
		{
			auto val = return_permutation<TIndex>(search_range, init_grid_number.size(), count);
			auto current_grid = grid_number;
			for(size_t i = 0; i != current_grid.size(); i++)
			{
				current_grid[i] += val[i];
			}
//			std::cin.get();

			all_unique = true;
			auto sample = std::make_shared<trie_based::TrieBased<trie_based::NodeCount<TIndex>,TIndex>>();
			sample->set_dimension(current_grid.size());
			for(size_t i = 0; i != in_sample.size(); ++i)
			{
				std::vector<TIndex> temp(in_sample[i].size());
				for(size_t j = 0; j != in_sample[i].size(); ++j)
				{
					temp[j] = get_the_closest_grid_node_to_the_value(lb[j], ub[j], current_grid[j], in_sample[i][j]);
				}
				if(!sample->search(temp))
					sample->insert(temp);
				else
				{
					all_unique = false;
					break;
				}
			}
			if(all_unique)
				grid_number = current_grid;
		}
		//std::cout << "new grid" << std::endl;
		//for(size_t i = 0; i != grid_number.size(); i++)
		//{
		//	std::cout << grid_number[i] << ' ';
		//}
		//std::cout << std::endl;

		dx.resize(grid_number.size());
		grid_ranges.resize(grid_number.size());
		for(size_t i = 0; i != grid_number.size(); i++)
		{
			grid_ranges[i] = ub[i] - lb[i];
			dx[i] = grid_ranges[i]/(TFloat(grid_number[i])*2);
		}
	}

	template <typename TIndex, typename TFloat>
	std::vector<size_t> Quantile<TIndex, TFloat>::get_grid_number() const
	{
		return grid_number;
	}

	template <typename TIndex, typename TFloat>
	inline TFloat Quantile<TIndex, TFloat>::get_grid_value(size_t current_dimension, size_t index) const
	{
		return lb[current_dimension] + index*grid_ranges[current_dimension]/TFloat(grid_number[current_dimension]);
	}


	template <typename TIndex, typename TFloat>
	size_t Quantile<TIndex, TFloat>::get_lb(TFloat lb, TFloat ub, size_t gridn, const TFloat &value) const
	{
		size_t it, first = 0;
		int count = gridn + 1, step;
		TFloat range_normalized = (ub - lb)/TFloat(gridn);
		while(count > 0)
		{
			it = first;
			step = count / 2;
			it += step;

			if(lb + range_normalized*(it + 0.5) < value)
			{
				first = ++it;
				count -= step + 1;
			}
			else
				count = step;
		}
		return first;
	}
	template <typename TIndex, typename TFloat>
	size_t Quantile<TIndex, TFloat>::get_the_closest_grid_node_to_the_value(TFloat lb, TFloat ub, size_t gridn, TFloat value) const
	{
		TFloat range_normalized = (ub - lb)/TFloat(gridn);
		size_t it = get_lb(lb, ub, gridn, value);
		if(it >= gridn)
			return gridn - 1;
		else if(it == 0)
			return 0;
		else
			return std::abs(lb + range_normalized*(it - 0.5) - value) < std::abs(lb + range_normalized*(it + 0.5) - value) ? it - 1 : it;
	}

	template <typename TIndex, typename TFloat>
	TFloat Quantile<TIndex, TFloat>::get_min_delta_from_grid_node(TFloat lb, TFloat ub, size_t gridn, TFloat value)
	{
		size_t ind = get_the_closest_grid_node_to_the_value(lb, ub, gridn, value);
		TFloat range_normalized = (ub - lb)/TFloat(gridn);
		TFloat node_value = lb + range_normalized*(ind + 0.5);
		return std::abs(node_value - value);
	}

	template <typename TIndex, typename TFloat>
	size_t Quantile<TIndex, TFloat>::get_optimal_gridn_linear(const TFloat lb, const TFloat ub, const TFloat value, const TFloat delta)
	{
		TIndex max_gridn = std::numeric_limits<TIndex>::max(), grid_size = 1;
		while(get_min_delta_from_grid_node(lb, ub, grid_size, value) > delta)
		{
			++grid_size;
		}
		return size_t(grid_size);
	}
}


#endif
