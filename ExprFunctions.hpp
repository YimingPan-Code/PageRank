/******************************************************************************
 * Copyright (c) 2015-2016, TigerGraph Inc.
 * All rights reserved.
 * Project: TigerGraph Query Language
 * udf.hpp: a library of user defined functions used in queries.
 *
 * - This library should only define functions that will be used in
 *   TigerGraph Query scripts. Other logics, such as structs and helper
 *   functions that will not be directly called in the GQuery scripts,
 *   must be put into "ExprUtil.hpp" under the same directory where
 *   this file is located.
 *
 * - Supported type of return value and parameters
 *     - int
 *     - float
 *     - double
 *     - bool
 *     - string (don't use std::string)
 *     - accumulators
 *
 * - Function names are case sensitive, unique, and can't be conflict with
 *   built-in math functions and reserve keywords.
 *
 * - Please don't remove necessary codes in this file
 *
 * - A backup of this file can be retrieved at
 *     <tigergraph_root_path>/dev_<backup_time>/gdk/gsql/src/QueryUdf/ExprFunctions.hpp
 *   after upgrading the system.
 *
 ******************************************************************************/

#ifndef EXPRFUNCTIONS_HPP_
#define EXPRFUNCTIONS_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <gle/engine/cpplib/headers.hpp>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <set>
#include <string.h>
#include <cmath>
#include <ctime>
#include <iostream>
#include <map>

/**     XXX Warning!! Put self-defined struct in ExprUtil.hpp **
 *  No user defined struct, helper functions (that will not be directly called
 *  in the GQuery scripts) etc. are allowed in this file. This file only
 *  contains user-defined expression function's signature and body.
 *  Please put user defined structs, helper functions etc. in ExprUtil.hpp
 */
#include "ExprUtil.hpp"
using std::vector;
using std::set;


namespace UDIMPL {
  typedef std::string string; //XXX DON'T REMOVE
  typedef unsigned long ulong;

  ulong *graph = NULL, begin_index[N + 1] = {}, count[N] = {}, old_count[N] = {};

  const double alpha = 0.15;

  std:map<Int,Float> PR;
  /****** BIULT-IN FUNCTIONS **************/
  /****** XXX DON'T REMOVE ****************/
  inline int64_t str_to_int (string str) {
    return atoll(str.c_str());
  }

  inline int64_t float_to_int (float val) {
    return (int64_t) val;
  }

  inline string to_string (double val) {
    char result[200];
    sprintf(result, "%g", val);
    return string(result);
  }

  template <typename tuple>
  inline uint64_t getDeltaQ (tuple tup) {
    return tup.deltaQ;
  }

  template<typename tup> 
  inline float getWeight(tup t) {
    return t.weight;
  }

template <class T>
double mean(const T X[])
{
	double sum = 0.0;
	for (int i = 0; i < N; ++i)
		sum += X[i];
	return sum / N;
}


template <class T>
double covariance(const T X[], const T Y[], const double x_mean, const double y_mean)
{
	double sum = 0.0;
	for (int i = 0; i < N; i++)
		sum += (X[i] - x_mean) * (Y[i] - y_mean);
	return sum / (N - 1);
}


template <class T>
double variance(const T X[], const double x_mean)
{
	return covariance(X, X, x_mean, x_mean);
}


//Returns the standard deviation of the sample
template <class T>
double standardDeviation(const T X[], const double x_mean)
{
	return sqrt(variance(X, x_mean));
}


//Calculate the correlation coefficient
template <class T>
double sampleCorrelationCoefficient(T X[], T Y[])
{
	double x_mean = mean(X);
	double y_mean = mean(Y);
	double x_variance = variance(X, x_mean);
	double y_variance = variance(Y, y_mean);
	double xy_covariance = covariance(X, Y, x_mean, y_mean);
	return xy_covariance / sqrt(x_variance * y_variance);
}
//init graph
void init_graph(MapAccum<INT,INT> map)
{
	srand(0);
	ulong size;
	begin_index[0] = 0;
	vector<ulong> g;
	set<ulong> out;
	ulong i=0;
	for (iter = map.begin(); iter != map.end(); iter++)
	{
		size = (map->second).size();
		begin_index[i + 1] = begin_index[i] + size;
		for (ulong j = 0; j < size; ++j)
		{
			out.emplace((map->second)[j]);
		}
		g.insert(g.end(), map->first);
		g.insert(g.end(), out.begin(), out.end()),++i;
	}
	graph = new ulong[g.size()];
	memcpy(graph, g.data(), g.size() * sizeof(ulong));
}
Map<INT,FLOAT> main(MapAccum<INT,SetAccum<INT> map)
{
	ulong current_vertex, out_size;
	clock_t start, end;
	init_graph();
	ulong total_iteration_time = 0;
	for (ulong source = 0; source < N; ++source)
	{
		if (begin_index[source] == begin_index[source + 1])
			continue;
		memset(count, 0, N * sizeof(ulong));
		memset(old_count, 0, N * sizeof(ulong));
		current_vertex = source;
		int iter;
		for (iter = 0; iter < MAX_ITER; ++iter)
		{
			out_size = begin_index[current_vertex + 1] - begin_index[current_vertex];
			if (out_size == 0 || rand() < alpha * RAND_MAX)
				current_vertex = source; //
			else
			{
				current_vertex = graph[begin_index[current_vertex] + rand() % out_size];
			}
			++count[current_vertex];
			if ((iter + 1) % STEP == 0)
			{
				if (iter + 1 > STEP)
				{
					double r = sampleCorrelationCoefficient(old_count, count);
					if (r > 0.999)
					{
						++iter;
						break;
					}
				}
				memcpy(old_count, count, N * sizeof(ulong));
			}
		}
		total_iteration_time += iter;
		for (ulong i = 0; i < N; ++i)
			PR.insert(g[begin_index[i]]double(count[i]) / (iter * N));
	}
	delete[] graph;
	return PR;
}
  
}
/****************************************/

#endif /* EXPRFUNCTIONS_HPP_ */
