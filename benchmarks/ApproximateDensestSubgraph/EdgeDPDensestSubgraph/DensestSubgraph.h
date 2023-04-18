// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include "gbbs/edge_map_reduce.h"
#include "gbbs/gbbs.h"

#include <set>

#define C_0 1. / 20.
#define C_1 1
#define C_2 1
#define C 1

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

namespace gbbs {

struct edge {
  size_t u;
  size_t v;
  size_t load;

  edge(size_t _u, size_t _v, size_t _load) : u(_u), v(_v), load(_load) {}
  edge(size_t _u, size_t _v) : u(_u), v(_v), load(0) {}

  void increase_load(size_t addtl_load) { load += addtl_load; }
};

/* 
Determines how many edges each vertex has which are connected to another vertex in the set of vertices, sums, then divides by the number of vertices. This does not seem to be the traditional definition of density, but other implementations in this repo seem to use it.
*/
template <class Graph>
double graphDensity(Graph &G, std::set<size_t> vertices){
  int relevant_edges = 0;
  for (size_t v : vertices){
    auto neighbors = G.get_vertex(v).out_neighbors();
    for (size_t i = 0; i < G.get_vertex(v).out_degree(); i++){
      if (vertices.find(std::get<0>(neighbors.get_ith_neighbor(i))) != vertices.end()){
        relevant_edges++;
      }
    }
  }
  return ((double)relevant_edges) / ((double)vertices.size());
}

template <class Graph>
double EdgeDPDensestSubgraph(Graph& G, double density, double nu, double epsilon){
  const size_t n = G.n;
  // const size_t m = G.m;
  const double T = C_0 * std::log(n) / (std::pow(nu, 3));

  std::cout << "T = " << T << std::endl;

  // geometric distributions
  std::default_random_engine generator;
  std::geometric_distribution<size_t> x_distribution(epsilon / (6 * T * std::log(n) / std::log(1 + nu)));
  std::geometric_distribution<size_t> y_distribution(epsilon / (3 * T * (4 * T + 1) * std::log(n) / std::log(1 + nu)));
  std::geometric_distribution<size_t> z_distribution(epsilon / (6 * T * (4 * T + 1) * std::log(n) / std::log(1 + nu)));

  // initialize a sequence of sequences to store pointers to edge objects for each vertex
  parlay::sequence<parlay::sequence<edge*>> edge_maps = parlay::tabulate(n, [&](size_t i) {
    return parlay::sequence<edge*>();
  });
  // initialize the edge objects
  for (size_t v = 0; v < n; v++){
    auto neighbors = G.get_vertex(v).out_neighbors();
    for (size_t i = 0; i < G.get_vertex(v).out_degree(); i++){
      size_t u = std::get<0>(neighbors.get_ith_neighbor(i));
      if (u < v){
        edge *e = new edge(u, v);
        edge_maps[v].push_back(e);
        edge_maps[u].push_back(e);
      }
    }
  }


  if (density == 0.) { // lines 6-7 of algorithm 7
    return 0.;
  }

  for(size_t t = 1; t <= T; t++){ 
    // std::cout << "t = " << t << std::endl;

    auto X = parlay::tabulate(n, [&](size_t i) { 
      return x_distribution(generator);
    });

    parlay::sequence<parlay::sequence<int>> alphas(n); // this might be a memory leak

    for (size_t v = 0; v < n; v++){
      auto edges = edge_maps[v];
      // sort edges by load, breaking ties by other vertex id
      edges = parlay::sort(edges, [&](edge* e1, edge* e2) {
        if (e1->load == e2->load){
          size_t a = e1->u == v ? e1->v : e1->u;
          size_t b = e2->u == v ? e2->v : e2->u;
          return a < b;
        } else {
          return e1->load < e2->load;
        }
      });
      size_t max_2_alpha = min(std::ceil(density / 2) - 1 + X[v], edges.size());
      auto alpha = parlay::tabulate(edges.size(), [&](size_t i) {
        return i < max_2_alpha ? 2 : 0;
      });
      alphas[v] = alpha;
    }
    
    for (size_t load = 0; load < 4*T; load++){
      std::set<size_t> V = std::set<size_t>();

      auto Z = parlay::tabulate(n, [&](size_t i) { 
        return z_distribution(generator);
      });

      for (size_t v = 0; v < n; v++){
        int qualifying_edges = 0;
        auto edges = edge_maps[v];
        for(auto e : edges){
          if (e->load <= load){
            qualifying_edges++;
          }
        }
        if (qualifying_edges >= Z[v] + std::ceil(density / 2) - (C_1 * std::pow(std::log(n), 4)) / epsilon){
          V.insert(v);
        }
      }

      auto Y = y_distribution(generator);

      // if the subgraph induced by the vertices in V is sufficiently dense, return it
      double current_density = graphDensity(G, V);

      if (current_density >= density + Y - (C_2 * std::pow(std::log(n), 4)) / epsilon){
        std::cout << "density = " << current_density << std::endl;
        return current_density;
      }
    }
    // update the load of each edge
    for (size_t v = 0; v < n; v++){
      auto edges = edge_maps[v];
      for (uint i = 0; i < edges.size(); i++){
        edges[i]->increase_load(X[v] + alphas[v][i]);
      }
    }
  }

  std::cout << "no subgraph found" << std::endl;
  return 0.;
}

}  // namespace gbbs
