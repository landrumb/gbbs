#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>

#include "pbbslib/random_shuffle.h"
#include "gbbs/gbbs.h"
#include "pbbslib/strings/string_basics.h"

namespace gbbs {
template <class Graph>
void randomReorder(Graph& GA, std::string& outfile) {
  using W = typename Graph::weight_type;
  size_t n = GA.n;
  size_t m = GA.m;
  auto perm = pbbs::random_permutation<uintE>(n);

  auto edges = pbbs::sequence<uintE>(m);
  auto offs = pbbs::sequence<uintT>(n);
  parallel_for(0, n, [&] (size_t i) {
    offs[perm[i]] = GA.get_vertex(i).getOutDegree();
  });
  size_t tot = pbbslib::scan_add_inplace(offs.slice());
  cout << "m = " << m << " tot = " << tot << endl;

  parallel_for(0, n, [&] (size_t i) {
    size_t off = offs[perm[i]];
    size_t next_off = (perm[i] == (n-1)) ? m : offs[perm[i]+1];
    auto map_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
      uintE ngh_perm = perm[ngh];
      edges[off++] = ngh_perm;
    };
    GA.get_vertex(i).mapOutNgh(i, map_f, false);
    std::sort(edges.begin()+off, edges.begin()+next_off, std::less<uintE>());
  });

  std::ofstream file (outfile, std::ios::out | std::ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << outfile << std::endl;
    exit(0);
  }
  file << "AdjacencyGraph" << endl;
  file << n << endl;
  file << m << endl;
  auto off_chars = pbbslib::sequence_to_string(offs);
  auto edges_chars = pbbslib::sequence_to_string(edges);
  file.write(off_chars.begin(), off_chars.size());
  file.write(edges_chars.begin(), edges_chars.size());

  file.close();
  cout << "Done" << endl;
}

template <class Graph>
double Reorderer(Graph& GA, commandLine P) {
  auto outfile = P.getOptionValue("-of", "/ssd0/graphs/bench_experiments/out.adj");
  randomReorder(GA, outfile);
  exit(0);
  return 1.0;
}
}  // namespace gbbs

generate_main(gbbs::Reorderer, false);
