#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <set>
#include <utility>
using namespace Rcpp;



// Internal tree structure: adjacency list built once and then called multiple times by sampleTaxonPairs for efficiency
struct TreeAdjacency {
  std::unordered_map<int, std::vector<int>> children;
  std::unordered_map<int, int> parent;
  std::unordered_map<int, int> depth;  // NEW: cached depth per node
};


//' Build a cached tree adjacency structure
//'
//' @param edge Integer matrix, the tree's edge matrix
//' @param root the node number of the root
//' @return An external pointer to the cached structure
//' @noRd
// [[Rcpp::export]]
SEXP rcpp_buildTree(IntegerMatrix edge, int root) {
  TreeAdjacency* tree = new TreeAdjacency();
  int nedge = edge.nrow();
  tree->children.reserve(nedge);
  tree->parent.reserve(nedge);

  for (int i = 0; i < nedge; i++) {
    int p = edge(i, 0);
    int c = edge(i, 1);
    tree->children[p].push_back(c);
    tree->parent[c] = p;
  }

  // Breadth first search from root to compute depth for every node once
  tree->depth[root] = 0;
  std::vector<int> queue;
  queue.push_back(root);
  size_t head = 0;
  while (head < queue.size()) {
    int node = queue[head++];
    auto it = tree->children.find(node);
    if (it != tree->children.end()) {
      for (int child : it->second) {
        tree->depth[child] = tree->depth[node] + 1;
        queue.push_back(child);
      }
    }
  }

  XPtr<TreeAdjacency> ptr(tree, true);
  return ptr;
}


//' Get leaf descendants using a cached tree structure
//'
//' @param tree_ptr External pointer from rcpp_buildTree()
//' @param node Integer node ID
//' @param nleaves Integer number of leaves
//' @return Integer vector of descendant node IDs
//' @noRd
// [[Rcpp::export]]
IntegerVector rcpp_getDescendantsFast(SEXP tree_ptr, int node, int nleaves) {
  XPtr<TreeAdjacency> tree(tree_ptr);

  std::vector<int> result;
  std::vector<int> stack;
  stack.push_back(node);

  while (!stack.empty()) {
    int current = stack.back();
    stack.pop_back();

    if (current <= nleaves){
      result.push_back(current);
    }

    auto it = tree->children.find(current);
    if (it != tree->children.end()) {
      for (int child : it->second) {
        stack.push_back(child);
      }
    }
  }

  return wrap(result);
}




//' Get leaves that are NOT descendants using a cached tree structure
//'
//' @param tree_ptr External pointer from rcpp_buildTree()
//' @param node Integer node ID
//' @param nleaves Integer number of leaves
//' @return Integer vector of descendant node IDs
//' @noRd
// [[Rcpp::export]]
IntegerVector rcpp_getNonDescendantsFast(SEXP tree_ptr, int node, int nleaves) {

  XPtr<TreeAdjacency> tree(tree_ptr);

  std::vector<bool> is_desc(nleaves + 1, false);  // 1-indexed, tips are 1..n

  std::vector<int> stack;
  stack.push_back(node);
  while (!stack.empty()) {
    int current = stack.back();
    stack.pop_back();
    if (current <= nleaves) is_desc[current] = true;

    auto it = tree->children.find(current);
    if (it != tree->children.end()) {
      for (int child : it->second) {
        stack.push_back(child);
      }
    }
  }

  std::vector<int> result;
  result.reserve(nleaves);
  for (int i = 1; i <= nleaves; i++) {
    if (!is_desc[i]) result.push_back(i);
  }

  return wrap(result);
}



//' Get the path between two nodes in a tree
//' @param tree_ptr External pointer from rcpp_buildTree()
//' @param from Node number at the start of the path
//' @param to Node number at the start of the path
//' @noRd
// [[Rcpp::export]]
IntegerVector rcpp_nodepathFast(SEXP tree_ptr, int from, int to) {
  XPtr<TreeAdjacency> tree(tree_ptr);

  auto depth_it_from = tree->depth.find(from);
  auto depth_it_to = tree->depth.find(to);
  if (depth_it_from == tree->depth.end()) {
    stop("Node 'from' (%d) not found in tree structure", from);
  }
  if (depth_it_to == tree->depth.end()) {
    stop("Node 'to' (%d) not found in tree structure", to);
  }

  std::vector<int> up_from;
  std::vector<int> up_to;

  int a = from, b = to;
  int da = depth_it_from->second;
  int db = depth_it_to->second;

  up_from.push_back(a);
  up_to.push_back(b);

  while (da > db) {
    auto it = tree->parent.find(a);
    if (it == tree->parent.end()) {
      stop("unexpected dev error 2541: Reached node with no parent while ascending 'from' (node %d)", a);
    }
    a = it->second;
    da--;
    up_from.push_back(a);
  }
  while (db > da) {
    auto it = tree->parent.find(b);
    if (it == tree->parent.end()) {
      stop("unexpected dev error 2542: Reached node with no parent while ascending 'to' (node %d)", b);
    }
    b = it->second;
    db--;
    up_to.push_back(b);
  }

  int safety = 0;
  int max_iter = tree->depth.size() + 5;  // cannot legitimately exceed tree depth
  while (a != b) {
    if (++safety > max_iter) {
      stop("unexpected dev error 2543: nodepath failed to converge");
    }
    auto ita = tree->parent.find(a);
    auto itb = tree->parent.find(b);
    if (ita == tree->parent.end() || itb == tree->parent.end()) {
      stop("unexpected dev error 2544: Reached root without convergence");
    }
    a = ita->second;
    b = itb->second;
    up_from.push_back(a);
    up_to.push_back(b);
  }

  std::vector<int> result(up_from);
  result.insert(result.end(), up_to.rbegin() + 1, up_to.rend());

  return wrap(result);
}



//' Fast MRCA using cached depth + parent lookup
//'
//' @param tree_ptr External pointer from rcpp_buildTree()
//' @param a Node number of one descendent
//' @param b Node number of the other descendent
//' @return Integer node ID of the most recent common ancestor
//' @noRd
// [[Rcpp::export]]
int rcpp_getMRCA(SEXP tree_ptr, int a, int b) {
  XPtr<TreeAdjacency> tree(tree_ptr);

  if (a == b) return a;

  int da = tree->depth[a];
  int db = tree->depth[b];

  // Bring deeper node up to same depth as the other
  while (da > db) {
    a = tree->parent[a];
    da--;
  }
  while (db > da) {
    b = tree->parent[b];
    db--;
  }

  // Walk both up together until they meet
  while (a != b) {
    a = tree->parent[a];
    b = tree->parent[b];
  }

  return a; 
}



