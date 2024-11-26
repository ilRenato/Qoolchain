//
// Created by W10 on 08/11/2022.
//

#include "GraphFunctions.h"
#include <stack>
#include <numeric>
#include <map>
#include <vector>
#include <list>

using namespace std;

/* Breadth-first search */
vector<bool> BFS(vector<vector<int> > &G, int node, bool reverse){
    vector<bool> visited;
    list<int> queue;
    int num_nodes = (int) G.size();
    // initialize visited array
    visited.resize(num_nodes, false);
    // mark as visited the initial node
    visited[node] = true;
    // insert the initial node in the queue
    queue.push_back(node);
    while (!queue.empty()){
        // remove an element from the front of the queue
        node = queue.front();
        queue.pop_front();
        // scan all the nodes reached by an arc from current node
        for(int i=0; i<num_nodes; i++){
            bool condition = (reverse && (G[i][node] != 0)) || (!reverse && (G[node][i] != 0));
            if (condition){
                // check if node already visited
                if (!visited[i]){
                    // mark the node as visited
                    visited[i] = true;
                    // put the node in the queue
                    queue.push_back(i);
                }
            }
        }
    }
    return visited;
}

/* Recursive function used by DFS */
static void DFS_recursion(const vector<vector<int> > &G, vector<bool> &visited, int &node){
    visited[node] = true;
    for(int i=0; i<G.size(); i++){
        if ((G[node][i] != 0) && (!visited[i]))
            DFS_recursion(G, visited, i);
    }
}

/* Recursive Depth-first search */
vector<bool> DFS(vector<vector<int> > &G, int node){
    int num_nodes = (int) G.size();
    vector<bool> visited(num_nodes, false);
    DFS_recursion(G, visited, node);
    return visited;
}

/* Function that detects the strongly connected components (SCCs) in the adjacency matrix of a graph passed as parameter.
 * SCCs are subgraphs in which each node is reachable from any other node. */
list<vector<int> > stronglyConnectedComponents(vector<vector<int> > &G){
    int numNodes = (int) G.size();
    int disc_counter = 0;
    vector<bool> visited(numNodes, false);
    vector<int> parent(numNodes);
    /* Initialize parent of each node as the node itself.
     * In this way when a DFS ends the next DFS will start on successive node. */
    iota(parent.begin(), parent.end(), 0);
    bool edge_found;
    vector<int> disc(numNodes,-1);
    vector<int> low(numNodes,numNodes);
    stack<int> scc_stack;
    vector<bool> in_stack(numNodes,false);
    list<vector<int> > scc;
    /* Find SSCs with a DFS starting from each node  */
    for (int node=0; node<numNodes; node++){
        if (disc[node] == -1){
            disc[node] = low[node] = disc_counter;
            disc_counter++;
            scc_stack.push(node);
            in_stack[node] = true;
            while (!visited[node]){
                edge_found = false;
                for (int i=0; (i<numNodes) && (!edge_found); i++){
                    if (G[node][i] != 0){
                        if (disc[i] == -1){
                            disc[i] = low[i] = disc_counter;
                            disc_counter++;
                            scc_stack.push(i);
                            in_stack[i] = true;
                            parent[i] = node;
                            node = i;
                            edge_found = true;
                        } else {
                            if (in_stack[i]){
                                low[node] = min(low[node], disc[i]);
                            }
                        }
                    }
                }
                if (!edge_found){
                    int scc_node = -1;
                    if (low[node] == disc[node]){
                        scc.emplace_back(vector<int>());
                        do {
                            scc_node = scc_stack.top();
                            scc.back().emplace_back(scc_node);
                            in_stack[scc_node] = false;
                            scc_stack.pop();
                        } while (scc_node != node);
                        /* Don't want to consider SCCs of just one element so delete them in the vector */
                        if (scc.back().size() == 1)
                            scc.pop_back();
                    }
                    visited[node] = true;
                    low[parent[node]] = min(low[parent[node]], low[node]);
                    node = parent[node];
                }
            }
        }
    }
    return scc;
}

/* function used to trace back to the oldest ancestor found for a given node */
int getAncestor(vector<int> &parent, int node){
    if (parent[node] == node)
        return node;
    return getAncestor(parent, parent[node]);
}

list<vector<int> > connectedComponents(vector<vector<int> > &G){
    int numNodes = (int) G.size();
    vector<int> parent(numNodes);
    /* Initialize parent of each node as the node itself */
    iota(parent.begin(), parent.end(), 0);
    /* navigate all edges so that going from any node to its parent and so on a common ancestor for each
     * subgraph is reached */
    for (int i=0; i<numNodes; i++){
        for (int j=0; j<numNodes; j++){
            if (G[i][j] != 0)
                parent[getAncestor(parent, i)] = getAncestor(parent, j);
        }
    }
    /* set as parent of each node the ancestor reachable from it. In this way all the nodes of a subgraph have
     * as parent a common ancestor */
    for (int i=0; i<numNodes; i++){
        parent[i] = getAncestor(parent, parent[i]);
    }
    /* group all nodes that have a common ancestor, that are the different subgraphs found */
    map<int, vector<int> > m;
    for (int i=0; i<numNodes; i++){
        m[parent[i]].emplace_back(i);
    }
    list<vector<int> > CC;
    for (auto & element : m){
        CC.push_back(element.second);
    }
    return CC;
}