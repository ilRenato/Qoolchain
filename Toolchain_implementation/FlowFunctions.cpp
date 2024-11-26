//
// Created by W10 on 13/11/2022.
//

#include "FlowFunctions.h"
#include <list>
#include <algorithm>
#include <vector>
#include <climits>

using namespace std;

/* Flow network matrix */
static vector<vector<int> > G;
static vector<list<int> > _G; //CONTROLLARE SE MEGLIO MATRICE O LISTA
/* Residual network matrix to be constructed */
static vector<vector<int> > G_res;
/* vector of excess function per vertex */
static vector<int> excess;
/* vector of distance function per vertex */
static vector<int> d;
/* array of lists of nodes with label equal to the index in which they are stored */
static vector<list<int> > B;
/* highest label among the active nodes */
static int b;
/* number of variables in the flow network */
static int numVars;
/* number of nodes in the flow network */
static int numNodes;

/* Function to check whether a vertex is active or not */
static bool isActive(int v){
    if ((excess[v] > 0) && (d[v] < numNodes))
        return true;
    else
        return false;
}

/* Function to move excess from a vertex v to a vertex w
 * An amount of flow equal to the minimum of the capacity and the residual of the edge is added to the flow (v,w),
 * subtracted from the flow (w,v), added to the excess of the vertex w and subtracted from the excess of the vertex v */
static void push(int v, int w){
    int delta = G_res[v][w] < excess[v] ? G_res[v][w] : excess[v];
    G[v][w] -= delta;
    G_res[v][w] -= delta;
    G_res[w][v] += delta;
    excess[v] -= delta;
    excess[w] += delta;
}

/* Function to increase the label of a vertex v: d(v) = d(w) + 1
 * where w is the vertex reached by an edge outgoing from v that has the minimum d(w) */
static void relabel(int v){
    int min = numNodes;
    bool min_found = false;
    for(int i=0; i<numNodes; i++){
        if (G_res[v][i] > 0){
            if (d[i] < min) {
                min = d[i];
                min_found = true;
            }
        }
    }
    if (min_found)
        d[v] = min + 1;
    else
        d[v] = numNodes;
}

/* Function that sets the order of push and relabel operations.
 * If applicable a push on the first edge available is performed,
 * otherwise the current edge is substituted by a new one,
 * or if the current edge is the last edge of v it is substituted with the first one and relabel is applied.
 * These operations are repeated until the excess of the v is null or if there is a relabeling operation */
static void discharge(int v){
    bool relabeled = false;
    int w = 0;
    while((excess[v] != 0) && !relabeled){
        /* Check if push is applicable */
        if ((G_res[v][w] > 0) && (d[v] == (d[w] + 1)) && (d[w] < numNodes)){
            push(v,w);
            if (isActive(w) && (w!=(numVars*2 + 1)) && (find(B[d[w]].begin(), B[d[w]].end(), w) == B[d[w]].end()))
                B[d[w]].push_back(w);
        } else {
            /* Check if (v,w) is the last edge */
            if (w == (numVars*2 + 1)){
                relabel(v);
                relabeled = true;
                /* Change current edge */
            } else {
                w++;
            }
        }
    }
    if (isActive(v)){
        B[d[v]].push_back(v);
        /* update highest label */
        if (d[v] > b)
            b = d[v];
    }
}

/* Function for initial labeling
 * BFS used to find distances of nodes first from the sink and then from the source */
static void labelingBFS(){
    vector<bool> visited;
    list<int> queue;
    /* First BFS starting from the source */
    int node = numVars;
    visited.resize(numNodes, false);
    visited[node] = true;
    queue.push_back(node);
    /* don't want to relabel the sink so mark it as visited */
    visited[(numVars*2) + 1] = true;
    while (!queue.empty()) {
        node = queue.front();
        queue.pop_front();
        // scan all the nodes reached by an arc from current node
        for (int i = 0; i < numNodes; i++) {
            if ((G[node][i] != 0) && (!visited[i])) {
                visited[i] = true;
                // Assign same distance of the source
                d[i] = numNodes;
                queue.push_back(i);
            }
        }
    }
    /* BFS starting from the sink */
    node = (numVars*2) + 1;
    queue.push_back(node);
    fill(visited.begin(), visited.end(), false);
    visited[node] = true;
    /* don't want to relabel the source so mark it as visited */
    visited[numVars] = true;
    while (!queue.empty()){
        node = queue.front();
        queue.pop_front();
        // scan all the nodes reached by an arc from current node
        for(int i=0; i<numNodes; i++){
            // Search backwards so go to node i if there is an edge from i to node
            if (G[i][node] != 0){
                if (!visited[i]){
                    visited[i] = true;
                    // Assign distance of node + 1
                    d[i] = d[node] + 1;
                    queue.push_back(i);
                }
            }
        }
    }
}

/* Postordering reverse depth-first search that gives as a result the topological ordering of the network.
 * Reverse means that DFS is performed in the transpose graph, i.e. edges have opposite directions.
 * Postordering means that each node is inserted in the array of visited nodes when it is last visited.
 * A topological ordering can be obtained only if there are no cycles in the graph. So if a cycle is detected,
 * flow in each edge of the cycle is reduced by an amount equal to the flow value of the edge with minimum flow,
 * in order to break the cycle. Then the nodes in the cycle are set as unvisited and the DFS is resumed. */
static vector<int> reverseDFS(vector<vector<int> > &G_par, int node){
    /* Color identifies the state of each node:
     * WHITE - unvisited
     * GREY - still under processing
     * BLACK - node and all the nodes reachable from it visited. */
    enum COLOR {WHITE, GREY, BLACK};
    vector<int> parent(numNodes);
    vector<int> nodes_colors(numNodes, WHITE);
    vector<int> topological_order;
    bool edge_found;
    /* DFS ends when the starting node has been completely processed so all its edges have been examined */
    while (nodes_colors[node] != BLACK){
        edge_found = false;
        /* Find an edge with positive flow */
        for (int i=0; (i<numNodes) && (!edge_found); i++){
            if ((G_par[i][node] - G_res[i][node]) > 0){
                /* If WHITE node is found visit it */
                if (nodes_colors[i] == WHITE){
                    nodes_colors[i] = GREY;
                    parent[i] = node;
                    node = i;
                    edge_found = true;
                }
                    /* If GREY node is found a cycle has been found */
                else if (nodes_colors[i] == GREY){
                    //int min = INT_MAX;
                    int cycle_start = node;
                    /* Find the minimum flow value among the cycle's edges */
                    int min = G_par[i][node] - G_res[i][node];
                    while (node != i) {
                        int flow = G_par[node][parent[node]] - G_res[node][parent[node]];
                        if (flow < min)
                            min = flow;
                        /* Set the nodes of the cycle as unvisited */
                        nodes_colors[node] = WHITE;
                        node = parent[node];
                    }
                    /* Decrease the flow of all the edges by the minimum found */
                    G_res[i][node] += min;
                    G_res[node][i] -= min;
                    node = cycle_start;
                    while (node != i){
                        G_res[node][parent[node]] += min;
                        G_res[parent[node]][node] -= min;
                        node = parent[node];
                    }
                    edge_found = true;
                }
            }
        }
        /* If no edge with positive flow has been found node is completely processed so add it in the oredering */
        if (!edge_found){
            nodes_colors[node] = BLACK;
            topological_order.emplace_back(node);
            node = parent[node];
        }
    }
    return topological_order;
}

/* Function used to transform the preflow in a flow, so obtaining the residual network.
 * For each node with a positive excess the topological ordering is found using a postordering DFS. Then nodes
 * are processed in the oppposite of this order, flow is reduced in any edge entering this node until the excess
 * becomes zero. Therefore, the excess is moved from each node to nodes coming after in the topological ordering
 * and this*/
static void preflowToFlow(vector<vector<int> > &G_par){
    /* Repeat the process for each node with positive excess */
    /* Nodes are examined in an inverted order because nodes represented by higher numbers are typically further
     * from the source node, so in this way is possible that the excess of nodes closer to the source becomes zero
     * before examining them. */
    for (int i=numNodes-2; i>=0; i--){
        if (excess[i] > 0){
            /* DFS used to obtain topological ordering */
            vector<int> topological_queue = reverseDFS(G_par, i);
            /* Topological ordering followed in the opposite way */
            for (int q_index=(int)topological_queue.size()-1; q_index>=0; q_index--){
                int node = topological_queue[q_index];
                int remaining = excess[node];
                /* Find every edge with positive flow in order to decrease it if needed */
                for (int from_node=0; (from_node<numNodes) && (remaining!=0); from_node++){
                    int flow = G_par[from_node][node] - G_res[from_node][node];
                    if (flow > 0){
                        /* If flow higher than the excess of the node all the excess can be moved on this edge
                         * otherwise move all the excess possible, so make the flow equal to 0, and move the remaining
                         * on the next edge. */
                        if (flow >= remaining){
                            excess[from_node] += remaining;
                            excess[node] = 0;
                            G_res[from_node][node] += remaining;
                            G_res[node][from_node] -= remaining;
                            remaining = 0;
                        } else {
                            remaining -= flow;
                            excess[from_node] += flow;
                            excess[node] -= flow;
                            G_res[from_node][node] += flow;
                            G_res[node][from_node] = 0;
                        }
                    }
                }
            }
        }
    }
}

/* Function that calculates the max flow and the residual network of a flow network.
 * Parameters:
 *     G_par should be the adjacency matrix of the flow network.
 *     This matrix is transformed by the function in the adjacency matrix of the residual network.
 * Return value:
 *     Returns an integer that is the max flow value.
 *
 * The algorithm is subdivided into two parts. In the first one the excess and labeling functions are initialized.
 * The excess for each node is defined as the difference between ingoing and outgoing flow.
 * The labeling is an estimation of the distance of each node from the sink.
 * Then for each active node if possible flow is moved to nodes closer to sink, if not its labeling is updated
 * and if still flow can't be moved it becomes inactive. These operations are performed for each node following
 * the order given by the HL strategy, that is the node with the highest label is processed first.
 * At the end of this stage the max flow value is found, but the flow is still a preflow so a second stage is needed
 * to obtain the residual network.
 * The second part of the algorithm transforms the preflow in a flow. For each node with positive excess flow is
 * reduced until every node has zero excess. */
int maxFlow(vector<vector<int> >& G_par){
    G = G_par;
    /* Initialize residual network to be equal to input network */
    G_res = G;
    numNodes = (int)G.size();
    numVars = (int)G.size()/2-1;
    /* First step: find maximum preflow */
    /* Initialize HL strategy */
    B = vector<list<int> >(numNodes,list<int>());
    b = 0;
    /* Initialize distances */
    d = vector<int>(numNodes, 0);
    d[numVars] = numNodes;
    labelingBFS();
    /* Initialize excess of edges out of the source with the value of those edges' capacities */
    excess = vector<int>(numNodes, 0);
    for(int i=0; i<numVars*2+1; i++){
        excess[i] += G[numVars][i];
        excess[numVars] -= G[numVars][i];
        G_res[i][numVars] += G[numVars][i];
        G[numVars][i] = 0;
        G_res[numVars][i] = 0;
        if (isActive(i)) {
            B[d[i]].push_back(i);
            /* update highest label */
            if (d[i] > b)
                b = d[i];
        }
    }
    /* Process all nodes following HL strategy */
    while (b >= 0){
        if (B[b].empty()){
            b--;
        } else {
            int v = B[b].front();
            B[b].pop_front();
            discharge(v);
        }
    }
    /* Second step: transform preflow in flow */
    preflowToFlow(G_par);
    /* Give as result the residual network */
    G_par = G_res;
    return excess[numVars*2+1];
}