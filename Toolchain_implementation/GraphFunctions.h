//
// Created by W10 on 08/11/2022.
//

#ifndef QUBO_PREPROCESSING_GRAPHFUNCTIONS_H
#define QUBO_PREPROCESSING_GRAPHFUNCTIONS_H

#include <vector>
#include <list>

using namespace std;

vector<bool> BFS(vector<vector<int> > &G, int node, bool reverse);
vector<bool> DFS(vector<vector<int> > &G, int node);
list<vector<int> > stronglyConnectedComponents(vector<vector<int> > &G);
list<vector<int> > connectedComponents(vector<vector<int> > &G);

#endif //QUBO_PREPROCESSING_GRAPHFUNCTIONS_H
