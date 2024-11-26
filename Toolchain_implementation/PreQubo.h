//
// Created by W10 on 23/12/2022.
//

#ifndef QUBO_PREPROCESSING_PREQUBO_H
#define QUBO_PREPROCESSING_PREQUBO_H

#include <utility>
#include <vector>
#include <map>
#include <list>
#include <string>

using namespace std;

pair<vector<pair<vector<vector<int>>, list<int> > >, map<int, bool> > PreQubo(vector<vector<int> >& f, int &offset, string &filename);
int getOffset();
float getLowerBound();
vector<pair<vector<pair<vector<vector<int>>, list<int>>>, map<int, bool>>> ShannonDecomposition(vector<vector<int>> &Q, int &n_iterations, int &flag);

#endif //QUBO_PREPROCESSING_PREQUBO_H
