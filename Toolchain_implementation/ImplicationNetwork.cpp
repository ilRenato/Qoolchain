//
// Created by W10 on 29/10/2022.
//

#include "ImplicationNetwork.h"
#include "FlowFunctions.h"
#include <iostream>
#include <vector>

using namespace std;

ImplicationNetwork::ImplicationNetwork() {

}

ImplicationNetwork::ImplicationNetwork(vector<vector<int> > func, int offset) {
    _numVariables = (int)func.size();
    _offset = offset;
    _m_adj = vector<vector<int> >(2*(_numVariables+1), vector<int>(2*(_numVariables+1), 0));
    // Set quadratic terms of the implication network
    for(int i=0; i<_numVariables; i++){
        for(int j=i+1; j<_numVariables; j++){
            if (func[i][j] >= 0){
                _m_adj[i][complement(j)] = func[i][j];
                _m_adj[j][complement(i)] = func[i][j];
            }
            else{
                _m_adj[i][j] = -func[i][j];
                _m_adj[complement(j)][complement(i)] = -func[i][j];
                func[i][i] += func[i][j];
            }
        }
    }
    // Set linear terms of the implication network
    for(int i=0; i<_numVariables; i++){
        if (func[i][i] >= 0){
            _m_adj[_numVariables][complement(i)] = func[i][i];
            _m_adj[i][complement(_numVariables)] = func[i][i];
        }
        else {
            _m_adj[_numVariables][i] = -func[i][i];
            _m_adj[complement(i)][complement(_numVariables)] = -func[i][i];
            _offset += func[i][i];
        }
    }
}

ImplicationNetwork::ImplicationNetwork(vector<vector<int>> &m_adj){
    _m_adj = m_adj;
    _offset = 0;
    _numVariables = ((int) _m_adj.size()/2) - 1;
}

ImplicationNetwork::~ImplicationNetwork() {

}

vector<vector<int>> ImplicationNetwork::computeMaxFlow() {
    vector<vector<int> > residual(_m_adj);
    _max_flow = ((float) maxFlow(residual)) / 2;
    makeResidualSymmetric(residual);
    /* Remove ingoing edges at the source and outgoing edges at the sink because not needed */
    for (int i=0; i<((2*_numVariables)+1); i++){
        if (residual[i][_numVariables] > 0)
            residual[i][_numVariables] = 0;
    }
    for (int i=0; i<((2*_numVariables)+1); i++){
        if (residual[(2*_numVariables)+1][i] > 0)
            residual[(2*_numVariables)+1][i] = 0;
    }
    return residual;
}

int ImplicationNetwork::getOffset() const {
    return _offset;
}

float ImplicationNetwork::getMaxFlow() const {
    return _max_flow;
}

float ImplicationNetwork::getRoofDual() const{
    return _max_flow+(float)_offset;
}

void ImplicationNetwork::print() {
    for(int i=0; i<(2*(_numVariables+1)); i++){
        for(int j=0; j<(2*(_numVariables+1)); j++){
            cout << _m_adj[i][j] << " ";
        }
        cout << endl;
    }
}

int ImplicationNetwork::complement(int variable) const {
    if (variable <= _numVariables)
        return variable + _numVariables + 1;
    else
        return variable - _numVariables - 1;
}

void ImplicationNetwork::makeResidualSymmetric(vector<vector<int> > &res) {
    for(int i=0; i<=_numVariables; i++){
        for(int j=0; j<=_numVariables; j++){
            int sum = res[i][j] + res[complement(j)][complement(i)];
            res[i][j] = sum;
            res[complement(j)][complement(i)] = sum;
            if (j>i){
                sum = res[i][complement(j)] + res[j][complement(i)];
                res[i][complement(j)] = sum;
                res[j][complement(i)] = sum;
                sum = res[complement(i)][j] + res[complement(j)][i];
                res[complement(i)][j] = sum;
                res[complement(j)][i] = sum;
            }
        }
    }
}
