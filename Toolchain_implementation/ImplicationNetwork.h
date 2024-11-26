//
// Created by W10 on 29/10/2022.
//

#ifndef QUBO_PREPROCESSING_IMPLICATIONNETWORK_H
#define QUBO_PREPROCESSING_IMPLICATIONNETWORK_H

#include <iostream>
#include <vector>

using namespace std;

class ImplicationNetwork {
private:
    vector<vector<int> > _m_adj;
    int _numVariables;
    int _offset;
    float _max_flow;
    int complement(int variable) const;
    void makeResidualSymmetric(vector<vector<int> > &res);

public:
    ImplicationNetwork();
    ImplicationNetwork(vector<vector<int> > func, int offset);
    ImplicationNetwork(vector<vector<int>> &m_adj);
    ~ImplicationNetwork();

    vector<vector<int>> computeMaxFlow();
    int getOffset() const;
    float getMaxFlow() const;
    float getRoofDual() const;
    void print();

};


#endif //QUBO_PREPROCESSING_IMPLICATIONNETWORK_H
