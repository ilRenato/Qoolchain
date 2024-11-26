//
// Created by W10 on 13/11/2022.
//

#ifndef QUBO_PREPROCESSING_FLOWFUNCTIONS_H
#define QUBO_PREPROCESSING_FLOWFUNCTIONS_H

#include <vector>

using namespace std;

/****************************************************
 * NOTE: THE MAX FLOW ALGORITHM WORKS ONLY IF IN A FLOW NETWORK ON N NODES
 * THE SOURCE IS NODE N/2-1
 * THE SINK NODE IS N-1
 ***************************************************/
int maxFlow(vector<vector<int> >& G_par);

#endif //QUBO_PREPROCESSING_FLOWFUNCTIONS_H
