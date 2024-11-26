#ifndef QUBO_PREPROCESSING_RESIDUALNETWORK_H
#define QUBO_PREPROCESSING_RESIDUALNETWORK_H

#include "ImplicationNetwork.h"
#include <utility>
#include <vector>
#include <list>
#include <iostream>
#include <map>

using namespace std;

class ResidualNetwork {
private:
    vector<vector<int> > _m_adj;
    int _numVariables;
    list<int> _var_names;
    float _max_flow;
    int _offset;
    float _lower_bound;
    int _numStrongP, _numStrongPscc, _numWeakP, _numProbing;
    int _num_completescc;

    /* Receives a variable as parameter and returns its complement, that is the negated version for a normal
     * variable and the normal version for a negated variable. */
    int complement(int variable);

    /* Starting from a variable with a persistent assignment to 1, it searches all the variables reachable
     * from it and marks them as persistencies to 1. */
    map<int, bool> fixWithBFS(int start);

    /* Receives as parameter all the persistencies fixed and removes them from the adjacency matrix
     * of the residual network. It updates the number of variables and the list of variables names as well. */
    void removeFixedVars(map<int, bool> &fixed);

    /* Receives as parameter all the persistencies fixed and substitutes their values in the adjacency matrix
     * of the residual network. Coefficients are not just removed, new terms (also linear) may appear.
     * It updates the number of variables and the list of variables names as well. The return value is
     * the constant term produced by variables assignments. */
    int removeFixedVarsInSCC(map<int, bool> &fixed);

    /* Receives a variable as parameter and removes it from the adjacency matrix passed as parameter.
     * If assigning the variable to the value val produces a constant term it returns it. */
    int removeVar(int var, bool val, vector<vector<int>> &M);

    /* Heuristic method to find upper bound on the minimum of the function. */
    int DDT_onepass();

public:
    ResidualNetwork();
    ResidualNetwork(const vector<vector<int> >& network, float max_flow, int off);
    ResidualNetwork(const vector<vector<int> >& network, list<int> &v_names);
    ~ResidualNetwork();

    /* fixes to 1 (or to 0 if negated) all the literals that can be reached from the source */
    map<int, bool> fixStrongVariables();

    /* fixes to 1 or 0 all the literals that will be present in at least one solution */
    map<int, bool> fixAllVariables();

    /* divides the network in its non-connected components.
     * Since the different components contain disjoint sets of variables they can be solved independently */
    vector<ResidualNetwork> trivialDecomposition();

    /* extracts the strongly connected subgraphs and modifies this network removing those subgraphs.
     * Note that the sum of the minima of the posiforms associated to SCCs is the minimum of the cost function */
    vector<ResidualNetwork> SCCDecomposition();

    /* Selects a variable according to a method defined by flag and generates two subgraphs in which this variable has
     * been removed. This process is applied recursively n_divisions times. In this way smaller subgraphs with reduced
     * connectivity are obtained, and it is possible to apply other persistency and decomposition techniques
     * to these new subgraphs */
    vector<pair<ResidualNetwork, map<int, bool>>> ShannonDecomposition(int n_divisions, int flag);

    /* tries fixing some variables to find better bounds estimation and to fix more persistencies */
    map<int, bool> probing();
    map<int, bool> probing_v2();

    /* Returns the lower bound for the minimum of the cost function associated with the residual network */
    float getLowerBound() const;

    /* Returns the residual network in the form of a cost function */
    vector<vector<int> > getCostFunction(int &offset);

    list<int> getVarNames() const;

    /* Function used to collect statistics about the operations performed
     *   - strongP contains the number of strong persistencies found exploring the graph from the source
     *   - strongPscc contains the number of strong persistencies found by means of strongly connected components
     *   - weakP contains the number of weak persistencies found by means of strongly connected components
     *   - complete_scc contains the number of complete strongly connected components found.
     *     With complete is meant components in which is present both a variable and its negated form */
    void getStats(int &strongP, int &strongPscc, int &weakP, int &numProbing, int &complete_scc) const;

    void print();
    /*************************/
    vector<vector<int> > getResidual();
    /*************************/

};


#endif //QUBO_PREPROCESSING_RESIDUALNETWORK_H
