#include "ResidualNetwork.h"
#include "ImplicationNetwork.h"
#include "GraphFunctions.h"
#include <numeric>
#include <algorithm>
#include <utility>
#include <vector>
#include <list>
#include <iostream>
#include <map>
#include <unordered_set>

using namespace std;

ResidualNetwork::ResidualNetwork() {

}

ResidualNetwork::ResidualNetwork(const vector<vector<int> >& network, float max_flow, int off) {
    _m_adj = network;
    _numVariables = (int) network.size()/2-1;
    _max_flow = max_flow;
    _offset = off;
    _lower_bound = max_flow + (float)off;
    _var_names = list<int>(_numVariables);
    iota(_var_names.begin(), _var_names.end(), 0);
    _numStrongP = _numStrongPscc = _numWeakP = 0;
    _num_completescc = 0;
}

ResidualNetwork::ResidualNetwork(const vector<vector<int> >& network, list<int> &v_names) {
    _m_adj = network;
    _numVariables = (int) network.size()/2-1;
    _var_names = v_names;
    _numStrongP = _numStrongPscc = _numWeakP = 0;
    _num_completescc = 0;
}

ResidualNetwork::~ResidualNetwork() {

}

/* Fixes all the variables reached from a starting node passed as parameter. */
map<int, bool> ResidualNetwork::fixWithBFS(int start) {
    map<int, bool> fixed;
    vector<bool> visited;
    /* find nodes connected to the source */
    visited = BFS(_m_adj, start, false);
    auto it = _var_names.begin();
    for(int i=0; i<2*_numVariables+1; i++){
        if (i == _numVariables){
            it = _var_names.begin();
        }
        else {
            if (visited[i] && i != (_numVariables * 2) + 1) {
                if (i < _numVariables) {
                    fixed.insert({*it, true});
                } else {
                    fixed.insert({*it, false});
                }
            }
            it++;
        }
    }
    return fixed;
}

/* Removes from the list of variables (_var_names) and from the residual network matrix (_m_adj)
 * the variables that have been fixed. */
void ResidualNetwork::removeFixedVars(map<int, bool> &fixed) {
    /* Search for persistencies in the list of variables and remove these variables */
    vector<int> old_indexes;
    int index = 0;
    auto it = _var_names.begin();
    while (it != _var_names.end()){
        if (fixed.find(*it) != fixed.end())
            it = _var_names.erase(it);
        else {
            it++;
            old_indexes.emplace_back(index);
        }
        index++;
    }
    /* update number of variables */
    int old_source = _numVariables;
    old_indexes.emplace_back(old_source);
    int old_sink = (_numVariables*2)+1;
    int new_num_var = (int) _var_names.size();
    /* construct new matrix removing the variables that have been fixed */
    vector<vector<int>> new_mat((new_num_var+1)*2, vector<int>((new_num_var+1)*2, 0));
    for (int i=0; i<(new_num_var+1)*2; i++){
        for (int j=0; j<(new_num_var+1)*2; j++) {
            if ((j <= new_num_var) && (i <= new_num_var))
                new_mat[i][j] = _m_adj[old_indexes[i]][old_indexes[j]];
            else if ((j > new_num_var) && (i <= new_num_var))
                new_mat[i][j] = _m_adj[old_indexes[i]][complement(old_indexes[j-new_num_var-1])];
            else if ((j <= new_num_var) && (i > new_num_var))
                new_mat[i][j] = _m_adj[complement(old_indexes[i-new_num_var-1])][old_indexes[j]];
            else if ((j > new_num_var) && (i > new_num_var))
                new_mat[i][j] = _m_adj[complement(old_indexes[i-new_num_var-1])][complement(old_indexes[j-new_num_var-1])];
        }
    }
    _numVariables = new_num_var;
    _m_adj = new_mat;
}

int ResidualNetwork::removeFixedVarsInSCC(map<int, bool> &fixed) {
    int offset = 0;
    list<int> old_var_list = _var_names;
    /* Remove fixed variables form the list of variables */
    auto it = _var_names.begin();
    list<int> old_indexes; // contains indexes of the variables that are not fixed
    int i = 0;
    while (it != _var_names.end()){
        if (fixed.find(*it) != fixed.end())
            it = _var_names.erase(it);
        else {
            it++;
            old_indexes.emplace_back(i);
        }
        i++;
    }
    old_indexes.emplace_back(_numVariables);
    int new_num_var = (int) _var_names.size();
    /* construct new matrix removing the variables that have been fixed */
    vector<vector<int>> new_mat((new_num_var+1)*2, vector<int>((new_num_var+1)*2, 0));
    auto it_i = old_indexes.begin();
    auto it_persi = fixed.begin();
    i = 0;
    for (int old_i=0; old_i<_numVariables+1; old_i++){
        auto it_j = old_indexes.begin();
        advance(it_j, i);
        auto it_persj = fixed.begin();
        if (distance(it_persj, fixed.end()) > old_i-i) // Boundary check
            advance(it_persj, old_i-i);
        int j = i;
        for (int old_j=old_i; old_j<_numVariables+1; old_j++){
            /* If old_i and old_j are part of old_indexes they are part of the new matrix */
            if ((old_i == *it_i) && (old_j == *it_j)){
                new_mat[i][j] += _m_adj[old_i][old_j];
                new_mat[i][j+new_num_var+1] += _m_adj[old_i][complement(old_j)];
                new_mat[i+new_num_var+1][j] += _m_adj[complement(old_i)][old_j];
                new_mat[i+new_num_var+1][j+new_num_var+1] += _m_adj[complement(old_i)][complement(old_j)];
                new_mat[j][i] += _m_adj[old_j][old_i];
                new_mat[j][i+new_num_var+1] += _m_adj[old_j][complement(old_i)];
                new_mat[j+new_num_var+1][i] += _m_adj[complement(old_j)][old_i];
                new_mat[j+new_num_var+1][i+new_num_var+1] += _m_adj[complement(old_j)][complement(old_i)];
            }
            /* If old_i and old_j are not part of old_indexes they generate an offset */
            else if ((old_i != *it_i) && (old_j != *it_j)){
                if ((_m_adj[old_i][old_j] != 0) && (it_persi->second) && (!it_persj->second)){
                    offset += _m_adj[old_i][old_j];
                }
                else if ((_m_adj[old_i][complement(old_j)] != 0) && (it_persi->second) && (it_persj->second)){
                    offset += _m_adj[old_i][complement(old_j)];
                }
                else if ((_m_adj[complement(old_i)][old_j] != 0) && (!it_persi->second) && (!it_persj->second)){
                    offset += _m_adj[complement(old_i)][old_j];
                }
                else if ((_m_adj[complement(old_i)][complement(old_j)] != 0) && (!it_persi->second) && (it_persj->second)){
                    offset += _m_adj[old_i][old_j];
                }
            }
            /* If only old_i is part of the new matrix linear terms are generated */
            else if (old_i == *it_i){
                if (it_persj->second){
                    if (_m_adj[old_i][complement(old_j)] != 0){
                        new_mat[new_num_var][i+new_num_var+1] += _m_adj[old_i][complement(old_j)];
                        new_mat[i][2*new_num_var+1] += _m_adj[old_i][complement(old_j)];
                    }
                    if (_m_adj[complement(old_i)][complement(old_j)] != 0){
                        new_mat[new_num_var][i] += _m_adj[complement(old_i)][complement(old_j)];
                        new_mat[i+new_num_var+1][2*new_num_var+1] += _m_adj[complement(old_i)][complement(old_j)];
                    }
                }
                else {
                    if (_m_adj[old_i][old_j] != 0){
                        new_mat[new_num_var][i+new_num_var+1] += _m_adj[old_i][old_j];
                        new_mat[i][2*new_num_var+1] += _m_adj[old_i][old_j];
                    }
                    if (_m_adj[complement(old_i)][old_j] != 0){
                        new_mat[new_num_var][i] += _m_adj[complement(old_i)][old_j];
                        new_mat[i+new_num_var+1][2*new_num_var+1] += _m_adj[complement(old_i)][old_j];
                    }
                }
            }
            /* If only old_j is part of the new matrix linear terms are generated */
            else if (old_j == *it_j){
                if (it_persi->second){
                    if (_m_adj[old_i][old_j] != 0){
                        new_mat[new_num_var][j] += _m_adj[old_i][old_j];
                        new_mat[j+new_num_var+1][2*new_num_var+1] += _m_adj[old_i][old_j];
                    }
                    if (_m_adj[old_i][complement(old_j)] != 0){
                        new_mat[new_num_var][j+new_num_var+1] += _m_adj[old_i][complement(old_j)];
                        new_mat[j][2*new_num_var+1] += _m_adj[old_i][complement(old_j)];
                    }
                }
                else {
                    if (_m_adj[complement(old_i)][old_j] != 0){
                        new_mat[new_num_var][j] += _m_adj[complement(old_i)][old_j];
                        new_mat[j+new_num_var+1][2*new_num_var+1] += _m_adj[complement(old_i)][old_j];
                    }
                    if (_m_adj[complement(old_i)][complement(old_j)] != 0){
                        new_mat[new_num_var][j+new_num_var+1] += _m_adj[complement(old_i)][complement(old_j)];
                        new_mat[j][2*new_num_var+1] += _m_adj[complement(old_i)][complement(old_j)];
                    }
                }
            }
            /* If old_j is part of the not fixed variables, then a new item in the new matrix has been created,
             * thus j and it_j can be increamented. if old_j is part of the fixed variables it_persj has been used,
             * and it is possible to pass onto the next fixed variables, so it can be incremented. */
            if (old_j == *it_j){
                it_j++;
                j++;
            }
            else{
                if (next(it_persj) != fixed.end()) // check that we are not at the end of the map
                    it_persj++;
            }
        }
        /* If old_i is part of the not fixed variables, then a new item in the new matrix has been created,
         * thus i and it_i can be increamented. if old_i is part of the fixed variables it_persi has been used,
         * and it is possible to pass onto the next fixed variables, so it can be incremented. */
        if (old_i == *it_i){
            it_i++;
            i++;
        }
        else{
            if (next(it_persi) != fixed.end()) // check that we are not at the end of the map
                it_persi++;
        }
    }
    _numVariables = new_num_var;
    _m_adj = new_mat;
    return offset;
}

/* Finds all literals that can be reached from the source with a BFS using the source as starting point.
 * All the nodes reached are fixed to true if not complemented, or fixed to false if complemented.
 * Then the residual network is updated removing the fixed variables. */
map<int, bool> ResidualNetwork::fixStrongVariables() {
    map<int, bool> fixed;
    /* Execute BFS starting from the source */
    fixed = fixWithBFS(_numVariables);
    /* Update variables and matrix */
    removeFixedVars(fixed);
    _numStrongP = (int) fixed.size();
    return fixed;
}

/* Finds all the strongly connected components and then fixes all the possible variables.
 * For each component if it is a strong persistency all the variables are fixed and then with a BFS
 * all the variables that can be reached from this component are fixed to true. If the component
 * is a weak persistency, either it has been already fixed by a previous BFS or it is still not fixed.
 * In the latter case the component either can be reached by another strongly connected component and
 * so it will be fixed by a BFS, or it is not preceded by any persistency and so it can be fixed
 * to true or false independently of other variables. Then the residual network is updated removing
 * the fixed variables. */
map<int, bool> ResidualNetwork::fixAllVariables() {
    list<vector<int> > scc = stronglyConnectedComponents(_m_adj);
    map<int, bool> fixed;
    vector<vector<int> > to_skip;
    for (auto & component : scc){
        /* if this component is a negated version of a previous processed one it has to be skipped */
        sort(component.begin(), component.end());
        if (find(to_skip.begin(), to_skip.end(), component) == to_skip.end()) {
            /* Check if the component has a twin with negated variables or if it contains both
             * variables and negated ones */
            if (find(component.begin(), component.end(), complement(component[0])) == component.end()){
                /* component with negated variables has to be skipped */
                vector<int> negated_twin;
                for (auto & element : component){
                    /* Store the negated version of the component to detect and skip it in the future */
                    negated_twin.emplace_back(complement(element));
                }
                /* Add the negated version to the components that will be skipped */
                sort(negated_twin.begin(), negated_twin.end());
                to_skip.emplace_back(negated_twin);
                /* Check if variables in component have already been fixed */
                if (fixed.find(component[0]) == fixed.end()){
                    /* check whether persistency is strong or weak, so if there is an edge between component
                     * and negated one */
                    bool strong = false;
                    bool value = true;
                    for (int i=0; (i<component.size()) && (!strong); i++){
                        int element = component[i];
                        for (int j=0; (j<component.size()) && (!strong); j++){
                            int c_element = complement(component[j]);
                            if (_m_adj[element][c_element] != 0) {
                                strong = true;
                                /* true vertices can not lead to false ones so assign them to false */
                                value = false;
                            }
                            if (_m_adj[c_element][element] != 0) {
                                strong = true;
                                /* false vertices can not be preceded by true vertices so assign them to true */
                                value = true;
                            }
                        }
                    }
                    /* If persistencies are strong fix them */
                    if (strong){
                        _numStrongPscc += (int) component.size();
                        map<int, bool> p_to_fix;
                        if (value)
                            p_to_fix = fixWithBFS(component[0]);
                        else
                            p_to_fix = fixWithBFS(complement(component[0]));
                        fixed.insert(p_to_fix.begin(), p_to_fix.end());
                    }
                    /* If persistencies are weak fix them only if not preceded by other persistencies */
                    else {
                        _numWeakP += (int) component.size();
                        /* with a reverse BFS look for other strongly connected components that may reach the
                         * current component. If there is at least one strongly connected component C that reaches
                         * it then the persistencies in the current component will be fixed by the BFS of C. */
                        vector<bool> reverse_visited = BFS(_m_adj, component[0], true);
                        bool preceded = false;
                        auto it_scc = scc.begin();
                        while ((it_scc != scc.end()) && (!preceded)){
                            if ((*it_scc != component) && (*it_scc != negated_twin)){
                                for (int i=0; (i<it_scc->size()) && (!preceded); i++){
                                    if (reverse_visited[(*it_scc)[i]])
                                        preceded = true;
                                }
                            }
                            it_scc++;
                        }
                        /* If this component is not preceded by any other component than both the assignment
                         * 1 and 0 are valid for the variables in the component. */
                        map<int, bool> p_to_fix;
                        if (!preceded)
                            p_to_fix = fixWithBFS(component[0]);
                        fixed.insert(p_to_fix.begin(), p_to_fix.end());
                    }
                }
            }
            else {
                _num_completescc++;
            }
        }
    }
    /* Update variables and matrix */
    removeFixedVars(fixed);
    return fixed;
}

vector<ResidualNetwork> ResidualNetwork::trivialDecomposition(){
    vector<ResidualNetwork> subnetworks;
    list<vector<int> > CC;
    vector<int> first_vertex;
    /* Remove edges from source and to sink */
    vector<vector<int> > g = _m_adj;
    for (int i=0; i<(2*_numVariables)+1; i++){
        g[_numVariables][i] = 0;
        g[i][(_numVariables*2)+1] = 0;
    }
    /* find all the subgraph not connected among each other */
    CC = connectedComponents(g);
    /* for each subgraph construct a new residual network */
    for (auto & component : CC){
        bool skip = false;
        int new_size = (int) component.size();
        /* Check if component is not only made by source or sink and skip it in this case */
        if ((component[0] == _numVariables) || (component[0] == (2*_numVariables)+1))
            skip = true;
        /* Check if current component is the negated counterpart of a previously processed one */
        for (auto & v : first_vertex){
            /* If an already processed vertex is present skip this component */
            if (find(component.begin(), component.end(), v) != component.end())
                skip = true;
        }
        if (!skip) {
            /* Check if component is symmetric, so if it has both a variable and its complemented counterpart */
            if (find(component.begin(), component.end(), complement(component[0])) == component.end()) {
                first_vertex.emplace_back(complement(component[0]));
                /* For each vertex add to the component the negated vertex */
                for (int i=0; i<new_size; i++) {
                    component.emplace_back(complement(component[i]));
                }
            }
            /* include source */
            component.emplace_back(_numVariables);
            /* include sink */
            component.emplace_back((_numVariables * 2) + 1);
            new_size = (int) component.size();
            sort(component.begin(), component.end());
            vector<vector<int> > subnetwork(new_size, vector<int>(new_size));
            /* construct new residual network adjacency matrix */
            for (int i=0; i<new_size; i++) {
                for (int j=0; j<new_size; j++) {
                    subnetwork[i][j] = _m_adj[component[i]][component[j]];
                }
            }
            /* create residual network */
            list<int> v_names;
            /* Iterate half of components nodes because of symmetry */
            for (int i=0; i<(component.size()/2)-1; i++){
                auto l = _var_names.begin();
                advance(l, component[i]);
                v_names.emplace_back(*l);
            }
            subnetworks.emplace_back(ResidualNetwork(subnetwork, v_names));
        }
    }
    return subnetworks;
}

vector<ResidualNetwork> ResidualNetwork::SCCDecomposition(){
    vector<ResidualNetwork> subnetworks;
    list<vector<int> > scc = stronglyConnectedComponents(_m_adj);
    list<int> no_scc_elements((_numVariables+1)*2);
    /* Initialize nodes list of subgraph with no SCC with all the nodes */
    iota(no_scc_elements.begin(), no_scc_elements.end(), 0);
    auto it = scc.begin();
    /* Remove all strongly connected components that don't have both a variable and its complement */
    while (it != scc.end()){
        /* Check if complement of first variable is not in the same strongly connected component */
        if (find(it->begin(), it->end(), complement((*it)[0])) == it->end())
            it = scc.erase(it);
        else
            it++;
    }
    /* Extract each strongly connected component and create a new residual network for it */
    for (auto & component : scc){
        /* Add sink and source to each subgraph */
        component.emplace_back(_numVariables);
        component.emplace_back((_numVariables*2)+1);
        int new_size = (int) component.size();
        vector<vector<int> > subnetwork(new_size, vector<int>(new_size));
        sort(component.begin(), component.end());
        /* Extract subgraph from original network */
        for (int i=0; i<new_size; i++){
            for (int j=0; j<new_size; j++){
                subnetwork[i][j] = _m_adj[component[i]][component[j]];
            }
            /* If there is no edge between a variable in an SCC and a variable outside
             * then the former variable is not part of the remaining network without SCCs */
            bool arc_found = false;
            if ((component[i] != _numVariables) && (component[i] != (2*_numVariables)+1)) {
                auto e_it = component.begin();
                for (int j = 0; (j < (2 * _numVariables) + 1) && !arc_found; j++) {
                    /* do not consider edges with other nodes of the same SCC,
                     * the presence of an edge is obvious in that case */
                    if (*e_it == j) {
                        e_it++;
                    } else {
                        if ((_m_adj[component[i]][j] != 0) || (_m_adj[j][component[i]]))
                            arc_found = true;
                    }
                }
                if (!arc_found)
                    no_scc_elements.remove(component[i]);
            }
        }
        /* create residual network */
        list<int> v_names;
        /* Iterate half of components nodes because of symmetry */
        for (int i=0; i<(component.size()/2)-1; i++){
            auto l = _var_names.begin();
            advance(l, component[i]);
            v_names.emplace_back(*l);
        }
        subnetworks.emplace_back(ResidualNetwork(subnetwork, v_names));
    }
    /* Create residual network without extracted strongly connected components */
    int new_size = (int) no_scc_elements.size();
    vector<vector<int> > new_madj(new_size, vector<int>(new_size));
    auto it1 = no_scc_elements.begin();
    for (int i=0; i<new_size; i++) {
        auto it2 = no_scc_elements.begin();
        for (int j = 0; j < new_size; j++) {
            new_madj[i][j] = _m_adj[*it1][*it2];
            it2++;
        }
        it1++;
    }
    /* Update residual matrix and number of variables */
    _m_adj = new_madj;
    list<int> v_names;
    /* Iterate half of components nodes because of symmetry */
    auto it_var = no_scc_elements.begin();
    for (int i=0; i<(no_scc_elements.size()/2)-1; i++){
        auto l = _var_names.begin();
        advance(l, *it_var);
        v_names.emplace_back(*l);
        it_var++;
    }
    _var_names = v_names;
    _numVariables = (new_size/2)-1;
    return subnetworks;
}

vector<pair<ResidualNetwork, map<int, bool>>> ResidualNetwork::ShannonDecomposition(int n_divisions, int flag){
    vector<pair<ResidualNetwork, map<int, bool>>> sub_Q;
    int selected_var = 0;
    if (n_divisions > 0){
        /* Simplest criteria of selection:
         * The node with the largest number of connections (degree) is selected, therefore all the off-diagonal terms
         * are checked computing the degree of each node. */
        if (flag == 0) {
            vector<int> degrees(_numVariables, 0);
            int max_degree = 0;
            for (int i = 0; i < _numVariables + 1; i++) {
                for (int j = i + 1; j < _numVariables + 1; j++) {
                    if (_m_adj[i][j] != 0) {
                        if (i != _numVariables){
                            degrees[i]++;
                            if (degrees[i] > max_degree) {
                                max_degree = degrees[i];
                                selected_var = i;
                            }
                        }
                        if (j != _numVariables){
                            degrees[j]++;
                            if (degrees[j] > max_degree) {
                                max_degree = degrees[j];
                                selected_var = j;
                            }
                        }
                    }
                    if (_m_adj[i][complement(j)] != 0) {
                        if (i != _numVariables){
                            degrees[i]++;
                            if (degrees[i] > max_degree) {
                                max_degree = degrees[i];
                                selected_var = i;
                            }
                        }
                        if (j != _numVariables){
                            degrees[j]++;
                            if (degrees[j] > max_degree) {
                                max_degree = degrees[j];
                                selected_var = j;
                            }
                        }
                    }
                    if (_m_adj[complement(j)][i] != 0) {
                        if (i != _numVariables){
                            degrees[i]++;
                            if (degrees[i] > max_degree) {
                                max_degree = degrees[i];
                                selected_var = i;
                            }
                        }
                        if (j != _numVariables){
                            degrees[j]++;
                            if (degrees[j] > max_degree) {
                                max_degree = degrees[j];
                                selected_var = j;
                            }
                        }
                    }
                    if (_m_adj[j][i] != 0) {
                        if (i != _numVariables){
                            degrees[i]++;
                            if (degrees[i] > max_degree) {
                                max_degree = degrees[i];
                                selected_var = i;
                            }
                        }
                        if (j != _numVariables){
                            degrees[j]++;
                            if (degrees[j] > max_degree) {
                                max_degree = degrees[j];
                                selected_var = j;
                            }
                        }
                    }
                }
            }
        }
        /* Store the name of the selected variable to successively save it in the map to return */
        list<int> tmp_var_names = _var_names;
        auto it = tmp_var_names.begin();
        advance(it, selected_var);
        int selected_var_name = *it;
        /* Map to substitute first the value 0 and then 1 to the selected variable */
        map<int, bool> var_to_del = {{selected_var_name, false}};

        /* Create the first subgraph by putting the selected variable equal to 0 */
        vector<vector<int>> tmp_mat = _m_adj;
        removeFixedVarsInSCC(var_to_del);
        ImplicationNetwork tmp_net(_m_adj);
        ResidualNetwork Q0(tmp_net.computeMaxFlow(), _var_names);

        /* Create the second subgraph by putting the selected variable equal to 1 */
        /* Restore the residual network as it was before the creation of Q0 */
        _m_adj = tmp_mat;
        _var_names = tmp_var_names;
        _numVariables++;
        var_to_del = {{selected_var_name, true}};
        removeFixedVarsInSCC(var_to_del);
        tmp_net = ImplicationNetwork(_m_adj);
        ResidualNetwork Q1(tmp_net.computeMaxFlow(), _var_names);

        /* Apply recursion on Q0 */
        sub_Q = Q0.ShannonDecomposition(n_divisions-1, flag);
        /* Add the selected variable as a persistency */
        for (auto & item : sub_Q){
            item.second.insert({selected_var_name, 0});
        }
        /* Apply recursion on Q1 */
        vector<pair<ResidualNetwork, map<int, bool>>> tmp = Q1.ShannonDecomposition(n_divisions-1, flag);
        /* Add the selected variable as a persistency */
        for (auto & item : tmp){
            item.second.insert({selected_var_name, 1});
        }
        sub_Q.insert(sub_Q.end(), tmp.begin(), tmp.end());
    }
    else {
        sub_Q.emplace_back(*this, map<int, bool>());
    }
    return sub_Q;
}

int ResidualNetwork::complement(int variable) {
    if (variable <= _numVariables)
        return variable + _numVariables + 1;
    else
        return variable - _numVariables - 1;
}

int ResidualNetwork::removeVar(int var, bool val, vector<vector<int>> &M){
    int const_term = 0;
    /* Check when the variable to be substituted with val appears in the quadratic terms of the matrix
     * and update all the coefficients. When the variable is true, and it appears in a quadratic term, then
     * the other variable forms a linear term with the same coefficient. If it is false the term disappears.
     * When the variable is false and its complement appears in a quadratic term, then the other variable
     * forms a linear term with the same coefficient. If it is true the term disappears. */
    for (int i=0; i<_numVariables-1; i++){
        for (int j=i+1; j<_numVariables; j++){
            if (i == var){
                if (val){
                    M[_numVariables][j] += M[i][j];
                    M[complement(j)][(2*_numVariables)+1] += M[i][j];
                    M[i][j] = 0;
                    M[complement(j)][complement(i)] = 0;
                    M[_numVariables][complement(j)] += M[i][complement(j)];
                    M[j][(2*_numVariables)+1] += M[i][complement(j)];
                    M[i][complement(j)] = 0;
                    M[j][complement(i)] = 0;
                    M[complement(i)][j] = 0;
                    M[complement(j)][i] = 0;
                    M[complement(i)][complement(j)] = 0;
                    M[j][i] = 0;
                } else {
                    M[i][j] = 0;
                    M[complement(j)][complement(i)] = 0;
                    M[i][complement(j)] = 0;
                    M[j][complement(i)] = 0;
                    M[_numVariables][j] += M[complement(i)][j];
                    M[complement(j)][(2*_numVariables)+1] += M[complement(i)][j];
                    M[complement(i)][j] = 0;
                    M[complement(j)][i] = 0;
                    M[_numVariables][complement(j)] += M[complement(i)][complement(j)];
                    M[j][(2*_numVariables)+1] += M[complement(i)][complement(j)];
                    M[complement(i)][complement(j)] = 0;
                    M[j][i] = 0;
                }
            } else if (j == var){
                if (val){
                    M[_numVariables][complement(i)] += M[i][complement(j)];
                    M[i][(2*_numVariables)+1] += M[i][complement(j)];
                    M[i][complement(j)] = 0;
                    M[j][complement(i)] = 0;
                    M[_numVariables][i] += M[complement(i)][complement(j)];
                    M[complement(i)][(2*_numVariables)+1] += M[complement(i)][complement(j)];
                    M[complement(i)][complement(j)] = 0;
                    M[j][i] = 0;
                    M[i][j] = 0;
                    M[complement(j)][complement(i)] = 0;
                    M[complement(i)][j] = 0;
                    M[complement(j)][i] = 0;
                } else {
                    M[_numVariables][complement(i)] += M[i][j];
                    M[i][(2*_numVariables)+1] += M[i][j];
                    M[i][j] = 0;
                    M[complement(j)][complement(i)] = 0;
                    M[_numVariables][i] += M[complement(i)][j];
                    M[complement(i)][(2*_numVariables)+1] += M[complement(i)][j];
                    M[complement(i)][j] = 0;
                    M[complement(j)][i] = 0;
                    M[complement(i)][complement(j)] = 0;
                    M[j][i] = 0;
                    M[i][complement(j)] = 0;
                    M[j][complement(i)] = 0;
                }
            }
        }
    }
    /* Check when the variable to be substituted with val appears in the linear terms of the matrix
     * and update all the coefficients. If the variable is true, and it appears in a linear term, then
     * a constant value equal to its coefficient is added. If it is false the linear term disappears.
     * if the variable is false and its complement appears in a linear term, then a constant value
     * equal to its coefficient is added. If it is true the term disappears. */
    if (val){
        const_term += M[_numVariables][complement(var)];
        M[_numVariables][complement(var)] = 0;
        M[var][complement(_numVariables)] = 0;
        M[_numVariables][var] = 0;
        M[complement(var)][complement(_numVariables)] = 0;
    }
    else {
        const_term += M[_numVariables][var];
        M[_numVariables][var] = 0;
        M[complement(var)][complement(_numVariables)] = 0;
        M[_numVariables][complement(var)] = 0;
        M[var][complement(_numVariables)] = 0;
    }
    return const_term;
}

int ResidualNetwork::DDT_onepass() {
    vector<vector<int>> M = _m_adj;
    int result = 0;
    /* Assign the value of one variable for each iteration */
    for (int iterations=0; iterations<_numVariables; iterations++){
        int max = 0;
        unordered_set<int> linear_terms;
        pair<int, int> couple(0, 0);
        /* Find the quadratic term with the largest coefficient */
        for (int i=0; i<_numVariables-1; i++){
            for (int j=i+1; j<_numVariables; j++){
                if (M[i][j] > max){
                    max = M[i][j];
                    couple = {i, complement(j)};
                }
                if (M[i][complement(j)] > max){
                    max = M[i][complement(j)];
                    couple = {i, j};
                }
                if (M[complement(i)][j] > max){
                    max = M[complement(i)][j];
                    couple = {complement(i), complement(j)};
                }
                if (M[complement(i)][complement(j)] > max){
                    max = M[complement(i)][complement(j)];
                    couple = {complement(i), j};
                }
            }
        }
        /* Find the linear term with the largest coefficient */
        for (int j=0; j<_numVariables; j++){
            if (M[_numVariables][j] > 0){
                linear_terms.insert(complement(j));
                if (M[_numVariables][j] > max){
                    max = M[_numVariables][j];
                    couple = {_numVariables, complement(j)};
                }
            }
            if (M[_numVariables][complement(j)] > 0){
                linear_terms.insert(j);
                if (M[_numVariables][complement(j)] > max){
                    max = M[_numVariables][complement(j)];
                    couple = {_numVariables, j};
                }
            }
        }
        /* If no term has been found the algorithm is completed */
        if (couple == pair<int, int>(0, 0)){
            break;
        }
        /* If the term with the largest coefficient is a linear one assign the variable to 0 */
        if (couple.first == _numVariables){
            if (couple.second < _numVariables)
                result += removeVar(couple.second, false, M);
            else
                result += removeVar(complement(couple.second), true, M);
        }
        /* If the term with the largest coefficient is a quadratic one find the quadratic term with the largest
         * coefficient that is formed by a variable that is also present in a linear term. */
        else {
            int var;
            bool flag = false;
            max = 0;
            for (int u=0; u<_numVariables; u++){
                for (int v=u+1; v<_numVariables; v++){
                    /* If the literal u is in a linear term search the max coefficient in its quadratic terms */
                    if (linear_terms.find(u) != linear_terms.end()) {
                        /* If u is a normal variable it appears complemented in the posiform, so if
                         * the linear term has to vanish u = 1, thus in a quadratic term uv -> v=0 */
                        if (M[complement(u)][v] > max) {
                            flag = true;
                            max = M[complement(u)][v];
                            var = complement(v);
                        } else if (M[complement(u)][complement(v)] > max) {
                            flag = true;
                            max = M[complement(u)][complement(v)];
                            var = v;
                        }
                    } else if (linear_terms.find(complement(u)) != linear_terms.end()) {
                        /* If u is a complemented variable it appears as a normal variable in the posiform, so if
                         * the linear term has to vanish u = 0, thus in a quadratic term (!u)v -> v=0 */
                        if (M[u][v] > max) {
                            flag = true;
                            max = M[u][v];
                            var = complement(v);
                        } else if (M[u][complement(v)] > max) {
                            flag = true;
                            max = M[u][complement(v)];
                            var = v;
                        }
                    }
                    /* If the literal v is in a linear term search the max coefficient in its quadratic terms */
                    else if (linear_terms.find(v) != linear_terms.end()) {
                        /* If v is a normal variable it appears complemented in the posiform, so if
                         * the linear term has to vanish v = 1, thus in a quadratic term uv -> u=0 */
                        if (M[u][v] > max){
                            flag = true;
                            max = M[u][v];
                            var = u;
                        } else if (M[complement(u)][v] > max){
                            flag = true;
                            max = M[complement(u)][v];
                            var = complement(u);
                        }
                    } else if (linear_terms.find(complement(v)) != linear_terms.end()) {
                        /* If v is a complemented variable it appears as a normal variable in the posiform, so if
                         * the linear term has to vanish v = 0, thus in a quadratic term u(!v) -> u=0 */
                        if (M[u][complement(v)] > max){
                            flag = true;
                            max = M[u][complement(v)];
                            var = u;
                        } else if (M[u][complement(v)] > max){
                            flag = true;
                            max = M[u][complement(v)];
                            var = complement(u);
                        }
                    }
                }
            }
            /* If a quadratic term composed of a variable that is also present but complemented in a linear term
             * has been found, then the other variable in the quadratic term is set to 0 */
            if (flag){
                if (var < _numVariables)
                    result += removeVar(var, false, M);
                else
                    result += removeVar(complement(var), true, M);
            }
            /* If no such quadratic term has been found, the first variable of the term with the largest coefficient
             * is set to 0 */
            else {
                if (couple.first < _numVariables)
                    result += removeVar(couple.first, false, M);
                else
                    result += removeVar(complement(couple.first), true, M);
            }
        }
    }
    return result;
}

map<int, bool> ResidualNetwork::probing() {
    map<int, bool> new_pers;
    int U = DDT_onepass() + 1;
    ImplicationNetwork net;
    map<int, bool> p_to_add;
    for (int i=0; i<_numVariables; i++){
        /* Add the penalty term for the complement of variable i */
        _m_adj[_numVariables][i] += U+1;
        _m_adj[complement(i)][complement(_numVariables)] += U+1;
        net = ImplicationNetwork(_m_adj);
        vector<vector<int>> res_net = net.computeMaxFlow();
        float new_flow1 = net.getMaxFlow()/2;
        ResidualNetwork res(res_net, new_flow1, _offset);
        map<int, bool> persistencies1 = res.fixStrongVariables();
        p_to_add = res.fixAllVariables();
        persistencies1.insert(p_to_add.begin(), p_to_add.end());

        /* Remove previous penalty term and add the penalty term for the variable i */
        _m_adj[_numVariables][i] -= U+1;
        _m_adj[complement(i)][complement(_numVariables)] -= U+1;
        _m_adj[_numVariables][complement(i)] += U+1;
        _m_adj[i][complement(_numVariables)] += U+1;
        net = ImplicationNetwork(_m_adj);
        res_net = net.computeMaxFlow();
        float new_flow2 = net.getMaxFlow()/2;
        res = ResidualNetwork(res_net, new_flow2, _offset);
        map<int, bool> persistencies2 = res.fixStrongVariables();
        p_to_add = res.fixAllVariables();
        persistencies2.insert(p_to_add.begin(), p_to_add.end());

        /* Remove the penalty term */
        _m_adj[_numVariables][complement(i)] -= U+1;
        _m_adj[i][complement(_numVariables)] -= U+1;

        /* Check if lower bound is improved */
        if ((new_flow1 + _max_flow + (float)_offset > _lower_bound) && (new_flow2 + _max_flow + (float)_offset > _lower_bound)){
            if (new_flow1 > new_flow2)
                _lower_bound = new_flow2 + _max_flow + (float)_offset;
            else
                _lower_bound = new_flow1 + _max_flow + (float)_offset;
        }

        /* Check if new persistencies have been found */
        auto it = _var_names.begin();
        if (new_flow1 > (float)U) {
            advance(it, i);
            new_pers.insert({*it, false});
        }
        else if (new_flow2 > (float)U) {
            advance(it, i);
            new_pers.insert({*it, true});
        }
        for (auto & item : persistencies1){
            if (persistencies2.find(item.first) != persistencies2.end()){
                if (item.second == persistencies2[item.first]) {
                    it = _var_names.begin();
                    advance(it, item.first);
                    new_pers.insert({*it, item.second});
                }
            }
        }
    }
    _offset += removeFixedVarsInSCC(new_pers);
    net = ImplicationNetwork(_m_adj);
    _m_adj = net.computeMaxFlow();
    p_to_add = fixStrongVariables();
    new_pers.insert(p_to_add.begin(), p_to_add.end());
    p_to_add = fixAllVariables();
    new_pers.insert(p_to_add.begin(), p_to_add.end());
    _numProbing = (int) new_pers.size();
    return new_pers;
}

map<int, bool> ResidualNetwork::probing_v2() {
    map<int, bool> new_pers;
    int old_num_var = _numVariables;
    ImplicationNetwork net;
    map<int, bool> p_to_add;
    int U = DDT_onepass() + 1;
    for (int i=0; i<_numVariables; i++){
        /* Add the penalty term for the complement of variable i */
        _m_adj[_numVariables][i] += U+1;
        _m_adj[complement(i)][complement(_numVariables)] += U+1;
        net = ImplicationNetwork(_m_adj);
        vector<vector<int>> res_net = net.computeMaxFlow();
        float new_flow1 = net.getMaxFlow()/2;
        ResidualNetwork res(res_net, new_flow1, _offset);
        map<int, bool> persistencies1 = res.fixStrongVariables();
        p_to_add = res.fixAllVariables();
        persistencies1.insert(p_to_add.begin(), p_to_add.end());

        /* Remove previous penalty term and add the penalty term for the variable i */
        _m_adj[_numVariables][i] -= U+1;
        _m_adj[complement(i)][complement(_numVariables)] -= U+1;
        _m_adj[_numVariables][complement(i)] += U+1;
        _m_adj[i][complement(_numVariables)] += U+1;
        net = ImplicationNetwork(_m_adj);
        res_net = net.computeMaxFlow();
        float new_flow2 = net.getMaxFlow()/2;
        res = ResidualNetwork(res_net, new_flow2, _offset);
        map<int, bool> persistencies2 = res.fixStrongVariables();
        p_to_add = res.fixAllVariables();
        persistencies2.insert(p_to_add.begin(), p_to_add.end());

        /* Remove the penalty term */
        _m_adj[_numVariables][complement(i)] -= U+1;
        _m_adj[i][complement(_numVariables)] -= U+1;

        /* Check if lower bound is improved */
        if ((new_flow1 + _max_flow + (float)_offset > _lower_bound) && (new_flow2 + _max_flow + (float)_offset > _lower_bound)){
            if (new_flow1 > new_flow2)
                _lower_bound = new_flow2 + _max_flow + (float)_offset;
            else
                _lower_bound = new_flow1 + _max_flow + (float)_offset;
        }

        /* Check if new persistencies have been found */
        map<int, bool> tmp_pers;
        auto it = _var_names.begin();
        if (new_flow1 > (float)U) {
            advance(it, i);
            tmp_pers.insert({*it, false});
        }
        else if (new_flow2 > (float)U) {
            advance(it, i);
            tmp_pers.insert({*it, true});
        }
        for (auto & item : persistencies1){
            if (persistencies2.find(item.first) != persistencies2.end()){
                if (item.second == persistencies2[item.first]) {
                    it = _var_names.begin();
                    advance(it, item.first);
                    tmp_pers.insert({*it, item.second});
                }
            }
        }
        _offset += removeFixedVarsInSCC(tmp_pers);
        new_pers.insert(tmp_pers.begin(), tmp_pers.end());
    }
    net = ImplicationNetwork(_m_adj);
    _m_adj = net.computeMaxFlow();
    p_to_add = fixStrongVariables();
    new_pers.insert(p_to_add.begin(), p_to_add.end());
    p_to_add = fixAllVariables();
    new_pers.insert(p_to_add.begin(), p_to_add.end());
    _numProbing = (int) new_pers.size();
    return new_pers;
}

float ResidualNetwork::getLowerBound() const {
    return _lower_bound;
}

vector<vector<int> > ResidualNetwork::getCostFunction(int &offset) {
    vector<vector<int> > f(_numVariables, vector<int>(_numVariables, 0));
    /* Convert quadratic terms */
    for (int i=0; i<_numVariables; i++){
        for (int j=i+1; j<_numVariables; j++){
            /* Outgoing edges */
            /* var i normal and j negated */
            if (_m_adj[i][complement(j)]){
                f[i][j] += _m_adj[i][complement(j)];
            }
            /* var i and var j normal */
            if (_m_adj[i][j] != 0){
                f[i][j] -= _m_adj[i][j];
                f[i][i] += _m_adj[i][j];
            }
            /* Ingoing edges */
            /* var i normal and j negated */
            if (_m_adj[complement(j)][i] != 0){
                offset += _m_adj[complement(j)][i];
                f[i][i] -= _m_adj[complement(j)][i];
                f[j][j] -= _m_adj[complement(j)][i];
                f[i][j] += _m_adj[complement(j)][i];
            }
            /* var i and var j normal */
            if (_m_adj[j][i] != 0){
                f[j][j] += _m_adj[j][i];
                f[i][j] -= _m_adj[j][i];
            }
        }
    }
    /* Convert linear terms */
    for (int i=0; i<(2*_numVariables)+1; i++){
        /* If edge from source to variable it has to be complemented so a negated term gives an offset
         * and a negative linear term */
        if (i < _numVariables){
            f[i][i] -= _m_adj[_numVariables][i];
            offset += _m_adj[_numVariables][i];
        }
        /* If edge from source to negated variable it only gives a positive linear term */
        if (i > _numVariables){
            f[complement(i)][complement(i)] += _m_adj[_numVariables][i];
        }
    }
    return f;
}

list<int> ResidualNetwork::getVarNames() const {
    return _var_names;
}

void ResidualNetwork::getStats(int &strongP, int &strongPscc, int &weakP, int &numProbing, int &complete_scc) const {
    strongP = _numStrongP;
    strongPscc = _numStrongPscc;
    weakP = _numWeakP;
    numProbing = _numProbing;
    complete_scc = _num_completescc;
}

void ResidualNetwork::print() {
    for(int i=0; i<(2*(_numVariables+1)); i++){
        for(int j=0; j<(2*(_numVariables+1)); j++){
            cout << _m_adj[i][j] << " ";
        }
        cout << endl;
    }
}

vector<vector<int> > ResidualNetwork::getResidual() {
    return _m_adj;
}
