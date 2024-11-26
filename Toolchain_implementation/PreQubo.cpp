//
// Created by W10 on 23/12/2022.
//

#include "PreQubo.h"
#include "ImplicationNetwork.h"
#include "ResidualNetwork.h"
#include <fstream>
#include <utility>
#include <vector>
#include <map>
#include <list>
#include <string>

using namespace std;

static int g_offset;
static float g_lower_bound;

pair<vector<pair<vector<vector<int>>, list<int> > >, map<int, bool> > PreQubo(vector<vector<int> >& f, int &offset, string &filename) {
    fstream statfile;
    map<int, bool> persistencies;
    g_offset = offset;
    if (!filename.empty())
        statfile.open(filename, ios::app);
    ImplicationNetwork net = ImplicationNetwork(f, g_offset);
    g_offset = net.getOffset();
    vector<vector<int>> res_net = net.computeMaxFlow();
    float m_flow = net.getMaxFlow();
    ResidualNetwork res(res_net, m_flow, g_offset);
    persistencies = res.fixStrongVariables();
    map<int, bool> p_to_add = res.fixAllVariables();
    persistencies.insert(p_to_add.begin(), p_to_add.end());
    int pers_tot = (int) persistencies.size();
    p_to_add = res.probing();
    int probing_vars = (int) p_to_add.size();
    persistencies.insert(p_to_add.begin(), p_to_add.end());
    g_lower_bound = res.getLowerBound();
    if (!filename.empty()) {
        int strongP, strongPscc, weakP, probing, complete_scc;
        res.getStats(strongP, strongPscc, weakP, probing, complete_scc);
        statfile << "Num_of_variables " << (int) f.size() << endl;
        statfile << "Strong_persistencies " << strongP << endl;
        //statfile << "Strong_persistencies_sccMethod " << strongPscc << endl;
        //statfile << "Weak_persistencies_sccMethod " << weakP << endl;
        statfile << "All_persistencies " << pers_tot << endl;
        statfile << "Complete_sccs " << complete_scc << endl;
        statfile << "Probing_persistencies " << probing_vars << endl;
        statfile << "Num_persistencies_found " << (int) persistencies.size() << endl;
        statfile << "Persistencies_percentage " << ((float) persistencies.size()) / ((float) f.size()) * 100 << "%"
                 << endl;
    }
    vector<ResidualNetwork> subgraphs = res.trivialDecomposition();
    int n_dec = (int) subgraphs.size();
    if (!filename.empty())
        statfile << "Num_trivial_Decomposition " << n_dec << endl;
    int num_sccd = 0;
    for (int i=0; i<n_dec; i++){
        vector<ResidualNetwork> tmp = subgraphs[i].SCCDecomposition();
        num_sccd += (int) tmp.size();
        /* Substitute graph if the whole graph is a strongly connected component */
        if (subgraphs[i].getVarNames().empty())
            subgraphs[i] = tmp[0];
        else
            subgraphs.insert(subgraphs.end(), tmp.begin(), tmp.end());
    }
    if (!filename.empty()) {
        statfile << "Num_SCC_Decomposition " << num_sccd << endl;
        statfile << "Total_Decompositions " << (int) subgraphs.size() << endl;
    }
    vector<pair<vector<vector<int> >, list<int>> > subqubos;
    int max_vars = 0;
    for (auto & subgraph : subgraphs){
        vector<vector<int>> func = subgraph.getCostFunction(g_offset);
        list<int> var_list;
        if (func.size() > max_vars)
            max_vars = (int) func.size();
        var_list = subgraph.getVarNames();
        subqubos.emplace_back(pair<vector<vector<int>>, list<int>>(func, var_list));
    }
    if (!filename.empty()) {
        statfile << "Largest_subfunction " << max_vars << endl;
        statfile << "Variable_reduction_percentage " << ((float) f.size() - (float) max_vars) / ((float) f.size()) * 100 << "%"
                 << endl;
        statfile << "-------------------------------------------" << endl;
        statfile.close();
    }
    return {subqubos, persistencies};
}

int getOffset() {
    return g_offset;
}

float getLowerBound() {
    return g_lower_bound;
}

vector<pair<vector<pair<vector<vector<int>>, list<int>>>, map<int, bool>>> ShannonDecomposition(vector<vector<int>> &Q, int &n_iterations, int &flag){
    vector<pair<vector<pair<vector<vector<int>>, list<int>>>, map<int, bool>>> data_to_ret;
    /* Transform QUBO to a residual network
     * (Ignore maximum flow since this decomposition method should be applied after the execution of PreQubo) */
    ImplicationNetwork net = ImplicationNetwork(Q, 0);
    vector<vector<int>> res_net = net.computeMaxFlow();
    ResidualNetwork res(res_net, 0, 0);

    /* Apply Shannon decomposition algorithm */
    vector<pair<ResidualNetwork, map<int, bool>>> decomposed_networks = res.ShannonDecomposition(n_iterations, flag);

    /* Find persistencies and apply decompositions to all the subnetworks obtained */
    for (auto & item : decomposed_networks){
        /* Find all persistencies */
        map<int, bool> persistencies;
        persistencies = item.first.fixStrongVariables();
        map<int, bool> p_to_add = item.first.fixAllVariables();
        persistencies.insert(p_to_add.begin(), p_to_add.end());

        /* Add persistencies to the set of selection variables extracted by the Shannon algorithm */
        item.second.insert(persistencies.begin(), persistencies.end());

        /* Apply trivial decomposition */
        vector<ResidualNetwork> new_subnets = item.first.trivialDecomposition();

        int n_dec = (int) new_subnets.size();
        for (int i=0; i<n_dec; i++) {
            /* Apply SCC decomposition and add store the newly found subnetworks */
            vector<ResidualNetwork> tmp = new_subnets[i].SCCDecomposition();
            /* Substitute graph if the whole graph is a strongly connected component */
            if (new_subnets[i].getVarNames().empty() && !tmp.empty())
                new_subnets.emplace_back(tmp[0]);
            else {
                new_subnets.insert(new_subnets.end(), tmp.begin(), tmp.end());
            }
        }

        /* Transform back networks to QUBO */
        vector<pair<vector<vector<int> >, list<int>> > subqubos;
        for (auto & subnet : new_subnets){
            int off = 0;
            vector<vector<int>> f = subnet.getCostFunction(off);
            list<int> var_list = subnet.getVarNames();
            subqubos.emplace_back(pair<vector<vector<int>>, list<int>>(f, var_list));
        }
        /* Append the newly found subQUBO matrices to the vector */
        data_to_ret.emplace_back(pair<vector<pair<vector<vector<int>>, list<int>>>, map<int, bool>>(subqubos, item.second));
    }
    return data_to_ret;
}