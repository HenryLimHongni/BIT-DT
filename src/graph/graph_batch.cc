#include "graph/graph_batch.h"
#include<set>
#include <fstream>
#include<unistd.h>
using namespace std;

void GraphBatch::construct_index() {
    double start = get_current_time();
    ui **e = new ui *[n];
    for (int u = 0; u < n; u++) {
        e[u] = new ui[degree[u]];
    }
    //double start = get_current_time();
    ui edgeIdx = 0;
    for (int u = 0; u < n; u++) {
        for (int i = 0; i < degree[u]; i++) {
            int v = nbr[u][i];
            if (u < v) {
                e[u][i] = edgeIdx++;
            } else {
                int l = 0, r = degree[v] - 1, mid = 0;
                while (l < r) {
                    mid = l + (r - l) / 2;
                    if (nbr[v][mid] < u) {
                        l = mid + 1;
                    } else {
                        r = mid;
                    }
                }
                e[u][i] = e[v][l];
            }
        }
    }


    std::vector<int> bloomNumber;
    int u, v, w;
    int *lastUseVertex = new int[n], *lastCount = new int[n],
            *visitedStatusForW = new int[n];
    std::memset(lastCount, 0, sizeof(int) * n);
    std::memset(visitedStatusForW, -1, sizeof(int) * n);

    int **butterflyCount = new int *[n];
    for (int i = 0; i < n; i++) {
        butterflyCount[i] = new int[degree[i]]();

    }
    //start = get_current_time();
    int lastUseIdx = 0;
    for (u = 0; u < n; u++) {
        for (int j = 0; j < lastUseIdx; j++) {
            lastCount[lastUseVertex[j]] = 0;
            visitedStatusForW[lastUseVertex[j]] = -1;
        }
        lastUseIdx = 0;

        for (int i = 0; i < degree[u]; i++) {
            v = nbr[u][i];
            //if (degree[v] > degree[u]) continue;
            for (int j = 0; j < degree[v]; j++) {
                w = nbr[v][j];
                if (w >= u || w >= v)
                    break;
                lastCount[w]++;
                if (lastCount[w] == 1)
                    lastUseVertex[lastUseIdx++] = w;
            }
        }
        for (int i = 0; i < degree[u]; i++) {
            v = nbr[u][i];
            for (int j = 0; j < degree[v]; j++) {
                w = nbr[v][j];
                if (w >= u || w >= v)
                    break;
                int lastCountNumber = lastCount[w];
                if (lastCountNumber > 1) {
                    butterflyCount[u][i] += lastCountNumber - 1;
                    butterflyCount[v][j] += lastCountNumber - 1;
                } else
                    continue;
                if (visitedStatusForW[w] == -1) {
                    ui indexV = edge[e[u][i]].add_host_bloom_and_twin_edge(
                            bloomCount, e[v][j]);
                    ui indexU = edge[e[v][j]].add_host_bloom_and_twin_edge(
                            bloomCount, e[u][i]);
                    edge[e[u][i]].add_host_bloom_index_of_twin_edge(indexU);
                    edge[e[v][j]].add_host_bloom_index_of_twin_edge(indexV);
                    bloomNumber.emplace_back(lastCountNumber);
                    visitedStatusForW[w] = bloomCount++;
                } else {

                    ui indexV = edge[e[u][i]].add_host_bloom_and_twin_edge(
                            visitedStatusForW[w], e[v][j]);
                    ui indexU = edge[e[v][j]].add_host_bloom_and_twin_edge(
                            visitedStatusForW[w], e[u][i]);
                    edge[e[u][i]].add_host_bloom_index_of_twin_edge(indexU);
                    edge[e[v][j]].add_host_bloom_index_of_twin_edge(indexV);
                }
            }
            //edge[e[u][i]].add_butterfly_support(butterflyCount[u][i]);

        }
    }
    delete[] lastUseVertex;
    lastUseVertex = nullptr;
    delete[] lastCount;
    lastCount = nullptr;
    delete[] visitedStatusForW;
    visitedStatusForW = nullptr;

    //start = get_current_time();
    bloom = new Bloom[bloomCount];
    for (int i = 0; i < bloomCount; i++) {
        bloom[i].id = i;
        bloom[i].bloomNumber = bloomNumber[i];
        bloom[i].initialize_space();
    }



    for (u = 0; u < n; u++) {
        for (int i = 0; i < degree[u]; i++) {
            edge[e[u][i]].add_butterfly_support(butterflyCount[u][i]);
        }
    }
    std::cout << std::fixed << std::setprecision(6)
              << " bloom construction time:\t" << get_current_time() - start
              << "sec\n";
    int start1 = get_current_time();





    extraBloom.id = bloomCount;
    extraBloom.bloomNumber = m;
    extraBloom.initialize_space_member_edge_only();
    for (ui i = 0; i < m; i++) {
        auto& currentEdge = edge[i];
        int butterflySupport = currentEdge.get_butterfly_support();

        
        if (butterflySupport == 0) {
            edgeToPeel--;
            continue;
        }

        
        currentEdge.compute_slack_value();

        
        for (int j = 0; j < currentEdge.get_host_bloom_number(); j++) {
            int bloomID = currentEdge.get_host_bloom_id_by_index(j);
            pair_t index = bloom[bloomID].add_member_edge(i, j, edge);

            currentEdge.add_reverse_index_in_host_bloom(index);
        }

        
        pair_t index = extraBloom.add_member_edge(i, edge);
        currentEdge.set_reverse_index_in_extra_bloom(index);
    }
    /*
    std::cout << std::fixed << std::setprecision(6)
              << "Index construction time:\t" << get_current_time() - start
              << "sec\n";
    */
    for (u = 0; u < n; u++) {
        delete[] e[u];
        delete[] butterflyCount[u];
    }
    delete[] e;
    delete[] butterflyCount;
    delete[] nbrAll;
    nbrAll = nullptr;

    delete[] nbr;
    nbr = nullptr;

    delete[] degree;
    degree = nullptr;
    e = nullptr;
    butterflyCount = nullptr;
    bloomNumber.clear();
    bloomNumber.shrink_to_fit(); 
    isAffected = new char[m];
    std::memset(isAffected, 0, sizeof(char) * m);
    bloomCounter = new int[bloomCount];
    std::memset(bloomCounter, 0, sizeof(int) * bloomCount);
    //affectedBloom.reserve(bloomCount / 10);
    //affectEdge.reserve(m / 100);
        std::cout << std::fixed << std::setprecision(6)
              << "Index construction time:\t" << get_current_time() - start1
              << "sec\n";
}

void GraphBatch::check_mature_edge(ui edgeID,
                              vector<ui> &peelList) {
    if (edge[edgeID].isPeel)
        return;

    int counterSum = collect_counter(edgeID);
    int requiredSupport = edge[edgeID].get_butterfly_support();
    int extraCounter = extraBloom.get_counter();

    if (counterSum + extraCounter >= requiredSupport) {
        edge[edgeID].isPeel = true;
        peelList.emplace_back(edgeID);
    } else {
        int trackValue = requiredSupport - counterSum - extraCounter;
        int temp = edge[edgeID].get_slack_value();
        edge[edgeID].compute_slack_value(trackValue);
        if(temp != edge[edgeID].get_slack_value()){

        // Process host bloom edges
        for (ui i = 0; i < edge[edgeID].get_host_bloom_number(); i++) {
            int bloomID = edge[edgeID].get_host_bloom_id_by_index(i);
            pair_t reverseIndex = edge[edgeID].get_reverse_index_in_host_bloom_by_index(i);
            remove_edge_from_bloom_by_index(bloomID, reverseIndex);
            pair_t index = bloom[bloomID].add_member_edge(edgeID, i, edge);
            edge[edgeID].set_reverse_index_by_index(i, index);
        }

        // Process extra bloom edges
        pair_t reverseIndex = edge[edgeID].get_reverse_index_in_extra_bloom();
        remove_edge_from_extra_bloom_by_index(reverseIndex);
        pair_t index = extraBloom.add_member_edge(edgeID, edge);
        edge[edgeID].set_reverse_index_in_extra_bloom(index);
        }
    }
}

void GraphBatch::bitruss_decomposition() {
        std::ifstream statm_file("/proc/self/statm");
    if (statm_file) {
        size_t size, resident, share, text, lib, data, dt;
        statm_file >> size >> resident >> share >> text >> lib >> data >> dt;
        std::cout << "Memory usage: " << resident * sysconf(_SC_PAGESIZE) / 1024 << " KB" << std::endl;
    } else {
        std::cerr << "Failed to open /proc/self/statm" << std::endl;
    }
    visitedEdge = 0;
    std::vector<ui> peelList;
    std::vector<ui> peelListTmp;
    //std::vector<ui> peelListTmp;
    std::vector<ui> matureList;
    /*
    std::ifstream statm_file("/proc/self/statm");
    if (statm_file) {
        size_t size, resident, share, text, lib, data, dt;
        statm_file >> size >> resident >> share >> text >> lib >> data >> dt;
        std::cout << "Memory usage: " << resident * sysconf(_SC_PAGESIZE) / 1024 << " KB" << std::endl;
    } else {
        std::cerr << "Failed to open /proc/self/statm" << std::endl;
    }
    */
    
    double start = get_current_time();
    std::cout << "bitruss decomposing..." << std::endl;
    //cout<<edge[0].get_butterfly_support()<<" "<<edge[0].get_host_bloom_number()<<" "<<edge[0].get_slack_value()<<endl;
    while (visitedEdge < edgeToPeel) {

        if (peelList.empty()) {

            extraBloom.send_value_to_member1( matureList, peelList,edge,isAffected);
            //cout<<edge[0].get_butterfly_support()<<" "<<extraBloom.get_counter()<<endl;

            //extraBloom.increse_counter(1);
            for (ui i = 0; i < matureList.size(); i++) {
                ui edgeID = matureList[i];
                check_mature_edge(edgeID, peelList);
                 isAffected[edgeID] = 0;
                
            }
            matureList.clear();
        } else {
            
            peel_edge(peelList, peelListTmp);
           //visitedEdge += (int)peelList.size();
            peelList.swap(peelListTmp);
            peelListTmp.clear();
        }
    }
    std::cout << std::fixed << std::setprecision(6)
              << "Bitruss decomposition time:\t" << get_current_time() - start
              << "sec\n";
}

void GraphBatch::peel_edge(std::vector<ui> edgeList,
                           std::vector<ui> &peelList) {
    
    affectedBloom.clear();
    affectEdge.clear();
    for (ui i = 0, size = edgeList.size(); i < size; i++) {
        ui edgeID = edgeList[i];
        //cout<<"edgeID:"<<edgeID<<endl;
        auto *peelEdge = &edge[edgeID];
        pair_t index = peelEdge->get_reverse_index_in_extra_bloom();

        remove_edge_from_extra_bloom_by_index(index);
        for (ui i = 0; i < peelEdge->get_host_bloom_number(); i++) {
            int bloomID = peelEdge->get_host_bloom_id_by_index(i);
            TwinInfo twinEdgeInfo = peelEdge->get_twin_edge_info_by_index(i);
            pair_t reverseIndex =
                    peelEdge->get_reverse_index_in_host_bloom_by_index(i);
            auto *currentBloom = &bloom[bloomID];
            int bloomNumber = currentBloom->bloomNumber;
            ui twinEdgeID = twinEdgeInfo.twinEdgeID;
            // remove peel edge from current bloom

            remove_edge_from_bloom_by_index(bloomID, reverseIndex);

            ui indexInTwinEdge = twinEdgeInfo.hostBloomIndex;
            pair_t indexInHostBloom =
                    edge[twinEdgeID].get_reverse_index_in_host_bloom_by_index(
                            indexInTwinEdge);

            if (bloomNumber <= 1) {
                remove_edge_from_bloom_by_index(bloomID, indexInHostBloom);
                remove_bloom_from_edge_by_index(twinEdgeID, indexInTwinEdge);
                edge[twinEdgeID].decrease_butterfly_support(
                        currentBloom->get_counter());
                currentBloom->bloomNumber--;
                continue;
            }
            currentBloom->send_value_to_member( bloomNumber - 1, indexInHostBloom,
                                            edge);
            remove_edge_from_bloom_by_index(bloomID, indexInHostBloom);
            remove_bloom_from_edge_by_index(twinEdgeID, indexInTwinEdge);
            edge[twinEdgeID].decrease_butterfly_support(
                    currentBloom->get_counter() + bloomNumber - 1);
            if (edge[twinEdgeID].check_maturity() && !isAffected[twinEdgeID]) {
                affectEdge.emplace_back(twinEdgeID);
                isAffected[twinEdgeID] = 1;
            }
            if (bloomCounter[bloomID] == 0) {
                    affectedBloom.emplace_back(bloomID);
            }
            bloomCounter[bloomID]++;
        }
        visitedEdge++;
        bitrussNumber[edgeID] = extraBloom.get_counter();
    }
    for (ui i = 0; i < affectedBloom.size(); i++) {
        int bloomID = affectedBloom[i];
        bloom[bloomID].bloomNumber -= bloomCounter[bloomID];
        if (bloom[bloomID].bloomNumber == 0) {
            bloomCounter[bloomID]--;
        }
        bloom[bloomID].send_value_to_member(bloomCounter[bloomID], affectEdge,peelList,
                                            edge, isAffected);
        //bloom[bloomID].increse_counter(bloomCounter[bloomID]);
        bloomCounter[bloomID] = 0;
    }
    

    for (ui i = 0; i < affectEdge.size(); i++) {
        ui currentEdgeID = affectEdge[i];
        check_mature_edge(currentEdgeID, peelList);
        isAffected[currentEdgeID] = 0;
    }



    
}