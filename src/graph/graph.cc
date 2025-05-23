#include "graph/graph.h"
#include<queue>
#include <fstream>
#include<unistd.h>

Graph::Graph(const std::string path) {
    std::ifstream fin;
    std::cout<<"path:"<<path<<std::endl;
    fin.open(path + "/graph-sort.bin", std::ios::binary | std::ios::in);
    fin.read((char *)&n, sizeof(int));
    fin.read((char *)&n, sizeof(int));
    fin.read((char *)&m, sizeof(long long));
    degree = new int[n];
    nbrAll = new int[m];
    nbr = new int *[n];
    fin.read((char *)degree, sizeof(int) * n);
    fin.read((char *)nbrAll, sizeof(int) * m);
    fin.close();

    std::cout << "n = " << n << ", m = " << m / 2 << std::endl;

    long long p = 0;
    for (int i = 0; i < n; i++) {
        nbr[i] = nbrAll + p;
        p += degree[i];
    }

    m /= 2;
    edge = new Edge[m];
    bitrussNumber = new int[m];
    //isInPeelList = new char[m];
    edgeToPeel = m;
    for (ui i = 0; i < m; i++) {
        edge[i].id = i;
    }
    std::memset(bitrussNumber, 0, sizeof(int) * m);
    //std::memset(isInPeelList, 0, sizeof(char) * m);
}

Graph::~Graph() {
    delete[] degree;
    degree = nullptr;
    delete[] nbrAll;
    nbrAll = nullptr;
    nbr = nullptr;
    delete[] edge;
    edge = nullptr;
    delete[] bloom;
    bloom = nullptr;
}

void Graph::construct_index() {
    double start = get_current_time();
    ui **e = new ui *[n];
    for (int u = 0; u < n; u++) {
        e[u] = new ui[degree[u]];
    }
    
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
    /*
    std::cout << std::fixed << std::setprecision(6)
              << "Edge ID assignment time:\t" << get_current_time() - start
              << "sec\n";
    */
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
    /*
    std::cout << std::fixed << std::setprecision(6)
              << "Butterfly counting time:\t" << get_current_time() - start
              << "sec\n";
    */
    delete[] lastUseVertex;
    lastUseVertex = nullptr;
    delete[] lastCount;
    lastCount = nullptr;
    delete[] visitedStatusForW;
    visitedStatusForW = nullptr;

    
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


    //int maxButterflySupport = 0;



    extraBloom.id = bloomCount;
    extraBloom.bloomNumber = m;
    //std::cout << "maxButterflySupport:" << maxButterflySupport << std::endl;
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
    bloomNumber.clear();
    bloomNumber.shrink_to_fit(); 

    e = nullptr;
    butterflyCount = nullptr;
    std::cout << std::fixed << std::setprecision(6)
              << "Index construction time:\t" << get_current_time() - start1
              << "sec\n";
}

void Graph::remove_edge_from_bloom_by_index(int bloomID, pair_t index) {
    affect_edge_t affectEdgeInfo = bloom[bloomID].remove_member_by_index(index);
    if (affectEdgeInfo.first == -1) {
        return;
    } else {
        ui affectEdgeID = affectEdgeInfo.first;
        ui affectIndex = affectEdgeInfo.second;
        edge[affectEdgeID].set_reverse_index_by_index(affectIndex, index);
    }
}

void Graph::remove_edge_from_extra_bloom_by_index(pair_t index) {
    ui affectEdgeID = extraBloom.remove_member_by_index_id_only(index);
    if (affectEdgeID == -1) {
        return;
    } else {
        edge[affectEdgeID].set_reverse_index_in_extra_bloom(index);
    }
}

void Graph::remove_bloom_from_edge_by_index(ui edgeID, ui index) {
    affect_bloom_t affectBloomInfo =
            edge[edgeID].remove_host_bloom_by_index(index);
    if (std::get<0>(affectBloomInfo) == -1) {
        return;
    } else {
        if (std::get<1>(affectBloomInfo).first != -1)
            bloom[std::get<0>(affectBloomInfo)].set_reverse_index_by_index(
                    std::get<1>(affectBloomInfo), index);
        edge[std::get<2>(affectBloomInfo)].set_twin_index_by_index(
                std::get<3>(affectBloomInfo), index);
    }
}
void Graph::bitruss_decomposition() {
    std::ifstream statm_file("/proc/self/statm");
    if (statm_file) {
        size_t size, resident, share, text, lib, data, dt;
        statm_file >> size >> resident >> share >> text >> lib >> data >> dt;
        std::cout << "Memory usage: " << resident * sysconf(_SC_PAGESIZE) / 1024 << " KB" << std::endl;
    } else {
        std::cerr << "Failed to open /proc/self/statm" << std::endl;
    }
    ui visitedEdge = 0;
    std::queue<ui> peelList;
    //std::vector<ui> peelListTmp;
    std::vector<ui> matureList;
    double start = get_current_time();
    std::cout << "bitruss decomposing..." << std::endl;
    while (visitedEdge < edgeToPeel) {
        if (peelList.empty()) {
            extraBloom.send_value_to_member( matureList, peelList,edge);
            //extraBloom.increse_counter(1);
            for (ui i = 0; i < matureList.size(); i++) {
                ui edgeID = matureList[i];
                check_mature_edge(edgeID, peelList);
            }
            matureList.clear();
        } else {
            while (!peelList.empty()) {
                ui edgeID = peelList.front(); 
                peelList.pop();                      
                
                if (bitrussNumber[edgeID] != 0)
                    continue;
                
                bitrussNumber[edgeID] = extraBloom.get_counter();
                peel_edge(edgeID, peelList);         
                visitedEdge++;
            }
        }
    }
    std::cout << std::fixed << std::setprecision(6)
              << "Bitruss decomposition time:\t" << get_current_time() - start
              << "sec\n";
}

void Graph::peel_edge(ui edgeID, std::queue<ui> &peelList) {
    auto *peelEdge = &edge[edgeID];
    pair_t index = peelEdge->get_reverse_index_in_extra_bloom();

    // remove peel edge from extra bloom
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

        // deal with twin edge
        currentBloom->send_value_to_member( bloomNumber - 1, indexInHostBloom,
                                           edge);
        remove_edge_from_bloom_by_index(bloomID, indexInHostBloom);
        remove_bloom_from_edge_by_index(twinEdgeID, indexInTwinEdge);
        edge[twinEdgeID].decrease_butterfly_support(
                currentBloom->get_counter() + bloomNumber - 1);
        if (edge[twinEdgeID].check_maturity()) {
            check_mature_edge(twinEdgeID, peelList);
        }

        currentBloom->bloomNumber--;
        std::vector<ui> matureList;
        currentBloom->send_value_to_member(matureList, peelList, edge);
        //currentBloom->increse_counter(1);
        for (ui j = 0; j < matureList.size(); j++) {
            ui currentEdgeID = matureList[j];
            check_mature_edge(currentEdgeID, peelList);
        }
    }
}

int Graph::collect_counter(ui edgeID) {
    int counter = 0;
    for (ui i = 0; i < edge[edgeID].get_host_bloom_number(); i++) {
        counter +=
                bloom[edge[edgeID].get_host_bloom_id_by_index(i)].get_counter();
    }
    return counter;
}

void Graph::check_mature_edge(ui edgeID,
                              std::queue<ui> &peelList) {
    if (edge[edgeID].isPeel)
        return;

    int counterSum = collect_counter(edgeID);
    int requiredSupport = edge[edgeID].get_butterfly_support();
    int extraCounter = extraBloom.get_counter();

    if (counterSum + extraCounter >= requiredSupport) {
        edge[edgeID].isPeel = true;
        peelList.push(edgeID);
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

void Graph::output_bitruss_number(std::string outputDir) {
    std::ofstream fout;
    fout.open(outputDir + "/bn.txt", std::ios::out);
    for (ui i = 0; i < m; i++) {
        fout << i << "\t" << bitrussNumber[i] << std::endl;
    }
    fout.close();
}
