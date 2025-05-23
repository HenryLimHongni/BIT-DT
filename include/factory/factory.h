#ifndef FACTORY_H
#define FACTORY_H

#include "graph/graph.h"
#include "graph/graph_batch.h"

class GraphFactory {
public:
    GraphFactory() {}

    Graph *generate_graph(const std::string filePath, const std::string method);
};

#endif