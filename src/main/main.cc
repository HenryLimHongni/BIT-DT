#include <iostream>
#include <string>

#include "factory/factory.h"
#include <gflags/gflags.h> 

DEFINE_string(filepath, "../dataset/github", "Filepath to the dataset.");
DEFINE_string(mode, "batch", "Mode of graph generation: batch or other modes.");

int main(int argc, char **argv) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    
    auto *factory = new GraphFactory();
    auto *graph = factory->generate_graph(FLAGS_filepath.c_str(), FLAGS_mode);
    graph->construct_index();
    graph->bitruss_decomposition();
    graph->output_bitruss_number(FLAGS_filepath.c_str());

    google::ShutDownCommandLineFlags();
    return 0;
}