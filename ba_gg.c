#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "pcg_basic.h"


double gettimenow() {
    return (double)clock() / CLOCKS_PER_SEC;
}

int main() {
    printf("Generating network...\n");

    const int m = 8;         // Initial number of nodes
    const int n = 5000000;    // Total number of nodes to be generated
    int i, j, k;
    int new_nodes_edge_targets[m];
    int *edge_ends = malloc((m + 2 * m * (n - m)) * sizeof(int));
    double t0 = gettimenow();
    FILE *fp;

    pcg32_srandom(42, 42);

    if (edge_ends == NULL) {
        perror("Error allocating memory");
        exit(EXIT_FAILURE);
    }

    // Open file for writing edges
    fp = fopen("network.txt", "w");
    if (fp == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Initialize the edge end list with eight unconnected nodes
    for (i = 0; i < m; ++i) {
        edge_ends[i] = i;
    }

    // Add nodes until we have n total
    for (i = m; i < n; ++i) {
        // For each new node, connect it with edges to m different nodes from the edge end list
        for (j = 0; j < m; ++j) {
            int is_duplicate;
            do {
                // draw random node from the edge end list (this provides preferential attachment)
                new_nodes_edge_targets[j] = edge_ends[(int)pcg32_boundedrand(m + 2 * m * (i - m))];
                // Check each new node against the previous ones and draw again if duplicate
                is_duplicate = 0;
                for (k = 0; k < j; ++k) {
                    if (new_nodes_edge_targets[k] == new_nodes_edge_targets[j]) {
                        is_duplicate = 1;
                        break;
                    }
                }
            } while (is_duplicate);
        }
        // Write edges for the new node to the file
        for (j = 0; j < m; ++j) {
            // Add the edge ends to the list
            edge_ends[m + 2 * m * (i - m) + 2 * j] = i;
            edge_ends[m + 2 * m * (i - m) + 2 * j + 1] = new_nodes_edge_targets[j];
            // Write the edge to the file
            fprintf(fp, "%d %d\n", i, new_nodes_edge_targets[j]);
        }
    }

    free(edge_ends);
    fclose(fp);
    printf("Graph generation complete.\n");

    // Print the time taken
    printf("Generation time cost: %f seconds\n", gettimenow() - t0);

    return 0;
}
