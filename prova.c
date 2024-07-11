#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#define PI 3.14159265358979323
#define TEMPERATURE_DECAY_RATE 0.98
#define MIN_MOVEMENT_THRESHOLD 250
#define MIN_DISTANCE 10.0

typedef struct {
    double x;
    double y;
} Force;

typedef struct {
    double x;
    double y;
} Node;

typedef struct {
    size_t start;
    size_t end;
} Edge;

typedef struct {
    Node* nodes;
    Edge* edges;
    size_t node_count;
    size_t edge_count;
} SimpleGraph;

int peso;
volatile sig_atomic_t stop = 0;
double temperature;
double k;
double scaleFactor;
double offsetX;
double offsetY;

SimpleGraph readGraphFile(char* file_name);
void calculateRepulsion(Force* net_forces, SimpleGraph* graph, size_t start, size_t end);
void calculateAttraction(Force* net_forces, SimpleGraph* graph, size_t start, size_t end);
Force* initializeForceVector(SimpleGraph graph);
void moveNodes(Force* net_forces, SimpleGraph* graph);
void getMaxNodeDimensions(SimpleGraph* graph, double* maxX, double* maxY);
void writeNodeData(FILE* output, SimpleGraph* graph, int start, int end);
void writeEdgeData(FILE* output, SimpleGraph* graph, int start, int end);

int main(int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0) {
            printf("Uso corretto: %s nome_file_grafo\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    char* file_name = argv[1];
    SimpleGraph graph;

    if (rank == 0) {
        graph = readGraphFile(file_name);
    }

    // Broadcast delle informazioni sui nodi e archi
    MPI_Bcast(&graph.node_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&graph.edge_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        graph.nodes = malloc(graph.node_count * sizeof(Node));
        graph.edges = malloc(graph.edge_count * sizeof(Edge));
    }

    // Broadcast dei dati dei nodi e archi
    MPI_Bcast(graph.nodes, graph.node_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(graph.edges, graph.edge_count, MPI_INT, 0, MPI_COMM_WORLD);

    int screenWidth = 2240 - 10;
    int screenHeight = 1400 - 10;
    double k = sqrt((screenWidth * screenHeight) / graph.node_count);

    double maX, maY;
    getMaxNodeDimensions(&graph, &maX, &maY);
    double scaleFactor = fmin(screenWidth / (2 * maX), screenHeight / (2 * maY));

    double temperature = 500.0;
    int it = 100;
    int i;
    int stable_count = 0;
    int stable_threshold = 10;
    int movement_detected = 0;

    for (i = 0; i < it; i++) {
        Force* net_forces = initializeForceVector(&graph);

        int nodes_per_process = graph.node_count / size;
        int start = rank * nodes_per_process;
        int end = (rank == size - 1) ? graph.node_count : (rank + 1) * nodes_per_process;

        calculateRepulsion(net_forces, &graph, start, end);

        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, net_forces, sizeof(Force) * graph.node_count, MPI_BYTE, MPI_COMM_WORLD);

        int edges_per_process = graph.edge_count / size;
        start = rank * edges_per_process;
        end = (rank == size - 1) ? graph.edge_count : (rank + 1) * edges_per_process;

        calculateAttraction(net_forces, &graph, start, end);

        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, net_forces, sizeof(Force) * graph.node_count, MPI_BYTE, MPI_COMM_WORLD);

        moveNodes(net_forces, &graph);
        free(net_forces);

        if (temperature > 1.0) {
            temperature *= TEMPERATURE_DECAY_RATE;
        }

        if (movement_detected < MIN_MOVEMENT_THRESHOLD) {
            stable_count++;
            if (stable_count >= stable_threshold) {
                break;
            }
        } else {
            stable_count = 0;
        }
    }
// Scrivi i dati su file
    char output_filename[100];
    sprintf(output_filename, "out_%d.txt", rank);
    FILE* output = fopen(output_filename, "w");
    if (!output) {
        printf("Processo %d: Impossibile aprire il file di output.\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    writeNodeData(output, &graph, local_start_node, local_end_node);
    writeEdgeData(output, &graph, local_start_edge, local_end_edge);

    fclose(output);

    // Codice per il controllo di interruzione e output finale
    if (rank == 0) {
        FILE* final_output = fopen("out_final.txt", "w");
        if (!final_output) {
            printf("Processo %d: Impossibile aprire il file di output finale.\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        // Codice per scrivere il file finale
        fclose(final_output);
    }

    MPI_Finalize();
    return 0;
}
// Funzione per scrivere i dati dei nodi su file
void writeNodeData(FILE* output, SimpleGraph* graph, int start, int end) {
    for (int i = start; i < end; i++) {
        fprintf(output, "%zu %f %f\n", i, graph->nodes[i].x, graph->nodes[i].y);
    }
}

// Funzione per scrivere i dati degli archi su file
void writeEdgeData(FILE* output, SimpleGraph* graph, int start, int end) {
    for (int i = start; i < end; i++) {
        fprintf(output, "%zu %zu\n", graph->edges[i].start, graph->edges[i].end);
    }
}

SimpleGraph readGraphFile(char* file_name) {
    FILE* input = fopen(file_name, "r");
    SimpleGraph graph;
    fscanf(input, "%zu", &graph.node_count);

    size_t temp_capacity = 100;
    size_t temp_edge_count = 0;
    Edge* temp_edges = malloc(temp_capacity * sizeof(Edge));

    size_t node1, node2;
    double pesett;

    if (!peso) {
        while (fscanf(input, "%zu %zu", &node1, &node2) == 2) {
            if (temp_edge_count >= temp_capacity) {
                temp_capacity *= 2;
                temp_edges = realloc(temp_edges, temp_capacity * sizeof(Edge));
            }
            temp_edges[temp_edge_count].start = node1;
            temp_edges[temp_edge_count].end = node2;
            temp_edge_count++;
        }
    } else {
        while (fscanf(input, "%zu %zu %lf", &node1, &node2, &pesett) == 3) {
            if (temp_edge_count >= temp_capacity) {
                temp_capacity *= 2;
                temp_edges = realloc(temp_edges, temp_capacity * sizeof(Edge));
            }
            temp_edges[temp_edge_count].start = node1;
            temp_edges[temp_edge_count].end = node2;
            temp_edge_count++;
        }
    }

    fclose(input);

    graph.edge_count = temp_edge_count;
    graph.edges = malloc(graph.edge_count * sizeof(Edge));
    memcpy(graph.edges, temp_edges, graph.edge_count * sizeof(Edge));

    graph.nodes = malloc(graph.node_count * sizeof(Node));
    for (size_t i = 0; i < graph.node_count; i++) {
        graph.nodes[i].x = cos((2 * PI * i) / graph.node_count);
        graph.nodes[i].y = sin((2 * PI * i) / graph.node_count);
    }

    free(temp_edges);
    return graph;
}

void calculateRepulsion(Force* net_forces, SimpleGraph* graph, size_t start, size_t end) {
    for (size_t i = start; i < end; i++) {
        for (size_t j = 0; j < graph->node_count; j++) {
            if (i == j) continue;

            double dx = graph->nodes[i].x - graph->nodes[j].x;
            double dy = graph->nodes[i].y - graph->nodes[j].y;
            double distance = sqrt(dx * dx + dy * dy);

            if (distance < MIN_DISTANCE) {
                distance = MIN_DISTANCE;
            }

            double repulsive_force = (k * k) / distance;
            net_forces[i].x += (dx / distance) * repulsive_force;
            net_forces[i].y += (dy / distance) * repulsive_force;
        }
    }
}

void calculateAttraction(Force* net_forces, SimpleGraph* graph, size_t start, size_t end) {
    for (size_t i = 0; i < graph->edge_count; i++) {
        size_t src = graph->edges[i].start;
        size_t dest = graph->edges[i].end;

        if (src < start || src >= end) continue;

        double dx = graph->nodes[src].x - graph->nodes[dest].x;
        double dy = graph->nodes[src].y - graph->nodes[dest].y;
        double distance = sqrt(dx * dx + dy * dy);

        if (distance < MIN_DISTANCE) {
            distance = MIN_DISTANCE;
        }

        double attractive_force = (distance * distance) / k;
        net_forces[src].x -= (dx / distance) * attractive_force;
        net_forces[src].y -= (dy / distance) * attractive_force;
        net_forces[dest].x += (dx / distance) * attractive_force;
        net_forces[dest].y += (dy / distance) * attractive_force;
    }
}

Force* initializeForceVector(SimpleGraph graph) {
    Force* net_forces = malloc(graph.node_count * sizeof(Force));
    for (size_t i = 0; i < graph.node_count; i++) {
        net_forces[i].x = 0.0;
        net_forces[i].y = 0.0;
    }
    return net_forces;
}

void moveNodes(Force* net_forces, SimpleGraph* graph) {
    for (size_t i = 0; i < graph->node_count; i++) {
        double dx = net_forces[i].x;
        double dy = net_forces[i].y;
        double distance = sqrt(dx * dx + dy * dy);

        if (distance > temperature) {
            dx = (dx / distance) * temperature;
            dy = (dy / distance) * temperature;
        }

        graph->nodes[i].x += dx;
        graph->nodes[i].y += dy;

        if (fabs(dx) < MIN_MOVEMENT_THRESHOLD && fabs(dy) < MIN_MOVEMENT_THRESHOLD) {
            temperature *= TEMPERATURE_DECAY_RATE;
        }
    }
}

void getMaxNodeDimensions(SimpleGraph* graph, double* maxX, double* maxY) {
    *maxX = *maxY = 0.0;
    for (size_t i = 0; i < graph->node_count; i++) {
        if (fabs(graph->nodes[i].x) > *maxX) {
            *maxX = fabs(graph->nodes[i].x);
        }
        if (fabs(graph->nodes[i].y) > *maxY) {
            *maxY = fabs(graph->nodes[i].y);
        }
    }
}
