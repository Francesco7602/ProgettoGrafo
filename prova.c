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

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (argc < 5) {
        if (world_rank == 0) {
            printf("Uso corretto: %s nome_file iterazioni temperatura peso\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    char* file_name = argv[1];
    int it = atoi(argv[2]);
    temperature = atof(argv[3]);
    peso = atoi(argv[4]);

    SimpleGraph graph;
    if (world_rank == 0) {
        graph = readGraphFile(file_name);
    }

    MPI_Bcast(&graph.node_count, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&graph.edge_count, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    if (world_rank != 0) {
        graph.nodes = malloc(graph.node_count * sizeof(Node));
        graph.edges = malloc(graph.edge_count * sizeof(Edge));
    }

    MPI_Bcast(graph.nodes, graph.node_count * sizeof(Node), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(graph.edges, graph.edge_count * sizeof(Edge), MPI_BYTE, 0, MPI_COMM_WORLD);

    int nodes_per_process = graph.node_count / world_size;
    int my_start = world_rank * nodes_per_process;
    int my_end = (world_rank == world_size - 1) ? graph.node_count : (world_rank + 1) * nodes_per_process;

    k = sqrt((2240 * 1400) / graph.node_count);

    double maX, maY;
    getMaxNodeDimensions(&graph, &maX, &maY);
    scaleFactor = fmin(2240 / (2 * maX), 1400 / (2 * maY));
    offsetX = 0.0;
    offsetY = 0.0;

    for (int i = 0; i < it; i++) {
        Force* net_forces = initializeForceVector(graph);

        calculateRepulsion(net_forces, &graph, my_start, my_end);
        calculateAttraction(net_forces, &graph, my_start, my_end);

        Force* global_forces = malloc(graph.node_count * sizeof(Force));
        MPI_Allgather(net_forces + my_start, nodes_per_process * sizeof(Force), MPI_BYTE,
                      global_forces, nodes_per_process * sizeof(Force), MPI_BYTE,
                      MPI_COMM_WORLD);

        if (world_rank == 0) {
            moveNodes(global_forces, &graph);
        }

        MPI_Bcast(graph.nodes, graph.node_count * sizeof(Node), MPI_BYTE, 0, MPI_COMM_WORLD);

        free(net_forces);
        free(global_forces);

        if (temperature > 1.0) {
            temperature *= TEMPERATURE_DECAY_RATE;
        }
    }

    if (world_rank == 0) {
        FILE* output = fopen("out.txt", "w");
        fprintf(output, "%zu %d\n", graph.node_count, peso);
        for (size_t i = 0; i < graph.node_count; i++) {
            fprintf(output, "%zu %f %f\n", i, graph.nodes[i].x, graph.nodes[i].y);
        }
        fprintf(output, "%zu\n", graph.edge_count);
        for (size_t i = 0; i < graph.edge_count; i++) {
            fprintf(output, "%zu %zu\n", graph.edges[i].start, graph.edges[i].end);
        }
        fclose(output);
    }

    if (world_rank != 0) {
        free(graph.nodes);
        free(graph.edges);
    }

    MPI_Finalize();
    return 0;
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
