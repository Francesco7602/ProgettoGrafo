#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

// Definizione di alcune costanti globali
#define PI 3.14159265358979323
#define TEMPERATURE_DECAY_RATE 0.98
#define MIN_DISTANCE 10.0  // Distanza minima tra i nodi

// Strutture
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

typedef struct {
    double x;
    double y;
    int i;
} ForceMomentanea;



// Prototipi delle funzioni
SimpleGraph readGraphFile(char* file_name);
void calculateRepulsion(Force* net_forces, SimpleGraph* graph, size_t start, size_t end);
void calculateAttraction(Force* net_forces, SimpleGraph* graph, size_t start, size_t end);
Force* initializeForceVector(SimpleGraph graph);
void moveNodes(Force* net_forces, SimpleGraph* graph, size_t start, size_t end);
void getMaxNodeDimensions(SimpleGraph* graph, double* maxX, double* maxY);


int peso;
double temperature;
double k; // Distanza ideale
double scaleFactor; // Fattore di scala iniziale
double offsetX; // Offeset iniziale
double offsetY;


int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 5) {
        if (rank == 0) {
            printf("Uso corretto: %s nome_file iterazioni temperatura peso\n", argv[0]);
        }
        MPI_Finalize();
        return 1; // Esce dal programma con codice di errore
    }

    char* file_name = argv[1];
    int it = atoi(argv[2]);
    temperature = atof(argv[3]);
    peso = atoi(argv[4]);

    SimpleGraph graph;
    if (rank == 0) {
        graph = readGraphFile(file_name);
    }

    // Broadcast del numero di nodi e archi a tutti i processi
    MPI_Bcast(&graph.node_count, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&graph.edge_count, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        // Allocazione della memoria per nodi e archi sugli altri processi
        graph.nodes = malloc(graph.node_count * sizeof(Node));
        graph.edges = malloc(graph.edge_count * sizeof(Edge));
    }


    if (rank == 0) {
        // Broadcast dei nodi e degli archi dal processo 0 a tutti i processi
        MPI_Bcast(graph.nodes, graph.node_count * sizeof(Node), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Bcast(graph.edges, graph.edge_count * sizeof(Edge), MPI_BYTE, 0, MPI_COMM_WORLD);
    } else {
        // Ricezione dei nodi e degli archi dagli altri processi
        MPI_Bcast(graph.nodes, graph.node_count * sizeof(Node), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Bcast(graph.edges, graph.edge_count * sizeof(Edge), MPI_BYTE, 0, MPI_COMM_WORLD);
    }



    double screenWidth = 2240 - 10;
    double screenHeight = 1400 - 10;
    k = sqrt((screenWidth * screenHeight) / graph.node_count);

    double maxX, maxY;
    getMaxNodeDimensions(&graph, &maxX, &maxY);
    double scaleFactor = fmin(screenWidth / (2 * maxX), screenHeight / (2 * maxY));
    double offsetX = 0.0;
    double offsetY = 0.0;
    Force* net_forces = initializeForceVector(graph);
    if (!net_forces) {
        fprintf(stderr, "Error allocating memory for forces\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (int i = 0; i < it; i++) {
        // Calcolo della repulsione
        calculateRepulsion(net_forces, &graph, rank * graph.node_count / size, (rank + 1) * graph.node_count / size);
        // Calcolo dell'attrazione
        calculateAttraction(net_forces, &graph, rank * graph.edge_count / size, (rank + 1) * graph.edge_count / size);

        MPI_Barrier(MPI_COMM_WORLD);

        if (rank == 0) {
            Force* dati = initializeForceVector(graph);

            // Buffer per ricevere tutte le forze dai processi
            Force* all_forces = malloc(size * graph.node_count * sizeof(Force));

            // Array per i counts e i displacements
            int* counts = malloc(size * sizeof(int));
            int* displacements = malloc(size * sizeof(int));

            // Impostazione dei counts e dei displacements
            for (int i = 0; i < size; i++) {
                counts[i] = graph.node_count * sizeof(Force);
                displacements[i] = i * graph.node_count * sizeof(Force);
            }

            // Raccogliere le forze da tutti i processi
            MPI_Gatherv(net_forces, graph.node_count * sizeof(Force), MPI_BYTE,
                        all_forces, counts, displacements, MPI_BYTE,
                        0, MPI_COMM_WORLD);

            // Sommare le forze ricevute
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < graph.node_count; j++) {
                    net_forces[j].x += all_forces[i * graph.node_count + j].x;
                    net_forces[j].y += all_forces[i * graph.node_count + j].y;
                }
            }

            moveNodes(net_forces, &graph, 0, graph.node_count);

            // Liberare la memoria allocata
            free(all_forces);
            free(counts);
            free(displacements);
        } else {
            // Invio delle forze al processo radice
            MPI_Gatherv(net_forces, graph.node_count * sizeof(Force), MPI_BYTE,
                        NULL, NULL, NULL, MPI_BYTE,
                        0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // Controllo della temperatura
        if (temperature > 1.0) {
            temperature *= TEMPERATURE_DECAY_RATE;
        }

        MPI_Bcast(graph.nodes, graph.node_count * sizeof(Node), MPI_BYTE, 0, MPI_COMM_WORLD);
        for (int i = 0; i< graph.node_count; i++){
            net_forces[i].x = 0.0;
            net_forces[i].y = 0.0;
        }

        // Sincronizzazione tra i processi
        MPI_Barrier(MPI_COMM_WORLD);

    }

    if (rank == 0) {
        FILE* output = fopen("out.txt", "w");
        if (!output) {
            printf("Impossibile aprire il file di output per la scrittura.\n");
            MPI_Finalize();
            return 0;
        }
        fprintf(output, "%zu %d\n", graph.node_count, peso);

        // Scrivi le posizioni dei nodi
        for (size_t i = 0; i < graph.node_count; i++) {
            fprintf(output, "%zu %f %f\n", i, graph.nodes[i].x, graph.nodes[i].y);
        }

        fprintf(output, "%zu\n", graph.edge_count);

        // Scrivi gli archi
        for (size_t i = 0; i < graph.edge_count; i++) {
            fprintf(output, "%zu %zu\n", graph.edges[i].start, graph.edges[i].end);
        }

        fclose(output);
    }

    free(graph.nodes);
    free(graph.edges);
    MPI_Finalize();
    return 0;
}

SimpleGraph readGraphFile(char* file_name) {
    FILE* input = fopen(file_name, "r");
    if (!input) {
        fprintf(stderr, "Error opening file\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    SimpleGraph graph;
    fscanf(input, "%zu", &graph.node_count); // Legge il numero totale di nodi

    // Creazione di un array temporaneo per memorizzare i dati degli archi
    size_t temp_capacity = 100; // Capacità iniziale, aumenta dinamicamente se necessario
    size_t temp_edge_count = 0;
    Edge* temp_edges = malloc(temp_capacity * sizeof(Edge));

    if (!temp_edges) {
        fprintf(stderr, "Error allocating memory for edges\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    size_t node1, node2;
    double pesett;

    if (!peso) { // Se non è presente il peso, legge due alla volta
        while (fscanf(input, "%zu %zu", &node1, &node2) == 2) {
            if (temp_edge_count >= temp_capacity) {
                // Aumenta la capacità dell'array temporaneo
                temp_capacity *= 2;
                temp_edges = realloc(temp_edges, temp_capacity * sizeof(Edge));
                if (!temp_edges) {
                    fprintf(stderr, "Error reallocating memory for edges\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
            temp_edges[temp_edge_count].start = node1;
            temp_edges[temp_edge_count].end = node2;
            temp_edge_count++;
        }
    } else { // Se è presente il peso, legge tre alla volta
        while (fscanf(input, "%zu %zu %lf", &node1, &node2, &pesett) == 3) {
            if (temp_edge_count >= temp_capacity) {
                // Aumenta la capacità dell'array temporaneo
                temp_capacity *= 2;
                temp_edges = realloc(temp_edges, temp_capacity * sizeof(Edge));
                if (!temp_edges) {
                    fprintf(stderr, "Error reallocating memory for edges\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
            temp_edges[temp_edge_count].start = node1;
            temp_edges[temp_edge_count].end = node2;
            temp_edge_count++;
        }
    }

    fclose(input);

    graph.nodes = malloc(graph.node_count * sizeof(Node));
    if (!graph.nodes) {
        fprintf(stderr, "Error allocating memory for nodes\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    graph.edges = malloc(temp_edge_count * sizeof(Edge));
    if (!graph.edges) {
        fprintf(stderr, "Error allocating memory for edges\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Inizializza le posizioni dei nodi
    srand(time(NULL));
    for (size_t i = 0; i < graph.node_count; i++) {
        graph.nodes[i].x = (double)rand() / RAND_MAX;
        graph.nodes[i].y = (double)rand() / RAND_MAX;
    }

    // Copia gli archi temporanei nell'array finale degli archi
    memcpy(graph.edges, temp_edges, temp_edge_count * sizeof(Edge));
    graph.edge_count = temp_edge_count;

    free(temp_edges);
    return graph;
}

void calculateRepulsion(Force* net_forces, SimpleGraph* graph, size_t start, size_t end) {
    for (size_t i = start; i < end; i++) {
        for (size_t j = 0; j < graph->node_count; j++) {
            if (i != j) {
                double dx = graph->nodes[i].x - graph->nodes[j].x;
                double dy = graph->nodes[i].y - graph->nodes[j].y;
                double distance = sqrt(dx * dx + dy * dy);

                if (distance < MIN_DISTANCE) {
                    distance = MIN_DISTANCE;
                }

                double force_magnitude = k * k / distance;
                double force_x = force_magnitude * (dx / distance);
                double force_y = force_magnitude * (dy / distance);

                net_forces[i].x += force_x;
                net_forces[i].y += force_y;

            }
        }
    }
}

void calculateAttraction(Force* net_forces, SimpleGraph* graph, size_t start, size_t end) {
    for (size_t i = start; i < end; i++) {
        size_t node1 = graph->edges[i].start;
        size_t node2 = graph->edges[i].end;
        double dx = graph->nodes[node2].x - graph->nodes[node1].x; //forse devo invertire l'ordine dei nodi, prima 1 e po 2'
        double dy = graph->nodes[node2].y - graph->nodes[node1].y;
        double distance = sqrt(dx * dx + dy * dy);
        double force = distance * distance / k;
        net_forces[node1].x += force * dx / distance; //forse devo fare (force * dx / distance) e invertire i+ con - e viceversa
        net_forces[node1].y += force * dy / distance;
        net_forces[node2].x -= force * dx / distance;
        net_forces[node2].y -= force * dy / distance;

    }
}

Force* initializeForceVector(SimpleGraph graph) {
    Force* net_forces = malloc(graph.node_count * sizeof(Force));
    if (!net_forces) {
        fprintf(stderr, "Error allocating memory for forces\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (size_t i = 0; i < graph.node_count; i++) {
        net_forces[i].x = 0;
        net_forces[i].y = 0;
    }
    return net_forces;
}

void moveNodes(Force* net_forces, SimpleGraph* graph, size_t start, size_t end) {
    for (size_t i = start; i < end; i++) {
        double dx = net_forces[i].x;
        double dy = net_forces[i].y;


        double displacement = sqrt(dx * dx + dy * dy);
        if (displacement > temperature) {
            dx = dx / displacement * temperature;
            dy = dy / displacement * temperature;
        }

        graph->nodes[i].x += dx;
        graph->nodes[i].y += dy;
    }
}

void getMaxNodeDimensions(SimpleGraph* graph, double* maxX, double* maxY) {
    *maxX = 0.0;
    *maxY = 0.0;

    for (size_t i = 0; i < graph->node_count; i++) {
        if (graph->nodes[i].x > *maxX) {
            *maxX = graph->nodes[i].x;
        }
        if (graph->nodes[i].y > *maxY) {
            *maxY = graph->nodes[i].y;
        }
    }
}
