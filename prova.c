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
    int i;
    double x;
    double y;
} Data;

// Prototipi delle funzioni
SimpleGraph readGraphFile(char* file_name);
void calculateRepulsion(Force* net_forces, SimpleGraph* graph, size_t start, size_t end);
void calculateAttraction(Force* net_forces, SimpleGraph* graph, size_t start, size_t end);
Force* initializeForceVector(SimpleGraph graph);
void moveNodes(Force* net_forces, SimpleGraph* graph, Data* dati, size_t start, size_t end);
void getMaxNodeDimensions(SimpleGraph* graph, double* maxX, double* maxY);
int peso;
double k, temperature;


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

    MPI_File input;
    MPI_File_open(MPI_COMM_WORLD, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &input);

    SimpleGraph graph;
    if (rank == 0) {
        graph = readGraphFile(file_name);
    }
        MPI_File_close(&input);


    // Distribuzione dei dati ai processi
    MPI_Bcast(&graph.node_count, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&graph.edge_count, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    graph.nodes = malloc(graph.node_count * sizeof(Node));
    graph.edges = malloc(graph.edge_count * sizeof(Edge));

    // Scatter dei nodi e degli archi
    MPI_Scatter(graph.nodes, graph.node_count * sizeof(Node), MPI_BYTE,
                graph.nodes, graph.node_count * sizeof(Node), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Scatter(graph.edges, graph.edge_count * sizeof(Edge), MPI_BYTE,
                graph.edges, graph.edge_count * sizeof(Edge), MPI_BYTE, 0, MPI_COMM_WORLD);


    double screenWidth = 2240 - 10;
    double screenHeight = 1400 - 10;
    k = sqrt((screenWidth * screenHeight) / graph.node_count);

    double maxX, maxY;
    getMaxNodeDimensions(&graph, &maxX, &maxY);
    double scaleFactor = fmin(screenWidth / (2 * maxX), screenHeight / (2 * maxY));
    double offsetX = 0.0;
    double offsetY = 0.0;

    int quit = 0;
    Data dati = malloc((graph.node_count / size) * sizeof(Data));
    
        
    


    for (int i = 0; i < it; i++) {
        // Scatter dei nodi e degli archi
    MPI_Scatter(graph.nodes, graph.node_count * sizeof(Node), MPI_BYTE,
                graph.nodes, graph.node_count * sizeof(Node), MPI_BYTE, 0, MPI_COMM_WORLD);
        Force* net_forces = initializeForceVector(graph);

        // Calcolo della repulsione
        calculateRepulsion(net_forces, &graph, rank * graph.node_count / size, (rank + 1) * graph.node_count / size);

        // Sincronizzazione tra i processi
        MPI_Barrier(MPI_COMM_WORLD);

        // Calcolo dell'attrazione
        calculateAttraction(net_forces, &graph, rank * graph.edge_count / size, (rank + 1) * graph.edge_count / size);

        

        moveNodes(net_forces, &graph, &dati, rank * graph.node_count / size, (rank + 1) * graph.node_count / size);
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0){
            //Data datiMomentanei = malloc((graph.node_count / size) * sizeof(Data));
            for(int i = 1; i < size; i++){
                MPI_Recv(dati, graph.node_count / size, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int j = 0, j <graph.node_count / size, j++){
                    graph.nodes[dati[j].i].x = dati[j].x;
                    graph.nodes[dati[j].i].y = dati[j].y;
                }
            }else{
                MPI_Ssend(dati, graph.node_count / size, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
            }
        }

        free(net_forces);
        // Sincronizzazione tra i processi
        MPI_Barrier(MPI_COMM_WORLD);

        // Controllo della temperatura
        if (temperature > 1.0) {
            temperature *= TEMPERATURE_DECAY_RATE;
        }

        if (quit) {
            break;
        }
    }

   /* Raccolta dei risultati
    MPI_Gather(rank == 0 ? MPI_IN_PLACE : graph.nodes, graph.node_count * sizeof(Node), MPI_BYTE,
               graph.nodes, graph.node_count * sizeof(Node), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Gather(rank == 0 ? MPI_IN_PLACE : graph.edges, graph.edge_count * sizeof(Edge), MPI_BYTE,
               graph.edges, graph.edge_count * sizeof(Edge), MPI_BYTE, 0, MPI_COMM_WORLD);*/

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
        //free(graph.nodes);
        //free(graph.edges);
    }

    MPI_Finalize();
    return 0;
}

SimpleGraph readGraphFile(char* file_name) {
    FILE* input = fopen(file_name, "r");
    SimpleGraph graph;
    fscanf(input, "%zu", &graph.node_count); // Legge il numero totale di nodi

    // Creazione di un array temporaneo per memorizzare i dati degli archi
    size_t temp_capacity = 100; // Capacità iniziale, aumenta dinamicamente se necessario
    size_t temp_edge_count = 0;
    Edge* temp_edges = malloc(temp_capacity * sizeof(Edge));

    size_t node1, node2;
    double pesett;

    if (!peso) { // Se non è presente il peso, legge due alla volta
        while (fscanf(input, "%zu %zu", &node1, &node2) == 2) {
            if (temp_edge_count >= temp_capacity) {
                // Aumenta la capacità dell'array temporaneo
                temp_capacity *= 2;
                temp_edges = realloc(temp_edges, temp_capacity * sizeof(Edge));
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
            }
            temp_edges[temp_edge_count].start = node1;
            temp_edges[temp_edge_count].end = node2;
            temp_edge_count++;
        }
    }

    fclose(input);

    // Alloca la memoria per gli archi nel grafico finale
    graph.edge_count = temp_edge_count;
    graph.edges = malloc(graph.edge_count * sizeof(Edge));
    memcpy(graph.edges, temp_edges, graph.edge_count * sizeof(Edge));

    // Alloca la memoria per i nodi e inizializza le posizioni
    graph.nodes = malloc(graph.node_count * sizeof(Node));
    for (size_t i = 0; i < graph.node_count; i++) {
        graph.nodes[i].x = cos((2 * PI * i) / graph.node_count);
        graph.nodes[i].y = sin((2 * PI * i) / graph.node_count);
    }

    free(temp_edges);
    return graph;
}

void getMaxNodeDimensions(SimpleGraph* graph, double* maxX, double* maxY) {
    *maxX = 0.0;
    *maxY = 0.0;
    for (size_t i = 0; i < graph->node_count; i++) {
        if (fabs(graph->nodes[i].x) > *maxX) {
            *maxX = fabs(graph->nodes[i].x);
        }
        if (fabs(graph->nodes[i].y) > *maxY) {
            *maxY = fabs(graph->nodes[i].y);
        }
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

void calculateRepulsion(Force* net_forces, SimpleGraph* graph, size_t start, size_t end) {
    for (size_t i = start; i < end; i++) {
        for (size_t j = 0; j < graph->node_count; j++) {
            if (i != j) {
                double dx = graph->nodes[i].x - graph->nodes[j].x;
                double dy = graph->nodes[i].y - graph->nodes[j].y;
                double distance = sqrt(dx * dx + dy * dy);
                if (distance > 0.0) {
                    double force = - (k * k) / distance;
                    net_forces[i].x -= force * dx / distance;
                    net_forces[i].y -= force * dy / distance;
                }
            }
        }
    }
}

void calculateAttraction(Force* net_forces, SimpleGraph* graph, size_t start, size_t end) {
    for (size_t i = start; i < end; i++) {
        size_t node1 = graph->edges[i].start;
        size_t node2 = graph->edges[i].end;
        double dx = graph->nodes[node2].x - graph->nodes[node1].x;
        double dy = graph->nodes[node2].y - graph->nodes[node1].y;
        double distance = sqrt(dx * dx + dy * dy);
        double force = distance * distance / k;
        net_forces[node1].x += force * dx / distance;
        net_forces[node1].y += force * dy / distance;
        net_forces[node2].x -= force * dx / distance;
        net_forces[node2].y -= force * dy / distance;
    }
}

void moveNodes(Force* net_forces, SimpleGraph* graph, Data* dati, size_t start, size_t end) {
    int j = 0;
    for (size_t i = start; i < end; i++) {
        double distance = sqrt(net_forces[i].x * net_forces[i].x + net_forces[i].y * net_forces[i].y);
        dati[j].i = i;
        if (distance > 0.0) {
            double dx = fmin(distance, temperature) * net_forces[i].x / distance;
            double dy = fmin(distance, temperature) * net_forces[i].y / distance;
            graph->nodes[i].x += dx;
            graph->nodes[i].y += dy;
          
        }
        dati[j].x = graph->nodes[i].x;
        dati[j].y = graph->nodes[i].x;
        j++;
    }
}
