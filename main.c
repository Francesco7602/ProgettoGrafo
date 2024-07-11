#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <SDL2/SDL.h>
#include <signal.h>
#include <unistd.h>

// Definizione di alcune costanti globali
#define PI 3.14159265358979323
#define TEMPERATURE_DECAY_RATE 0.98
#define MIN_MOVEMENT_THRESHOLD 250 // Movimento minimo desiderato
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
    SimpleGraph* graph;
    Force* net_forces;
    size_t start;
    size_t end;
} ThreadData;

// Varibili globali
int peso;
volatile sig_atomic_t stop = 0; // Variabile volatile per indicare l'interruzione
double temperature;
double k; // Distanza ideale
int movement_detected = 0;
int NUM_THREADS;
double scaleFactor; // Fattore di scala iniziale
double offsetX; // Offeset iniziale
double offsetY;

void handle_signal(int sig) {
    if (sig == SIGINT) {
        stop = 1;
    }
}

// Prototipi delle funzioni
void drawGraph(SDL_Renderer* renderer, SimpleGraph* graph, int screenWidth, int screenHeight);
SimpleGraph readGraphFile(char* file_name);
void calculateRepulsion(Force* net_forces, SimpleGraph* graph, size_t start, size_t end);
void calculateAttraction(Force* net_forces, SimpleGraph* graph, size_t start, size_t end);
void* threadRepulsion(void* arg);
void* threadAttraction(void* arg);
Force* initializeForceVector(SimpleGraph graph);
void moveNodes(Force* net_forces, SimpleGraph* graph);
void repeat(int* response);
void getMaxNodeDimensions(SimpleGraph* graph, double* maxX, double* maxY);

// Funzione principale
int main(int argc, char* argv[]) {
    if (argc < 5) {
        printf("Uso corretto: %s nome_file iterazioni(consiglio di iniziare con 100) temperatura(consiglio di iniziare con 500) peso(1 si 0 no)\n", argv[0]);
        return 1; // Esce dal programma con codice di errore
    }
    char* file_name = argv[1];
    int it = atoi(argv[2]);
    temperature = atoi(argv[3]);
    peso = atoi(argv[4]);
    FILE* input = fopen(file_name, "r");
    if (!input) {
        printf("Il file non esiste.\n");
        return 1;
    }
    fclose(input);
    NUM_THREADS = sysconf(_SC_NPROCESSORS_ONLN) *8;
    SimpleGraph graph = readGraphFile(file_name);

    // Initialize SDL
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL_Init Error: %s\n", SDL_GetError());
        exit(EXIT_FAILURE);
    }
    SDL_DisplayMode display_mode;
    if (SDL_GetDesktopDisplayMode(0, &display_mode) != 0) {
        fprintf(stderr, "SDL_GetDesktopDisplayMode Error: %s\n", SDL_GetError());
        SDL_Quit();
        exit(EXIT_FAILURE);
    }
    int screenWidth = display_mode.w - 10;
    int screenHeight = display_mode.h - 10;
    k = sqrt((screenWidth * screenHeight) / graph.node_count);

    SDL_Window* window = SDL_CreateWindow("Graph Visualization",
                                          SDL_WINDOWPOS_UNDEFINED,
                                          SDL_WINDOWPOS_UNDEFINED,
                                          screenWidth, screenHeight,
                                          SDL_WINDOW_FULLSCREEN);
    if (window == NULL) {
        printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (renderer == NULL) {
        printf("Renderer could not be created! SDL_Error: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    double maX, maY;
    getMaxNodeDimensions(&graph, &maX, &maY);// Calcola le dimensioni massime dei nodi per adattare le dimensioni della finestra
    scaleFactor = fmin(screenWidth / (2 * maX), screenHeight / (2 * maY));
    offsetX = 0.0;
    offsetY = 0.0;
    //drawGraph(renderer, &graph, screenWidth, screenHeight); // Disegno il grafico
    signal(SIGINT, handle_signal); // Gestore di segnale, chiama handle_signal quando riceve SIGINT
    SDL_Event e; // Variabile per ricevere un evento
    int quit = 0;
    int i = 0;
    int stable_count = 0;
    int stable_threshold = 10; // Numero di iterazioni consecutive stabili necessarie per interrompere il ciclo

   
        
        while (i < it) {
            Force* net_forces = initializeForceVector(graph); // Inizializza il vettore delle forze
            getMaxNodeDimensions(&graph, &maX, &maY);
            scaleFactor = fmin(screenWidth / (2 * maX), screenHeight / (2 * maY));
            // Creazione dei thread per il calcolo della repulsione
            pthread_t repulsion_threads[NUM_THREADS];
            ThreadData repulsion_data[NUM_THREADS];
            size_t nodes_per_thread = graph.node_count / NUM_THREADS;
            for (size_t t = 0; t < NUM_THREADS; t++) { // Calcola la repulsione
                repulsion_data[t].graph = &graph;
                repulsion_data[t].net_forces = net_forces;
                repulsion_data[t].start = t * nodes_per_thread;
                repulsion_data[t].end = (t == NUM_THREADS - 1) ? graph.node_count : (t + 1) * nodes_per_thread;
                pthread_create(&repulsion_threads[t], NULL, threadRepulsion, &repulsion_data[t]);
            }
            for (size_t t = 0; t < NUM_THREADS; t++) {
                pthread_join(repulsion_threads[t], NULL);
            }

            // Creazione dei thread per il calcolo dell'attrazione
            pthread_t attraction_threads[NUM_THREADS];
            ThreadData attraction_data[NUM_THREADS];
            size_t edges_per_thread = graph.edge_count / NUM_THREADS;
            for (size_t t = 0; t < NUM_THREADS; t++) { // Calcolal'attrazione
                attraction_data[t].graph = &graph;
                attraction_data[t].net_forces = net_forces;
                attraction_data[t].start = t * edges_per_thread;
                attraction_data[t].end = (t == NUM_THREADS - 1) ? graph.edge_count : (t + 1) * edges_per_thread;
                pthread_create(&attraction_threads[t], NULL, threadAttraction, &attraction_data[t]);
            }
            for (size_t t = 0; t < NUM_THREADS; t++) {
                pthread_join(attraction_threads[t], NULL);
            }

            moveNodes(net_forces, &graph); // Aggiorna la posizione dei nodi
            //drawGraph(renderer, &graph, screenWidth, screenHeight);
            free(net_forces);

            if (quit || stop) {
                break;
            }

            // Controllo temperatura ed eventuale decadimento
            if (temperature > 1.0) {
                temperature *= TEMPERATURE_DECAY_RATE;
            }
            if (i == it- 1) { // Controlla se ha finito il numero di iterazioni
                if (movement_detected != 0) { // Se il grafo ancora si muove
                    i = it -10; // Fa 10 iterazioni in più
                }
            }
            i += 1;
            //SDL_Delay(5);
        }
        //drawGraph(renderer, &graph, screenWidth, screenHeight);
    
    writeFinalPositions("final_positions.txt", &graph);
    // Cleanup SDL
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    free(graph.nodes);
    free(graph.edges);
    return 0;
}



void writeFinalPositions(char* file_name, SimpleGraph* graph) {
    FILE* output = fopen(file_name, "w");
    if (!output) {
        fprintf(stderr, "Errore nell'apertura del file %s.\n", file_name);
        return;
    }

    for (size_t i = 0; i < graph->node_count; i++) {
        fprintf(output, "Nodo %zu: x = %.2f, y = %.2f\n", i, graph->nodes[i].x, graph->nodes[i].y);
    }

    fclose(output);
}


void drawGraph(SDL_Renderer* renderer, SimpleGraph* graph, int screenWidth, int screenHeight) {
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255); // Colore bianco
    SDL_RenderClear(renderer); // Pulisce lo schermo
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); // Colore nero per gli archi
    for (size_t i = 0; i < graph->edge_count; i++) {
        size_t node1 = graph->edges[i].start;
        size_t node2 = graph->edges[i].end;
        SDL_RenderDrawLine(renderer,
                           (int)((graph->nodes[node1].x + offsetX) * scaleFactor + screenWidth / 2),
                           (int)((graph->nodes[node1].y + offsetY) * scaleFactor + screenHeight / 2),
                           (int)((graph->nodes[node2].x + offsetX) * scaleFactor + screenWidth / 2),
                           (int)((graph->nodes[node2].y + offsetY) * scaleFactor + screenHeight / 2)); //Disegna la linea tra due punti
    }

    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255); // Colore rosso per i nodi
    for (size_t i = 0; i < graph->node_count; i++) {
        int x = (int)((graph->nodes[i].x + offsetX) * scaleFactor + screenWidth / 2);
        int y = (int)((graph->nodes[i].y + offsetY) * scaleFactor + screenHeight / 2);
        SDL_Rect rect = {x - (2* scaleFactor), y - (2* scaleFactor), 4* scaleFactor, 4* scaleFactor};
        SDL_RenderFillRect(renderer, &rect); //Disegno i nodi
    }
    SDL_RenderPresent(renderer); // Disegna su schermo
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
    memcpy(graph.edges, temp_edges, graph.edge_count * sizeof(Edge)); // copio i blocchi di memoria

    // Alloca la memoria per i nodi e inizializza le posizioni
    graph.nodes = malloc(graph.node_count * sizeof(Node));
    for (size_t i = 0; i < graph.node_count; i++) {
        graph.nodes[i].x = cos((2 * PI * i) / graph.node_count);
        graph.nodes[i].y = sin((2 * PI * i) / graph.node_count);
    }

    free(temp_edges);
    return graph;
}


void getMaxNodeDimensions(SimpleGraph* graph, double* maxX, double* maxY) { // Trova le coordinate massime dei nodi, in valore assoluto
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

Force* initializeForceVector(SimpleGraph graph) { // inizializza i vettori forza a 0
    Force* net_forces = malloc(graph.node_count * sizeof(Force));
    for (size_t i = 0; i < graph.node_count; i++) {
        net_forces[i].x = 0.0;
        net_forces[i].y = 0.0;
    }
    return net_forces;
}

void calculateRepulsion(Force* net_forces, SimpleGraph* graph, size_t start, size_t end) { // Forza repulsiva
    for (size_t i = start; i < end; i++) {
        for (size_t j = 0; j < graph->node_count; j++) {
            if (i != j) { // Se non si tratta dello stesso nodo
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

void calculateAttraction(Force* net_forces, SimpleGraph* graph, size_t start, size_t end) { // Forza attrativa
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

void* threadRepulsion(void* arg) { //
    ThreadData* data = (ThreadData*)arg;
    calculateRepulsion(data->net_forces, data->graph, data->start, data->end);
    return NULL;
}

void* threadAttraction(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    calculateAttraction(data->net_forces, data->graph, data->start, data->end);
    return NULL;
}

void moveNodes(Force* net_forces, SimpleGraph* graph) { // Aggiorna le posizioni dei nodi
    for (size_t i = 0; i < graph->node_count; i++) {
        double distance = sqrt(net_forces[i].x * net_forces[i].x + net_forces[i].y * net_forces[i].y);
        if (distance > 0.0) {
            double dx = fmin(distance, temperature) * net_forces[i].x / distance;
            double dy = fmin(distance, temperature) * net_forces[i].y / distance;
            graph->nodes[i].x += dx;
            graph->nodes[i].y += dy;
            if (distance > MIN_MOVEMENT_THRESHOLD) {
                movement_detected = 1;
            }
            else {
                movement_detected = 0;
            }
        }
    }
}
