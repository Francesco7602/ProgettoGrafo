#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <SDL2/SDL.h>

// Definizione di alcune costanti globali
#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 600
#define SCALE_FACTOR 30.0
#define OFFSET_X SCREEN_WIDTH / 2
#define OFFSET_Y SCREEN_HEIGHT / 2

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

void drawGraph(SDL_Renderer* renderer, SimpleGraph* graph, int screenWidth, int screenHeight) {
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255); // Colore bianco
    SDL_RenderClear(renderer); // Pulisce lo schermo
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); // Colore nero per gli archi

    for (size_t i = 0; i < graph->edge_count; i++) {
        size_t node1 = graph->edges[i].start;
        size_t node2 = graph->edges[i].end;
        SDL_RenderDrawLine(renderer,
                           (int)((graph->nodes[node1].x + OFFSET_X) * SCALE_FACTOR),
                           (int)((graph->nodes[node1].y + OFFSET_Y) * SCALE_FACTOR),
                           (int)((graph->nodes[node2].x + OFFSET_X) * SCALE_FACTOR),
                           (int)((graph->nodes[node2].y + OFFSET_Y) * SCALE_FACTOR)); //Disegna la linea tra due punti
    }

    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255); // Colore rosso per i nodi
    for (size_t i = 0; i < graph->node_count; i++) {
        int x = (int)((graph->nodes[i].x + OFFSET_X) * SCALE_FACTOR);
        int y = (int)((graph->nodes[i].y + OFFSET_Y) * SCALE_FACTOR);
        SDL_Rect rect = {x - 2, y - 2, 4, 4};
        SDL_RenderFillRect(renderer, &rect); //Disegno i nodi
    }

    SDL_RenderPresent(renderer); // Disegna su schermo
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input_file\n", argv[0]);
        return 1;
    }

    char* file_name = argv[1];
    FILE* input = fopen(file_name, "r");
    if (!input) {
        fprintf(stderr, "Error: Cannot open input file %s\n", file_name);
        return 1;
    }

    // Lettura del file e creazione del grafo
    SimpleGraph graph;
    char line[100]; // Assumiamo che ogni riga del file non superi i 100 caratteri

    // Contiamo il numero di nodi (righe nel file)
    graph.node_count = 0;
    while (fgets(line, sizeof(line), input)) {
        graph.node_count++;
    }

    // Allocazione dei nodi
    graph.nodes = malloc(graph.node_count * sizeof(Node));

    // Riavvolgiamo il file per iniziare da capo
    rewind(input);

    // Lettura delle coordinate dei nodi
    size_t node_index = 0;
    while (fgets(line, sizeof(line), input)) {
        double x, y;
        sscanf(line, "Nodo %*zu: x = %lf, y = %lf", &x, &y);
        graph.nodes[node_index].x = x;
        graph.nodes[node_index].y = y;
        node_index++;
    }

    fclose(input);

    // Creazione di una finestra SDL per visualizzare il grafo
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL_Init Error: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Window* window = SDL_CreateWindow("Graph Visualization", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                                          SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
    if (window == NULL) {
        fprintf(stderr, "SDL_CreateWindow Error: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (renderer == NULL) {
        fprintf(stderr, "SDL_CreateRenderer Error: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    // Ciclo principale di gestione degli eventi e disegno
    SDL_Event event;
    int quit = 0;
    while (!quit) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = 1;
            }
        }

        // Disegna il grafo sulla finestra
        drawGraph(renderer, &graph, SCREEN_WIDTH, SCREEN_HEIGHT);

        SDL_Delay(10); // Delay per il controllo degli eventi
    }

    // Cleanup SDL
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    // Deallocazione della memoria
    free(graph.nodes);

    return 0;
}
