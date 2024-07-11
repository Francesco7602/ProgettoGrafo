#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <SDL2/SDL.h>
#include <signal.h>



typedef struct {
    double x;
    double y;
} Node;

typedef struct {
    size_t start;
    size_t end;
} Edge;


// Varibili globali
int peso;
int node_count = 0;
int edge_count= 0;
volatile sig_atomic_t stop = 0; // Variabile volatile per indicare l'interruzione
double scaleFactor; // Fattore di scala iniziale
double offsetX; // Offeset iniziale
double offsetY;



void handle_signal(int sig) {
    if (sig == SIGINT) {
        stop = 1;
    }
}

// Prototipi delle funzioni
void handle_signal(int sig);
void drawGraph(SDL_Renderer* renderer, Node *nodes, Edge *edges, int screenWidth, int screenHeight);
void getMaxNodeDimensions(double* maxX, double* maxY,  Node *nodes);


// Funzione principale
int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Uso corretto: %s nome_file\n", argv[0]);
        return 1; // Esce dal programma con codice di errore
    }
    char* file_name = argv[1];
    FILE* input = fopen(file_name, "r");
    if (!input) {
        printf("Il file non esiste.\n");
        return 1;
    }

    // Lettura del numero di nodi e peso (che non viene utilizzato qui)
    fscanf(input, "%d %*d", &node_count);

    // Allocazione dinamica degli array di nodi
    Node *nodes = malloc(node_count * sizeof(Node));

    // Lettura delle coordinate dei nodi
    for (int i = 0; i < node_count; i++) {
        fscanf(input, "%*d %lf %lf\n", &nodes[i].x, &nodes[i].y);
    }

    // Lettura del numero di archi
    fscanf(input, "%d", &edge_count);

    // Allocazione dinamica degli array di archi
    Edge *edges = malloc(edge_count * sizeof(Edge));

    // Lettura degli indici degli archi
    for (int i = 0; i < edge_count; i++) {
        fscanf(input, "%zu %zu\n", &edges[i].start, &edges[i].end);
    }

    fclose(input);




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
    getMaxNodeDimensions(&maX, &maY, nodes);// Calcola le dimensioni massime dei nodi per adattare le dimensioni della finestra
    scaleFactor = fmin(screenWidth / (2 * maX), screenHeight / (2 * maY));
    offsetX = 0.0;
    offsetY = 0.0;
    //drawGraph(renderer, &graph, screenWidth, screenHeight); // Disegno il grafico
    signal(SIGINT, handle_signal); // Gestore di segnale, chiama handle_signal quando riceve SIGINT
    SDL_Event e; // Variabile per ricevere un evento
    int quit = 0;
    int i = 0;
    while (!quit && !stop) {
        while (SDL_PollEvent(&e) != 0) { // Controllo se ci sono eventi
            if (e.type == SDL_QUIT) {
                quit = 1;
            } else if (e.type == SDL_MOUSEWHEEL) { // Scorrere con la rotella del mouse
                // Zoom in/out
                if (e.wheel.y > 0) { // Scroll up
                    scaleFactor *= 1.1;
                } else if (e.wheel.y < 0) { // Scroll down
                    scaleFactor /= 1.1;
                }
                //drawGraph(renderer, screenWidth, screenHeight); //Aggiorno il grafico
            }else if (e.type == SDL_FINGERMOTION){ // Scroll da touchpad
                if (e.tfinger.fingerId == 0 && e.type == SDL_FINGERMOTION) {
                    // Calcola la differenza di posizione tra due dita
                    double scrollDistance = sqrt(pow(e.tfinger.dx, 2) + pow(e.tfinger.dy, 2));
                    // Aggiorna scaleFactor in base a scrollDistance e alla direzione del movimento
                    scaleFactor *= scrollDistance > 0 ? (1.0 + scrollDistance * 0.01) : 1.0; // Aumento o diminuzione dello scale factor in base alla direzione
                    //drawGraph(renderer, screenWidth, screenHeight);
                }

            }
            else if (e.type == SDL_KEYDOWN) {
                // Gestione del movimento con i tasti freccia
                switch (e.key.keysym.sym) {
                    case SDLK_LEFT:
                        offsetX += 500.0* (1/scaleFactor); // Spostamento a sinistra di 10 pixel
                        //drawGraph(renderer, screenWidth, screenHeight); //Aggiorno il grafico
                        break;
                    case SDLK_RIGHT:
                        offsetX -= 500.0* (1/scaleFactor); // Spostamento a destra di 10 pixel
                        //drawGraph(renderer, screenWidth, screenHeight); //Aggiorno il grafico
                        break;
                    case SDLK_UP:
                        offsetY += 500.0* (1/scaleFactor); // Spostamento verso l'alto di 10 pixel
                        //drawGraph(renderer, screenWidth, screenHeight); //Aggiorno il grafico
                        break;
                    case SDLK_DOWN:
                        offsetY -= 500.0* (1/scaleFactor); // Spostamento verso il basso di 10 pixel
                        //drawGraph(renderer, screenWidth, screenHeight); //Aggiorno il grafico
                        break;
                    default:
                        break;
                }
            }
        }
        drawGraph(renderer,nodes,edges,screenWidth, screenHeight);
    }

    // Cleanup SDL
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
void drawGraph(SDL_Renderer* renderer, Node *nodes, Edge *edges, int screenWidth, int screenHeight) {
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255); // Colore bianco
    SDL_RenderClear(renderer); // Pulisce lo schermo
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); // Colore nero per gli archi
    for (size_t i = 0; i < edge_count; i++) {
        size_t node1 = edges[i].start;
        size_t node2 = edges[i].end;
        SDL_RenderDrawLine(renderer,
                           (int)((nodes[node1].x + offsetX) * scaleFactor + screenWidth / 2),
                           (int)((nodes[node1].y + offsetY) * scaleFactor + screenHeight / 2),
                           (int)((nodes[node2].x + offsetX) * scaleFactor + screenWidth / 2),
                           (int)((nodes[node2].y + offsetY) * scaleFactor + screenHeight / 2)); //Disegna la linea tra due punti
    }

    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255); // Colore rosso per i nodi
    for (size_t i = 0; i < node_count; i++) {
        int x = (int)((nodes[i].x + offsetX) * scaleFactor + screenWidth / 2);
        int y = (int)((nodes[i].y + offsetY) * scaleFactor + screenHeight / 2);
        SDL_Rect rect = {x - (2* scaleFactor), y - (2* scaleFactor), 4* scaleFactor, 4* scaleFactor};
        SDL_RenderFillRect(renderer, &rect); //Disegno i nodi
    }
    SDL_RenderPresent(renderer); // Disegna su schermo
}


void getMaxNodeDimensions(double* maxX, double* maxY, Node *nodes) { // Trova le coordinate massime dei nodi, in valore assoluto
    *maxX = 0.0;
    *maxY = 0.0;
    for (size_t i = 0; i < node_count; i++) {
        if (fabs(nodes[i].x) > *maxX) {
            *maxX = fabs(nodes[i].x);
        }
        if (fabs(nodes[i].y) > *maxY) {
            *maxY = fabs(nodes[i].y);
        }
    }
}
