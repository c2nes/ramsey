
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Size of cliques to find */
#define CLIQUE_N 5

/* File to load matrix to check, and its degree */
#define ADJ_MATRIX_FILE "g55.42"
#define ADJ_MATRIX_ORDER 42

/* Debug flags */
#define DUMP_CLIQUES 1

/* Color of each edge (0 -> red, 1 -> blue) */
typedef uint8_t color;

static color** load_matrix(void) {
    FILE* f;
    color** adj;
    color* adj_flat;
    int c, i;

    adj = malloc(sizeof(color*) * ADJ_MATRIX_ORDER);
    adj_flat = malloc(sizeof(color) * ADJ_MATRIX_ORDER * ADJ_MATRIX_ORDER);
    if(adj == NULL || adj_flat == NULL) {
        perror("Could not alloc");
        exit(EXIT_FAILURE);
    }

    f = fopen(ADJ_MATRIX_FILE, "r");
    if(f == NULL) {
        perror("Could not open adjacency matrix");
        exit(EXIT_FAILURE);
    }

    /* Populate matrix with rows */
    for(i = 0; i < ADJ_MATRIX_ORDER; i++) {
        adj[i] = adj_flat + (i * ADJ_MATRIX_ORDER);
    }

    /* Populate adjacency matrix from file */
    c = 0;
    i = 0;
    while(c != EOF) {
        c = fgetc(f);
        if(c == '0') {
            adj_flat[i++] = 0;
        } else if(c == '1') {
            adj_flat[i++] = 1;
        }
    }

    /* If we exit the loop and haven't read the correct number of values
       "something bad" has happened */
    if(i != (ADJ_MATRIX_ORDER * ADJ_MATRIX_ORDER)) {
        fprintf(stderr, "Error: matrix not fully populated\n");
        exit(EXIT_FAILURE);
    }

    return adj;
}

static inline void swap_rows(color** matrix, int order, int n, int m) {
    color* row_n = NULL;
    color* row_m = NULL;

    if(n == m) {
        return;
    }

    row_n = malloc(order * sizeof(color));
    row_m = malloc(order * sizeof(color));
    memcpy(row_n, matrix[n], order * sizeof(color));
    memcpy(row_m, matrix[m], order * sizeof(color));
    memcpy(matrix[n], row_m, order * sizeof(color));
    memcpy(matrix[m], row_n, order * sizeof(color));

    /*
    for(int i = 0; i < order; i++) {
        matrix[i][m] = row_n[i];
        matrix[i][n] = row_m[i];
        }*/

    free(row_n);
    free(row_m);
}

/* Increment the given clique. If size = CLIQUE_N then next_clique should be
   used because it is much faster */
static inline bool next_n_clique(uint16_t* clique, uint16_t size) {
    int i, j;

    /* Think about AJD_MATRIX_ORDER bins with CLIQUE_N markers. We move the
       right most marker by one to the right which can be moved without change
       the relative order of the markers or causing a marker to leave the last
       bin. When a marker is moved, all the markers to its right are moved into
       the bins to its immediate right */

    for(i = size - 1; i >= 0; i--) {
        /* We found a marker we can move */
        if(clique[i] < clique[i + 1] - 1) {
            /* Move the marker */
            clique[i]++;
            
            /* Move all markers further right to be after the one we just
               moved */
            for(j = i + 1; j < size; j++) {
                clique[j] = clique[j - 1] + 1;
            }

            /* The clique has changed, return true */
            return true;
        }
    }

    /* Nothing could be moved (last combination) */
    return false;
}

/* Is the given n-clique monochromatic in matrix */
static inline bool is_n_monochromatic(uint16_t* clique, int n, color** matrix) {
    color cc = matrix[clique[0]][clique[1]];
    color* row;
    int i, j;

    for(i = 0; i < n; i++) {
        row = matrix[clique[i]];
        for(j = i + 1; j < n; j++) {
            if(row[clique[j]] ^ cc) {
                return false;
            }
        }
    }

    return true;
}

/* Find the monochrome n-cliques in matix. Return the list of n-cliques and
   store the number of cliques found in cliques_found */
static uint16_t** find_monochromatic_n_cliques(color** matrix, int order, int n, int* cliques_found) {
    uint16_t* current_clique = malloc(sizeof(uint16_t) * (n + 1));
    uint16_t** cliques = NULL;
    uint16_t count = 0;
    int i;
    
    current_clique[n] = order;
    for(i = 0; i < n; i++) {
        current_clique[i] = i;
    }

    do {
        if(is_n_monochromatic(current_clique, n, matrix)) {
            count++;
            cliques = realloc(cliques, sizeof(uint16_t*) * count);
            cliques[count - 1] = malloc(sizeof(uint16_t) * n);
            memcpy(cliques[count - 1], current_clique, sizeof(uint16_t) * n);
        }
    } while(next_n_clique(current_clique, n));
    
    *cliques_found = count;
    return cliques;
}

static void dump_graph(color** matrix, int order) {
    for(int i = 0; i < order; i++) {
        for(int j = 0; j < order; j++) {
            printf("%d", matrix[i][j]);
        }
        printf("\n");
    }
}

int main(void) {
    /* The adjacency matrix being inspected for mono-chromatic cliques */
    color** matrix;

    /* Clique count */
    uint16_t** cliques = NULL;
    int count = 0;

    matrix = load_matrix();
    printf("Successfully loaded matrix\n");


    for(int i = 0; i < ADJ_MATRIX_ORDER; i++) {
        for(int j = i; j < ADJ_MATRIX_ORDER; j++) {
            swap_rows(matrix, ADJ_MATRIX_ORDER, i, j);
            dump_graph(matrix, ADJ_MATRIX_ORDER);
            printf("\n");
            cliques = find_monochromatic_n_cliques(matrix, ADJ_MATRIX_ORDER, CLIQUE_N, &count);
            if(count > 0) {
                printf("Found %d %d-cliques\n", count, CLIQUE_N);
            }
            swap_rows(matrix, ADJ_MATRIX_ORDER, i, j);
        }
    }

#if DUMP_CLIQUES
    for(int i = 0; i < count; i++) {
        for(int j = 0; j < CLIQUE_N; j++) {
            printf("%2d ", cliques[i][j]);
        }
        printf("\n");
    }
#endif
    
    return 0;
}
