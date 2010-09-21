
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
#define SHOW_PERMUTATIONS 0

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
        if(i < ADJ_MATRIX_ORDER * ADJ_MATRIX_ORDER) {
            if(c == '0') {
                adj_flat[i++] = 0;
            } else if(c == '1') {
                adj_flat[i++] = 1;
            }
        } else if(c == '0' || c == '1') {
            i++;
            break;
        }
    }

    /* If we exit the loop and haven't read the correct number of values
       "something bad" has happened */
    if(c != EOF || i != (ADJ_MATRIX_ORDER * ADJ_MATRIX_ORDER)) {
        fprintf(stderr, "Error: invalid matrix size\n");
        exit(EXIT_FAILURE);
    }

    return adj;
}

static void dump_graph(color** matrix, int order) {
    for(int i = 0; i < order; i++) {
        for(int j = 0; j < order; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

/* Expand the given n by n matrix to be a n + 1 by n + 1 matrix with the
   original matrix embedded at 0, 0. New values are all initilized to 0 */
static color** expand(color** matrix, int n) {
    color** new_matrix = malloc((n + 1) * sizeof(color*));
    new_matrix[0] = malloc((n + 1) * (n + 1) * sizeof(color));

    for(int i = 0; i < n + 1; i++) {
        /* Reconstruct the pointer array */
        new_matrix[i] = new_matrix[0] + (i * (n + 1));

        /* Zero the last new column value */
        new_matrix[i][n] = 0;
    }
    
    for(int i = n - 1; i >= 0; i--) {
        /* Move each row into its new place */
        memcpy(new_matrix[i], matrix[i], n * sizeof(color));
    }

    /* Set the bottom row to all 0's */
    memset(new_matrix[n], 0, n * sizeof(color));

    /* Free old matrix */
    free(matrix[0]);
    free(matrix);

    return new_matrix;
}

/* Construct the next permutation of the edges of the new node.
 * O(n), n = order of the matrix
 */
static inline bool next_graph(color** matrix, int order) {
    int i = 0;
    color* row = matrix[order - 1];

    while(row[i] == 1 && i < order) {
        row[i] = 0;
        matrix[i][order - 1] = 0;
        
        i++;
    }

    if(i == order) {
        return false;
    }
    
    row[i] = 1;
    matrix[i][order - 1] = 1;

    return true;
}

/* O(n), n = CLIQUE_N */
static inline bool next_clique(uint16_t* state) {
    int i, j;

    /* Think about AJD_MATRIX_ORDER bins with CLIQUE_N markers. We move the
       right most marker by one to the right which can be moved without change
       the relative order of the markers or causing a marker to leave the last
       bin. When a marker is moved, all the markers to its right are moved into
       the bins to its immediate right */

    for(i = CLIQUE_N - 1; i >= 0; i--) {
        /* We found a marker we can move */
        if(state[i] < state[i + 1] - 1) {
            /* Move the marker */
            state[i]++;
            
            /* Move all markers further right to be after the one we just
               moved */
            for(j = i + 1; j < CLIQUE_N; j++) {
                state[j] = state[j - 1] + 1;
            }

            /* The state has changed, return true */
            return true;
        }
    }

    /* Nothing could be moved (last combination) */
    return false;
}

/* Worse case O(n^2), n = CLIQUE_N */
static inline bool is_monochromatic(uint16_t* clique, color** matrix) {
    color cc = matrix[clique[0]][clique[1]];
    color* row;
    int i, j;

    for(i = 0; i < CLIQUE_N; i++) {
        row = matrix[clique[i]];

        /* This only needs to be j = i + 1, but taking j = 0, and using
           -funroll-loops this will be much faster */
        for(j = 0; j < CLIQUE_N; j++) {
            if(row[clique[j]] ^ cc) {
                return false;
            }
        }
    }

    return true;
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

int main(void) {
    /* The adjacency matrix being inspected for mono-chromatic cliques */
    color** matrix;

    /* Order of the matix (i.e. the order of the complete graph) */
    int order = ADJ_MATRIX_ORDER;

    /* The current clique being inspected */
    uint16_t clique[CLIQUE_N + 1];

    /* Edge check */
    bool monochromatic = true;

    /* Iterators */
    int i;

    /* 4-cliques */
    uint16_t** four_cliques;
    int four_clique_count = 0;

    /* Possible 5-cliques */
    uint16_t** five_cliques;

    /* Clique count */
    int count;

    matrix = load_matrix();
    printf("Successfully loaded matrix\n");

    /* Find all four cliques in the existing graph */
    four_cliques = find_monochromatic_n_cliques(matrix, order, CLIQUE_N - 1, &four_clique_count);
    printf("Found %d 4-cliques\n", four_clique_count);

    /* Construct a list of potential five cliques using the four cliques and the
       new node */
    five_cliques = malloc(sizeof(uint16_t*) * four_clique_count);
    five_cliques[0] = malloc(sizeof(uint16_t) * 5 * four_clique_count);
    for(i = 0; i < four_clique_count; i++) {
        five_cliques[i] = five_cliques[0] + (i * 5);
        memcpy(five_cliques[i], four_cliques[i], sizeof(uint16_t) * 4);
        five_cliques[i][4] = 42;
    }

    /* Expand the matrix */
    matrix = expand(matrix, order);
    order++;

    /* This "extra" element in the clique gives the order of the graph which
       next_clique uses so that its algorithm is more consistent */
    clique[CLIQUE_N] = order;

    count = 0;

    int max = 0;
    while(true) {
#if SHOW_PERMUTATIONS
        for(i = order; i > 0; i--) {
            printf("%d", matrix[order - 1][i - 1]);
        }
        printf("\n");
#endif

        i = 0;
        while(i < four_clique_count && !(monochromatic = is_monochromatic(five_cliques[i], matrix))) {
            i++;
        }
        if(i > max) {
            max = i;
            for(i = order; i > 0; i--) {
                printf("%d", matrix[order - 1][i - 1]);
            }
            printf("\n");
        }

        /* Successfully found a graph with 0 monochromatic cliques */
        if(!monochromatic) {
            break;
        }

        if(!next_graph(matrix, order) || count == (1 << 22)) {
            printf("Exhausted possibilities! No such extension of the current graph\n");
            printf("%d\n", max);
            exit(1);
        }

        count++;
    }

    /* Found a 0-clique matix! Woot! */
    printf("Found clique-less extension: \n\n");
    dump_graph(matrix, order);
    
    return 0;
}
