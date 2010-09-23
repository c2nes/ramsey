
/**
 * Purpose: Given a 42-vertex complete graph with red/blue coloring such that no
 *  monocromatic 5-clique exists, attempt to construct a 43-vertex complete graph
 *  also containing no monocromatic 5-clique.
 *
 * Theory of operation: First the 42-vertex graph is loaded and all
 *  monochromatic 4-cliques are located. Any monochromatic 5-clique in the
 *  43-vetex graph will have one of these 4-cliques as a subgraph. Then the
 *  possible colorings of the 43rd node are permuted over, checking each
 *  coloring against each of the 4-cliques. If a graph is found which produces
 *  no monochromatic 5-cliques, this graph is printed. If every graph is
 *  attepmted without success, the program will simply exit.
 */

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

/* Size of block that permutations of the edge colors are checked in (should be less than 32) */
#define PERM_BLOCK_SIZE 26
#define PERM_BLOCK_HIGH_MASK ((uint32_t)(~((1 << PERM_BLOCK_SIZE) - 1)))

/* Debug flags */
#define SHOW_PERMUTATIONS 0
#define TRACK_MAX_SUCCESS 0

/* Color of each edge (0 -> red, 1 -> blue) */
typedef uint8_t color;

struct Perm_s {
    uint32_t perm;
    struct Perm_s* prev;
    struct Perm_s* next;
    uint32_t __pad;
};
typedef struct Perm_s Perm;

static color** load_matrix(void);
static void dump_graph(color** matrix, int order);
static void print_bin(uint32_t n, uint8_t width);
static color** expand(color** matrix, int n);
static inline bool next_graph(color** matrix, int order);
static inline bool is_monochromatic(uint16_t* clique, color** matrix);

static inline bool next_n_clique(uint16_t* clique, uint16_t size);
static inline bool is_n_monochromatic(uint16_t* clique, int n, color** matrix);
static uint16_t** find_monochromatic_n_cliques(color** matrix, int order, int n, int* cliques_found);

static void perm_alloc(void);
static void perm_init(void);
static inline void perm_remove(Perm* p);
static inline Perm* next_perm_with_mask(Perm* p, uint32_t mask, color cc);
static inline bool perm_mask(uint16_t* clique);
static inline uint32_t perm_next(void);

/* Permutation generator state */
static Perm* perm_block = NULL;
static Perm* perm_head = NULL;
static Perm* perm_reset_head = NULL;
static uint32_t perm_count = 0;
static uint32_t* perm_filtered = NULL;
static uint32_t* perm_filtered_start = NULL;
static uint32_t* perm_filtered_end = NULL;

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

    fclose(f);

    return adj;
}

static void dump_graph(color** matrix, int order) {
    for(int i = 0; i < order; i++) {
        for(int j = 0; j < order; j++) {
            printf("%d", matrix[i][j]);
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

/* Construct the next permutation of the edges of the new node. This only
   permutes the row values of the new node (i.e. the column value shouldn't be
   used). This is an optimization

   O(n), n = order of the matrix
 */
static inline bool next_graph(color** matrix, int order) {
    static color* row = NULL;
    uint32_t p;

    if(row == NULL) {
        row = matrix[order - 1];
    }

    //if(perm_head == NULL) {
    if(perm_filtered == perm_filtered_end) {
        int i = PERM_BLOCK_SIZE;

        while(row[i] && i != order - 1) {
            row[i] = 0;
            i++;
        }

        if(i == order - 1) {
            return false;
        }

        row[i] = 1;

        //perm_init();
        //perm_head = perm_reset_head;
        perm_filtered = perm_filtered_start;
    }

    //p = perm_next();
    //p = perm_head->perm;
    //perm_head = perm_head->next;
    p = *perm_filtered;
    perm_filtered++;
    for(uint32_t i = 0; i < PERM_BLOCK_SIZE; i++) {
        row[i] = (p >> i) & 1;
    }

    return true;
}

/* Check if a clique is monotone. Assumes the clique is presented in least to
   greatest order, and that the first four nodes in the clique form a
   monochromatic clique. Additionally, the color of the clique should be stored
   after the clique values in clique[CLIQUE_N]

   Worse case O(n^2), n = CLIQUE_N 
*/
static inline bool is_monochromatic(uint16_t* clique, color** matrix) {
    static color* row = NULL;
    color cc = (color) clique[CLIQUE_N];

    if(row == NULL) {
        row = matrix[clique[CLIQUE_N - 1]];
    }

    /* Only check the colors of the new edges */
    for(int i = 0; i < CLIQUE_N - 1; i++) {
        if(row[clique[i]] ^ cc) {
            return false;
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
    free(current_clique);

    return cliques;
}

static void print_bin(uint32_t n, uint8_t width) {
    for(int i = width - 1; i >= 0; i--) {
        printf("%d", (n >> i) & 1);
    }
    printf("\n");
}

static void perm_alloc(void) {
    perm_block = malloc(sizeof(Perm) * (1 << PERM_BLOCK_SIZE));

    for(uint32_t i = 0; i < 1 << PERM_BLOCK_SIZE; i++) {
        perm_block[i].perm = i;
    }
}

static void perm_init(void) {
    for(uint32_t i = 0; i < (1 << PERM_BLOCK_SIZE) - 1; i++) {
        perm_block[i].next = &(perm_block[i + 1]);
    }

    for(uint32_t i = 1; i < 1 << PERM_BLOCK_SIZE; i++) {
        perm_block[i].prev = &(perm_block[i - 1]);
    }

    perm_head = perm_block;
    perm_count = 1 << PERM_BLOCK_SIZE;
    perm_block[(1 << PERM_BLOCK_SIZE) - 1].next = NULL;
    perm_block[0].prev = NULL;
}

static inline void perm_remove(Perm* p) {
    if(p->prev) {
        p->prev->next = p->next;
    } else {
        perm_head = p->next;
    }

    if(p->next) {
        p->next->prev = p->prev;
    }

    p->next = p->prev = NULL;
    perm_count--;
}

static inline Perm* next_perm_with_mask(Perm* p, uint32_t mask, color cc) {
    uint32_t x_mask = mask;

    if(cc) {
        x_mask = 0;
    }

    p = p->next;
    while(p && ((p->perm ^ x_mask) & mask) != mask) {
        p = p->next;
    }

    return p;
}

static inline bool perm_mask(uint16_t* clique) {
    color cc = (color) clique[CLIQUE_N];
    uint32_t mask = 0;
    uint32_t x_mask = 0;
    Perm* p;
    Perm* pp;

    for(int i = 0; i < CLIQUE_N - 1; i++) {
        if(clique[i] > PERM_BLOCK_SIZE) {
            return false;
        }

        mask |= (((uint32_t)1) << clique[i]);
    }
    
    if(!cc) {
        x_mask = mask;
    }
    
    p = perm_head;

    if(p == NULL) {
        return false;
    }

    if(((p->perm ^ x_mask) & mask) != mask) {
        p = next_perm_with_mask(p, mask, cc);
    }

    while(p) {
        if(p->next || p->prev || p == perm_head) {
            pp = p->next;
            perm_remove(p);
            p = pp;
        }

        if(p && ((p->perm ^ x_mask) & mask) != mask) {
            p = next_perm_with_mask(p, mask, cc);
        }
    }

    return true;
}

static inline void perm_regroup(void) {
    Perm* p = perm_head;
    uint32_t j = 0;
    
    while(p) {
        memcpy(perm_block + j, p, sizeof(Perm));
        
        if(perm_block[j].next) {
            perm_block[j].next->prev = &perm_block[j];
        }
        if(perm_block[j].prev) {
            perm_block[j].prev->next = &perm_block[j];
        } else {
            perm_head = &perm_block[j];
        }
        
        p = perm_block[j].next;
        j++;
    }
}

static inline uint32_t perm_next(void) {
    Perm* p = perm_head;
    uint32_t v = p->perm;

    perm_head = p->next;
    if(perm_head) {
        perm_head->prev = NULL;
    }

    p->prev = p->next = NULL;

    perm_count--;

    return v;
}

static inline uint32_t clique_to_mask(uint16_t* clique) {
    uint32_t mask = 0;

    for(int i = 0; i < CLIQUE_N - 1; i++) {
        if(clique[i] < 32) {
            mask |= (1 << clique[i]);
        } else {
            return 0;
        }
    }

    return mask;
}

static inline uint16_t overlap_factor(uint32_t mask1, uint32_t mask2) {
    uint32_t mask_overlap;
    uint16_t overlap = 0;
    
    mask_overlap = mask1 & mask2;
    
    while(mask_overlap) {
        overlap += mask_overlap & 1;
        mask_overlap >>= 1;
    }

    return overlap;
}

static inline uint16_t** sort_cliques(uint16_t** cliques, int clique_count) {
    uint16_t** new_cliques = malloc(clique_count * sizeof(uint16_t*));
    uint16_t copied = 0;
    uint32_t mask, temp_mask;
    uint16_t min_over, min_i, over, i, last_color;
    new_cliques[0] = malloc(clique_count * 8 * sizeof(uint16_t));
    
    for(i = 0; i < clique_count; i++) {
        new_cliques[i] = new_cliques[0] + (i * 8);
    }

    memcpy(new_cliques[0], cliques[0], sizeof(uint16_t) * 8);
    cliques[0][6] = 1;
    last_color = cliques[0][5];
    mask = clique_to_mask(cliques[0]);
    copied++;

    while(copied < clique_count) {
        min_i = 0;
        min_over = 1024;

        for(i = 0; i < clique_count; i++) {
            if(cliques[i][6] == 0 && cliques[i][5] != last_color) {
                temp_mask = clique_to_mask(cliques[i]);
                if(temp_mask) {
                    over = overlap_factor(mask, temp_mask);
                    if(over < min_over) {
                        min_over = over;
                        min_i = i;
                    }
                }
            }
        }

        if(min_over < 1024) {
            memcpy(new_cliques[copied], cliques[min_i], sizeof(uint16_t) * 8);
            cliques[min_i][6] = 1;
            mask |= clique_to_mask(cliques[min_i]);
            last_color = cliques[min_i][5];
            copied++;
        } else {
            break;
        }
    }

    for(i = 0; i < clique_count; i++) {
        if(cliques[i][6] == 0) {
            memcpy(new_cliques[copied], cliques[i], sizeof(uint16_t) * 8);
            cliques[i][6] = 1;
            copied++;
        }
    }

    if(copied != clique_count) {
        printf(":(\n");
    }

    return new_cliques;
}

int main(void) {
    /* The adjacency matrix being inspected for mono-chromatic cliques */
    color** matrix;

    /* Order of the matix (i.e. the order of the complete graph) */
    int order = ADJ_MATRIX_ORDER;

    /* Edge check */
    bool monochromatic = true;

    /* Iterator */
    int i;

    /* 4-cliques */
    uint16_t** four_cliques;
    int four_clique_count = 0;

    /* Possible 5-cliques */
    uint16_t** five_cliques;

    /* Number of graphs checked */
    int count = 0;

#if TRACK_MAX_SUCCESS
    /* A permutation which fails, will fail after a set number of cliques being
       checked. This keeps track of the largest number of cliques which had to
       be checked to invalidate a permutation */
    int max = 0;
#endif

    /* Allocate memory to the permutation generator */
    perm_alloc();
    perm_init();

    matrix = load_matrix();
    printf("Successfully loaded matrix\n");

    /* Find all four cliques in the existing graph */
    four_cliques = find_monochromatic_n_cliques(matrix, order, CLIQUE_N - 1, &four_clique_count);
    printf("Found %d 4-cliques\n", four_clique_count);

    /* Construct a list of potential five cliques using the four cliques and the
       new node */
    five_cliques = malloc(sizeof(uint16_t*) * four_clique_count);
    five_cliques[0] = malloc(sizeof(uint16_t) * 8 * four_clique_count);
    for(i = 0; i < four_clique_count; i++) {
        five_cliques[i] = five_cliques[0] + (i * 8);
        memcpy(five_cliques[i], four_cliques[i], sizeof(uint16_t) * 4);
        five_cliques[i][4] = order;
        
        /* Store the color of the clique as well */
        five_cliques[i][5] = matrix[five_cliques[i][0]][five_cliques[i][1]];
        five_cliques[i][6] = 0;
        
        free(four_cliques[i]);
    }
    free(four_cliques);

    /* Expand the matrix */
    matrix = expand(matrix, order);
    order++;

    /* Sort, lol */
#if 0
    five_cliques = sort_cliques(five_cliques, four_clique_count);

    for(i = 0; i < 100; i++) {
        print_bin(((1 << five_cliques[i][0]) |
                   (1 << five_cliques[i][1]) |
                   (1 << five_cliques[i][2]) |
                   (1 << five_cliques[i][3])), 32);
            
    }
#endif

    int s;
    int last_resize = perm_count;
    Perm* p;
    printf("Filtering... ");fflush(stdout);
    s = perm_count;
#if 0
    for(i = 0; i < four_clique_count; i++) {
        if(perm_mask(five_cliques[i])) {
            printf("Removed %.2f%% of permutations\n", 100 * ((float)s - perm_count) / s);
        }

        if(perm_count < (last_resize >> 1) && perm_count > 1024) {
            printf("Regrouping\n");
            last_resize = perm_count;
            perm_regroup();
        }
    }
#else
#define BLOCKY 32
    for(int j = 0; j < BLOCKY; j++) {
        for(i = j; i < four_clique_count; i += BLOCKY) {
            if(perm_mask(five_cliques[i])) {
                printf("Removed %.2f%% of permutations\n", 100 * ((float)s - perm_count) / s);
            }
        }
        printf("Next group\n");
        perm_regroup();
    }
#endif
    printf("Done! Removed %.2f%% of permutations (%u/%u)\n", 100 * ((float)s - perm_count) / s, s - perm_count, s);
    printf("Permutation space: %u\n", (1 << (order - PERM_BLOCK_SIZE)) * perm_count);
    printf("Load factor: %.4f\n",   ((float)((1 << ((order - PERM_BLOCK_SIZE) - 16)) * perm_count)) / (1 << 20) );
    perm_reset_head = perm_head;
    exit(0);

    perm_filtered = malloc(perm_count * sizeof(uint32_t));

    p = perm_head;
    i = 0;
    while(p) {
        perm_filtered[i] = p->perm;
        p = p->next;
        i++;
    }
    perm_filtered_start = perm_filtered;
    perm_filtered_end = perm_filtered + perm_count;

    int large_blocks = 0;

    while(true) {
        /* Attempt to move to the next graph */
        if(next_graph(matrix, order) == false) {
            printf("Exhausted possibilities! No such extension of the current graph\n");
            break;
        }

#if 1
        if(perm_filtered == perm_filtered_end) {
            for(i = order - 1; i > 0; i--) {
                printf("%d", matrix[order - 1][i - 1]);
            }
            printf("\n");
            large_blocks++;
        }
#endif

#if SHOW_PERMUTATIONS
        for(i = order - 1; i > 0; i--) {
            printf("%d", matrix[order - 1][i - 1]);
        }
        printf("\n");
#endif

        /* Check all cliques under this permutation */
        i = 0;
        while(i < four_clique_count && !(monochromatic = is_monochromatic(five_cliques[i], matrix))) {
            i++;
        }

#if TRACK_MAX_SUCCESS
        if(i > max) {
            max = i;
            for(int j = order - 1; j > 0; j--) {
                printf("%d", matrix[order - 1][j - 1]);
            }
            printf(" (%d) \n", i);
        }
#endif

        /* Successfully found a graph with 0 monochromatic cliques */
        if(monochromatic == false) {
            printf("Found clique-less extension: \n\n");
            dump_graph(matrix, order);
            break;
        }

        count++;
    }

    free(matrix[0]);
    free(matrix);
    free(five_cliques[0]);
    free(five_cliques);
    
    return 0;
}
