/**
 * @file   sudoku.c
 * @author Merlin Nimier-David <merlin.nimier-david@epfl.ch>
 * @author Jean-Cédric Chappelier <jean-cedric.chappelier@epfl.ch>
 *
 * @copyright EPFL 2021
**/
/**
 * @section DESCRIPTION
 *
 * Template du second homework du cours CS-207, année 2021.
**/

// This file requires at least C99 to compile

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h> // isspace()

// ======================================================================
#define INPUT_NO_VALUE '.'
#define INITIAL_LINE_READ_SIZE 32

// ======================================================================

typedef struct Cell Cell;

typedef const Cell** Region;

struct Cell {
    Region constraints;
    char* candidates;
    size_t nb_candidates;
};

typedef struct {
    Cell* content;
    const size_t n;
}Sudoku;


// ======================================================================
// Some function prototypes

void release(Sudoku* s);
void fprint(FILE* stream, const Sudoku* s);
#define print(X) fprint(stdout, X)

void copy_cell(const Sudoku* from, Sudoku* to, size_t ci, size_t cj);

// ======================================================================
Cell* get_cell(const Sudoku* s, size_t i, size_t j)
{
    if (s == NULL) return NULL;

    if(i < 0 && i >= s->n && j < 0 && j >= s->n) return NULL;
    return &s->content[(i*s->n)+j];
}

/* ======================================================================
 * Upper bound on the number of stored contrainsts:
 * 3 x (n-1) :
 *   a. Row    except the cell itself (thus n-1)
 *   b. Column except the cell itself (thus n-1)
 *   c. Region except the cell itself (thus n-1)
 *      The region will most likely overlap with the rows and columns,
 *       but for the sake of simplicity, we allow duplicate constraints.
 */
#define number_of_constraints(s) (3 * ((s)->n - 1))

// ======================================================================
/**
 * Prints all possible values for each cell, one cell per line.
 * This function is provided for debuging purposes.
 */
void print_states(const Sudoku* s)
{
    for (size_t i = 0; i < s->n; ++i) {
        printf("Row %zu:\n", i);
        for (size_t j = 0; j < s->n; ++j) {
            const Cell * const cell = get_cell(s, i, j);
            printf("    Cell (%zu, %zu): %zu candidates: \"%s\"\n",
                   i + 1, j + 1, cell->nb_candidates, cell->candidates);
        }
    }
}

// ======================================================================
/**
 * Tool function to check candidate value k of cell (i,j) against its contraints.
 * i and j are provided only for information.
 */
int check_candidates_vs_constraints(const Sudoku* s,
                                    size_t i, size_t j, size_t k)
{
    const size_t nc = number_of_constraints(s);
    const Cell * const cell = get_cell(s, i, j);
    const char c = cell->candidates[k];

    for (size_t l = 0; l < nc; ++l) {
        const Cell * const other = cell->constraints[l];
        if (other->nb_candidates == 1 && other->candidates[0] == c) {
            fprintf(stderr, "Candidate %zu of cell (%zu, %zu) has value %c, "
                            "which conflicts with its %zu-th neighbor.\n",
                    k + 1, i + 1, j + 1, c, l + 1);
            return 0;
        }
    }
    return 1;
}

// ======================================================================
/**
 * Checks all the values (and candidate values) in the cells.
 * Returns 0 (and prints error message) if some values are wrong.
 */
int check_values(const Sudoku* s)
{
    if (s == NULL) return 1; // If no Sudoku, then all the values are correct ;-)

    const char last = '1' + (char) s->n;

    for (size_t i = 0; i < s->n; ++i) {
        for (size_t j = 0; j < s->n; ++j) {
            const Cell * const cell = get_cell(s, i, j);

            for (size_t k = 0; k < cell->nb_candidates; ++k) {
                const char c = cell->candidates[k];
                // 1. All candidates should be in the alphabet
                if (c < '1' || c > last) {
                    fprintf(stderr, "Candidate %zu of cell (%zu, %zu) has "
                                    "invalid candidate: %c\n",
                            k + 1, i + 1, j + 1, c);
                    fprint(stderr, s);
                    return 0;
                }

                // 2. None of the candidates should conflict
                if (!check_candidates_vs_constraints(s, i, j, k)) {
                    fprintf(stderr, "Constraints for cell (%zu, %zu) candidate %zu"
                                    "not satisfied!\n",
                            i + 1, j + 1, k + 1);
                    fprint(stderr, s);
                    return 0;
                }
            }
        }
    }
    return 1;
}

// ======================================================================
/**
 * Test whether a Sudoku is solved or not.
 */
int is_solved(const Sudoku* s)
{
    if (s == NULL) return 1; // no Sudoku is a solved Sudoku ;-)

    for (size_t i = 0; i < s->n; ++i) {
        for (size_t j = 0; j < s->n; ++j) {
            if (get_cell(s, i, j)->nb_candidates != 1)
                return 0;
        }
    }
    return 1;
}

// ======================================================================
/**
 * Print a Sudoku (on stream).
 */
void fprint(FILE* stream, const Sudoku* s)
{
    if (s == NULL || s->content == NULL)
        return;

    const size_t n = s->n;
    fprintf(stream, "╔");
    for (size_t i = 0; i < 2 * n + 1; ++i)
        fprintf(stream, "═");
    fprintf(stream, "╗\n");

    // N'affiche que les cells dont le nb_values est à 1
    for (size_t i = 0; i < n; ++i) {
        fprintf(stream, "║ ");
        for (size_t j = 0; j < n; ++j) {
            const Cell * const cell = get_cell(s, i, j);
            char c = 'x'; // this will be printed in case of lack of candidates (error)
            if (cell->nb_candidates > 1)
                c = INPUT_NO_VALUE;
            else if (cell->nb_candidates == 1)
                c = cell->candidates[0];
            fprintf(stream, "%c ", c);
        }
        fprintf(stream, "║\n");
    }
    fprintf(stream, "╚");
    for (size_t i = 0; i < 2 * n + 1; ++i)
        fprintf(stream, "═");
    fprintf(stream, "╝\n");
}

// ======================================================================

/**
 * round: for a given index and grid size, computes the index of the
 * corresponding region's top left corner in the case of a regular Sudoku.
 * For instance, round(
 */
#define round(x, y) ((y) * ((x)/(y)))


// ======================================================================
int set_cell_constraints(Sudoku* s, size_t ci, size_t cj,
                         const size_t *custom_regions)
{


    if(s == NULL) return 0;
    if(ci < 0 && ci >= s->n && cj < 0 && cj >= s->n) return 0;
    Region reg = calloc(number_of_constraints(s), sizeof(Cell*));
    if(reg == NULL){
        free(reg);
        return 0;
    }

    size_t z = 0;

    for (int i = 0; i < s->n; ++i) {
        if(i != cj)
            reg[z++] = get_cell(s, ci, i);
    }
    for (int i = 0; i < s->n; ++i) {
        if(i != ci)
            reg[z++] = get_cell(s, i, cj);
    }

    for (size_t i = 0; i < s->n; ++i) {
        for (size_t j = 0; j < s->n; ++j) {
            if(custom_regions == NULL){
                if(round( i, (size_t)sqrt(s->n)) == round( ci, (size_t)sqrt(s->n)) &&
                   round( j, (size_t)sqrt(s->n)) == round( cj, (size_t)sqrt(s->n)))
                    if(i != ci || j != cj)
                        reg[z++] = get_cell(s, i, j);
            }else{
                if(custom_regions[(i*s->n)+j] == custom_regions[(ci*s->n)+cj])
                    if(i != ci || j != cj)
                        reg[z++] = get_cell(s, i, j);
            }
        }
    }

    get_cell(s, ci, cj)->constraints = reg;

    return 1;
}

// ======================================================================
/**
 * Check whether c is a valid value for the Sudoku of size n.
 */
int is_valid_char(char c, size_t n)
{
    return (c >= '1') && (c - '1' < (int) n);
}

// ======================================================================
/**
 * Creates a special grid with arbitrary size
 * and (optionally) custom regions.
 * Example: https://fr.wikipedia.org/wiki/Sudoku#/media/Fichier:Sudoku-irregulier-by-skratt.png
 */
int create(Sudoku* s, size_t n,
           const char** initial_values, const size_t* regions)
{
    if (s == NULL) return 0;


    if(initial_values == NULL) return 0;
    s->content = calloc(n*n, sizeof(Cell));
    if(s->content == NULL){
        release(s);
        return 0;
    }

    size_t* ptr = &s->n;
    *ptr = n;

    size_t z;
    for (size_t i = 0; i < n; ++i) {
        z = 0;
        for (size_t j = 0; j < n; ++j) {
            while(isspace(initial_values[i][z])){
                z++;
            }
            get_cell(s, i, j)->candidates = calloc(n + 1, sizeof(char));
            if(get_cell(s, i, j)->candidates == NULL){
                release(s);
                return 0;
            }
            get_cell(s, i, j)->candidates[n] = '\0';
            if(is_valid_char(initial_values[i][z], n) == 1){
                get_cell(s, i, j)->candidates[0] = initial_values[i][z];
                get_cell(s, i, j)->nb_candidates = 1;
            }else{
                for (size_t k = 0; k < n; ++k) {
                    get_cell(s, i, j)->candidates[k] = '0' + (char) (k+1);
                }
                get_cell(s, i, j)->nb_candidates = n;
            }
            set_cell_constraints(s, i, j, regions);
            z++;
        }
    }
    return 1;
}

// ======================================================================
/**
 * Backs up the candidates of all cells, but does
 * not do a deep copy of the constraints (pointers) since they
 * don't need to be mutated during solve.
 */
int create_shallow_backup(const Sudoku* s, Sudoku* out)
{

    if (s == NULL)   return 0;

    out->content = calloc(s->n*s->n, sizeof(Cell));
    if(out->content == NULL){
        release(out);
        return 0;
    }

    size_t *ptr = &out->n;
    *ptr = s->n;

    for (size_t i = 0; i < s->n; ++i) {
        for (size_t j = 0; j < s->n; ++j) {
            get_cell(out, i, j)->candidates = calloc(out->n + 1, sizeof(char));
            if(get_cell(out, i, j)->candidates == NULL){
                release(out);
                return 0;
            }
            copy_cell(s, out, i, j);
            get_cell(out, i, j)->candidates[out->n] = '\0';
        }
    }
    return 1;
}

// ======================================================================
void restore_shallow_backup(const Sudoku* from, Sudoku* to)
{

    if (from == NULL)   return;
    if (to == NULL)   return;

    for (size_t i = 0; i < from->n; ++i) {
        for (size_t j = 0; j < from->n; ++j) {
            copy_cell(from, to, i, j);
        }
    }
}
/**
 * Copies the content of the from cell to the to cell
 */
void copy_cell(const Sudoku* from, Sudoku* to, size_t ci, size_t cj){
    if (from == NULL) return;
    if(to == NULL) return;

    strcpy(get_cell(to, ci, cj)->candidates, get_cell(from, ci, cj)->candidates);
    get_cell(to, ci, cj)->nb_candidates = get_cell(from, ci, cj)->nb_candidates;
}

// ======================================================================
void release(Sudoku* s)
{
    if (s == NULL || s->content == NULL)
        return;

    for (size_t i = 0; i < s->n; ++i) {
        for (size_t j = 0; j < s->n; ++j) {
            Cell* const cell = get_cell(s, i, j);
            free(cell->constraints); // free is protected against NULL
            free(cell->candidates);
        }
    }
    free(s->content);
    s->content = NULL;
}

// ======================================================================
int update_candidates(Cell* cell, const Cell* constraint, int* has_changed)
{

    if(cell == NULL || constraint == NULL || has_changed == NULL) return 0;
    if(constraint->nb_candidates == 1){
        int j = 0;
        int found = 0;
        for (int i = 0; i < cell->nb_candidates; ++i) {
            if(cell->candidates[i] == constraint->candidates[0]){
                found = 1;
                *has_changed = 1;
            }else{
                cell->candidates[j++] = cell->candidates[i];
            }
        }
        if(found == 1){
            cell->nb_candidates --;
        }
        cell->candidates[cell->nb_candidates] = '\0';
        if(cell->nb_candidates == 0){
            return 0;
        }
    }

    return 1;
}

// ======================================================================
Cell* cell_with_fewest_candidates(const Sudoku* s)
{
    Cell* best_c = NULL;
    size_t best_n = (size_t) -1;
    for (size_t i = 0; i < s->n; ++i) {
        for (size_t j = 0; j < s->n; ++j) {
            Cell* const cell = get_cell(s, i, j);
            const size_t nc = cell->nb_candidates;
            if (nc < best_n && nc > 1) {
                best_c = cell;
                best_n = cell->nb_candidates;
            }
        }
    }
    return best_c;
}

// ======================================================================
/* Simple backtracking solving algorithm
 */
int backtracking_solver(Sudoku* s, size_t depth)
{
    if (s == NULL)    return 1;
    if (is_solved(s)) return 1;

    if (depth >= s->n * s->n + 1) {
        fprintf(stderr, "backtracking_solver: exceeded max depth: %zu\n", depth);
        return 0;
    }

    const size_t n = s->n;
    const size_t nc = number_of_constraints(s);

    // --- Phase 1: apply constraints until the grid no longer evolves

    int has_changed = 0;
    do {
        has_changed = 0;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                Cell* const cell = get_cell(s, i, j);
                for (size_t k = 0; k < nc; ++k) {
                    const int success = update_candidates(cell, cell->constraints[k], &has_changed);
                    if (!success) {
                        // Failure to update candidates, thus indicates an unsolvable grid
                        return 0;
                    }
                }
            }
        }
    } while (has_changed);

    if (is_solved(s)) {
        return check_values(s);
    }

    /* --- Phase 2: we should now be able to start solving
     *     the grid, possibly by fixing a value arbitrarily. */

    Cell* const best_cell = cell_with_fewest_candidates(s);
    if (best_cell == NULL || best_cell->nb_candidates == 0)
        return 0;

    Sudoku backup;
    if (!create_shallow_backup(s, &backup)) {
        fprintf(stderr, "Failed to create shallow backup at depth %zu\n", depth);
        return 1;
    }

    /* Try each candidate in turn. Each time, attempt solving
     * recursively and backtrack if the solve failed.          */
    int success = 0;
    for (size_t k = 0; !success && k < best_cell->nb_candidates; ++k) {
        /* Any mutation we do here will be undone by a backup
         * restore in case of failure.                        */
        best_cell->candidates[0] = best_cell->candidates[k];
        best_cell->candidates[1] = '\0';
        best_cell->nb_candidates = 1;

        success |= backtracking_solver(s, depth + 1);
        if (success) {
            success = check_values(s);
        } else {
            restore_shallow_backup(&backup, s);
        }
    }

    release(&backup);
    return success;
}

// ======================================================================
char* realloc_or_free(char* p, size_t s)
{
    char* const new = realloc(p, s);
    if (new == NULL) {
        free(p);
    }
    return new;
}

// ======================================================================
// Allocates `where` and reads a line into it (from stdin).
// Reads the newline character but does not add it to `where` content.
// `how_many` is set to the number of character read, without newline,
// i.e. *how_many == strlen(*where).
int get_line(char** where, size_t* how_many)
{

    size_t nb = INITIAL_LINE_READ_SIZE;
    char* buffer = calloc(nb + 1, sizeof(char));
    if(buffer == NULL){
        free(buffer);
        return 0;
    }
    buffer[nb] = '\0';
    char in;
    size_t count = 0;
    fflush(stdin);
    do{
        if(count >= nb){
            nb *= 2;
            buffer = realloc_or_free(buffer, nb + 1);
            buffer[nb] = '\0';
        }
        in = (char) getc(stdin);
        buffer[count] = in;
        count++;
    }while(in != '\n');

    *how_many = count - 1;
    *where = buffer;

    return 1;
}

// ======================================================================
// Tool function
void free_tmp_inputs(size_t n, char** initial_values,
                     size_t* regions, char* tmp)
{
    if (initial_values != NULL) {
        for (size_t i = 0; i < n; ++i)
            free(initial_values[i]); // Okay even if already NULL
        free(initial_values);
    }
    if (regions != NULL)
        free(regions);
    free(tmp);
}

// ======================================================================
/**
 * Reads a grid to solve from standard input.
 * Expected format:
 *   <grid width n> <0|1: whether the grid has custom regions>
 *   <n lines with the initial values>
 * And if there are custom regions:
 *   <n lines, one region index per grid cell
 */
int read_grid(Sudoku* s)
{
    if (s == NULL) return 0;

    size_t n = 0;
    if (scanf("%zu", &n) != 1 || n < 2 || n > 9) {
        fprintf(stderr, "Could not read a valid grid size, expected n between 2 and 9\n");
        return 0;
    }

    int has_custom_regions = 0;
    if (scanf("%d", &has_custom_regions) != 1) {
        fprintf(stderr, "Could not read `has_custom_regions`\n");
        return 0;
    }
    if (!has_custom_regions && n != 4 && n != 9) {
        fprintf(stderr, "It's impossible to have standard regions without a square number\n");
        return 0;
    }

    char** initial_values = calloc(n, sizeof(char *));
    if (!initial_values) {
        fprintf(stderr, "Could not allocate `initial_values`\n");
        return 0;
    }

    size_t* regions = NULL;
    size_t i = 0;
    while (i < n) {
        char* dest = NULL;
        size_t line_size = 0;
        if (!get_line(&dest, &line_size)) {
            fprintf(stderr, "Could not read line %zu of `initial_values`\n", i + 1);
            free_tmp_inputs(n, initial_values, regions, dest);
            return 0;
        }
        if (line_size == 0) {
            // Skip empty lines
            free(dest);
        } else {
            if (line_size < n) {
                fprintf(stderr, "Line %zu of `initial_values` has invalid"
                                " length (expected >= %lu, found %lu): \"%s\"\n",
                        i + 1, n, line_size, dest);
                free_tmp_inputs(n, initial_values, regions, dest);
                return 0;
            }
            initial_values[i++] = dest;
        }
    }

    int success = 0;
    if (has_custom_regions) {
        const size_t n2 = n * n;
        regions = calloc(n2, sizeof(size_t));
        if (regions == NULL) {
            fprintf(stderr, "Could not allocate `regions`\n");
            free_tmp_inputs(n, initial_values, regions, NULL);
            return 0;
        }
        for (size_t j = 0; j < n2; ++j) {
            if (scanf("%zu", &regions[j]) != 1) {
                fprintf(stderr, "Could not read %zu-th region\n", j + 1);
                free_tmp_inputs(n, initial_values, regions, NULL);
                return 0;
            }
        }
    }
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wincompatible-pointer-types"

    success = create(s, n, initial_values, regions);

#pragma GCC diagnostic pop

    free_tmp_inputs(n, initial_values, regions, NULL);
    return success;
}

// ======================================================================
int main(void)
{
    Sudoku s;
    if (!read_grid(&s)) {
        fprintf(stderr, "Failed to read or create grid!\n");
        return EXIT_FAILURE;
    }

    printf("~~~~~~ Initial ~~~~~~\n");
    print(&s);
    if (!backtracking_solver(&s, 0))
        fprintf(stderr, "The grid could not be solved!\n");
    else
        printf("~~~*~~ Solved! ~~*~~~\n");
    print(&s);
    release(&s);

    return EXIT_SUCCESS;
}
