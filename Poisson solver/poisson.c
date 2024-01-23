#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>

/**
 * poisson.c
 * Implementation of a Poisson solver with Dirichlet boundary conditions.
 *
 * BUILDING:
 * gcc -o poisson poisson.c -lpthread
 * 
 * to enable profiling
 * g++ -pg -g -o poisson poisson.c -lpthread
 * [note: linking pthread isn't strictly needed until you add your
 *        multithreading code]
 */

// Define to enable timer based profiling
#define PROFILE

// Set to true when operating in debug mode to enable verbose logging
static bool debug = false;

typedef struct
{
    int thread_id;      // Unique id of the worker thread
    int start;          // Start index of the worker thread
    int n;              // End index of the worker thread
    float delta;
    double *curr;
    double *next;
    double *source;
} ThreadArgs;

void* TopLayer(void* pargs)
{
    ThreadArgs* args = (ThreadArgs*)pargs;
    int n = args->n;
    int n2 = n * n;
    int k = args->start;
    double d2 = args->delta * args->delta;

    double *curr = args->curr;
    double *next = args->next;
    double *source = args->source;

    int idx = n2 * k - 1;

    int idx_n;
    int idx_s;
    int idx_e;
    int idx_w;
    int idx_d;

    double ce;
    double cw;
    double cn;
    double cs;
    double cd;

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            // Calculate common indices (may make it easier for compilers to optimize the code)
            idx++;
            idx_n = idx + n;
            idx_s = idx - n;
            idx_e = idx + 1;
            idx_w = idx - 1;
            idx_d = idx - n2;

            // Calculate common values (may make it easier for compilers to optimize the code)
            ce = curr[idx_e];
            cw = curr[idx_w];
            cn = curr[idx_n];
            cs = curr[idx_s];
            cd = curr[idx_d];

            
            if (j == 0) {
                if (i == 0) {
                    next[idx] = (2 * ce + 2 * cn + 2 * cd - d2 * source[idx]) / 6.0;
                } else if (i == n - 1) {
                    next[idx] = (2 * cw + 2 * cn + 2 * cd - d2 * source[idx]) / 6.0;
                } else {
                    next[idx] = (ce + cw + 2 * cn + 2 * cd - d2 * source[idx]) / 6.0;
                }
            } else if (j == n - 1) {
                if (i == 0) {
                    next[idx] = (2 * ce + 2 * cs + 2 * cd - d2 * source[idx]) / 6.0;
                } else if (i == n - 1) {
                    next[idx] = (2 * cw + 2 * cs + 2 * cd - d2 * source[idx]) / 6.0;
                } else {
                    next[idx] = (ce + cw + 2 * cs + 2 * cd - d2 * source[idx]) / 6.0;
                }
            } else {
                if (i == 0) {
                    next[idx] = (2 * ce + cs + cn + 2 * cd - d2 * source[idx]) / 6.0;
                } else if (i == n - 1) {
                    next[idx] = (2 * cw + cs + cn + 2 * cd - d2 * source[idx]) / 6.0;
                } else {
                    next[idx] = (ce + cw + cs + cn + 2 * cd - d2 * source[idx]) / 6.0;
                }
            }

        }
    }
    return NULL;
}

void* MiddleLayer(void* pargs)
{
    ThreadArgs* args = (ThreadArgs*)pargs;
    int n = args->n;
    int n2 = n * n;
    int k = args->start;
    double d2 = args->delta * args->delta;

    double *curr = args->curr;
    double *next = args->next;
    double *source = args->source;

    int idx = n2 * k - 1;

    double ce;
    double cw;
    double cn;
    double cs;
    double cu;
    double cd;

    int idx_n;
    int idx_s;
    int idx_e;
    int idx_w;
    int idx_u;
    int idx_d;

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            // Calculate common indices (may make it easier for compilers to optimize the code)
            idx++;
            idx_n = idx + n;
            idx_s = idx - n;
            idx_e = idx + 1;
            idx_w = idx - 1;
            idx_u = idx + n2;
            idx_d = idx - n2;

            // Calculate common values (may make it easier for compilers to optimize the code)
            ce = curr[idx_e];
            cw = curr[idx_w];
            cn = curr[idx_n];
            cs = curr[idx_s];
            cu = curr[idx_u];
            cd = curr[idx_d];

            if (j == 0) {
                if (i == 0) {
                    next[idx] = (2 * ce + 2 * cn + cu + cd - d2 * source[idx]) / 6.0;
                } else if (i == n - 1) {
                    next[idx] = (2 * cw + 2 * cn + cu + cd - d2 * source[idx]) / 6.0;
                } else {
                    next[idx] = (ce + cw + 2 * cn + cu + cd - d2 * source[idx]) / 6.0;
                }
            } else if (j == n - 1) {
                if (i == 0) {
                    next[idx] = (2 * ce + 2 * cs + cu + cd - d2 * source[idx]) / 6.0;
                } else if (i == n - 1) {
                    next[idx] = (2 * cw + 2 * cs + cu + cd - d2 * source[idx]) / 6.0;
                } else {
                    next[idx] = (ce + cw + 2 * cs + cu + cd - d2 * source[idx]) / 6.0;
                }
            } else {
                if (i == 0) {
                    next[idx] = (2 * ce + cs + cn + cu + cd - d2 * source[idx]) / 6.0;
                } else if (i == n - 1) {
                    next[idx] = (2 * cw + cs + cn + cu + cd - d2 * source[idx]) / 6.0;
                } else {
                    next[idx] = (ce + cw + cs + cn + cu + cd - d2 * source[idx]) / 6.0;
                }
            }
        }
    }
    return NULL;
}

/**
 * @brief Solve Poissons equation for a given cube with Dirichlet boundary
 * conditions on all sides.
 *
 * @param n             The edge length of the cube. n^3 number of elements.
 * @param source        Pointer to the source term cube, a.k.a. forcing function.
 * @param iterations    Number of iterations to perform.
 * @param threads       Number of thread to use for solving.
 * @param delta         Grid spacing.
 * @return double*      Solution to Poissons equation.  Caller must free.
 */
double* poisson_dirichlet (int n, double *source, int iterations, int threads, float delta)
{

    threads = n-1;

    if (debug)
    {
        printf ("Starting solver with:\n"
               "n = %i\n"
               "iterations = %i\n"
               "threads = %i\n"
               "delta = %f\n",
               n, iterations, threads, delta);
    }
    
    // Allocate some buffers to calculate the solution in
    double *curr = (double*)calloc (n * n * n, sizeof (double));
    double *next = (double*)calloc (n * n * n, sizeof (double));

    // Ensure we haven't runPROFILE out of memory
    if (curr == NULL || next == NULL)
    {
        fprintf (stderr, "Error: ran out of memory when trying to allocate %i sized cube\n", n);
        exit (EXIT_FAILURE);
    }

    pthread_t threadLS[threads];
    ThreadArgs argsLS[threads];


    for (int i = 0; i < threads; i++)
    {
        // Fill in the arguments to the workeroid* layer (void* pargs)
        argsLS[i].thread_id = i;
        argsLS[i].start = i+1;
        argsLS[i].n = n;
        argsLS[i].delta = delta;
        argsLS[i].next = next;
        argsLS[i].curr = curr;          
        argsLS[i].source = source;
    }
    
    // Solve Poisson's equation for the given inputs

    for (int iter = 0 ; iter < iterations; iter++) 
    {
        for (int i = 0 ; i < threads-1; i++)
        { 
            if (pthread_create (&threadLS[i], NULL, &MiddleLayer, &argsLS[i]) != 0)
            {
                fprintf (stderr, "Error creating worker thread!\n");
            }           
        }
        if (pthread_create (&threadLS[threads-1], NULL, &TopLayer, &argsLS[threads-1]) != 0)
        {
            fprintf (stderr, "Error creating worker thread!\n");
        }    

        for (int i = 0; i < threads; i++)
        {
            pthread_join (threadLS[i], NULL);
        }

        for (int i = 0 ; i < threads ; i++)
        {
            double* temp = argsLS[i].curr;
            argsLS[i].curr = argsLS[i].next;
            argsLS[i].next = temp; 
        }
    }

    // Free one of the buffers and return the correct answer in the other.
    // ThePROFILE caller is nmemcpy(curr, next, n3 * sizeof(double));ow responsible for free'ing the returned pointer.
    free (next);

    if (debug)
    {
        printf ("Finished solving.\n"); 
    }

    return curr;
}


int main (int argc, char **argv)
{
    // Default settings for solver
    int iterations = 100;
    int n = 901;
    int threads = n - 1;
    float delta = 1;

    // Parse the command line arguments
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp (argv[i], "-h") == 0 || strcmp (argv[i], "--help") == 0)
        {
            printf ("Usage: poisson [-n size] [-i iterations] [-t threads] [--debug]\n");
            return EXIT_SUCCESS;
        }

        if (strcmp (argv[i], "-n") == 0)
        {
            if (i == argc - 1)
            {
                fprintf (stderr, "Error: expected size after -n!\n");
                return EXIT_FAILURE;
            }

            n = atoi (argv[++i]);
        }

        if (strcmp (argv[i], "-i") == 0)
        {
            if (i == argc - 1)
            {
                fprintf (stderr, "Error: expected iterations after -i!\n");
                return EXIT_FAILURE;
            }

            iterations = atoi (argv[++i]);
        }

        if (strcmp (argv[i], "-t") == 0)
        {
            if (i == argc - 1)
            {
                fprintf (stderr, "Error: expected threads after -t!\n");
                return EXIT_FAILURE;
            }

            threads = atoi (argv[++i]);
        }

        if (strcmp (argv[i], "--debug") == 0)
        {
            debug = true;
        }
    }

    // Ensure we have an odd sized cube
    if (n % 2 == 0)
    {
        fprintf (stderr, "Error: n should be an odd number!\n");
        return EXIT_FAILURE;
    }

    // Create a source term with a single point in the centre
    double *source = (double*)calloc (n * n * n, sizeof (double));
    if (source == NULL)
    {
        fprintf (stderr, "Error: failed to allocated source term (n=%i)\n", n);
        return EXIT_FAILURE;
    }

    source[(n * n * n) / 2] = 1;

    #ifdef PROFILE
    clock_t time_start = clock();
    #endif

    // Calculate the resulting field with Dirichlet conditions
    double *result = poisson_dirichlet (n, source, iterations, threads, delta);

    #ifdef PROFILE
    clock_t time_finish = clock();
    double time_spent = (double)(time_finish - time_start) / CLOCKS_PER_SEC;
    printf("Time taken: %f\n",time_spent);
    #endif

    // Print out the middle slice of the cube for validation
    for (int x = 0; x < n; ++x)
    {
        for (int y = 0; y < n; ++y)
        {
            printf ("%0.5f ", result[((n / 2) * n + y) * n + x]);
        }
        printf ("\n");
    }


    free (source);
    free (result);

    return EXIT_SUCCESS;
}
