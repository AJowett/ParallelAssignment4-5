/***************************************************************************/
/* Template for Asssignment 4/5 ********************************************/
/* Team Names Here              **(*****************************************/
/***************************************************************************/

/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>

#include<clcg4.h>

#include<mpi.h>
#include<pthread.h>

// #define BGQ 1 // when running BG/Q, comment out when testing on mastiff

#ifdef BGQ
#include<hwi/include/bqc/A2_inlines.h>
#else
#define GetTimeBase MPI_Wtime            
#endif

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define ALIVE 1
#define DEAD  0

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
unsigned long long g_start_cycles=0;
unsigned long long g_end_cycles=0;

int total_ticks = 256;

double dim = 32768;

// You define these
typedef struct arg_t{
    int** chunk;
    int start_row;
    int end_row;
    int mpi_myrank;
    int thread_num;
    double threshold;
    int* recv_row;
    int* send_row;
} arg_t;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these
void * process_row(void* arg);

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
//    int i = 0;
    int mpi_myrank;
    int mpi_commsize;

    unsigned long long start_cycle = 0;
    unsigned long long end_cycle = 0;

    int num_threads= atoi(argv[1]);
    double threshold = atof(argv[2]);
    int exper_type = atoi(argv[3]);

// Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    
// Init 32,768 RNG streams - each rank has an independent stream
    InitDefault();
    
// Note, used the mpi_myrank to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.    
    printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
       mpi_myrank, mpi_commsize, GenVal(mpi_myrank));

    if(mpi_myrank == 0){
        start_time =  GetTimeBase();
    }

    int** my_rank_chunk = (int **) calloc(dim/mpi_commsize, sizeof(int *);
    for(unsigned int i = 0; i < dim/mpi_commsize){
        my_rank_chunk[i] = (int *) calloc(dim, sizeof(int));
        for(unsigned int j = 0; j < dim; j++){
            my_rank_chunk[j] = ALIVE;
        }
    }

    int* recv_row = (int *) calloc(dim, sizeof(int));
    int* send_row = (int *) calloc(dim, sizeof(int));

    pthread_t tid[num_threads];
      
    MPI_Request send_request, recv_request;
    MPI_Status status;

    int* alive_cells
    if(mpi_myrank == 0){
        alive_cells = (int *) calloc(total_ticks, sizeof(int));
    }

    for(int i = 0; i < total_ticks; i++){
        //Communicate ghost ranks
        if(mpi_myrank != mpi_commsize - 1){
            MPI_Isend(recv_row, dim, MPI_INT, 0, 0 MPI_COMM_WORLD, &send_request);    
        }
        else{
            MPI_Isend(recv_row, dim, MPI_INT, mpi_myrank + 1, 0, MPI_COMM_WORLD, &send_request);
        }

        if(mpi_myrank == 0){
            MPI_Irecv(recv_row, dim, MPI_INT, mpi_commsize - 1, 0, MPI_COMM_WORLD, &recv_request);
        }
        else{
            MPI_Irecv(recv_row, dim, MPI_INT, mpi_myrank - 1, 0, MPI_COMM_WORLD, &recv_request);
        }
        
        MPI_Wait(&recv_request, &status);
    
        int rows_per_thread = dim/mpi_commsize/num_threads;

        unsigned long long rank_alive_cells = 0;
        unsigned long long global_alive_cells = 0;
        for(int j = 0; j < num_threads; j++){
            int* thread_alive_cells;

            arg_t thread_arg;
            thread_arg.start_row = j * rows_per_thread;
            thread_arg.end_row = (j + 1) * rows_per_thread - 1;
            thread_arg.chunk = my_rank_chunk;
            pthread_create(tid[j], process_row, (void *) thread_arg);
        }

        for(int j = 0; j < num_threads; j++){
            pthread_join(tid[j], (void **) &thread_alive_cells);
            rank_alive_cells += *thread_alive_cells;
            free(thread_alive_cells);
            thread_alive_cells = NULL;
        }

        MPI_Reduce(&rank_alive_cells, &global_alive_cells, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        if(mpi_myrank == 0){
            alive_cells[i] = global_alive_cells;
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }   

    MPI_Barrier( MPI_COMM_WORLD );
    if(mpi_myrank == 0){
        end_cycle = GetTimeBase();
    }

    if(exper_type){

    }

    if(exper_type){

    }
// Insert your code
    

// END -Perform a barrier and then leave MPI
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/
