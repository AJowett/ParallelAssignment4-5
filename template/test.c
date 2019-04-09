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

#define dim 1024
/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
unsigned long long g_start_cycles=0;
unsigned long long g_end_cycles=0;

// You define these
int ghost_row_above[dim] = {1};
int ghost_row_below[dim] = {1};
unsigned long long num_alive_rank;

pthread_barrier_t recv_barrier1;
pthread_barrier_t recv_barrier2;
pthread_barrier_t update_barrier;
pthread_barrier_t add_barrier;

pthread_mutex_t rank_alive_lock;

typedef struct arg_t {
    int start_row;
    int end_row;
    int num_ticks;
    int thread_num;
    int num_threads;
    int mpi_myrank;
    int mpi_commsize;
    int threshold;
    int** chunk;
    unsigned long long* rank_alive_cells;
} arg_t;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

int count_neighbs(int row, int col, int** chunk){
    int num_neighbs = 0;
    if(row == 0){
        if(col == 0){
            num_neighbs += chunk[row][dim - 1]; //Neighbor to the left
            num_neighbs += ghost_row_above[dim - 1]; //Diagonally upper left neighbor
            num_neighbs += chunk[row + 1][dim - 1]; //Diagonally lower left neighbor

            num_neighbs += chunk[row][col + 1]; //Neighbor to the right
            num_neighbs += ghost_row_above[col + 1];   //Diagonally upper right neighbor
            num_neighbs += chunk[row + 1][col + 1]; //Diagonally lower right neighbor
        }

        else if(col == dim - 1){
            num_neighbs += chunk[row][0]; //Neighbor to the right
            num_neighbs += ghost_row_above[0];   //Diagonally upper right neighbor
            num_neighbs += chunk[row + 1][0]; //Diagonally lower right neighbor

            num_neighbs += chunk[row][col - 1]; //Neighbor to the left
            num_neighbs += ghost_row_above[col - 1]; //Diagonally upper left neighbor
            num_neighbs += chunk[row + 1][col- 1]; //Diagonally lower left neighbor
        }

        else{
            num_neighbs += chunk[row][col - 1]; //Neighbor to the left
            num_neighbs += ghost_row_above[col - 1]; //Diagonally upper left neighbor
            num_neighbs += chunk[row + 1][col- 1]; //Diagonally lower left neighbor

            num_neighbs += chunk[row][col + 1]; //Neighbor to the right
            num_neighbs += ghost_row_above[col + 1];   //Diagonally upper right neighbor
            num_neighbs += chunk[row + 1][col + 1] ;//Diagonally lower right neighbor
        }

        num_neighbs += ghost_row_above[col]; //Neighbor above us
        num_neighbs += chunk[row + 1][col]; //Neighbor below us
    }

    else if(row == dim - 1){
        if(col == 0){
            num_neighbs += chunk[row][dim - 1]; //left neighbor
            num_neighbs += ghost_row_below[dim - 1]; //diagonally lower left neighbor
            num_neighbs += chunk[row - 1][dim - 1]; //diagonally upper left neighbor

            num_neighbs += chunk[row][col + 1]; //right neighbor
            num_neighbs += ghost_row_below[col + 1]; //diagonally lower right neighbor
            num_neighbs += chunk[row - 1][col + 1]; //diagonally upper right neighbor
         }
        else if(col == dim - 1){
            num_neighbs += chunk[row][0]; //right neighbor
            num_neighbs += ghost_row_below[0]; //diagonally lower right neighbor
            num_neighbs += chunk[row - 1][0]; //diagonally upper right neighbor

            num_neighbs += chunk[row][col - 1]; //left neighbor
            num_neighbs += ghost_row_below[col - 1]; //diagonally Lower left neighbor
            num_neighbs += chunk[row - 1][col - 1]; //diagonally upper left neighbor
        }
        else{
            num_neighbs += chunk[row][col - 1]; //left neighbor
            num_neighbs += ghost_row_below[col - 1]; //diagonally Lower left neighbor
            num_neighbs += chunk[row - 1][col - 1]; //diagonally upper left neighbor

            num_neighbs += chunk[row][col + 1]; //right neighbor
            num_neighbs += ghost_row_below[col + 1]; //diagonally lower right neighbor
            num_neighbs += chunk[row - 1][col + 1]; //diagonally upper right neighbor
        }

        num_neighbs += chunk[row-1][col]; //neighbor above us
        num_neighbs += ghost_row_below[col]; //neighbor below us
    }
    else{
        if(col == 0){
            num_neighbs += chunk[row][dim - 1]; //left neighbor
            num_neighbs += chunk[row + 1][dim -1]; //diagonally lower left
            num_neighbs += chunk[row - 1][dim -1]; //diagonally upper left
        
            num_neighbs += chunk[row][col + 1]; //right neighbor
            num_neighbs += chunk[row + 1][col + 1]; //diagonally lower right
            num_neighbs += chunk[row - 1][col + 1]; //diagonally upper right
        }

        else if(col == dim - 1){
            num_neighbs += chunk[row][0]; //right neighbor
            num_neighbs += chunk[row + 1][0]; //diagonally lower right
            num_neighbs += chunk[row - 1][0]; //diagonally upper right
        
            num_neighbs += chunk[row][col - 1]; //left neighbor
            num_neighbs += chunk[row + 1][col - 1]; //diagonally lower left
            num_neighbs += chunk[row - 1][col - 1]; //diagonally uppper left
        }

        else{
            num_neighbs += chunk[row][col - 1]; //left neighbor
            num_neighbs += chunk[row + 1][col - 1]; //diagonally lower left
            num_neighbs += chunk[row - 1][col - 1]; //diagonally uppper left

            num_neighbs += chunk[row][col + 1]; //right neighbor
            num_neighbs += chunk[row + 1][col + 1]; //diagonally lower right
            num_neighbs += chunk[row - 1][col + 1]; //diagonally upper right
        }

        num_neighbs += chunk[row - 1][col]; //above us
        num_neighbs += chunk[row + 1][col]; //below us
    }

    return num_neighbs;
}


// You define these
void * process_rows(void * arg){
    arg_t thread_arg = *(arg_t *) arg;
    int start_row = thread_arg.start_row;
    int end_row = thread_arg.end_row;
    int num_ticks = thread_arg.num_ticks;
    int thread_num = thread_arg.thread_num;
    int num_threads = thread_arg.num_threads;
    int mpi_myrank = thread_arg.mpi_myrank;
    int mpi_commsize = thread_arg.mpi_commsize;
    int threshold = thread_arg.threshold;
    int** chunk = thread_arg.chunk;
    unsigned long long* rank_alive_cells = thread_arg.rank_alive_cells;

    int* send_buff = (int *) calloc(dim, sizeof(int));

    MPI_Request send_request1, send_request2;
    MPI_Request recv_request1, recv_request2;
    MPI_Status status;

    int num_alive_thread = (end_row - start_row) * dim;

    printf("rank: %d, thread: %d, start_row: %d, end_row: %d\n", mpi_myrank, thread_num, start_row, end_row);

    InitDefault();

    for(unsigned int tick = 0; tick < num_ticks; tick++){
        if(thread_num == 0){
            for(unsigned int i = 0; i < dim; i++){
                send_buff[i] = chunk[dim - 1][i];
            }
            if(mpi_myrank == mpi_commsize - 1)
                MPI_Isend(&send_buff, dim, MPI_INT, 0, 0, MPI_COMM_WORLD, &send_request1);
            else
                MPI_Isend(&send_buff, dim, MPI_INT, mpi_myrank + 1, 0, MPI_COMM_WORLD, &send_request1);
        
            if(mpi_myrank == 0)
                MPI_Irecv(ghost_row_above, dim, MPI_INT, mpi_commsize - 1, 0, MPI_COMM_WORLD, &recv_request1);
            else
                MPI_Irecv(ghost_row_above, dim, MPI_INT, mpi_myrank - 1, 0, MPI_COMM_WORLD, &recv_request1);
        
            MPI_Wait(&recv_request1, &status);
        }
        
        pthread_barrier_wait(&recv_barrier1);
        
        if(thread_num == 0){
            for(unsigned int i = 0; i < dim; i++){
                send_buff[i] = chunk[0][i];
            }
            if(mpi_myrank == mpi_commsize - 1)
                MPI_Isend(&send_buff, dim, MPI_INT, 0, 0, MPI_COMM_WORLD, &send_request2);
            else
                MPI_Isend(&send_buff, dim, MPI_INT, mpi_myrank + 1, 0, MPI_COMM_WORLD, &send_request2);
        
            if(mpi_myrank == 0)
                MPI_Irecv(ghost_row_below, dim, MPI_INT, mpi_commsize - 1, 0, MPI_COMM_WORLD, &recv_request2);
            else
                MPI_Irecv(ghost_row_below, dim, MPI_INT, mpi_myrank - 1, 0, MPI_COMM_WORLD, &recv_request2);
        
            MPI_Wait(&recv_request2, &status);
        }
        
        pthread_barrier_wait(&recv_barrier2);

        for(unsigned int i = start_row; i < end_row; i++){
            for(unsigned int j = 0; j < dim; j++){

                double val = GenVal(i);
                if(val > threshold){
                    int num_neighbs = count_neighbs(i, j, chunk);
                    if(chunk[i][j] == ALIVE){
                        if(num_neighbs < 2 || num_neighbs > 3){
                            chunk[i][j] = DEAD;
                            num_alive_thread--;
                        }
                    }

                    else{
                        if(num_neighbs > 3){
                            chunk[i][j] = ALIVE;
                            num_alive_thread++;
                        }
                    }
                }

                else{
                    val = GenVal(i + (dim/mpi_commsize * num_threads));
                    if(val < 0.5 && chunk[i][j] == ALIVE){
                        chunk[i][j] = DEAD;
                        num_alive_thread--;
                    }

                    else if(val < 0.5 && chunk[i][j] == DEAD){
                        chunk[i][j] = ALIVE;
                        num_alive_thread++;
                    }
                }
            }
        }

        pthread_barrier_wait(&update_barrier);

        if(thread_num == 0){
            pthread_mutex_lock(&rank_alive_lock);
            num_alive_rank += num_alive_thread;
            pthread_mutex_unlock(&rank_alive_lock);
        }

        pthread_barrier_wait(&add_barrier);

        if(thread_num == 0){
            printf("rank %d: %llu alive at tick %d\n", mpi_myrank, );
            pthread_mutex_lock(&rank_alive_lock);
            rank_alive_cells[tick] = num_alive_rank;
            pthread_mutex_unlock(&rank_alive_lock);
        }
    }

    return NULL;
}

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
//    int i = 0;
    int mpi_myrank;
    int mpi_commsize;
// Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    
    if(argc != 5){
        fprintf(stderr, "Usage: ./a.out <num_threads> <threshold> <exper_type> <num_ticks>\n");
        return EXIT_FAILURE;
    }

    int num_threads= atoi(argv[1]);
    double threshold = atof(argv[2]);
    int exper_type = atoi(argv[3]);
    int num_ticks = atoi(argv[4]);
// Init 32,768 RNG streams - each rank has an independent stream
    InitDefault();
    
// Note, used the mpi_myrank to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.    
    printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
	   mpi_myrank, mpi_commsize, GenVal(mpi_myrank));
    
    MPI_Barrier( MPI_COMM_WORLD );
    
// Insert your code
    pthread_barrier_init(&recv_barrier1, NULL, num_threads);
    pthread_barrier_init(&recv_barrier2, NULL, num_threads);
    pthread_barrier_init(&update_barrier, NULL, num_threads);
    pthread_barrier_init(&add_barrier, NULL, num_threads);

    pthread_mutex_init(&rank_alive_lock, NULL);

    int** chunk = (int **) calloc(dim, sizeof(int *));
    for(unsigned int i = 0; i < dim; i++){
        chunk[i] = (int *) calloc(dim, sizeof(int));
        for(unsigned int j = 0; j < dim; j++){
            chunk[i][j] = ALIVE;
        }
    }

    int rows_per_thread = dim/mpi_commsize/num_threads;
    unsigned int rows_per_rank = dim/mpi_commsize;

    unsigned long long* rank_alive_cells = (unsigned long long*) calloc(num_ticks, sizeof(unsigned long long *));
    rank_alive_cells[0] = rows_per_rank * dim;
    unsigned long long* global_alive_cells = (unsigned long long*) calloc(num_ticks, sizeof(unsigned long long *));
    global_alive_cells[0] = mpi_commsize * rank_alive_cells[0];

    if(mpi_myrank == 0){
        g_start_cycles = GetTimeBase();
    }
    pthread_t tid[num_threads];
    for(unsigned int i = 0; i < num_threads; i++){
        arg_t* arg = (arg_t *) malloc(sizeof(arg_t));
        
        arg->start_row = i * rows_per_thread;
        arg->end_row = (i + 1) * rows_per_thread;
        arg->num_ticks = num_ticks;
        arg->thread_num = i;
        arg->num_threads = num_threads;
        arg->mpi_myrank = mpi_myrank;
        arg->mpi_commsize = mpi_commsize;
        arg->threshold = threshold;
        arg->chunk = chunk;
        arg->rank_alive_cells = rank_alive_cells;

        pthread_create(&tid[i], NULL, process_rows, arg);
    }

    for(unsigned int i = 0; i < num_threads; i++){
        pthread_join(tid[i], NULL);
    }
    if(mpi_myrank == 0){
        g_end_cycles = GetTimeBase();
    }

    MPI_Reduce(rank_alive_cells, global_alive_cells, num_ticks, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(mpi_myrank == 0){
        for(unsigned int tick = 0; tick < num_ticks; tick++){
            printf("Num alive at tick %d: %llu \n", tick, global_alive_cells[tick]);
        }
    }

// END -Perform a barrier and then leave MPI
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/