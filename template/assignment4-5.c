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

#include"clcg4.h"

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
double processor_frequency = 1600000000.0; // processing speed for BG/Q

pthread_barrier_t send_barrier;
pthread_barrier_t recv_barrier;
pthread_barrier_t mpi_thread_barrier;
pthread_barrier_t rank_sum_barrier;

pthread_mutex_t rank_alive_lock;

//int total_ticks = 256;

int dim = 1024;

// You define these
typedef struct arg_t{
    int** chunk;
    int start_row;
    int end_row;

    int mpi_myrank;
    int mpi_commsize;
    int rows_per_thread;
    int rows_per_rank;
    int thread_num;
    int num_ticks;
    unsigned int* rank_alive_cells;

    double threshold;
    int* recv_row_start;
    int* recv_row_end;
} arg_t;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these
int check_neighbors(int row, int col, int** chunk, int rows_per_rank, int* recv_row_start, int* recv_row_end, int mpi_myrank, int thread_num){
    int num_neighbs = 0;
    
    if(row == 0){
        if(col == 0){
            num_neighbs += chunk[row][dim - 1]; //Neighbor to the left
            num_neighbs += recv_row_start[dim - 1]; //Diagonally upper left neighbor
            num_neighbs += chunk[row + 1][dim - 1]; //Diagonally lower left neighbor

            num_neighbs += chunk[row][col + 1]; //Neighbor to the right
            num_neighbs += recv_row_start[col + 1];   //Diagonally upper right neighbor
            num_neighbs += chunk[row + 1][col + 1]; //Diagonally lower right neighbor
        }

        else if(col == dim - 1){
            num_neighbs += chunk[row][0]; //Neighbor to the right
            num_neighbs += recv_row_start[0];   //Diagonally upper right neighbor
            num_neighbs += chunk[row + 1][0]; //Diagonally lower right neighbor

            num_neighbs += chunk[row][col - 1]; //Neighbor to the left
            num_neighbs += recv_row_start[col - 1]; //Diagonally upper left neighbor
            num_neighbs += chunk[row + 1][col- 1]; //Diagonally lower left neighbor
        }

        else{
            num_neighbs += chunk[row][col - 1]; //Neighbor to the left
            num_neighbs += recv_row_start[col - 1]; //Diagonally upper left neighbor
            num_neighbs += chunk[row + 1][col- 1]; //Diagonally lower left neighbor

            num_neighbs += chunk[row][col + 1]; //Neighbor to the right
            num_neighbs += recv_row_start[col + 1];   //Diagonally upper right neighbor
            num_neighbs += chunk[row + 1][col + 1] ;//Diagonally lower right neighbor
        }

        num_neighbs += recv_row_start[col]; //Neighbor above us
        num_neighbs += chunk[row + 1][col]; //Neighbor below us
    }

    else if(row == dim - 1){
        if(col == 0){
            num_neighbs += chunk[row][dim - 1]; //left neighbor
            num_neighbs += recv_row_end[dim - 1]; //diagonally lower left neighbor
            num_neighbs += chunk[row - 1][dim - 1]; //diagonally upper left neighbor

            num_neighbs += chunk[row][col + 1]; //right neighbor
            num_neighbs += recv_row_end[col + 1]; //diagonally lower right neighbor
            num_neighbs += chunk[row - 1][col + 1]; //diagonally upper right neighbor
         }
        else if(col == dim - 1){
            num_neighbs += chunk[row][0]; //right neighbor
            num_neighbs += recv_row_end[0]; //diagonally lower right neighbor
            num_neighbs += chunk[row - 1][0]; //diagonally upper right neighbor

            num_neighbs += chunk[row][col - 1]; //left neighbor
            num_neighbs += recv_row_end[col - 1]; //diagonally Lower left neighbor
            num_neighbs += chunk[row - 1][col - 1]; //diagonally upper left neighbor
        }
        else{
            num_neighbs += chunk[row][col - 1]; //left neighbor
            num_neighbs += recv_row_end[col - 1]; //diagonally Lower left neighbor
            num_neighbs += chunk[row - 1][col - 1]; //diagonally upper left neighbor

            num_neighbs += chunk[row][col + 1]; //right neighbor
            num_neighbs += recv_row_end[col + 1]; //diagonally lower right neighbor
            num_neighbs += chunk[row - 1][col + 1]; //diagonally upper right neighbor
        }

        num_neighbs += chunk[row-1][col]; //neighbor above us
        num_neighbs += recv_row_end[col]; //neighbor below us
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

    printf("rank: %d, thread: %d, row: %d, col: %d num_neighbs: %d\n", mpi_myrank, thread_num, row, col, num_neighbs);
    return num_neighbs;
}

void * process_rows(void * arg){
    arg_t thread_arg = *(arg_t *) arg;
    int mpi_myrank = thread_arg.mpi_myrank;
    int mpi_commsize = thread_arg.mpi_commsize;
    unsigned int start_row = thread_arg.start_row;
    unsigned int end_row  = thread_arg.end_row;
    int thread_num = thread_arg.thread_num;
    int rows_per_thread = thread_arg.rows_per_thread;
    int rows_per_rank = thread_arg.rows_per_rank;
    int num_ticks = thread_arg.num_ticks;
    int** chunk = thread_arg.chunk;
    int* recv_row_start = thread_arg.recv_row_start;
    int* recv_row_end = thread_arg.recv_row_end;
    unsigned int* rank_alive_cells = thread_arg.rank_alive_cells;
    
    printf("rank: %d, thread number %d with start_row: %d started, rows_per_thread: %d num_ticks: %d\n", mpi_myrank, thread_num, start_row, rows_per_thread, num_ticks);

    MPI_Request send_request, recv_request;
    MPI_Status status;
    //Initialize our RNG stream
    InitDefault();

    //All cells are alive initially
    int thread_alive_cells = dim * (thread_arg.end_row - thread_arg.start_row);
    for(int i = 0; i < num_ticks; i++){
        //Communicate ghost ranks
        printf("rank %d, thread %d, about to do communication\n", mpi_myrank, thread_num);
        if(thread_num == 0 && mpi_commsize > 1){
            printf("rank: %d, mpi_commsize: %d\n", mpi_myrank, mpi_commsize);
            if(mpi_myrank == mpi_commsize - 1){
                MPI_Isend(chunk[dim - 1], dim, MPI_INT, 0, 0, MPI_COMM_WORLD, &send_request);    
            }

            else{
                MPI_Isend(chunk[dim - 1], dim, MPI_INT, mpi_myrank + 1, 0, MPI_COMM_WORLD, &send_request);
            }

            if(mpi_myrank == 0){
                MPI_Irecv(recv_row_start, dim, MPI_INT, mpi_commsize - 1, 0, MPI_COMM_WORLD, &recv_request);
            }

            else{
                MPI_Irecv(recv_row_start, dim, MPI_INT, mpi_myrank - 1, 0, MPI_COMM_WORLD, &recv_request);
            }
        }

        //Wait for all recieves to go through
        if(thread_num == 0 && mpi_commsize > 1){
            //printf("About to wait\n");
            MPI_Wait(&recv_request, &status);
            printf("rank: %d Done with waiting\n", mpi_myrank);
        }

        if(mpi_myrank == 0 && thread_num == 0){
            //printf("Tick: %d\n", i);
        }
        pthread_barrier_wait(&recv_barrier);
        
        if(thread_num == 0 && mpi_commsize > 1){
            printf("part 2 rank: %d, mpi_commsize: %d\n", mpi_myrank, mpi_commsize);
            if(mpi_myrank == mpi_commsize - 1){
                MPI_Isend(chunk[0], dim, MPI_INT, 0, 0, MPI_COMM_WORLD, &send_request);    
            }

            else{
                MPI_Isend(chunk[0], dim, MPI_INT, mpi_myrank + 1, 0, MPI_COMM_WORLD, &send_request);
            }

            if(mpi_myrank == 0){
                MPI_Irecv(recv_row_end, dim, MPI_INT, mpi_commsize - 1, 0, MPI_COMM_WORLD, &recv_request);
            }

            else{
                MPI_Irecv(recv_row_end, dim, MPI_INT, mpi_myrank - 1, 0, MPI_COMM_WORLD, &recv_request);
            }
        }

        //Wait for all recieves to go through
        if(thread_num == 0 && mpi_commsize > 1){
            //printf("About to wait\n");
            MPI_Wait(&recv_request, &status);
            printf("rank: %d Done with waiting\n", mpi_myrank);
        }

        pthread_barrier_wait(&recv_barrier);

        printf("rank: %d thread: %d about to update rows\n", mpi_myrank, thread_num);
        //Go through all the assigned rows for this thread
        for(unsigned int j = start_row; j < end_row; j++){
            for(unsigned int k = 0; k < dim; k++){
                //calculate the number of alive neighbors 
                int num_neighbs = 0;
                printf("rank: %d, thread: %d, Checking neighbors for cell: (%d, %d)\n", mpi_myrank, thread_num, j, k);
                num_neighbs = check_neighbors(j, k, chunk, rows_per_rank, recv_row_start, recv_row_end, mpi_myrank, thread_num);
                //RNG + threshold
                double val = GenVal(j + (thread_arg.rows_per_rank * mpi_myrank));
                if(val > thread_arg.threshold){
                    if(chunk[j][k] == ALIVE){
                        if(num_neighbs < 2){
                            chunk[j][k] = DEAD;
                            thread_alive_cells--;
                        }
                        else if(num_neighbs > 3){
                            chunk[j][k] = DEAD;
                            thread_alive_cells--;
                        }
                    }

                    else{
                        if(num_neighbs == 3){
                            chunk[j][k] = ALIVE;
                            thread_alive_cells++;
                        }
                    }
                }
                
                else{
                    val = GenVal(j + (thread_arg.rows_per_rank * mpi_myrank));
                    //Cell becomes alive
                    if(val < 0.5){
                        if(chunk[j][k] == DEAD){
                            chunk[j][k] = ALIVE;
                            thread_alive_cells++;
                        }
                    }
                    //Cell becomes dead
                    else{
                        if(chunk[j][k] == ALIVE){
                            chunk[j][k] = DEAD;
                            thread_alive_cells--;
                        }
                    }
                }
            }
        }

        printf("rank %d, thread %d: done with row update for current tick\n", mpi_myrank, thread_num);

        //Have all threads add the number alive in their rows to the rank total
        pthread_mutex_lock(&rank_alive_lock);
        rank_alive_cells[i] += thread_alive_cells;
        pthread_mutex_unlock(&rank_alive_lock);

        printf("rank %d, thread %d, added alive to rank_alive\n", mpi_myrank, thread_num);

        //Wait for everyone to get here
        MPI_Barrier(MPI_COMM_WORLD);
        pthread_barrier_wait(&mpi_thread_barrier);

        printf("rank %d, thread %d, Done waiting\n", mpi_myrank, thread_num);
    }  

    printf("Done with main loop\n");
    //Make sure all threads have finished updating rows
    pthread_barrier_wait(&mpi_thread_barrier);


    //Return the array with alive cells per tick
    if(thread_num == 0){
        pthread_exit(&thread_alive_cells);
    }
    else{
        return NULL;
    }
}

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

    if(argc != 5){
        fprintf(stderr, "Usage: ./a.out <num_threads> <threshold> <exper_type> <num_ticks>\n");
        return EXIT_FAILURE;
    }

    int num_threads= atoi(argv[1]);
    double threshold = atof(argv[2]);
    int exper_type = atoi(argv[3]);
    int num_ticks = atoi(argv[4]);

// Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    
// Init 32,768 RNG streams - each rank has an independent stream
    InitDefault();
    
// Note, used the mpi_myrank to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.    
    printf("Rank %d of %d has been started and a first Random Value of %lf\n", mpi_myrank, mpi_commsize, GenVal(mpi_myrank));

    unsigned int* rank_alive_cells = (unsigned int *) calloc(num_ticks, sizeof(unsigned int));
    unsigned int* global_alive_cells = NULL;

    if(mpi_myrank == 0){
        start_cycle =  GetTimeBase();
        global_alive_cells = (unsigned int*) calloc(num_ticks, sizeof(unsigned int));
    }

    int** my_rank_chunk = (int **) calloc(dim/mpi_commsize, sizeof(int *));
    for(unsigned int i = 0; i < dim/mpi_commsize; i++){
        my_rank_chunk[i] = (int *) calloc(dim, sizeof(int));
        for(unsigned int j = 0; j < dim; j++){
            my_rank_chunk[i][j] = ALIVE;
        }
    }

    int* recv_row_start = (int *) calloc(dim, sizeof(int));
    int* recv_row_end = (int *) calloc(dim, sizeof(int));
    for(unsigned int i = 0; i < dim; i++){
        recv_row_start[i] = ALIVE;
        recv_row_end[i] = ALIVE;
    }

    printf("Hi\n");
    pthread_t tid[num_threads];
    
    pthread_barrier_init(&send_barrier, NULL, num_threads);
    pthread_barrier_init(&recv_barrier, NULL, num_threads);
    pthread_barrier_init(&mpi_thread_barrier, NULL, num_threads);
    pthread_barrier_init(&rank_sum_barrier, NULL, num_threads);

    pthread_mutex_init(&rank_alive_lock, NULL);

    MPI_Status status;

    int rows_per_thread = dim/mpi_commsize/num_threads;
    unsigned int rows_per_rank = dim/mpi_commsize;

    for(int j = 0; j < num_threads; j++){

        arg_t* thread_arg = (arg_t *) malloc(sizeof(arg_t));

        thread_arg->chunk = my_rank_chunk;
        thread_arg->start_row = j * rows_per_thread;
        thread_arg->end_row = (j + 1) * rows_per_thread - 1;

        thread_arg->mpi_myrank = mpi_myrank;
        thread_arg->mpi_commsize = mpi_commsize;
        thread_arg->rows_per_thread = rows_per_thread;
        thread_arg->rows_per_rank = rows_per_rank;
        thread_arg->thread_num = j;
        thread_arg->num_ticks = num_ticks;
        thread_arg->rank_alive_cells = rank_alive_cells;

        thread_arg->threshold = threshold;
        thread_arg->recv_row_start = recv_row_start;
        thread_arg->recv_row_end = recv_row_end;

        printf("THread number %d with start row %d about to start\n", thread_arg->thread_num, j * rows_per_thread);
        pthread_create(&tid[j], NULL, process_rows, thread_arg);
    }
    printf("Done creating threads\n");

    for(int j = 0; j < num_threads; j++){
        pthread_join(tid[j], NULL);
    } 
    printf("Rank: %d finished simulation ticks\n", mpi_myrank);

    //Sum the arrays for every tick
    MPI_Reduce(rank_alive_cells, global_alive_cells, num_ticks, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    //MPI_Barrier( MPI_COMM_WORLD );
    if(mpi_myrank == 0){
        end_cycle = GetTimeBase();
    }

    unsigned long long io_start_cycle = 0;
    unsigned long long io_end_cycle = 0;
    if(exper_type == 1){
        if(mpi_myrank == 0)
            io_start_cycle = GetTimeBase();

        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, "universe.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        
        for(int i = 0; i < dim/mpi_commsize; i++){
            MPI_File_write_at(fh, i + (rows_per_rank * mpi_myrank), my_rank_chunk[i], dim, MPI_INT, &status);
        }
        MPI_File_close(&fh);

        if(mpi_myrank == 0)
            io_end_cycle = GetTimeBase();

        MPI_Barrier(MPI_COMM_WORLD);
    }

    if(exper_type == 2){
        unsigned int heatmap_rows = rows_per_rank / 32;
        int** heatmap = (int **) calloc(rows_per_rank/32, sizeof(int *));
        for(unsigned int i = 0; i < heatmap_rows; i++){
            heatmap[i] = (int *) calloc(1024, sizeof(int));
        }

        for(unsigned int i = 0; i < rows_per_rank; i++){
            for(unsigned int j = 0; j < dim; j++){
                heatmap[i/32][j/32] += my_rank_chunk[i][j];
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, "heatmap.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

        for(unsigned int i = 0; i < heatmap_rows; i++){
            MPI_File_write_at(fh, i + (heatmap_rows * mpi_myrank), heatmap[i], 1024, MPI_INT, &status);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    if(mpi_myrank == 0){
        for(int i = 0; i < num_ticks; i++){
            printf("Tick %d, %du cells alive\n", i, global_alive_cells[i]);
        }

        printf("Simulation complete in %e seconds\n", ((double) (end_cycle - start_cycle)) / processor_frequency);

        if(exper_type == 1){
            printf("IO complete in %e seconds\n", ((double) (io_end_cycle - io_start_cycle)) / processor_frequency);
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
