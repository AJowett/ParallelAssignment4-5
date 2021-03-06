/***************************************************************************/
/* Template for Asssignment 4/5 ********************************************/
/* Andrew Jowett            **(*****************************************/
/* Min Park                 **(*****************************************/
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

//Dimensions of the array
#define dim 1024
/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
unsigned long long g_start_cycles=0;
unsigned long long g_end_cycles=0;

//Our rank's ghost row and number of alive cells
int ghost_row_above[dim] = {1};
int ghost_row_below[dim] = {1};
unsigned long long num_alive_rank = 0;

//Our rank's chunk of the array
int** chunk;

pthread_barrier_t recv_barrier1;
pthread_barrier_t recv_barrier2;
pthread_barrier_t update_barrier;
pthread_barrier_t add_barrier;

pthread_mutex_t rank_alive_lock;

//Arugment struct for our threads
typedef struct arg_t {
    int start_row;
    int end_row;
    int num_ticks;
    int thread_num;
    int num_threads;
    int mpi_myrank;
    int mpi_commsize;
    double threshold;
    int** chunk;
    unsigned long long* rank_alive_cells;
} arg_t;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

//Helper function to count how many of a cell's neighbors are alive
int count_neighbs(int row, int col){
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

//Function invoked by our threads, resonsible for processing 
//the rows assigned to a given thread
void * process_rows(void * arg){
    //Unpacking our thread_arg
    arg_t thread_arg = *(arg_t *) arg;
    int start_row = thread_arg.start_row;
    int end_row = thread_arg.end_row;
    int num_ticks = thread_arg.num_ticks;
    int thread_num = thread_arg.thread_num;
    int num_threads = thread_arg.num_threads;
    int mpi_myrank = thread_arg.mpi_myrank;
    int mpi_commsize = thread_arg.mpi_commsize;
    double threshold = thread_arg.threshold;
    unsigned long long* rank_alive_cells = thread_arg.rank_alive_cells;

    //Inter rank communication variables
    int* send_buff = (int *) calloc(dim, sizeof(int));
    int* recv_buff = (int *) calloc(dim, sizeof(int));
    MPI_Request send_request1, send_request2;
    MPI_Request recv_request1, recv_request2;
    MPI_Status status;

    int num_alive_thread = (end_row - start_row) * dim;

    //Setup thread's RNG stream
    InitDefault();

    for(unsigned int tick = 0; tick < num_ticks; tick++){
        if(thread_num == 0){
            //Send our rank's bottom row to the rank beneath us
            for(unsigned int i = 0; i < dim; i++){
                send_buff[i] = chunk[dim - 1][i];
            }
        
            if(mpi_myrank == mpi_commsize - 1)
                MPI_Isend(send_buff, dim, MPI_INT, 0, 0, MPI_COMM_WORLD, &send_request1);
            else
                MPI_Isend(send_buff, dim, MPI_INT, mpi_myrank + 1, 0, MPI_COMM_WORLD, &send_request1);
            
            //Recieve the ghost row that goes above us
            if(mpi_myrank == 0)
                MPI_Irecv(recv_buff, dim, MPI_INT, mpi_commsize - 1, 0, MPI_COMM_WORLD, &recv_request1);
            else
                MPI_Irecv(recv_buff, dim, MPI_INT, mpi_myrank - 1, 0, MPI_COMM_WORLD, &recv_request1);
        
            MPI_Wait(&recv_request1, &status);
            
                    
            for(unsigned int i = 0; i < dim; i++){
                ghost_row_above[i] = recv_buff[i];
            }            
        }
        
        pthread_barrier_wait(&recv_barrier1);
        
        //Send our thread's top row
        if(thread_num == 0){
            for(unsigned int i = 0; i < dim; i++){
                send_buff[i] = chunk[0][i];
            }
            if(mpi_myrank == mpi_commsize - 1)
                MPI_Isend(send_buff, dim, MPI_INT, 0, 0, MPI_COMM_WORLD, &send_request2);
            else
                MPI_Isend(send_buff, dim, MPI_INT, mpi_myrank + 1, 0, MPI_COMM_WORLD, &send_request2);
            
            //Recieve the ghost row that goes below us
            if(mpi_myrank == 0)
                MPI_Irecv(recv_buff, dim, MPI_INT, mpi_commsize - 1, 0, MPI_COMM_WORLD, &recv_request2);
            else
                MPI_Irecv(recv_buff, dim, MPI_INT, mpi_myrank - 1, 0, MPI_COMM_WORLD, &recv_request2);
        
            MPI_Wait(&recv_request2, &status);
                        
            for(unsigned int i = 0; i < dim; i++){
                ghost_row_below[i] = recv_buff[i];
            }            
        }
        
        pthread_barrier_wait(&recv_barrier2);

        //Each thread goes through all the rows its assigned and updates each entry
        for(unsigned int i = start_row; i < end_row; i++){
            for(unsigned int j = 0; j < dim; j++){

                //Calculate if a cell is alive or dead
                double val = GenVal(i);
                if(val > threshold){
                    int num_neighbs = count_neighbs(i, j);
                    if(chunk[i][j] == ALIVE){
                        if(num_neighbs < 2 || num_neighbs > 3){
                            chunk[i][j] = DEAD;
                            num_alive_thread--;
                        }
                    }

                    else{
                        if(num_neighbs == 3){
                            chunk[i][j] = ALIVE;
                            num_alive_thread++;
                        }
                    }
                }

                else{
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

        //Update how many cells are alive in the rank based on our thread
        if(thread_num == 0){
            pthread_mutex_lock(&rank_alive_lock);
            num_alive_rank += num_alive_thread;
            pthread_mutex_unlock(&rank_alive_lock);
        }

        pthread_barrier_wait(&add_barrier);

        //Update the alive array fat the current tick
        if(thread_num == 0){
            pthread_mutex_lock(&rank_alive_lock);
            rank_alive_cells[tick] = num_alive_rank;
            num_alive_rank = 0;
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
    int mpi_myrank;
    int mpi_commsize;

    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    
    if(argc < 5){
        fprintf(stderr, "Usage: ./a.out <num_threads> <threshold> <exper_type> <num_ticks>\n");
        return EXIT_FAILURE;
    }

    int num_threads= atoi(argv[1]);
    //Current rank acts as a thread
    if(num_threads == 0)
        num_threads = 1;
    double threshold = atof(argv[2]);
    int exper_type = atoi(argv[3]);
    int num_ticks = atoi(argv[4]);

    //File to which we print the end state of the universe
    char* universeFilename;
    if(exper_type == 1 && argc != 6){
        fprintf(stderr, "Usage: ./a.out <num_threads> <threshold> <exper_type> <num_ticks> <universe file>\n");
        return EXIT_FAILURE;
    }
    else{
        universeFilename = argv[5];
    }

    //File to which we print the heat map of the end state of the universe
    char* heatmapFilename;
    if(exper_type == 2 && argc != 7){
        fprintf(stderr, "Usage: ./a.out <num_threads> <threshold> <exper_type> <num_ticks> <universe file> <heatmap file>\n");
        return EXIT_FAILURE;
    }
    else{
        heatmapFilename = argv[6];
    }
// Init 32,768 RNG streams - each rank has an independent stream
    InitDefault();
    
// Note, used the mpi_myrank to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.    
    printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
       mpi_myrank, mpi_commsize, GenVal(mpi_myrank));
    
    MPI_Barrier( MPI_COMM_WORLD );
    
    //Initialize our barriers and other variables
    pthread_barrier_init(&recv_barrier1, NULL, num_threads);
    pthread_barrier_init(&recv_barrier2, NULL, num_threads);
    pthread_barrier_init(&update_barrier, NULL, num_threads);
    pthread_barrier_init(&add_barrier, NULL, num_threads);

    pthread_mutex_init(&rank_alive_lock, NULL);

    chunk = (int **) calloc(dim, sizeof(int *));
    for(unsigned int i = 0; i < dim; i++){
        chunk[i] = (int *) calloc(dim, sizeof(int));
        for(unsigned int j = 0; j < dim; j++){
            chunk[i][j] = ALIVE;
        }
    }

    int rows_per_thread = dim/mpi_commsize/num_threads;
    unsigned int rows_per_rank = dim/mpi_commsize;

    unsigned long long* rank_alive_cells = (unsigned long long*) calloc(num_ticks, sizeof(unsigned long long *));
    //rank_alive_cells[0] = rows_per_rank * dim;
    unsigned long long* global_alive_cells = (unsigned long long*) calloc(num_ticks, sizeof(unsigned long long *));
    //global_alive_cells[0] = mpi_commsize * rank_alive_cells[0];

    if(mpi_myrank == 0){
        g_start_cycles = GetTimeBase();
    }

    //Start all our threads
    pthread_t tid[num_threads];
    for(unsigned int i = 0; i < num_threads; i++){

        //Create a new thread arg and pass it into the new thread
        arg_t* arg = (arg_t *) malloc(sizeof(arg_t));
        
        arg->start_row = i * rows_per_thread;
        arg->end_row = (i + 1) * rows_per_thread;
        arg->num_ticks = num_ticks;
        arg->thread_num = i;
        arg->num_threads = num_threads;
        arg->mpi_myrank = mpi_myrank;
        arg->mpi_commsize = mpi_commsize;
        arg->threshold = threshold;
        arg->rank_alive_cells = rank_alive_cells;

        pthread_create(&tid[i], NULL, process_rows, arg);
    }

    for(unsigned int i = 0; i < num_threads; i++){
        pthread_join(tid[i], NULL);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(mpi_myrank == 0){
        g_end_cycles = GetTimeBase();
    }

    //Fill in the global alive array
    MPI_Reduce(rank_alive_cells, global_alive_cells, num_ticks, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if(mpi_myrank == 0){
        for(unsigned int tick = 0; tick < num_ticks; tick++){
            printf("Num alive at tick %d: %llu \n", tick, global_alive_cells[tick]);
        }
        double time_in_secs = ((double) (g_end_cycles - g_start_cycles)) / g_processor_frequency;
        printf("Simulation complete in %e seconds\n", time_in_secs);
    }

    //We want to print out the final game universe
    if(exper_type == 1){
        //Reset our start cycles
        g_start_cycles = 0;
        g_end_cycles = 0;
        if(mpi_myrank == 0){
            g_start_cycles = GetTimeBase();
        }

        MPI_Status status;
        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, universeFilename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        
        //Set our view so that we are writing ints to the file
        MPI_File_set_view(fh, 0, MPI_INT, MPI_INT, (char *) NULL, MPI_INFO_NULL);


        //Go through the array and calculate the offset for each cell
        unsigned long long rank_offset;
        unsigned long long row_offset;
        for(unsigned int i = 0; i < rows_per_rank; i++){
            rank_offset = mpi_myrank * rows_per_rank * dim;
            row_offset = i * dim;
            MPI_File_write_at(fh, rank_offset + row_offset, chunk[i], dim, MPI_INT, &status);
        }

        MPI_File_close(&fh);

        MPI_Barrier(MPI_COMM_WORLD);

        if(mpi_myrank == 0){
            g_end_cycles = GetTimeBase();

            double time_in_secs = ((double) (g_end_cycles - g_start_cycles)) / g_processor_frequency;
            printf("Parallel io complete in %e seconds\n", time_in_secs);

        }
    }

    //We want to write the heatmap to a file
    if(exper_type == 2){
        if(mpi_myrank == 0){
            g_start_cycles = GetTimeBase();
        }

        int** heatmap;

        //the number of wide blocks that can be made out of the rank rows
        unsigned int heatmap_rows = rows_per_rank / 32;
        if(mpi_myrank != 0){
            heatmap = (int **) calloc(heatmap_rows, sizeof(int *));
            for(unsigned int i = 0; i < heatmap_rows; i++){
                heatmap[i] = (int *) calloc(1024, sizeof(int));
            }

            //Add the chunk items to the 32 wide block they belong to
            for(unsigned int i = 0; i < rows_per_rank; i++){
                for(unsigned int j = 0; j < dim; j++){
                    heatmap[i/32][j/32] += chunk[i][j];
                }
            }

            //Send each row of the heatmap to the 0th MPI rank
            MPI_Request send_request;
            for(unsigned int i = 0; i < rows_per_rank/32; i++){
                MPI_Isend(heatmap[i], 1024, MPI_INT, 0, 0, MPI_COMM_WORLD, &send_request);
            }
        }

        else{
            //Fill in the rows that we're responsible for in the heatmap
            heatmap = (int **) calloc(1024, sizeof(int *));
            for(unsigned int i = 0; i < 1024; i++){
                heatmap[i] = (int *) calloc(1024, sizeof(int));
            }

            for(unsigned int i = 0; i < rows_per_rank; i++){
                for(unsigned int j = 0; j < dim; j++){
                    heatmap[i/32][j/32] += chunk[i][j];
                }
            }

            //Receive the rows from all other threads
            MPI_Request recv_request;
            MPI_Status status;
            int* recv_buff = (int *) calloc(1024, sizeof(int));
            for(unsigned int i = 1; i < mpi_commsize; i++){
                for(unsigned int j = 1; j < rows_per_rank/32; j++){
                    MPI_Irecv(recv_buff, 1024, MPI_INT, i, 0, MPI_COMM_WORLD, &recv_request);
                    MPI_Wait(&recv_request, &status);
                    for(unsigned int k = 0; k < 1024; k++){
                        heatmap[i][k] = recv_buff[k];
                    }
                }
            }

            //Print each row to the heatmap file
            FILE* heatmapFile = fopen(heatmapFilename, "w");
            for(unsigned int  i = 0; i < 1024; i++){
                for(unsigned int j = 0; j < 1024; j++){
                    fprintf(heatmapFile, "%d ", heatmap[i][j]);
                }

                fprintf(heatmapFile, "\n");

            }

            fclose(heatmapFile);
        }


        MPI_Barrier(MPI_COMM_WORLD);

        if(mpi_myrank == 0){
            g_end_cycles = GetTimeBase();
            printf("heatmap complete in %e seconds\n", ((double) (g_end_cycles - g_start_cycles)) / g_processor_frequency);
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