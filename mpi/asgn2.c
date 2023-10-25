#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <memory.h>
#include <string.h>

#define DEBUG 0
#define MASTER 0
#define BASE 0
#define GRID 1
#define EXTRA 2
#define TABLE_SIZE 5
#define CHARGING_PORTS 5
#define LOW_PORTS 0
#define CYCLE_DURATION_SECONDS 5
#define USE_DURATION 2
#define TERMINATE 10

int base_station_io(MPI_Comm world_comm, MPI_Comm comm, int m, int n);
int charging_node_io(MPI_Comm world_comm, MPI_Comm comm, int m, int n);
int extra_node_io(MPI_Comm world_comm, MPI_Comm comm);

int main(int argc, char *argv[])
{
    int my_rank, size, provided;

    // MPI_Init(&argc, &argv);
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int m, n;

    MPI_Wtime();

    if (my_rank == MASTER)
    {
        if (DEBUG)
        {
            printf("DEBUG is on\n");
            fflush(stdout);
        }
        // receive simulation grid size input
        printf("Specify the simulation EV charging nodes grid size M * N (e.g. 4 3): ");
        fflush(stdout);
        scanf("%d%d", &m, &n);
    }

    // broadcast grid size to every MPI process
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // splitting base and grid communicator (as well as excess processes)
    int color;
    if (my_rank == MASTER)
    {
        color = BASE;
    }
    else if (my_rank <= m * n)
    {
        color = GRID;
    }
    else
    {
        color = EXTRA;
    }

    MPI_Comm new_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, my_rank, &new_comm);

    int my_new_comm_rank;
    MPI_Comm_rank(new_comm, &my_new_comm_rank);

    if (DEBUG)
    {
        printf("[MPI process %d] I am now MPI process %d in group %d.\n", my_rank, my_new_comm_rank, color);
        fflush(stdout);
    }

    // run specific routines for each type of processes
    if (my_rank == MASTER)
        base_station_io(MPI_COMM_WORLD, new_comm, m, n);
    else if (my_rank <= m * n)
        charging_node_io(MPI_COMM_WORLD, new_comm, m, n);
    else
    {
        extra_node_io(MPI_COMM_WORLD, new_comm);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}

int base_station_io(MPI_Comm world_comm, MPI_Comm comm, int m, int n)
{
    int iter = 0;
    int number_of_nodes = m * n;
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // getting neighbor and coords of each node
    int *number_of_neighbors = (int *)calloc(number_of_nodes, sizeof(int));
    int **neighbors = (int **)malloc(number_of_nodes * sizeof(int *));
    for (int i = 0; i < number_of_nodes; i++)
    {
        neighbors[i] = (int *)malloc(4 * sizeof(int));
        for (int j = 0; j < 4; j++)
        {
            neighbors[i][j] = -2;
        }
    }
    MPI_Status receive_neighbor_status[number_of_nodes];
#pragma omp parallel for num_threads(number_of_nodes)
    for (int i = 0; i < number_of_nodes; i++)
    {
        int temp[4];
        MPI_Recv(temp, 4, MPI_INT, i + 1, 4, world_comm, &receive_neighbor_status[i]);
        if (DEBUG)
        {
            printf("[BASE] received neighbors from grid node %d.\n", i);
            fflush(stdout);
        }
        for (int j = 0; j < 4; j++)
        {
            neighbors[i][j] = temp[j];
            if (temp[j] >= 0)
            {
                number_of_neighbors[i] += 1;
            }
        }
    }
    int **coords = (int **)malloc(number_of_nodes * sizeof(int *));
    for (int i = 0; i < number_of_nodes; i++)
    {
        coords[i] = (int *)malloc(2 * sizeof(int));
    }
    MPI_Status receive_coord_status[number_of_nodes];
#pragma omp parallel for num_threads(number_of_nodes)
    for (int i = 0; i < number_of_nodes; i++)
    {
        int temp[2];
        MPI_Recv(temp, 2, MPI_INT, i + 1, 5, world_comm, &receive_coord_status[i]);
        if (DEBUG)
        {
            printf("[BASE] received coords from grid node %d.\n", i);
            fflush(stdout);
        }
        for (int j = 0; j < 2; j++)
        {
            coords[i][j] = temp[j];
        }
    }

    // compute neighbors of neighbors
    int **neighbor_of_neighbors = (int **)malloc(number_of_nodes * sizeof(int *));
    for (int i = 0; i < number_of_nodes; i++)
    {
        neighbor_of_neighbors[i] = (int *)malloc(8 * sizeof(int));
        for (int j = 0; j < 8; j++)
        {
            neighbor_of_neighbors[i][j] = -2;
        }
        int counter = 0;
        for (int j = 0; j < 4; j++)
        {
            int n = neighbors[i][j];
            if (n >= 0)
            {
                for (int k = 0; k < 4; k++)
                {
                    int non = neighbors[n][k];
                    if (non >= 0 && non != i)
                    {
                        neighbor_of_neighbors[i][counter] = non;
                        counter += 1;
                    }
                }
            }
        }
    }

    int start_time;
    int current_time = (int)time(NULL);

    // creating log file
    FILE *file;
    file = fopen("base_log.txt", "w");

    if (file == NULL)
    {
        perror("Error opening the file");
        return 0;
    }

    int term = 0;
    while (!term)
    {
        start_time = current_time;
        printf("\n\nITERATION: %d \n\n", iter + 1);
        fflush(stdout);

        MPI_Barrier(world_comm);

        MPI_Request receive_request[number_of_nodes];
        MPI_Status receive_status[number_of_nodes];
        int *alert_time = (int *)calloc(number_of_nodes, sizeof(int));
        double *report_time = (double *)calloc(number_of_nodes, sizeof(double));
        double *elapsed_time = (double *)calloc(number_of_nodes, sizeof(double));
        int *recv_completed = (int *)calloc(number_of_nodes, sizeof(int));
        char **all_reports = (char **)malloc(number_of_nodes * sizeof(char *));

        for (int i = 0; i < number_of_nodes; i++)
        {
            all_reports[i] = (char *)malloc(256 * sizeof(char));
        }

// receiving reports from grid nodes
#pragma omp parallel for num_threads(number_of_nodes)
        for (int i = 0; i < number_of_nodes; i++)
        {
            MPI_Irecv(all_reports[i], 256, MPI_CHAR, i + 1, 2, world_comm, &receive_request[i]);
            if (DEBUG)
            {
                printf("[BASE] [ITER %d] waiting to receive grid node %d.\n", iter, i);
                fflush(stdout);
            }
        }

        // wait cycle duration
        while (current_time - start_time < CYCLE_DURATION_SECONDS)
        {
#pragma omp parallel for num_threads(number_of_nodes)
            for (int i = 0; i < number_of_nodes; i++)
            {
                // testing if report is received
                if (!recv_completed[i])
                {
                    MPI_Test(&receive_request[i], &recv_completed[i], &receive_status[i]);
                    if (recv_completed[i])
                    {
                        report_time[i] = MPI_Wtime();
                    }
                }
            }
            usleep(500000);
            current_time = (int)time(NULL);
        }

        // store report
        int *availability_table = (int *)malloc(number_of_nodes * sizeof(int));
        for (int i = 0; i < number_of_nodes; i++)
        {
            availability_table[i] = CHARGING_PORTS;
        }

        for (int i = 0; i < number_of_nodes; i++)
        {
            // cancel receive requests that are did not receive any message
            if (!recv_completed[i])
            {
                MPI_Cancel(&receive_request[i]);
            }
            // reading report from nodes
            else
            {
                char *time_str;
                char *availability_str;
                char *neighbors_availability_str;
                char *send_time_str;
                char delimiter[] = "/";
                char array_delimiter[] = ",";

                time_str = strtok(all_reports[i], delimiter);
                availability_str = strtok(NULL, delimiter);
                neighbors_availability_str = strtok(NULL, delimiter);
                send_time_str = strtok(NULL, delimiter);

                if (DEBUG)
                {
                    printf("[BASE] report received from grid node rank: %d\n[REPORT] node: %d // time: %s, availability: %s // neighbors availability: %s \n", i, i, time_str, availability_str, neighbors_availability_str);
                    // printf("[BASE] report received from grid node rank: %d\n[REPORT] node: %d // %s\n", i, i, all_reports[i]);
                    fflush(stdout);
                }

                alert_time[i] = atoi(time_str);
                availability_table[i] = atoi(availability_str);
                double send_time;
                sscanf(send_time_str, "%lf", &send_time);
                elapsed_time[i] = report_time[i] - send_time;
                int neighbors_availability[4];
                for (int j = 0; j < 4; j++)
                {
                    if (j == 0)
                    {
                        neighbors_availability[j] = atoi(strtok(neighbors_availability_str, array_delimiter));
                    }
                    else
                    {
                        neighbors_availability[j] = atoi(strtok(NULL, array_delimiter));
                    }
                    if (neighbors[i][j] >= 0)
                    {
                        availability_table[neighbors[i][j]] = neighbors_availability[j];
                    }
                }
            }
        }

        // computing if available neighbor of neighbors exist
        int **available_nons = (int **)malloc(number_of_nodes * sizeof(int *));
        for (int i = 0; i < number_of_nodes; i++)
        {
            available_nons[i] = (int *)malloc(8 * sizeof(int));
            for (int j = 0; j < 8; j++)
            {
                available_nons[i][j] = -2;
            }
        }

#pragma omp parallel for num_threads(number_of_nodes)
        for (int i = 0; i < number_of_nodes; i++)
        {
            // if (availability_table[i] <= LOW_PORTS)
            if (recv_completed[i])
            {
                int counter = 0;
                for (int j = 0; j < 8; j++)
                {
                    if (neighbor_of_neighbors[i][j] >= 0 && availability_table[neighbor_of_neighbors[i][j]] > LOW_PORTS)
                    {
                        int exists = 0;
                        for (int k = 0; k < 8; k++)
                        {
                            if (available_nons[i][k] == neighbor_of_neighbors[i][j])
                            {
                                exists = 1;
                            }
                        }
                        if (!exists)
                        {
                            available_nons[i][counter] = neighbor_of_neighbors[i][j];
                            counter++;
                        }
                    }
                }
                int buf[8] = {-1};
                for (int j = 0; j < counter; j++)
                {
                    buf[j] = available_nons[i][j];
                }
                MPI_Request send_request_node;
                MPI_Status send_status_node;
                int reply_time = MPI_Wtime();
                MPI_Isend(buf, 8, MPI_INT, i + 1, 3, world_comm, &send_request_node);
                MPI_Wait(&send_request_node, &send_status_node);
                elapsed_time[i] += MPI_Wtime() - reply_time;
            }
        }

        // logging
        for (int i = 0; i < number_of_nodes; i++)
        {
            if (recv_completed[i])
            {
                fprintf(file, "------------------------------------------------------------------------------------------------------------\n");
                fprintf(file, "Iteration: %d\n", iter + 1);
                char log_time_str[50];
                time_t current_time;
                time(&current_time);
                strftime(log_time_str, sizeof(log_time_str), "%A, %Y-%m-%d %H:%M:%S", gmtime(&current_time));
                char alert_time_str[50];
                time_t alert = (time_t)alert_time[i];
                strftime(alert_time_str, sizeof(alert_time_str), "%A, %Y-%m-%d %H:%M:%S", gmtime(&alert));
                fprintf(file, "Logged time: %s\nAlert reported time (Local time): %s\n", log_time_str, alert_time_str);
                fprintf(file, "Number of adjacent nodes: %d\nAvailability to be considered full: %d\n\n", number_of_neighbors[i], LOW_PORTS);
                fprintf(file, "Reporting Node \tCoord \tPort Value \tAvailable Port \n%d \t\t\t\t(%d, %d) \t%d \t\t\t%d\n\n", i, coords[i][0], coords[i][1], CHARGING_PORTS, availability_table[i]);
                fprintf(file, "Adjacent Nodes \tCoord \tPort Value \tAvailable Port \n");
                for (int j = 0; j < 4; j++)
                {
                    if (neighbors[i][j] >= 0)
                    {
                        fprintf(file, "%d \t\t\t\t(%d, %d) \t%d \t\t\t%d\n", neighbors[i][j], coords[neighbors[i][j]][0], coords[neighbors[i][j]][1], CHARGING_PORTS, availability_table[neighbors[i][j]]);
                    }
                }
                fprintf(file, "\nNearby Nodes \tCoord \n");
                for (int j = 0; j < 8; j++)
                {
                    if (neighbor_of_neighbors[i][j] >= 0)
                    {
                        fprintf(file, "%d \t\t\t\t(%d, %d)\n", neighbor_of_neighbors[i][j], coords[neighbor_of_neighbors[i][j]][0], coords[neighbor_of_neighbors[i][j]][1]);
                    }
                }
                fprintf(file, "\nAvailable station nearby: ");
                for (int j = 0; j < 8; j++)
                {
                    if (available_nons[i][j] >= 0)
                    {
                        if (j == 0)
                        {
                            fprintf(file, "%d", available_nons[i][j]);
                        }
                        else
                        {
                            fprintf(file, ", %d", available_nons[i][j]);
                        }
                    }
                }
                fprintf(file, "\nCommunication Time (seconds): %.10f \nTotal Messages send between reporting node and base station: 2 \n\n", elapsed_time[i]);
                fflush(file);
            }
        }

        // free allocated memory
        free(recv_completed);
        for (int i = 0; i < number_of_nodes; i++)
        {
            free(all_reports[i]);
        }
        free(all_reports);
        free(alert_time);
        free(report_time);
        free(elapsed_time);
        free(availability_table);
        for (int i = 0; i < number_of_nodes; i++)
        {
            free(available_nons[i]);
        }
        free(available_nons);

        iter += 1;
        MPI_Barrier(world_comm);

        // check for termination requirements and broadcast
        if (iter == TERMINATE)
        {
            term = 1;
            printf("[BASE] TERMINATE\n");
            fflush(stdout);
        }
        MPI_Bcast(&term, 1, MPI_INT, 0, world_comm);
    }

    // free allocated memory
    for (int i = 0; i < number_of_nodes; i++)
    {
        free(neighbors[i]);
    }
    free(neighbors);
    free(number_of_neighbors);
    for (int i = 0; i < number_of_nodes; i++)
    {
        free(coords[i]);
    }
    free(coords);
    for (int i = 0; i < number_of_nodes; i++)
    {
        free(neighbor_of_neighbors[i]);
    }
    free(neighbor_of_neighbors);

    fclose(file);
    return 0;
}

int charging_node_io(MPI_Comm world_comm, MPI_Comm comm, int m, int n)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    // creating 2d grid topology
    int dims[2] = {m, n};    // (3x3)
    int periods[2] = {0, 0}; // non-periodic
    int grid_rank;

    MPI_Comm GRID_COMM, ROW_COMM, COL_COMM;
    MPI_Cart_create(comm, 2, dims, periods, 1, &GRID_COMM);
    MPI_Comm_rank(GRID_COMM, &grid_rank);

    // getting coordinates in grid
    int coords[2];
    MPI_Cart_coords(GRID_COMM, grid_rank, 2, coords);

    if (DEBUG)
    {
        printf("[MPI process %d] I am now MPI process %d in GRID. (%d,%d)\n", my_rank, grid_rank, coords[0], coords[1]);
        fflush(stdout);
    }

    // creating row and col comm
    int free_coords[2];

    free_coords[0] = 0;
    free_coords[1] = 1;
    MPI_Cart_sub(GRID_COMM, free_coords, &ROW_COMM);

    free_coords[0] = 1;
    free_coords[1] = 0;
    MPI_Cart_sub(GRID_COMM, free_coords, &COL_COMM);

    // initialize avalability table and variables
    srand((int)time(NULL) * my_rank);
    int k = CHARGING_PORTS;
    int table_size = TABLE_SIZE;
    int **availability_table = (int **)malloc(table_size * sizeof(int *));
    for (int i = 0; i < table_size; i++)
    {
        availability_table[i] = (int *)calloc(2, sizeof(int));
    }
    int iter = 0;
    int start_time;
    int current_time = (int)time(NULL);
    int *ports_in_use = (int *)calloc(k, sizeof(int));

    // find neighbors
    int nbr_i_lo, nbr_i_hi;
    int nbr_j_lo, nbr_j_hi;

    MPI_Cart_shift(GRID_COMM, 0, 1, &nbr_i_lo, &nbr_i_hi);
    MPI_Cart_shift(GRID_COMM, 1, 1, &nbr_j_lo, &nbr_j_hi);
    int neighbors[4] = {nbr_i_lo, nbr_i_hi, nbr_j_lo, nbr_j_hi};

    // send neighbors and coords to base
    MPI_Send(neighbors, 4, MPI_INT, 0, 4, world_comm);
    if (DEBUG)
    {
        printf("[GRID %d] send neighbors to BASE.\n", grid_rank);
        fflush(stdout);
    }
    MPI_Send(coords, 2, MPI_INT, 0, 5, world_comm);
    if (DEBUG)
    {
        printf("[GRID %d] send coords to BASE.\n", grid_rank);
        fflush(stdout);
    }

    int term = 0;
    while (!term)
    {
        start_time = current_time;

        // update availability table
        availability_table[iter][0] = current_time;
        availability_table[iter][1] = k;

        int use[k];

// using omp threads to simulate each port
#pragma omp parallel for num_threads(k)
        for (int i = 0; i < k; i++)
        {

            if (ports_in_use[i] == 0)
            {
                int dice = rand() % 10;
                use[i] = (dice < 7) ? 1 : 0;
                // 70 chance to use available ports
                if (use[i] == 1)
                {
                    ports_in_use[i] = USE_DURATION;
                }
            }
            else
            {
                use[i] = 1;
                ports_in_use[i] -= 1;
            }
#pragma omp critical
            {
                availability_table[iter][1] -= use[i];
            }
        }

        if (DEBUG)
        {
            // printing table
            printf("\n");
            fflush(stdout);
            printf("Availability Table: %d\n", grid_rank);
            fflush(stdout);
            for (int i = 0; i < 5; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    printf("%d\t", availability_table[i][j]);
                    fflush(stdout);
                }
                printf("\n");
                fflush(stdout);
            }
            printf("Ports in use:\t");
            fflush(stdout);
            for (int i = 0; i < k; i++)
            {
                printf("%d\t", ports_in_use[i]);
                fflush(stdout);
            }
            printf("\n");
            fflush(stdout);
        }

        MPI_Barrier(world_comm);

        // variables to request data from neighbor
        MPI_Request send_request[4];
        MPI_Request receive_request[4];
        MPI_Status send_status[4];
        MPI_Status receive_status[4];

        // variables to listen request from neighbor
        MPI_Request receive_request_nbr[4];
        MPI_Status receive_status_nbr[4];
        MPI_Request send_request_nbr[4];
        MPI_Status send_status_nbr[4];
        char recvValues_nbr[4] = {-1, -1, -1, -1};
        int *recvValues = (int *)malloc(4 * sizeof(int));
        for (int i = 0; i < 4; i++)
        {
            recvValues[i] = -1;
        }
        int recv_completed_nbr[4] = {0, 0, 0, 0};
        int recv_completed[4] = {0, 0, 0, 0};

        // listening to message from neighbors
        for (int i = 0; i < 4; i++)
        {
            if (neighbors[i] != MPI_PROC_NULL)
            {
                MPI_Irecv(&recvValues_nbr[i], 1, MPI_CHAR, neighbors[i], 0, GRID_COMM, &receive_request_nbr[i]);
            }
        }

        for (int i = 0; i < 4; i++)
        {
            if (neighbors[i] != MPI_PROC_NULL)
            {
                char msg;
                // prompt neighbors if all ports in use
                if (availability_table[iter % table_size][1] <= LOW_PORTS)
                {
                    msg = 'L';
                    MPI_Isend(&msg, 1, MPI_CHAR, neighbors[i], 0, GRID_COMM, &send_request[i]);
                    MPI_Wait(&send_request[i], &send_status[i]);
                    // if (DEBUG)
                    // {
                    //     printf("[GRID %d] request sent to neighbor rank %d.\n", grid_rank, neighbors[i]);
                    //     fflush(stdout);
                    // }
                }
                // else send a message to be ignored
                else
                {
                    msg = 'N';
                    MPI_Isend(&msg, 1, MPI_CHAR, neighbors[i], 0, GRID_COMM, &send_request[i]);
                    MPI_Wait(&send_request[i], &send_status[i]);
                }
            }
        }

        // cycle duration first half
        while (current_time - start_time < (CYCLE_DURATION_SECONDS / 2))
        {
            for (int i = 0; i < 4; i++)
            {
                if (neighbors[i] != MPI_PROC_NULL && !recv_completed_nbr[i])
                {
                    // checking if data request received from neighbor
                    MPI_Test(&receive_request_nbr[i], &recv_completed_nbr[i], &receive_status_nbr[i]);
                    // sending data if data request received
                    if (recv_completed_nbr[i])
                    {
                        if (recvValues_nbr[i] == 'L')
                        {
                            // if (DEBUG)
                            // {
                            //     printf("[GRID %d] message received from neighbor rank: %d ||| %c\n", grid_rank, receive_status_nbr[i].MPI_SOURCE, recvValues_nbr[i]);
                            //     fflush(stdout);
                            // }
                            MPI_Isend(&availability_table[iter % table_size][1], 1, MPI_INT, neighbors[i], 1, GRID_COMM, &send_request_nbr[i]);
                            MPI_Wait(&send_request_nbr[i], &send_status_nbr[i]);
                            // if (DEBUG)
                            // {
                            //     printf("[GRID %d] replied to neighbor rank: %d ||| %d\n", grid_rank, neighbors[i], availability_table[iter % table_size][1]);
                            //     fflush(stdout);
                            // }
                        }
                    }
                }
            }
            // usleep(500000);
            current_time = (int)time(NULL);
        }

        // syncing all nodes
        MPI_Barrier(GRID_COMM);

        for (int i = 0; i < 4; i++)
        {
            if (neighbors[i] != MPI_PROC_NULL)
            {
                // receiving neighbor data
                if (availability_table[iter % table_size][1] <= LOW_PORTS)
                {
                    MPI_Irecv(&recvValues[i], 1, MPI_INT, neighbors[i], 1, GRID_COMM, &receive_request[i]);
                    MPI_Wait(&receive_request[i], &receive_status[i]);
                }
            }
        }

        // cycle duration second half
        int report_sent = 0;
        int nearest_availability = LOW_PORTS;
        int nearest_available_node = -1;
        while (current_time - start_time < CYCLE_DURATION_SECONDS)
        {
            for (int i = 0; i < 4; i++)
            {
                if (neighbors[i] != MPI_PROC_NULL && !recv_completed[i] && availability_table[iter % table_size][1] <= LOW_PORTS)
                {
                    MPI_Test(&receive_request[i], &recv_completed[i], &receive_status[i]);
                }
                if (availability_table[iter % table_size][1] <= LOW_PORTS && !report_sent)
                {
                    int all_received = 1;
                    // check if all neighbor replied with their data
                    if (neighbors[i] != MPI_PROC_NULL && !recv_completed[i])
                    {
                        all_received = 0;
                    }
                    // send report to base station
                    if (all_received)
                    {
                        // verify neighbor availability
                        for (int j = 0; j < 4; j++)
                        {
                            if (neighbors[j] != MPI_PROC_NULL)
                            {
                                // if (DEBUG)
                                // {
                                //     printf("[GRID %d] Comparison neighbor rank %d: %d > %d nearest availability\n", grid_rank, neighbors[j], recvValues[j], nearest_availability);
                                //     fflush(stdout);
                                // }
                                if (recvValues[j] > nearest_availability)
                                {
                                    nearest_availability = recvValues[j];
                                    nearest_available_node = neighbors[j];
                                    // if (DEBUG)
                                    // {
                                    //     printf("[GRID %d] Comparison success, nearest available node: %d\n", grid_rank, nearest_available_node);
                                    //     fflush(stdout);
                                    // }
                                }
                            }
                        }
                        if (nearest_availability <= LOW_PORTS)
                        {
                            MPI_Request send_request_report;
                            MPI_Status send_status_report;
                            // creating report in string
                            char report_buf[256];
                            char temp[50];
                            report_buf[0] = '\0';
                            sprintf(temp, "%d", availability_table[iter % table_size][0]);
                            strcat(report_buf, temp);
                            strcat(report_buf, "/");
                            sprintf(temp, "%d", availability_table[iter % table_size][1]);
                            strcat(report_buf, temp);
                            strcat(report_buf, "/");
                            for (int j = 0; j < 4; j++)
                            {
                                sprintf(temp, "%d", recvValues[j]);
                                strcat(report_buf, temp);
                                strcat(report_buf, ",");
                            }
                            strcat(report_buf, "/");
                            sprintf(temp, "%.10f", MPI_Wtime());
                            strcat(report_buf, temp);
                            strcat(report_buf, "///");

                            if (DEBUG)
                            {
                                printf("[GRID %d] sent report: %s\n", grid_rank, report_buf);
                                fflush(stdout);
                            }

                            // sending report
                            MPI_Isend(report_buf, strlen(report_buf + 1), MPI_CHAR, BASE, 2, world_comm, &send_request_report);
                            MPI_Wait(&send_request_report, &send_status_report);
                            report_sent = 1;
                        }
                    }
                }
            }
            // usleep(500000);
            current_time = (int)time(NULL);
        }

        // printing nearest available node
        printf("[GRID %d] currently have %d availability(s).\n", grid_rank, availability_table[iter % table_size][1]);
        fflush(stdout);
        if (availability_table[iter % table_size][1] <= LOW_PORTS)
        {
            // if neighbors all unavailable, receive neighbor of neighbor from base station
            if (nearest_availability <= 0)
            {
                printf("[GRID %d] CURRENTLY LOW ||| neighbors all unavailable.\n", grid_rank);
                fflush(stdout);
                MPI_Request receive_request_base;
                MPI_Status receive_status_base;
                int buf[8];
                MPI_Irecv(buf, 8, MPI_INT, 0, 3, world_comm, &receive_request_base);
                MPI_Wait(&receive_request_base, &receive_status_base);
                nearest_available_node = buf[0];
                if (nearest_available_node >= 0)
                {
                    nearest_availability = CHARGING_PORTS;
                }
                if (nearest_availability <= LOW_PORTS)
                {
                    printf("[GRID %d] CURRENTLY LOW ||| no nearest available nodes (including neighbor of neighbors).\n", grid_rank);
                    fflush(stdout);
                }
                else
                {
                    printf("[GRID %d] CURRENTLY LOW ||| nearest available node (neighbor of neighbor): %d.\n", grid_rank, nearest_available_node);
                    fflush(stdout);
                }
            }
            else
            {
                printf("[GRID %d] CURRENTLY LOW ||| nearest available node: %d.\n", grid_rank, nearest_available_node);
                fflush(stdout);
            }
        }

        // increment cycle iteration
        iter = (iter + 1) % table_size;

        // free allocated memory
        free(recvValues);

        // sync base stations and all nodes
        MPI_Barrier(world_comm);

        // check for termination
        MPI_Bcast(&term, 1, MPI_INT, 0, world_comm);
        if (term)
        {
            printf("[GRID %d] TERMINATE\n", grid_rank);
            fflush(stdout);
        }
    }

    // free allocated memory
    for (int i = 0; i < table_size; i++)
    {
        free(availability_table[i]);
    }
    free(availability_table);
    free(ports_in_use);

    MPI_Comm_free(&GRID_COMM);
    return 0;
}

int extra_node_io(MPI_Comm world_comm, MPI_Comm comm)
{
    int iter = 0;
    int term = 0;
    while (!term)
    {
        MPI_Barrier(world_comm);
        iter += 1;
        MPI_Barrier(world_comm);
        // check for termination
        MPI_Bcast(&term, 1, MPI_INT, 0, world_comm);
    }
    return 0;
}