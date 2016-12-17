#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


double* init_jacobi(int N)
{
	double* J = new double[N*N];
	std::cout << "Initial matrix\n";
	for (int i=0; i<N*N; i++) {
			J[i] = i / N;
			std::cout << J[i] << " ";
	}

	return J;
}

void exchangeBorders()
{

}

void calculateJacobi()
{

}

//TODO: REPLACE WITH BOOL
int isJcobiSteady(int steady_tester_value, int steady_tester_process) {
	return 0; // 0 means steady
}


int main()
{
	int N = 5;
	int max_iter = 1;
	int epsila = 1;

	int initial_parameters[3] = {N, max_iter, epsila};

	int rows_per_process;
	int elems_per_process;

	int* sendcounts;
	int* displs;
	int sum = 0;
	int last_rows;

	double* subJ;
	double* exchange_buffer_one;
	double* exchange_buffer_two;
	double* exchange_buffer_three;
	double* exchange_buffer_four;


	MPI_Init(NULL, NULL);
	// MPI_Barrier(MPI_COMM_WORLD);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


	MPI_Bcast(&initial_parameters, 3, MPI_INT, 0, MPI_COMM_WORLD);
	N = initial_parameters[0];
	max_iter = initial_parameters[1];
	epsila = initial_parameters[2];

	rows_per_process = N / world_size;
	last_rows = N % world_size;

	sendcounts = new int[world_size];
	displs = new int[world_size];

	//calculate send counts and displacements for the main matrix splitting
	for (int i = 0; i < world_size; i++) {
			sendcounts[i] = N*rows_per_process;
			if (last_rows > 0) {
					sendcounts[i] += N;
					last_rows--;
			}

			displs[i] = sum;
			sum += sendcounts[i];
	}

	if (0 == world_rank) {
			for (int i = 0; i < world_size; i++) {
					printf("sendcounts[%d] = %d\trows[%d] = %d\n", i, sendcounts[i], i, sendcounts[i]/N);
			}
	}

	subJ = new double[sendcounts[world_rank]];
	exchange_buffer_one = new double[N];
	exchange_buffer_two = new double[N];
	exchange_buffer_three = new double[N];
	exchange_buffer_four = new double[N];

	double* J = NULL;
	int success_steady_root_check;
	if (world_rank == 0) {
			J = init_jacobi(N);
			success_steady_root_check = 0;
	}

	MPI_Scatterv(J, sendcounts, displs, MPI_DOUBLE, subJ, sendcounts[world_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	int must_continue = 1;
	int is_steady;
	int steady_tester = 0;
	int skip = sendcounts[world_rank] - N; // for borders exchange

	do {
		// exchange borders
		if (world_rank % 2 == 0) {
			for (int i=0; i<N; i++) {
					exchange_buffer_one[i] = subJ[i+skip];
					exchange_buffer_two[i] = subJ[i+skip];
					exchange_buffer_three[i] = subJ[i];
					exchange_buffer_four[i] = subJ[i];
			}
			if (world_rank == 0) { //send and receive once for first subJ
				MPI_Send(exchange_buffer_one, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD);
				MPI_Recv(exchange_buffer_two, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (int i=0; i<N; i++) {
						subJ[i+skip] = exchange_buffer_two[i];
				}
			}
			else if (world_rank == world_size-1) {//send and receive once for last subJ
				MPI_Send(exchange_buffer_four, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD);
				MPI_Recv(exchange_buffer_three, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (int i=0; i<N; i++) {
						subJ[i] = exchange_buffer_three[i];
				}
			}
			else { //send and receive twice
				MPI_Send(exchange_buffer_four, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD);
				MPI_Recv(exchange_buffer_three, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(exchange_buffer_one, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD);
				MPI_Recv(exchange_buffer_two, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (int i=0; i<N; i++) {
						subJ[i] = exchange_buffer_three[i];
						subJ[i+skip] = exchange_buffer_two[i];
				}
			}
		}
		else if (world_rank % 2 != 0 ) {
			for (int i=0; i<N; i++) {
					exchange_buffer_one[i] = subJ[i];
					exchange_buffer_two[i] = subJ[i];
					exchange_buffer_three[i] = subJ[i+skip];
					exchange_buffer_four[i] = subJ[i+skip];
			}

			MPI_Send(exchange_buffer_two, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD);
			MPI_Recv(exchange_buffer_one, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (int i=0; i<N; i++) {
					subJ[i] = exchange_buffer_one[i];
			}
			if (world_rank != world_size-1) {
				MPI_Send(exchange_buffer_three, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD);
				MPI_Recv(exchange_buffer_four, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (int i=0; i<N; i++) {
						subJ[i+skip] = exchange_buffer_four[i];
				}
			}
		}

		calculateJacobi();
		is_steady = isJcobiSteady(steady_tester, world_rank);
		steady_tester++;
		MPI_Reduce(&is_steady, &success_steady_root_check, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (world_rank == 0 && success_steady_root_check == 0) {
			must_continue = 0;
		}
		MPI_Bcast(&must_continue, 1, MPI_INT, 0, MPI_COMM_WORLD);
	} while (must_continue);


	double* final_J = NULL;
	if (world_rank == 0) {
			final_J = new double[N*N];
	}
	MPI_Gatherv(subJ, sendcounts[world_rank], MPI_DOUBLE, final_J, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	std::cout << "\n";
	if (world_rank == 0) {
		std::cout << "Final matrix after exhcange\n";
		for (int i=0; i<N*N; i++) {
			std::cout << final_J[i] << " ";
		}
		std::cout << "\n";
	}

	if (world_rank == 0) {
		delete[] J;
		delete[] final_J;
	}
	delete[] sendcounts;	delete[] displs; delete[] subJ;
	delete[] exchange_buffer_one;  delete[] exchange_buffer_two;
	delete[] exchange_buffer_three; delete[] exchange_buffer_four;

	MPI_Finalize();

	return 0;
}
