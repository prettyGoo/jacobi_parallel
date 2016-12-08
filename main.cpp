#include <iostream>
#include <stdlib.h>
#include <mpi.h>


double* init_jacobi(int N)
{
	double* J = new double[N*N];
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
	// 0 means steady
	if (steady_tester_process == 1 && steady_tester_value > 2500000) {
		return 0;
	}
	else if (steady_tester_process == 1) {
		return 1;
	}
	return 0;
}


int main()
{
	int N = 6;
	int rows_per_process;
	int elems_per_process;

	int* sendcounts;
	int* displs;
	int sum = 0;
	int last_rows;
	double* subJ;


	MPI_Init(NULL, NULL);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	elems_per_process = (N * N) / world_size;
	rows_per_process = N / world_size;
	last_rows = N % world_size;

	sendcounts = new int[world_size];
	displs = new int[world_size];

	//calculate send counts and displacements
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
					printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
			}
	}

	subJ = new double[sendcounts[world_size]];



	double* J = NULL;
	int success_steady_root_check;
	if (world_rank == 0) {
			J = init_jacobi(N);
			success_steady_root_check = 0;
	}

	MPI_Scatterv(J, sendcounts, displs, MPI_DOUBLE, subJ, sendcounts[world_size], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	int must_continue = 1;
	int is_steady;
	int steady_tester = 0;

	do {
		exchangeBorders();
		calculateJacobi();
		is_steady = isJcobiSteady(steady_tester, world_rank);
		// MPI_Barrier(MPI_COMM_WORLD);
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

	MPI_Finalize();

	return 0;
}
