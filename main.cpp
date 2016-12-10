#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <mpi.h>


#define MAX_FILE_OUTPUT 10

using namespace std;

char filename[11] = "output.txt";
int iterations = 0;


double* create_1d_matrix(int size)
{
  double* J = new double[size^2];

  return J;
}
double** create_2d_matrix(int row_size, int column_size)
{
  double** J = new double*[row_size];
  for (int i=0; i<row_size; i++) {
    J[i] = new double[column_size];
  }

  return J;
}
void delete_1d_matrix(double* J)
{
  delete[] J;

  return;
}
void delete_2d_matrix(double** J, int row_size)
{
  for (int i=0; i<row_size; i++)
    delete[] J[i];

  return;
}
void transform_1d_to_2d(double* matrix_1d, double** matrix_2d, int row_size, int column_size)
{
  int column_index = 0;
  for (int i=0; i<row_size*column_size; i++) {
    matrix_2d[i/column_size][column_index] = matrix_1d[i];
    column_index++;

    if (column_index == column_size)
      column_index = 0;
  }

  return;
}


double* createAndInit(int row_size, int column_size, double init_value)
{
  //This must be 1d-array
  //TODO: replace with a 1d initialization
  double* J = create_1d_matrix(row_size*column_size);
  for (int i=0; i<row_size*column_size; i++) {
    J[i] = init_value;
  }

  return  J;
}
double** createAndInitSub(double* sub_J_1d, int row_size, int column_size)
{
  double** sub_J_2d = create_2d_matrix(row_size, column_size);
  transform_1d_to_2d(sub_J_1d, sub_J_2d, row_size, column_size);

  return sub_J_2d;
}
void calculateJacobi(double** current_J, double** next_J, int row_size, int column_size) //TODO next_J and curren_J: chagne places???
{
  for (int row = 1; row < row_size-1; row++)
	{
		for (int column = 1; column < column_size-1; column++) {
			next_J[row][column] = (current_J[row - 1][column] + current_J[row][column - 1] + current_J[row + 1][column] + current_J[row][column + 1]) / 4;
		}
	}

  return;
}
int jacobiIsSteady(double** current_J, double** next_J, int row_size, int column_size, double epsila)
{
  for (int row = 1; row < row_size - 1; row++)
	{
		for (int column = 1; column < column_size - 1; column++) {
			double A = ceil(10000 * current_J[row][column]) / 10000.0;
			double B = ceil(10000 * next_J[row][column]) / 10000.0;
			if (abs(A - B) > epsila)
				return 1;
		}
	}

	return 0;
}


void outputJacobi(double** J, int row_size, int column_size) {

	for (int row = 0; row < row_size; row++)
	{
		for (int column = 0; column < column_size; column++)
			if (column == 0)
				cout << left << setw(10) << setprecision(5) << J[row][column];
			else if (column == column_size - 1)
				cout << left << setprecision(5) << J[row][column];
			else
				cout << left << setw(10) << setprecision(5) << J[row][column];
		cout << endl;
	}
	cout << endl;
}
// void fileOutputJacobi(char* name, double** J, int N) {
// 	ofstream outputFile;
// 	outputFile.open(filename, ios::app);
//
// 	N += 2;
//
// 	outputFile << name << endl;
//
// 	for (int row = 0; row < N; row++)
// 	{
// 		for (int column = 0; column < N; column++)
// 			if (column == 0)
// 				outputFile << left << setw(12) << setprecision(5) << J[row][column];
// 			else if (column == N - 1)
// 				outputFile << left << setprecision(5) << J[row][column];
// 			else
// 				outputFile << left << setw(12) << setprecision(5) << J[row][column];
// 		outputFile << endl;
// 	}
// 	outputFile << endl;
//
// 	outputFile.close();
// }


int main()
{
	//remove(filename);
	// int N;
	// std::cout << "Input the number of rows and columns: ";
	// std::cin >> N;
	// int n_iters;
	// std::cout << "Input the max value of iterations: ";
	// std::cin >> n_iters;
	// double epsila;
	// std::cout << "Input precision of the calculation: ";
	// std::cin >> epsila;
	// double init_value;
	// std::cout << "Input an ititial value: ";
	// std::cin >> init_value;
	//auto t_1 = std::chrono::high_resolution_clock::now();

  int N = 4;
	int n_iters = 10;
  double init_value = 1.0;
	double epsila = 0.001;
	int init_int_params[2] = {N, n_iters};
  double init_double_params[2] = {init_value, epsila};

	int rows_per_process;

	int* sendcounts;
	int* displs;
	int sum = 0;
	int last_rows;

	double* sub_J;
	double* exchange_buffer_one;
  double* exchange_buffer_two;
  double* exchange_buffer_three;
	double* exchange_buffer_four;


  MPI_Init(NULL, NULL);

  int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  MPI_Bcast(&init_int_params, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&init_double_params, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //TODO separate function for init params
  N = init_int_params[0];
  n_iters = init_int_params[1];
  init_value = init_double_params[0];
  epsila = init_double_params[1];

  rows_per_process = N / world_size;
  last_rows = N % world_size;


  //TODO separate function
  sendcounts = new int[world_size];
  displs = new int[world_size];
  // calculate send counts and displacements for the main matrix splitting
  for (int i = 0; i < world_size; i++) {
      sendcounts[i] = (rows_per_process+1)*(N+2);
      if (last_rows > 0) {
          sendcounts[i] += N;
          last_rows--;
      }

      displs[i] = sum;
      sum += sendcounts[i];
  }


  sub_J = new double[sendcounts[world_rank]];
  exchange_buffer_one = new double[N];
  exchange_buffer_two = new double[N];
  exchange_buffer_three = new double[N];
  exchange_buffer_four = new double[N];


  double* J = NULL;
  int success_steady_root_check; // 1 for true, 0 for false
  int current_iter = 0;

  if (world_rank == 0) {
      J = createAndInit(N+2, N+2, 1);
      success_steady_root_check = 0;
  }

  MPI_Scatterv(J, sendcounts, displs, MPI_DOUBLE, sub_J, sendcounts[world_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  int must_continue = 1;
  int is_steady; // 0 means steady
  int skip = sendcounts[world_rank] - N; // for borders exchange
  int row_size = sendcounts[world_rank] / (N+2);
  int column_size = N+2;

  double** current_J = createAndInitSub(sub_J, row_size, column_size);
  double** next_J = createAndInitSub(sub_J, row_size, column_size);

  do {
    // exchange borders
    // if (world_rank % 2 == 0) {
    //   for (int i=0; i<N; i++) {
    //       exchange_buffer_one[i] = sub_J[i+skip];
    //       exchange_buffer_two[i] = sub_J[i+skip];
    //       exchange_buffer_three[i] = sub_J[i];
    //       exchange_buffer_four[i] = sub_J[i];
    //   }
    //   if (world_rank == 0) { //send and receive once for first sub_J
    //     MPI_Send(exchange_buffer_one, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD);
    //     MPI_Recv(exchange_buffer_two, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     for (int i=0; i<N; i++) {
    //         sub_J[i+skip] = exchange_buffer_two[i];
    //     }
    //   }
    //   else if (world_rank == world_size-1) {//send and receive once for last sub_J
    //     MPI_Send(exchange_buffer_four, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD);
    //     MPI_Recv(exchange_buffer_three, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     for (int i=0; i<N; i++) {
    //         sub_J[i] = exchange_buffer_three[i];
    //     }
    //   }
    //   else { //send and receive twice
    //     MPI_Send(exchange_buffer_four, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD);
    //     MPI_Recv(exchange_buffer_three, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     MPI_Send(exchange_buffer_one, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD);
    //     MPI_Recv(exchange_buffer_two, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     for (int i=0; i<N; i++) {
    //         sub_J[i] = exchange_buffer_three[i];
    //         sub_J[i+skip] = exchange_buffer_two[i];
    //     }
    //   }
    // }
    // else if (world_rank % 2 != 0 ) {
    //   for (int i=0; i<N; i++) {
    //       exchange_buffer_one[i] = sub_J[i];
    //       exchange_buffer_two[i] = sub_J[i];
    //       exchange_buffer_three[i] = sub_J[i+skip];
    //       exchange_buffer_four[i] = sub_J[i+skip];
    //   }
    //
    //   MPI_Send(exchange_buffer_two, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD);
    //   MPI_Recv(exchange_buffer_one, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //   for (int i=0; i<N; i++) {
    //       sub_J[i] = exchange_buffer_one[i];
    //   }
    //   if (world_rank != world_size-1) {
    //     MPI_Send(exchange_buffer_three, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD);
    //     MPI_Recv(exchange_buffer_four, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     for (int i=0; i<N; i++) {
    //         sub_J[i+skip] = exchange_buffer_four[i];
    //     }
    //   }
    // }

    if (current_iter % 2 == 0) {//TODO current_iter needs replacement???
      calculateJacobi(current_J, next_J, row_size, column_size);
      outputJacobi(next_J,row_size, column_size );
    }
		else {
      calculateJacobi(next_J, current_J, row_size, column_size);
      outputJacobi(current_J,row_size, column_size );
    }


		is_steady = jacobiIsSteady(current_J, next_J, row_size, column_size, epsila);
    std::cout << is_steady << "\n";


    MPI_Reduce(&is_steady, &success_steady_root_check, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    current_iter++;
    if (world_rank == 0 && (success_steady_root_check == 0 || current_iter == n_iters)) {
      must_continue = 0;
    }
    MPI_Bcast(&must_continue, 1, MPI_INT, 0, MPI_COMM_WORLD);
  } while (must_continue);

  if (world_rank == 0)
    std::cout << "Number of iterations: " << current_iter << "\n";

  MPI_Finalize();
	//cout << "Number of iterations: " << iterations << endl;
	//auto t_2 = std::chrono::high_resolution_clock::now();
	//std::cout << "Total time: " << std::chrono::duration<double, std::milli>(t_2 - t_1).count() << " ms " << std::endl;

	return 0;
}
