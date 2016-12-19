#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <mpi.h>


#define MAX_FILE_OUTPUT 10

using namespace std;

char filename[11] = "output.txt";
int iterations = 0;


double* create_1d_matrix(int size)
{
	double* J = new double[size];

	return J;
}

double** create_2d_matrix(int row_size, int column_size)
{
	double** J = new double*[row_size];
	for (int i = 0; i < row_size; i++) {
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
	for (int i = 0; i < row_size; i++)
		delete[] J[i];

	return;
}

void transform_1d_to_2d(double* matrix_1d, double** matrix_2d, int row_size, int column_size)
{
	int row_index = 0, column_index = 0;
	//cout << "Size: " << row_size * column_size << endl;
	for (int i = 0; i < row_size * column_size; i++) {
		matrix_2d[row_index][column_index] = matrix_1d[i];
		column_index++;

		if (column_index == column_size)
		{
			row_index++;
			column_index = 0;
		}
	}

	return;
}

void transform_2d_to_1d(double** matrix_2d, double* matrix_1d, int row_size, int column_size)
{
	int index = 0;
	for(int i = 0; i < row_size; i++)
		for(int j = 0; j < column_size; j++)
		{
			matrix_1d[index] = matrix_2d[i][j];
			index++;
		}

	return;
}

double* createAndInit(int row_size, int column_size, double init_value)
{
  //This must be 1d-array
	double* J = create_1d_matrix(row_size * column_size);


	int index = 0;

	for(int i = 0; i < row_size; i++) {
		if (i == 0) { //the upper border
			for(int j = 0; j < row_size; j++) {
				J[index] = init_value;
				index++;
			}
		}
		else if(i == row_size - 1) { //the lower border
			for(int j = 0; j < row_size; j++) {
				J[index] = init_value;
				index++;
			}
		}
		else {
			J[index] = init_value;
			index++;
			for(int j = 0; j < row_size - 2; j++) {
				J[index] = 0;
				index++;
			}
			J[index] = init_value;
			index++;
		}
	}

	return  J;
}

double** createAndInitSub(double* sub_J_1d, int row_size, int column_size, int extra_row_param)
{
  int rows;
	double** sub_J;
  if (extra_row_param == -1 || extra_row_param == 1) { //need one extra border
    sub_J = create_2d_matrix(row_size+1, column_size);
    rows = row_size + 1;
  }
  else if (extra_row_param == 0){ //need two extra borders
    sub_J = create_2d_matrix(row_size+2, column_size);
    rows = row_size + 2;
  }


  double** sub_J_2d = create_2d_matrix(row_size, column_size);
	transform_1d_to_2d(sub_J_1d, sub_J_2d, row_size, column_size); //1d to 2d without borders

  int sub_i=0; int sub_j=0;
  for (int i=0; i<rows; i++) {
    if (sub_j == column_size) {
      sub_i++;
      sub_j = 0;
    }
    for (int j=0; j<column_size; j++) {
      if ((i==0 && (extra_row_param==1 || extra_row_param==0)) || (i==rows-1 && (extra_row_param==-1 || extra_row_param==0))) {
				if (j==0 || j==column_size-1)
        	sub_J[i][j] = 1; // just nothins for now
				else
					sub_J[i][j] = 0;
      }
      else {
        sub_J[i][j] = sub_J_2d[sub_i][sub_j];
        sub_j++;
      }
    }
  }


	return sub_J; //subJ with borders
}

double* getFinalSub(double** J, int row_size, int column_size, int extra_row_param)
{//transforms 2d matrix to 1d and remove borders

	int rows;
	double* final_J_1d = create_1d_matrix(row_size*column_size);
	double** final_J_2d = create_2d_matrix(row_size, column_size);

	int index_tuner;
  if (extra_row_param == -1) // the bottom border must be removed
		index_tuner = 0;
	else if (extra_row_param == 1 || extra_row_param == 0) //remove top border
		index_tuner = 1;


	for (int i=0; i<row_size; i++) {
		for (int j=0; j<column_size; j++) {
			final_J_2d[i][j] = J[i+index_tuner][j];
		}
	}

	transform_2d_to_1d(final_J_2d, final_J_1d, row_size, column_size);
	return final_J_1d;
}

void calculateJacobi(double** current_J, double** next_J, int row_size, int column_size, int extra_row_param)
{
  if (extra_row_param == -1 || extra_row_param == 1)
    row_size +=1;
  else if (extra_row_param == 0)
    row_size +=2;

	for (int row = 1; row < row_size-1; row++)
	{
		for (int column = 1; column < column_size-1; column++) {
			next_J[row][column] = (current_J[row - 1][column] + current_J[row][column - 1] + current_J[row + 1][column] + current_J[row][column + 1]) / 4;
		}
	}

	return;
}

int jacobiIsSteady(double** current_J, double** next_J, int row_size, int column_size, int extra_row_param, double epsila)
{
  if (extra_row_param == -1 || extra_row_param == 1)
    row_size +=1;
  else if (extra_row_param == 0)
    row_size +=2;

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

void outputJacobi(double** J, int row_size, int column_size)
{
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

void fileOutputJacobi(char* name, double** J, int N) {
	ofstream outputFile;
	outputFile.open(filename, ios::app);

	N += 2;

	outputFile << name << endl;

	for (int row = 0; row < N; row++)
	{
		for (int column = 0; column < N; column++)
			if (column == 0)
				outputFile << left << setw(12) << setprecision(5) << J[row][column];
			else if (column == N - 1)
				outputFile << left << setprecision(5) << J[row][column];
			else
				outputFile << left << setw(12) << setprecision(5) << J[row][column];
		outputFile << endl;
	}
	outputFile << endl;

	outputFile.close();
}

int main()
{
	//remove(filename);
	int N;
	int n_iters;
	double epsila;
	double init_value;
	int init_int_params[2];
	double init_double_params[2];

	std::chrono::high_resolution_clock::time_point t_1;
	std::chrono::high_resolution_clock::time_point t_2;

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

	if (world_rank == 0) {
		std::cout << "Input the number of rows and columns: ";
		std::cin >> N;
		std::cout << "Input the max value of iterations: ";
		std::cin >> n_iters;
		init_int_params[0] = N; init_int_params[1] = n_iters;

		std::cout << "Input precision of the calculation: ";
		std::cin >> epsila;
		std::cout << "Input an ititial value: ";
		std::cin >> init_value;
		init_double_params[0] = init_value; init_double_params[1] = epsila;

		t_1 =  std::chrono::high_resolution_clock::now();
	}
	MPI_Bcast(&init_int_params, 2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&init_double_params, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	N = init_int_params[0];
	n_iters = init_int_params[1];
	init_value = init_double_params[0];
	epsila = init_double_params[1];

	rows_per_process = N / world_size;
	last_rows = N % world_size;

	sendcounts = new int[world_size];
	displs = new int[world_size];

  // calculate send counts and displacements for the main matrix splitting
	for (int i = 0; i < world_size; i++) {
		if (world_rank == 0 || world_rank == world_size-1)
			sendcounts[i] = (rows_per_process + 1) * (N + 2); //TODO: check if it works with more than two processes
		else
			sendcounts[i] = rows_per_process * (N + 2);

		if (last_rows > 0) {
			sendcounts[i] += (N+2);
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
		J = createAndInit(N + 2, N + 2, 1); //1 = init value
		success_steady_root_check = 0;
	}

  	MPI_Scatterv(J, sendcounts, displs, MPI_DOUBLE, sub_J, sendcounts[world_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  	int must_continue = 1;
  	int is_steady; // 0 means steady
  	int skip = sendcounts[world_rank] - N; // for borders exchange
		int column_size = N + 2;
  	int row_size = sendcounts[world_rank] / column_size;
    int extra_row_param; //determine where border will added

    if (world_rank == 0) //extra row will be added below
      extra_row_param = -1;
    else if (world_rank = world_size-1) //extra row will be added above
      extra_row_param = 1;
    else // both above and below extra rows are needed
      extra_row_param = 0;

  	double** current_J = createAndInitSub(sub_J, row_size, column_size, extra_row_param); // it has one or two additional rows
  	double** next_J = createAndInitSub(sub_J, row_size, column_size, extra_row_param); // it has one or two additional rows

  	do {

			if (world_rank == 0) {
				if (current_iter % 2 == 0) {
					for (int i=1; i<column_size-1; i++)
							exchange_buffer_one[i-1] = current_J[row_size-1][i];
				}
				else if (current_iter % 2 != 0) {
					for (int i=1; i<column_size-1; i++)
							exchange_buffer_one[i-1] = next_J[row_size-1][i];
				}
				MPI_Send(exchange_buffer_one, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD);
				MPI_Recv(exchange_buffer_two, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				if (current_iter % 2 == 0) {
					for (int i=1; i<column_size-1; i++)
							current_J[row_size][i] = exchange_buffer_two[i-1];
				}
				else if (current_iter % 2 != 0) {
					for (int i=1; i<column_size-1; i++)
							next_J[row_size][i] = exchange_buffer_two[i-1];
				}
			}
			else if (world_rank == world_size-1) {//send and receive once for last sub_J
				if (current_iter % 2 == 0) {
					for (int i=1; i<column_size-1; i++)
							exchange_buffer_three[i-1] = current_J[1][i];
				}
				else if (current_iter % 2 != 0) {
					for (int i=1; i<column_size-1; i++)
							exchange_buffer_three[i-1] = next_J[1][i];
				}
				MPI_Send(exchange_buffer_three, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD);
				MPI_Recv(exchange_buffer_four, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				if (current_iter % 2 == 0) {
					for (int i=1; i<column_size-1; i++)
							current_J[0][i] = exchange_buffer_four[i-1];
				}
				else if (current_iter % 2 != 0) {
					for (int i=1; i<column_size-1; i++)
							next_J[0][i] = exchange_buffer_four[i-1];
				}
			}
			else {
				if (current_iter % 2 == 0) {
					for (int i=1; i<column_size-1; i++) {
						exchange_buffer_one[i-1] = current_J[row_size-1][i];
						exchange_buffer_three[i-1] = current_J[1][i];
					}
				}
				else if (current_iter % 2 != 0) {
					for (int i=1; i<column_size-1; i++) {
						exchange_buffer_one[i-1] = next_J[row_size-1][i];
						exchange_buffer_three[i-1] = next_J[1][i];
					}
				}

				MPI_Send(exchange_buffer_one, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD);
				MPI_Recv(exchange_buffer_two, N, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(exchange_buffer_three, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD);
				MPI_Recv(exchange_buffer_four, N, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				if (current_iter % 2 == 0) {
					for (int i=1; i<column_size-1; i++) {
						current_J[row_size][i] = exchange_buffer_two[i-1];
						current_J[0][i] = exchange_buffer_four[i-1];
					}
				}
				else if (current_iter % 2 != 0) {
					for (int i=1; i<column_size-1; i++) {
						next_J[row_size][i] = exchange_buffer_two[i-1];
						next_J[0][i] = exchange_buffer_four[i-1];
					}
				}
			}

      //TODO Deal with row size
    	if (current_iter % 2 == 0) {
    		calculateJacobi(current_J, next_J, row_size, column_size, extra_row_param);
    		// outputJacobi(next_J,row_size, column_size );
    	}
    	else {
    		calculateJacobi(next_J, current_J, row_size, column_size, extra_row_param);
    		// outputJacobi(current_J,row_size, column_size);
    	}

    	is_steady = jacobiIsSteady(current_J, next_J, row_size, column_size, extra_row_param, epsila);

    	MPI_Reduce(&is_steady, &success_steady_root_check, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    	current_iter++;
    	if (world_rank == 0 && (success_steady_root_check == 0 || current_iter == n_iters)) {
    		must_continue = 0;
    	}

    	MPI_Bcast(&must_continue, 1, MPI_INT, 0, MPI_COMM_WORLD);
	} while (must_continue);

	double* final_sub;
  if (current_iter % 2 == 0)
    final_sub = getFinalSub(next_J, row_size, column_size, extra_row_param);
  else
    final_sub = getFinalSub(current_J, row_size, column_size, extra_row_param);

  double* final_J = NULL;
  if (world_rank == 0) {
    final_J = create_1d_matrix((N+2)*(N+2));
  }

  MPI_Gatherv(final_sub, sendcounts[world_rank], MPI_DOUBLE, final_J, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double** result = NULL;
  if (world_rank == 0) {
    result = create_2d_matrix(N+2, N+2);
    transform_1d_to_2d(final_J, result, N+2, N+2);
    outputJacobi(result, N+2, N+2);
    std::cout << "Iterations: " << current_iter << "\n";

		t_2 = std::chrono::high_resolution_clock::now();
		std::cout << "Total time: " << std::chrono::duration<double, std::milli>(t_2 - t_1).count() << " ms " << std::endl;
  }

	MPI_Finalize();


	return 0;
}
