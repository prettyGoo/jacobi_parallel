#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
// #include <chrono>
#include <mpi.h>

// int world_size;
// int world_rank;
#define MAX_FILE_OUTPUT 10

using namespace std;

char filename[11] = "output.txt";
int iterations = 0;



double** _createJacobi(int N) {

	N += 2;
	double** J = new double*[N];

	for (int index = 0; index < N; index++) {
		J[index] = new double[N];
	}

	return J;
}

void _initJacobi(double** J, int N, double init_value)
{
	N += 2;

	for (int column = 1; column < N-1; column++) { //init top and bottom borders
		J[0][column] = init_value;
		J[N - 1][column] = init_value;
	}

	J[0][0] = 0; J[0][N-1] = 0; J[N-1][0] = 0; J[N-1][N-1] = 0; //all corners to 0

	for (int row = 1; row < N - 1; row++) // init left and right borders and J itself
	{
		for (int column = 0; column < N; column++) {
			if (column == 0 || column == N - 1) {
				J[row][column] = init_value;
			}
			else {
				J[row][column] = 0;
			}
		}
	}
}

void _outputJacobi(double** J, int N) {
	N += 2;

	for (int row = 0; row < N; row++)
	{
		if (row == 0 || row == N-1) {
			for (int column = 0; column < N; column++)
					cout << left << setw(10) << setprecision(5) << 1;
			cout << endl;
		}
		else {
			for (int column = 0; column < N; column++)
				if (column == 0)
					cout << left << setw(10) << setprecision(5) << 1;
				else if (column == N - 1)
					cout << left << setprecision(5) << 1;
				else
					cout << left << setw(10) << setprecision(5) << J[row][column];
			cout << endl;
		}
	}
	cout << endl;
}


void _deleteJacobi(double** J, int N) {
	N += 2;

	for (int count = 0; count < N; count++)
		delete[] J[count];

}


bool _jacobiIsSteady(double** J, double** nextJ, int N, double e)
{
	for (int row = 1; row < N - 1; row++)
	{
		for (int column = 1; column < N - 1; column++) {
			double A = ceil(10000 * J[row][column]) / 10000.0;
			double B = ceil(10000 * nextJ[row][column]) / 10000.0;
			if (abs(A - B) > e)
			{
				return false;
			}
		}
	}

	return true;
}


void _calculateJacobi(int N, double** current_J, double** final_J)
{
	for (int row = 1; row < N - 1; row++)
	{
		for (int column = 1; column < N - 1; column++) {
			current_J[row][column] = (final_J[row - 1][column] + final_J[row][column - 1] + final_J[row + 1][column] + final_J[row][column + 1]) / 4;
		}
	}
}


void _startJacobi(double** J, double** nextJ, int N, int n_iters, double e)
{
	N += 2;
	bool is_steady = false;
	int current_iter = 1;

	double** current_J = J;
	double** final_J = nextJ;

	while (current_iter <= n_iters && !is_steady)
	{
		iterations++;
		current_iter++;

		if (current_iter % 2 == 0)
		{
			_calculateJacobi(N, current_J, final_J);

		}
		else
		{
			_calculateJacobi(N, final_J, current_J);
		}

		is_steady = _jacobiIsSteady(current_J, final_J, N, e);
	}
	std::cout << "Iterations: " << current_iter-1 << "\n";
	_deleteJacobi(current_J, N - 2);
	_deleteJacobi(final_J, N - 2);
}
//===============================

int main(int argc, char* argv[])
{
	//remove(filename);
	int N;
	int n_iters;
	double epsila;
	double init_value;

	double t_1, t_2;


	MPI_Init(&argc, &argv);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	if (world_rank == 0) {
		std::cout << "Input the number of rows and columns: ";
		std::cin >> N;
		std::cout << "Input the max value of iterations: ";
		std::cin >> n_iters;

		std::cout << "Input precision of the calculation: ";
		std::cin >> epsila;
		std::cout << "Input an ititial value: ";
		std::cin >> init_value;


		t_1 =  MPI_Wtime();

  double** J;
	double** nextJ;

	J = _createJacobi(N);
	nextJ = _createJacobi(N);

	_initJacobi(J, N, init_value);
	_initJacobi(nextJ, N, init_value);

	_startJacobi(J, nextJ, N, n_iters, epsila);

	_outputJacobi(J, N);

	t_2 = MPI_Wtime();
	std::cout << "Total time: " << t_2 - t_1 << " ms" << std::endl;
  }

	MPI_Finalize();


	return 0;
}
