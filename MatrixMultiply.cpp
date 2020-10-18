// HelloMPI.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>   
#include<math.h>
#include<string>
#include<vector>
#include<fstream>
#include<mpi.h>
using namespace std;

void getInput(string& filename, int& N, int *&A, int *&B, int *&C, int rank2D, int& np, MPI_Comm comm) {
	ifstream source;
	string path = "G:/github/mpi_codes/HelloMPI/";
	path += filename;
	source.open(path);
	string line;

	getline(source, line);

	int Nodes = stoi(line);
	int p = sqrt(np);
	N = Nodes / p;

	int t = N * N;
	A = (int*)malloc(t * sizeof(int));
	B = (int*)malloc(t * sizeof(int));

	if (!rank2D) {

		vector<vector<int>> AdjMat;
		while (getline(source, line)) {
			vector<int> row;
			for (auto to : line) row.push_back(to - '0');
			AdjMat.push_back(row);
		}
	
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				A[i * N + j] = B[i * N + j] = AdjMat[i][j];
			}
		}


		for (int rank = 1; rank < np; rank++) {
			int startI = (rank / p) * N;
			int startJ = (rank % p) * N;
			int endI = startI + N;
			int endJ = startJ + N;

			int* BlockMatrix;
			BlockMatrix = (int*)malloc(t * sizeof(int));

			for (int i = startI; i < endI; i++) {
				for (int j = startJ; j < endJ; j++) BlockMatrix[(i - startI) * N + (j - startJ)] = AdjMat[i][j];
			}
			
			MPI_Send(BlockMatrix, N * N, MPI_INT, rank, 0, comm);
			MPI_Send(BlockMatrix, N * N, MPI_INT, rank, 1, comm);
		}
	}
	else {
		MPI_Recv(A, N * N, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
		MPI_Recv(B, N * N, MPI_INT, 0, 1, comm, MPI_STATUS_IGNORE);
	}
}

void MatrixMultiply(int* A, int* B, int* C, int& N) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				C[i * N + j] += A[i * N + k] * B[k * N + j];
			}
		}
	}
}

int main(int argc, char **argv) {

	/* Parameters */
	int rank, np;  // Rank of process and Number of process
	int dims[2], periods[2]; // dimension and period for 2D Cartesian Topology
	MPI_Status status;
	MPI_Comm comm; // new communicator for 2D
	int rank2D; // new rank in 2D
	int source, dest; // Stores block rank for shifting
	int rankLeft, rankRight, rankUp, rankDown; // Rank of left, right, up, down block (Required to shift the block by one)
	int coords[2]; // Coordinates of the block

	// Matrices
	string filename(argv[1]);
	
	/* Set up*/
	MPI_Init(NULL, NULL);

	/* Start Timer */
	MPI_Barrier(MPI_COMM_WORLD);
	double elapsed_time = -MPI_Wtime();

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	/* Dimension and Period for Cartesian topology */
	dims[0] = dims[1] = sqrt(np);
	periods[0] = periods[1] = 1;  // 1 means true, 0 means false (Periodic)

	/* Creating MPI Cartesian topology */
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm);

	/* Getting Rank and Coordinates of the block */
	MPI_Comm_rank(comm, &rank2D);
	MPI_Cart_coords(comm, rank2D, 2, coords);

	// Matrices
	int N; // Block Size
	int *A, *B, *C;
	getInput(filename, N, A, B, C, rank2D, np, comm);
	C = (int*)malloc((N * N) * sizeof(int));
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) C[i * N + j] = 0;
	}
	MPI_Barrier(comm);
	/* Compute the Rank of Left block and Up block and store it 
	So, we dont need to compute them again and again while doing shifts*/
	MPI_Cart_shift(comm, 1, -1, &rankRight, &rankLeft); // 2nd argument is 0 (x-axis), rankRight is source and rankLeft is destination
	MPI_Cart_shift(comm, 0, -1, &rankDown, &rankUp); // 2nd argument is 1 (y-axis), rankDown is source and rankUp is destination

	/*if (rank2D == 5) {
		cout << "Left : " << rankLeft << '\n';
		cout << "Right : " << rankRight << '\n';

		cout << "Up : " << rankUp << '\n';
		cout << "Down : " << rankDown << '\n';
	}*/

	/* Cannon's Algorithm*/

	/* Initial Alignment for A*/
	MPI_Cart_shift(comm, 1, -coords[0], &source, &dest); // Get source and dest rank for A
	MPI_Sendrecv_replace(A, N * N, MPI_INT, dest, 1, source, 1, comm, &status);

	/* Initial Alignment for B*/
	MPI_Cart_shift(comm, 0, -coords[1], &source, &dest); // Get source and dest rank for B
	MPI_Sendrecv_replace(B, N * N, MPI_INT, dest, 1, source, 1, comm, &status);

	for (int shift = 0; shift < dims[0]; shift++) {
		MatrixMultiply(A, B, C, N);

		/* Shift block A by 1 to the left */
		MPI_Sendrecv_replace(A, N * N, MPI_INT, rankLeft, 1, rankRight, 1, comm, &status);
		/* Shift block B by 1 to up */
		MPI_Sendrecv_replace(B, N * N, MPI_INT, rankUp, 1, rankDown, 1, comm, &status);
	}

	/* Calculate Time */
	elapsed_time += MPI_Wtime();

	if (rank2D == 3) {
		cout << setprecision(10) << fixed << elapsed_time << '\n';
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) cout << C[i * N + j];
			cout << '\n';
		}
	}
	MPI_Finalize();
	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
