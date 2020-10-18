// HelloMPI.cpp : This file contains the 'main' function. Program execution begins and ends there.
// -

#include <iostream>
#include <iomanip>   
#include<math.h>
#include<string>
#include<vector>
#include<fstream>
#include<mpi.h>
using namespace std;

void getInput(string& fileA, string& fileG, int& N, int *&A, int *&B, int *&C, int rank2D, int& np, MPI_Comm comm) {
	ifstream sourceA, sourceG;

	/* Path for fileA */
	string pathA = "G:/github/PC-Assignment1/";
	pathA += fileA;
  // pathA = fileA;
	sourceA.open(pathA);

	/* Path for fileB */
	string pathG = "G:/github/PC-Assignment1/";
	pathG += fileG;
  // pathG = fileG;
	sourceG.open(pathG);

	string line;
	getline(sourceG, line);

	int Nodes = stoi(line);
	int p = sqrt(np);
	N = Nodes / p;

	int t = N * N;
	A = (int*)malloc(t * sizeof(int));
	B = (int*)malloc(t * sizeof(int));

	if (!rank2D) {

		vector<vector<int>> AMat, GMat;
		
		for (int i = 0; i < Nodes; i++) {
			vector<int> row;
			int edge;
			for (int j = 0; j < Nodes; j++) sourceG >> edge, row.push_back(edge);
			GMat.push_back(row);
		}

		getline(sourceA, line);
		for (int i = 0; i < Nodes; i++) {
			vector<int> row;
			int edge;
			for (int j = 0; j < Nodes; j++) sourceA >> edge, row.push_back(edge);
			AMat.push_back(row);
		}
	
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				A[i * N + j] = AMat[i][j];
				B[i * N + j] = GMat[i][j];
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
				for (int j = startJ; j < endJ; j++) BlockMatrix[(i - startI) * N + (j - startJ)] = AMat[i][j];
			}
			
			MPI_Send(BlockMatrix, N * N, MPI_INT, rank, 0, comm);

			for (int i = startI; i < endI; i++) {
				for (int j = startJ; j < endJ; j++) BlockMatrix[(i - startI) * N + (j - startJ)] = GMat[i][j];
			}

			MPI_Send(BlockMatrix, N * N, MPI_INT, rank, 1, comm);
		}
	}
	else {
		MPI_Recv(A, N * N, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
		MPI_Recv(B, N * N, MPI_INT, 0, 1, comm, MPI_STATUS_IGNORE);
	}
}

void printOutput(int *&C, int &N, int rank2D, int np, MPI_Comm comm) {
	if (!rank2D) {
		int *Block;
		int p = sqrt(np);
		vector<vector<int>> Result(N * p, vector<int>(N * p, 0));
		Block = (int*)malloc((N * N) * sizeof(int));
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) Result[i][j] = C[i * N + j];
		}
		for (int rank = 1; rank < np; rank++) {
			MPI_Recv(Block, N * N, MPI_INT, rank, 0, comm, MPI_STATUS_IGNORE);

			int startI = (rank / p) * N;
			int startJ = (rank % p) * N;
			int endI = startI + N;
			int endJ = startJ + N;

			for (int i = startI; i < endI; i++) {
				for (int j = startJ; j < endJ; j++) {
					Result[i][j] = Block[(i - startI) * N + j - startJ];
				}
			}
		}
		
		fstream sourceA;
		string pathA = "G:/github/PC-Assignment1/";
		pathA += "A.txt";
    // pathA = "A.txt";
		sourceA.open(pathA, ofstream::out | ofstream::trunc);
		sourceA << N * p;
		sourceA << "\n";
		for (int i = 0; i < N * p; i++) {
			for (int j = 0; j < N * p; j++) sourceA << Result[i][j] << ' ';
			sourceA << '\n';
		}
		sourceA.close();
	}
	else {
		MPI_Send(C, N * N, MPI_INT, 0, 0, comm);
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

	// Files
	string fileA(argv[1]);
	string fileG(argv[2]);

	// Matrices
	int N; // Block Size
	int *A, *B, *C;
	getInput(fileA, fileG, N, A, B, C, rank2D, np, comm);
	C = (int*)malloc((N * N) * sizeof(int));
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) C[i * N + j] = 0;
	}

	MPI_Barrier(comm);

	/* Compute the Rank of Left block and Up block and store it 
	So, we dont need to compute them again and again while doing shifts*/
	MPI_Cart_shift(comm, 1, -1, &rankRight, &rankLeft); // 2nd argument is 0 (x-axis), rankRight is source and rankLeft is destination
	MPI_Cart_shift(comm, 0, -1, &rankDown, &rankUp); // 2nd argument is 1 (y-axis), rankDown is source and rankUp is destination

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

	if (!rank2D) {
		cout << setprecision(10) << fixed << "Elapsed time by Rank 0 : " << elapsed_time << '\n';
	}

	printOutput(C, N, rank2D, np, comm);

	MPI_Finalize();
	return 0;
}
