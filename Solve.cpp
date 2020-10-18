#include<iostream>
#include <stdlib.h>   
#include<math.h>
#include<string>
#include<vector>
#include<fstream>
using namespace std;

void generateGraph(int Nodes) {
    vector<vector<int>> Graph(Nodes, vector<int>(Nodes, 0));
    for (int i = 0; i < Nodes; i++) {
      for (int j = i + 1; j < Nodes; j++) {
        Graph[i][j] = rand() % 2; // Binary Graph
        Graph[j][i] = Graph[i][j];
      } 
    }

    fstream sourceG, sourceA;
    sourceG.open("G.txt", ofstream::out | ofstream::trunc);
    sourceA.open("A.txt", ofstream::out | ofstream::trunc);
    sourceG << Nodes;
		sourceG << "\n";
    sourceA << Nodes;
		sourceA << "\n";
		for (int i = 0; i < Nodes; i++) {
			for (int j = 0; j < Nodes; j++) {
        sourceG << Graph[i][j] << ' ';
        sourceA << Graph[i][j] << ' ';
      }
			sourceG << '\n';
			sourceA << '\n';
		}
		sourceG.close();
		sourceA.close();
}

int getA(int S, int D) {
  ifstream sourceA;
  sourceA.open("A.txt");
  int Nodes;
  sourceA >> Nodes;
  vector<vector<int>> A(Nodes, vector<int>(Nodes, 0));
  for (int i = 0; i < Nodes; i++) {
    for (int j = 0; j < Nodes; j++) sourceA >> A[i][j];
  }
  
  return (A[S][D] ? 1 : 0);
}

int main() {
  int Nodes, S, D;
  cout << "Enter number of Nodes : \n";
  cin >> Nodes;
  cout << "Enter Source vertex: \n";
  cin >> S;
  if (S < 0 || S >= Nodes) {
    cout << "Sorce Vertex must be in between " << 0 << ' ' << Nodes - 1;
    exit(0);
  }
  cout << "Enter Destination vertex: \n";
  cin >> D;
  if (D < 0 || D >= Nodes) {
    cout << "Destination Vertex must be in between " << 0 << ' ' << Nodes - 1;
    exit(0);
  }
  if (S == D) {
    cout << "Path length is 0\n"; 
  }
  else {
    generateGraph(Nodes);
    int np;
    cout << "Enter number of process you want to run :";
    cin >> np;
    int p = sqrt(np);
    if (p * p != np) {
      cout << "Number of processes must be perfect square\n";
      exit(0);
    }
    if ((Nodes * Nodes) % np != 0) {
      cout << "Nodes must be divisible by sqrt(p)\n";
      exit(0);
    }

    string str = "mpiexec -np " + to_string(np) +" HelloMPI.exe A.txt G.txt";
    
    int length = 1;
    while (!getA(S, D) && length < Nodes) {
      const char *command = str.c_str(); 
      system(command);
      length++;
    } 

    if (length == Nodes) cout << "Path does not exist\n";
    else cout << "Path length is " << length << '\n';
  }
}