#include <iostream>
#include <vector>
#include "mpi.h"

using namespace std;

int MPI_BcastAnneau(
    void *buf,
    int count,
    MPI_Datatype type,
    int root,
    MPI_Comm comm)
{
  // A FAIRE ...
  return MPI_SUCCESS;
} 

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int root = 1; // Processus root qui possede le tableau entier
  int N = 16; // La taille du tableau a broadcaster
  std::vector<int> A(N); // Grand tableau, a initialiser que dans le processus root
  if (rank == root) { 
    for (int i = 0; i < N; i++) { A[i] = i; }
  }
  
  // Appeler MPI_BcastAnneau(...) sur A, et avec la root. Tester contre MPI_Bcast(...)
  // A FAIRE ...
//  MPI_BcastAnneau(...);
//  MPI_Bcast(...);

  // Afficher A pour chaque processus pour tester
  for (auto i : A) { cout << i << ":" << rank << endl; }

  MPI_Finalize();

  return 0;
}
