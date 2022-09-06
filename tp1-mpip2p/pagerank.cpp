#include <cassert>
#include <iostream>
#include <random>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>
#include <cmath>
#include "mpi.h"

using namespace std;

int N = 1024, numElems = 5000;
const int maxIters = 1000;

vector<int> colBegin, rowIdx;
int procRank, numProcs;

void generateAdjMatrix()
{
  random_device rd;
  mt19937 gen(1);
  uniform_int_distribution<> dist(0, N - 1);
  set<pair<int, int>> edgeSet;
  for (int i = 0; i < numElems; i++) {
    std::pair<int, int> newEdge;
    do {
      newEdge.first = dist(gen) % N;
      newEdge.second = dist(gen) % N;
    } while (newEdge.second == newEdge.first || edgeSet.find(newEdge) != edgeSet.end());
    edgeSet.insert(newEdge);
  }
  colBegin = std::vector<int>(N + 1);
  for (auto &i : edgeSet) { colBegin[i.second]++; }
  for (int i = 1; i <= N; i++) { colBegin[i] += colBegin[i - 1]; }
  rowIdx = std::vector<int>(colBegin[N]);
  for (auto &i : edgeSet) { rowIdx[--colBegin[i.second]] = i.first; }
}

void printAdjMatrix()
{
  cout << N << " " << numElems << endl;
  for (int colIdx = 0; colIdx < N; colIdx++) {
    for (int rowPtr = colBegin[colIdx]; rowPtr < colBegin[colIdx + 1]; rowPtr++) {
      cout << rowIdx[rowPtr] << " " << colIdx << endl;
    }
  }
}

void calculatePageRankSeq()
{
  std::vector<double> xold(N);
  // Initialiser xold à 1/N
  for (int i = 0; i < N; i++) { xold[i] = 1.0 / N; }
  // Calculer le nombre d'elements par ligne (ce qui correspond a |A(i, :)| dans le calcul de PageRank)
  std::vector<int> rowCount(N);
  for (int j = 0; j < N; j++) {
    for (int rowPtr = colBegin[j]; rowPtr < colBegin[j + 1]; rowPtr++) { rowCount[rowIdx[rowPtr]]++; }
  }
  // Iterer jusqu'à la convergence ou l'atteinte du maxIters
  std::vector<double> x(N);
  for (int iter = 0; iter < maxIters; iter++) {
    cout << "iter = " << iter << endl;
    // Calculer le pageprocRank pour chaque x[j] à partir du vecteur xold
    for (int j = 0; j < N; j++) {
      x[j] = 0.0;
      for (int rowPtr = colBegin[j]; rowPtr < colBegin[j + 1]; rowPtr++) {
        int i = rowIdx[rowPtr];
        x[j] += xold[i] * (1.0 / rowCount[i]);
      }
    }
    // Calculer la norme de l'x
    double normX = 0.0;
    for (int j = 0; j < N; j++) { normX += x[j] * x[j]; }
    normX = sqrt(normX);
    cout << "norm(x) = " << normX << endl;
    // Normaliser x
    for (int j = 0; j < N; j++) { x[j] /= normX; }
    double normDiff = 0.0;
    // Calculer la norme du (x - xold)
    for (int j = 0; j < N; j++) { normDiff += (x[j] - xold[j]) * (x[j] - xold[j]); }
    normDiff = sqrt(normDiff);
    cout << "norm(x - xold) = " << normDiff << endl;
    // Verifier la convergence
    if (normDiff < 1e-6) { break; }
    // Remplacer xold avec x pour l'iteration suivante
    for (int i = 0; i < N; i++) { xold[i] = x[i]; }
  }
}

void calculatePageRankPar()
{
  // Initialiser xold à 1/N en parallele, chaque processus initialise N / P elements consecutives du xold. En suite, 
  // les processus communiquent xold avec MPI_Allgather.
  // A FAIRE ...
  std::vector<double> xold(N);

  // Calculer le nombre d'elements par ligne (ce qui correspond a |A(i, :)| dans le calcul de PageRank). Chaque
  // processus fait le calcul pour N / P lignes consecutives, puis les processus communiquent rowCount avec
  // MPI_Allreduce.
  // A FAIRE ...
  std::vector<int> rowCount(N);

  // Iterer jusqu'à la convergence ou l'atteinte du maxIters. Chaque processus calcul N / P elements du vecteur x.
  std::vector<double> x(N);
  for (int iter = 0; iter < maxIters; iter++) {
    if (procRank == 0) { cout << "iter = " << iter << endl; }
    // Calculer le pageprocRank en parallele pour chaque x[j] à partir du vecteur xold
    // A FAIRE ...

    // Calculer la norme de l'x en parallele
    double normX = 0.0;
    // A FAIRE ...
    if (procRank == 0) { cout << "norm(x) = " << normX << endl; }

    // Normaliser x en parallele
    // A FAIRE ...
    
    // Calculer la norme du (x - xold) en parallele
    double normDiff = 0.0;
    // A FAIRE ...
    if (procRank == 0) { cout << "norm(x - xold) = " << normDiff << endl; }

    // Verifier la convergence
    if (normDiff < 1e-6) { break; }

    // Remplacer xold avec x en parallele pour l'iteration suivante 
    // A FAIRE ...
  }
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  generateAdjMatrix();
//  printAdjMatrix();
  if (procRank == 0) { calculatePageRankSeq(); }
  MPI_Barrier(MPI_COMM_WORLD);
  calculatePageRankPar();
  MPI_Finalize();
  return 0;
