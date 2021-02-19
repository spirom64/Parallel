#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"

#define N 8192

int Override(int proc, int nb, int mb)
{
	int flag = ((proc + 1) * mb % nb);
	return ((flag == 0) && (mb < nb)) || ((flag != 0) && (flag > mb));
}

int main(int argc, char **argv)
{
	int i, j, k;
	int sz, rank;
	int nprow, npcol, myrow, mycol, mb, nb;
	int dims[2];
	int gli, glj;
	double *A, *B, *C, *D, *E;
	double time0;
	MPI_Status *st;
	MPI_Request *req;

	MPI_Init(&argc, &argv);
	time0 = MPI_Wtime();
	MPI_Comm_size(MPI_COMM_WORLD, &sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	dims[0] = dims[1] = 0;
	MPI_Dims_create(sz, 2, dims);

	nprow = dims[0];
	npcol = dims[1];
	mb = N / nprow;
	nb = N / npcol;
	myrow = rank / nprow;
	mycol = rank % nprow;

	A = (double *) malloc(mb * nb * sizeof(double));
	B = (double *) malloc(mb * nb * sizeof(double));
	C = (double *) malloc(nb * nb * sizeof(double));

	for (i = 0; i < nb; i++)
		for (j = 0; j < mb; j++)
		{
			gli = myrow * nb + i;
			glj = mycol * mb + j;
			A[i * mb + j] =  ((gli == glj)? (double)(N - 2) / N : (double) -2 / N);
			gli = mycol * mb + j;
			glj = myrow * nb + i;
			B[j * nb + i] = ((gli == glj)? (double) (N - 2) / N : (double) -2 / N);
		}

	if (!Override(mycol, nb, mb))
	{
		D = (double *) malloc(N * nb * sizeof(double));
		E = (double *) malloc(N * nb * sizeof(double));
		// Инициализируем принадлежащие нашему процессу фрагменты матриц
		for (i = 0; i < mb; i++)
			for (j = 0; j < nb; j++)
			{ 
				if (((mycol + 1) * mb - 1) / nb == myrow) E[mycol * nb * mb + i * nb + j] = B[i * nb + j];
				D[mycol * nb * mb + j * mb + i] = A[j * mb + i];
			}
		for (i = 0; i < nb; i++)
			for (j = 0; j < nb; j++) C[i * nb + j] = 0;
		st = (MPI_Status *) malloc(sizeof(MPI_Status) * nprow * 2);
		req = (MPI_Request *) malloc(sizeof(MPI_Request) * nprow * 2);
		// Асинхронно собираем оставшиеся необходимые для умножения фрагменты от остальных процессов
		for (i = 0; i < nprow; i++)
		{
			if (i != mycol) MPI_Irecv((double *) D + i * nb * mb, nb * mb, MPI_DOUBLE, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, &req[(i > mycol)? i-1 : i]);
			if (((i == mycol) && (((mycol + 1) * mb - 1) / nb != myrow)) || (i != mycol)) MPI_Irecv((double *) E + i * nb * mb, nb * mb, MPI_DOUBLE, MPI_ANY_SOURCE, nprow + i, MPI_COMM_WORLD, &req[(i > mycol)? nprow + i - 2 : nprow + i - 1]);
		}
	}

	// Высылаем данные другим процессам
	for (i = 0; i < nprow; i++)
	{
		if ((i != mycol) && !Override(i, nb, mb)) MPI_Send(A, mb * nb, MPI_DOUBLE, nprow * myrow + i, mycol, MPI_COMM_WORLD);
	}
	for (i = 0; i < npcol; i++)
	{
		gli = (myrow * nb) / mb;
		if (i * nprow + gli != rank) MPI_Send(B, nb * mb, MPI_DOUBLE, i * nprow + gli, nprow + mycol, MPI_COMM_WORLD);
	}

	if (!Override(mycol , nb, mb))
	{
		// Дожидаемся, пока все процессы соберут данные
		MPI_Waitall(2 * (nprow - 1), req, st);
		glj = (mycol + 1) * mb / nb - (mb == nb);
		glj *= nb;
		gli = myrow * nb;
		// и собственно перемножаем фрагменты, принадлежащие нашему процессу
		for (j = 0; j < nb; j++)
			for (i = 0; i < nb; i++)
				for (k = 0; k < N; k++)
				{
					C[i * nb + j] += D[(k / mb) * nb * mb + i * mb + (k % mb)] * E[k * nb + j];
				}
	}
 
	MPI_Barrier(MPI_COMM_WORLD);
	time0 = MPI_Wtime() - time0;
	if (rank == 0) printf("Problem size: %d. Number of processes: %d. Time: %lf sec\n", N, sz, time0);
	MPI_Finalize();

	return 0;
}
