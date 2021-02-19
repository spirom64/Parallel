#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "mkl_pblas.h"

#define N 8192

int main(int argc, char **argv)
{
	int i, j;
	int size = N;
	double al = 1;
	double bt = 0;
	int cr = 0; int one = 1;
	int desca[9], descb[9], descc[9];
	int icntxt;
	int sz, rank;
	int nprow, npcol, myrow, mycol, ma, na, info, mb, nb;
	int dims[2];
	int gli, glj;
	double *A, *B, *C;
	double time0;
	char *endptr = NULL;
	int clmn = 0;
	char top[16];
	int lld;

	MPI_Init(&argc, &argv);
	time0 = MPI_Wtime();
	MPI_Comm_size(MPI_COMM_WORLD, &sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	dims[0] = dims[1] = 0;
	MPI_Dims_create(sz, 2, dims);
	nb = 0;
	i = 1;
	strcpy(top, "Row-column");
	if (argc > 1) {
		if (!(argc % 2)) {
			if (!rank) printf("Invalid command line: both option and value expected\n");
			MPI_Finalize();
			return 1;
		}
		while (argc > 1) {
			if (!strcmp(argv[i], "-b")) {
				nb = strtol(argv[i + 1], &endptr, 10);
				i += 2;
				argc -= 2;
				if (!nb) {
					if (!rank) printf("Invalid block size!\n");
					MPI_Finalize();
					return 1;
				}
			} else if (!strcmp(argv[i], "-t")) {
				if (!strcmp(argv[i + 1], "r")) {
					dims[0] = sz;
					dims[1] = 1;
					strcpy(top, "Row");
				} else if (!strcmp(argv[i + 1], "c")) {
					dims[0] = 1;
					dims[1] = sz;
					strcpy(top, "Column");
					clmn = 1
				} else if (strcmp(argv[i + 1], "b") && !rank)
					printf("Unknown topology. Defaulting to row-column\n");
				i += 2;
				argc -= 2;
			} else {
				if (!rank) printf("Unknown option %s = %s. Skipping\n", argv[i] + 1, argv[i + 1]);
				argc -= 2;
				i += 2;
			}
		}
	}

	Cblacs_pinfo(&rank, &sz);
	Cblacs_get(-1, 0, &icntxt);
	Cblacs_gridinit(&icntxt, (clmn)? "C" : "R", dims[0], dims[1]);
	nprow = dims[0];
	npcol = dims[1];
	Cblacs_gridinfo(icntxt, &nprow, &npcol, &myrow, &mycol);
	mb = N / nprow;
	if (!nb) nb = N / sz;
	ma = numroc_(&size, &nb, &myrow, &cr, &nprow);
	na = numroc_(&size, &nb, &mycol, &cr, &npcol);
	if (clmn) lld = na; else lld = ma;
	if (!lld) lld = 1;

	info = 0;
	descinit_(desca, &size, &size, &nb, &nb, &cr, &cr, &icntxt, &lld, &info);
	descinit_(descb, &size, &size, &nb, &nb, &cr, &cr, &icntxt, &lld, &info);
	descinit_(descc, &size, &size, &nb, &nb, &cr, &cr, &icntxt, &lld, &info);
	A = (double *) malloc(ma * na * sizeof(double));
	B = (double *) malloc(ma * na * sizeof(double));
	C = (double *) malloc(ma * na * sizeof(double));

	for (i = 0; i < ma; i++)
		for (j = 0; j < na; j++)
		{
			gli = ((i / nb) + myrow) * nb + (i % nb);
			glj = ((j / nb) + mycol) * nb + (j % nb);
			A[i * na + j] =  ((gli == glj)? (double)(N - 2) / N : (double) -2 / N);
			B[i * na + j] = A[i * na + j];
		}
 
	pdgemm_("N", "N", &size, &size, &size, &al, A, &one, &one, desca, B, &one, &one, descb, &bt, C, &one, &one, descc);

	Cblacs_barrier(icntxt, "A");
	Cblacs_gridexit(0);

	time0 = MPI_Wtime() - time0;
	if (rank == 0) printf("Problem size: %d. Block size: %d. Number of processes: %d. Topology: %s. Time: %lf sec\n", N, nb, sz, top, time0);
	MPI_Finalize();

	return 0;
}
