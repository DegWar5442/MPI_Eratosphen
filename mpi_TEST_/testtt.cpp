#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

#define MIN(a,b) ((a)<(b)?(a):(b))


inline int BLOCK_LOW(int id, int p, int n) { return (id * n) / p; }
inline int BLOCK_HIGH(int id, int p, int n) { return BLOCK_LOW(id + 1, p, n) - 1; }
inline int BLOCK_SIZE(int id, int p, int n) { return BLOCK_HIGH(id, p, n) - BLOCK_LOW(id, p, n) + 1; }

int main(int argc, char* argv[]) {
    int id, p, n, low_value, high_value, size, proc0_size;
    int count = 0, global_count = 0;
    double elapsed_time = 0.0;

    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (id == 0) {
        std::cout << "input end number (n): ";
        std::cin >> n;
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    low_value = 2 + BLOCK_LOW(id, p, n - 1);
    high_value = 2 + BLOCK_HIGH(id, p, n - 1);
    size = BLOCK_SIZE(id, p, n - 1);
    proc0_size = (n - 1) / p;

    if ((2 + proc0_size) < static_cast<int>(std::sqrt(static_cast<double>(n)))) {
        if (id == 0) std::cerr << "Too many processes\n";
        MPI_Finalize();
        return 1;
    }

    std::vector<char> marked(size, 0);

    int index = 0;
    int prime = 2;

    do {
        int first;
        if (prime * prime > low_value) {
            first = prime * prime - low_value;
        }
        else {
            if (low_value % prime == 0)
                first = 0;
            else
                first = prime - (low_value % prime);
        }

        for (int i = first; i < size; i += prime)
            marked[i] = 1;

        if (id == 0) {
            while (marked[++index]);
            prime = index + 2;
        }

        MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } while (prime * prime <= n);

    for (int i = 0; i < size; ++i)
        if (!marked[i]) count++;

    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    elapsed_time += MPI_Wtime();

    if (id == 0) {
        std::cout << global_count << " primes are less than or equal to " << n << "\n";
        std::cout << "Total elapsed time: " << elapsed_time << " seconds\n";
    }

    MPI_Finalize();
    return 0;
}
