#include <iostream>
#include "mpi.h"
#include "lab1.h"

#define A 100
#define N 4

struct complex_number {
    int real;
    int imag;

    complex_number(int a, int b) {
        this->real = a;
        this->imag = b;
    }

    complex_number() {
       
    }
};

void fill_matrix(complex_number m[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            m[i][j].real = rand() % 21 - 10;
            m[i][j].imag = rand() % 21 - 10;
        }
    }
}

void print_matrix(complex_number m[N][N]) {
    std::cout << "-----Print matrix-----" << std::endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << m[i][j].real << "\+\i" << m[i][j].imag << "   ";
        }
        std::cout << std::endl;
    }
    std::cout << "-----Matrix printed-----" << std::endl;
}

complex_number complex_mul(complex_number a, complex_number b) {
    complex_number result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}


void count(complex_number matrix_array[A][N][N]) {
    complex_number matrix_result[N][N];
    complex_number matrix_tmp[N][N];

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix_result[i][j] = complex_number(0, 0);
            for (int k = 0; k < N; k++) {
                complex_number tmp = complex_mul(matrix_array[0][i][k], matrix_array[1][k][j]);
                matrix_result[i][j].real += tmp.real;
                matrix_result[i][j].imag += tmp.imag;
            }
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix_tmp[i][j] = matrix_result[i][j];
        }
    }

    for (int matrix = 2; matrix < A; matrix++) {

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                matrix_result[i][j] = complex_number(0, 0);
                for (int k = 0; k < N; k++) {
                    complex_number tmp = complex_mul(matrix_tmp[i][k], matrix_array[matrix][k][j]);
                    matrix_result[i][j].real += tmp.real;
                    matrix_result[i][j].imag += tmp.imag;
                }
            }
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                matrix_tmp[i][j] = matrix_result[i][j];
            }
        }
    }
    //std::cout << "Single proccess result: " << std::endl;
    //print_matrix(matrix_result);
}



int main(int argc, char *argv[])
{
    complex_number matrix_array[A][N][N];
    int ProcNum, ProcRank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    MPI_Datatype COMPLEX_NUMBER;
    MPI_Type_contiguous(2, MPI_INT, &COMPLEX_NUMBER);
    MPI_Type_commit(&COMPLEX_NUMBER);

    MPI_Datatype COMPLEX_ARRAY;
    MPI_Type_contiguous(N, COMPLEX_NUMBER, &COMPLEX_ARRAY);
    MPI_Type_commit(&COMPLEX_ARRAY);

    int from = ProcRank * N / ProcNum;
    int to = (ProcRank+1) * N / ProcNum;
    double time_from;
    if (ProcRank == 0) {
        for (int i = 0; i < A; i++) {
            fill_matrix(matrix_array[i]);
            //print_matrix(matrix_array[i]);
        }
        time_from = MPI_Wtime();
        count(matrix_array);
        double time = MPI_Wtime() - time_from;
        std::cout << "Single proccess runtime: " << time << std::endl;
    }
    /*for (int matrix = 0; matrix < A - 1; matrix++) {
        MPI_Bcast(matrix_array[i])
    }*/
    time_from = MPI_Wtime();
    

    complex_number matrix_A[N][N];
    complex_number matrix_C[N][N];

    MPI_Bcast(matrix_array[1], N * N, COMPLEX_NUMBER, 0, MPI_COMM_WORLD);
    MPI_Scatter(matrix_array[0], N * N / ProcNum, COMPLEX_NUMBER, matrix_A[from], N * N / ProcNum, COMPLEX_NUMBER, 0, MPI_COMM_WORLD);



    complex_number matrix_result[N][N];
    complex_number matrix_tmp[N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix_result[i][j] = complex_number(0, 0);
            for (int k = 0; k < N; k++) {
                complex_number tmp = complex_mul(matrix_A[i][k], matrix_array[1][k][j]);
                matrix_result[i][j].real += tmp.real;
                matrix_result[i][j].imag += tmp.imag;
            }
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix_tmp[i][j] = matrix_result[i][j];
        }
    }

    for (int matrix = 2; matrix < A; matrix++) {
        MPI_Bcast(matrix_array[matrix], N * N, COMPLEX_NUMBER, 0, MPI_COMM_WORLD);
        //MPI_Scatter(matrix_result, N * N / ProcNum, COMPLEX_NUMBER, matrix_A[from], N * N / ProcNum, COMPLEX_NUMBER, 0, MPI_COMM_WORLD);

        for (int matrix = 2; matrix < A; matrix++) {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    matrix_result[i][j] = complex_number(0, 0);
                    for (int k = 0; k < N; k++) {
                        complex_number tmp = complex_mul(matrix_tmp[i][k], matrix_array[matrix][k][j]);
                        matrix_result[i][j].real += tmp.real;
                        matrix_result[i][j].imag += tmp.imag;
                    }
                }
            }

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    matrix_tmp[i][j] = matrix_result[i][j];
                }
            }
        }

       /* for (int i = from; i < to; i++) {
            for (int j = 0; j < N; j++) {
                matrix_C[i][j] = complex_number(0, 0);
                for (int k = 0; k < N; k++) {
                    complex_number tmp = complex_mul(matrix_A[i][k], matrix_array[1][k][j]);
                    matrix_C[i][j].real += tmp.real;
                    matrix_C[i][j].imag += tmp.imag;
                }
            }
        }*/
    }
    MPI_Gather(matrix_result[from], N * N / ProcNum, COMPLEX_NUMBER, matrix_C[from], N * N / ProcNum, COMPLEX_NUMBER, 0, MPI_COMM_WORLD);
    

    if (ProcRank == 0) {
        // matrix_array[2][0] = matrix_C[0];
        //std::cout << "MPI result: " << std::endl;
        //print_matrix(matrix_C);
        double time_to = MPI_Wtime();
        double runtime = time_to - time_from;
        std::cout << "MPI runtime: " << runtime << std::endl;
    }

    
    

    MPI_Type_free(&COMPLEX_ARRAY);
    MPI_Type_free(&COMPLEX_NUMBER);
    MPI_Finalize();
   

    return 0;
}
