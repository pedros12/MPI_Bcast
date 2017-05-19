#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void Get_n(int *n, int my_rank);
void Read_Matrix(int* matriz, int my_rank,int n, int * block_lenght, int * displacement);
void clear_matrix(int *matriz, int n);

int main(){
    int comm_sz, my_rank;
    MPI_Comm comm;
    int* matriz;
    int * block_lenght, *displacement;
    int n;

    //matriz = NULL;
    MPI_Datatype triangulo;

    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    Get_n(&n,my_rank);

        matriz = malloc(n*n*sizeof(int));

        block_lenght = malloc(n * sizeof(int));
        displacement = malloc(n * sizeof(int));

    clear_matrix(matriz,n);
    Read_Matrix(matriz , my_rank , n , block_lenght , displacement);
    MPI_Type_indexed(n, block_lenght, displacement, MPI_INT, &triangulo);
    MPI_Type_commit(&triangulo);

    if(my_rank == 0)
        MPI_Send(matriz,1,triangulo, 1, 0, MPI_COMM_WORLD);

    if(my_rank == 1){
      int i,j;

      MPI_Recv(matriz, 1 , triangulo, 0 , 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for(i = 0 ; i < n ; i++){
          for(j = 0 ; j < n; j++){
            printf("%d ",matriz[i*n + j] );
          }
          printf("\n");
      }
    }

    free(matriz);
    free(block_lenght);
    free(displacement);
    MPI_Type_free(&triangulo);
    MPI_Finalize();
}

void Get_n(int *n, int my_rank){
    if(my_rank == 0){
        printf("Enter the Matrix Order\n" );
        scanf("%d",n );
    }
    MPI_Bcast(n , 1, MPI_INT, 0 , MPI_COMM_WORLD);
}

void Read_Matrix(int* matriz, int my_rank,int n, int * block_lenght, int * displacement){
    int i;

    if(my_rank == 0){
        for(i = 0; i < n*n; i++)
            scanf("%d", &matriz[i] );
        // 4 4 displacement = 0, 5, 12,
        for(i = 0 ; i < n ; i++){
            block_lenght[i] = n - i;
            displacement[i] = i * (n + 1) ;
        }
        int j;
      }

    MPI_Bcast(displacement, n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(block_lenght, n, MPI_INT, 0, MPI_COMM_WORLD);
}

void clear_matrix(int *matriz, int n){
    int i;
    for(i = 0; i < n*n ; i ++)
        matriz[i] = 0;
}
