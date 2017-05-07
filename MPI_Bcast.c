/* File:     mpi_trap2.c
 * Purpose:  Use MPI to implement a parallel version of the trapezoidal
 *           rule.  This version accepts input of the endpoints of the
 *           interval and the number of trapezoids.
 *
 * Input:    The endpoints of the interval of integration and the number
 *           of trapezoids
 * Output:   Estimate of the integral from a to b of f(x)
 *           using the trapezoidal rule and n trapezoids.
 *
 * Compile:  mpicc -g -Wall -o mpi_trap2 mpi_trap2.c
 * Run:      mpiexec -n <number of processes> ./mpi_trap2
 *
 * Algorithm:
 *    1.  Each process calculates "its" interval of
 *        integration.
 *    2.  Each process estimates the integral of f(x)
 *        over its interval using the trapezoidal rule.
 *    3a. Each process != 0 sends its integral to 0.
 *    3b. Process 0 sums the calculations received from
 *        the individual processes and prints the result.
 *
 * Note:  f(x) is all hardwired.
 *
 * IPP:   Section 3.3.2  (pp. 100 and ff.)
 */
#include <stdio.h>
#include <math.h>
/* We'll be using MPI routines, definitions, etc. */
#include <mpi.h>

/* Get the input values */
void Get_input(int my_rank, int comm_sz, double* a_p, double* b_p,
      int* n_p);

/* Calculate local integral  */
double Trap(double left_endpt, double right_endpt, int trap_count,
   double base_len);

/* Function we're integrating */
double f(double x);

int main(void) {
   int my_rank, comm_sz, n, local_n;
   double a, b, h, local_a, local_b;
   double local_int, total_int;
   int source;
   double starttime, endtime;

   /* Let the system do what it needs to start up MPI */
   MPI_Init(NULL, NULL);
   // iniciar o tempo de execução
   starttime = MPI_Wtime();

   /* Get my process rank */
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   /* Find out how many processes are being used */
   MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

   Get_input(my_rank, comm_sz, &a, &b, &n);

   h = (b-a)/n;          /* h is the same for all processes */
   local_n = n/comm_sz;  /* So is the number of trapezoids  */

   /* Length of each process' interval of
    * integration = local_n*h.  So my interval
    * starts at: */
   local_a = a + my_rank*local_n*h;
   local_b = local_a + local_n*h;
   local_int = Trap(local_a, local_b, local_n, h);

   /* Add up the integrals calculated by each process */
   if (my_rank != 0)
      MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0,
            MPI_COMM_WORLD);
   else {
      total_int = local_int;
      for (source = 1; source < comm_sz; source++) {
         MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         total_int += local_int;
      }
   }

   /* Print the result */
   if (my_rank == 0) {
      printf("With n = %d trapezoids, our estimate\n", n);
      printf("of the integral from %f to %f = %.15e\n",
          a, b, total_int);
      // exibir tempo total do procedimento
      endtime = MPI_Wtime();
      printf("Tempo total: %f\n",endtime-starttime);
   }

   /* Shut down MPI */
   MPI_Finalize();

   return 0;
} /*  main  */

/*------------------------------------------------------------------
 * Function:     Get_input
 * Purpose:      Get the user input:  the left and right endpoints
 *               and the number of trapezoids
 * Input args:   my_rank:  process rank in MPI_COMM_WORLD
 *               comm_sz:  number of processes in MPI_COMM_WORLD
 * Output args:  a_p:  pointer to left endpoint
 *               b_p:  pointer to right endpoint
 *               n_p:  pointer to number of trapezoids
 */
void Get_input(int my_rank, int comm_sz, double* a_p, double* b_p,
      int* n_p) {
   int dest;
   int i, j;
   int flag = 0;

    if (my_rank == 0) {

      FILE *myFile;
	    myFile = fopen("input.txt", "r");
	    fscanf(myFile, "%lf\n", a_p);
	    fscanf(myFile, "%lf\n", b_p);
	    fscanf(myFile, "%i\n", n_p);
	    fclose(myFile);

      printf("a_p %lf\n",*a_p );
      printf("b_p %lf\n",*b_p );
      printf("n_p %d\n",*n_p );
        // cenário base: comm_sz = 8

        for (i = 1; i <= log2(comm_sz); i++) { // executa por 3 vezes
            // distribuir para os outros processos
            // a partir do processo 0: 4 / 2 / 1
        	dest = my_rank + (comm_sz/pow(2,i));
            MPI_Send(a_p, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
            MPI_Send(b_p, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
            MPI_Send(n_p, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
	      }
    } else { /* my_rank != 0 */
          if(my_rank % 2 == 1){
              MPI_Recv(a_p, 1, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE);
              MPI_Recv(b_p, 1, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE);
              MPI_Recv(n_p, 1, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE);
              //printf("My_rank %d receving from %d = %d\n", my_rank , my_rank - 1, *n_p );;
          }
          else{
              int div = comm_sz/2;
              for(i = log2(comm_sz)-1 ; i > 0  ; i--){
                  if (my_rank % div == 0 && flag != 1) {
                      MPI_Recv(a_p, 1, MPI_DOUBLE, my_rank - div, 0, MPI_COMM_WORLD,
                                MPI_STATUS_IGNORE);
                      MPI_Recv(b_p, 1, MPI_DOUBLE, my_rank - div, 0, MPI_COMM_WORLD,
                                MPI_STATUS_IGNORE);
                      MPI_Recv(n_p, 1, MPI_INT, my_rank - div, 0, MPI_COMM_WORLD,
                                MPI_STATUS_IGNORE);
                      //printf("My_rank %d receving from %d = %d\n",my_rank , my_rank - div, *n_p );
                      //se o rank pertence ao nivel 2 dá árvore ele irá enviar para 2 cores
                      for (j = 0; j < i; j++) {
                          div /= 2;
                          MPI_Send(a_p, 1, MPI_DOUBLE, my_rank + div, 0, MPI_COMM_WORLD);
                          MPI_Send(b_p, 1, MPI_DOUBLE, my_rank + div, 0, MPI_COMM_WORLD);
                          MPI_Send(n_p, 1, MPI_INT, my_rank + div, 0, MPI_COMM_WORLD);

                      }
                      flag = 1; // Core já recebeu e enviou
                  }
                  //evitar que um core com flag 1 faça divisão por 0
                  else if (flag != 1){
                      div /=2;
                  }
              }
          }
        }
        printf("My_rank %d a_n = %lf b_n = %lf n_p = %d\n",my_rank, *a_p, *b_p, *n_p );
        //se o rank pertence ao nivel 2 dá árvore ele irá enviar para 2 co
    } // fim else
  /* Get_input */

/*------------------------------------------------------------------
 * Function:     Trap
 * Purpose:      Serial function for estimating a definite integral
 *               using the trapezoidal rule
 * Input args:   left_endpt
 *               right_endpt
 *               trap_count
 *               base_len
 * Return val:   Trapezoidal rule estimate of integral from
 *               left_endpt to right_endpt using trap_count
 *               trapezoids
 */
double Trap(
      double left_endpt  /* in */,
      double right_endpt /* in */,
      int    trap_count  /* in */,
      double base_len    /* in */) {
   double estimate, x;
   int i;

   estimate = (f(left_endpt) + f(right_endpt))/2.0;
   for (i = 1; i <= trap_count-1; i++) {
      x = left_endpt + i*base_len;
      estimate += f(x);
   }
   estimate = estimate*base_len;

   return estimate;
} /*  Trap  */


/*------------------------------------------------------------------
 * Function:    f
 * Purpose:     Compute value of function to be integrated
 * Input args:  x
 */
double f(double x) {
   return x*x;
} /* f */
