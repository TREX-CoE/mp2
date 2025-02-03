/* Includes */


#include <stdio.h>
#include <trexio.h>
#include <stdint.h>
#include <err.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

/* Reordering function */

/*   There are permutation symmetries in the indices of integrals, which */
/*   are exploited to reduce the storage of the integrals in the TREXIO file. */

/*   \[ */
/*   \langle i j | k l \rangle =  */
/*   \langle k j | i l \rangle =  */
/*   \langle k l | i j \rangle =  */
/*   \langle i l | k j \rangle =  */
/*   \langle j i | l k \rangle =  */
/*   \langle j k | l i \rangle =  */
/*   \langle l k | j i \rangle =  */
/*   \langle l i | j k \rangle  */
/*   \] */

/*   The following function swaps the four indices such that any of the 8 */
/*   possible combinations gives the same quartet. */


void reorder(int *i_, int *j_, int *k_, int *l_) {
  {
    assert (*i_ >= 0); assert (*j_ >= 0);
    assert (*k_ >= 0); assert (*l_ >= 0);

    const int i = *i_; const int j = *j_;
    const int k = *k_; const int l = *l_;

    if (k<i) {
      *k_ = i ; *i_ = k;
    }
    if (l<j) {
      *l_ = j ; *j_ = l;
    }
  }
  {
    const int i = *i_; const int j = *j_;
    const int k = *k_; const int l = *l_;

    if (j<i) {
      *i_ = j ; *j_ = i;
      *k_ = l ; *l_ = k;
    }
  }
}

/* MP2 program */


int main(int argc, char** argv)
{

/* Open the TREXIO file */

/*    The name of the TREXIO file should be given as a command-line argument. */


  if (argc < 2) {
    fprintf(stderr, "usage: mp2 trexio_file.hdf5\n");
    exit(1);
  }

  trexio_exit_code rc = TREXIO_SUCCESS;
  trexio_t* trexio_file = trexio_open(argv[1], 'r', TREXIO_HDF5, &rc);

  if (rc != TREXIO_SUCCESS) {
    fprintf(stderr, "Error opening file %s", argv[1]);
    exit(1);
  }
  assert (trexio_file != NULL);

/* Read parameters from TREXIO */

/*    We need to read small scalar variables. The first ones are the */
/*    number of up-spin and down-spin electrons to define the number of */
/*    occupied orbitals ~n_occ~, and check that we are in a closed-shell system. */


  int n_up = 0;
  rc = trexio_read_electron_up_num(trexio_file, &n_up);
  if (rc != TREXIO_SUCCESS) {
    fprintf(stderr, "Error reading n_up");
    exit(1);
  }
  assert(n_up > 0);

  int n_dn = 0;
  rc = trexio_read_electron_dn_num(trexio_file, &n_dn);
  if (rc != TREXIO_SUCCESS) {
    fprintf(stderr, "Error reading n_dn");
    exit(1);
  }
  assert(n_up > 0);

  if (n_up != n_dn) {
    fprintf(stderr, "This code is works only for n_up = n_dn");
    exit(1);
  }

  const int n_occ = n_up;



/* We also need to read the total number of molecular orbitals to */
/* compute the number of virtual orbitals ~n_virt~. */


  int mo_num = 0;
  rc = trexio_read_mo_num(trexio_file, &mo_num);
  if (rc != TREXIO_SUCCESS) {
    fprintf(stderr, "Error reading mo_num");
    exit(1);
  }
  assert(mo_num > 0);

  const int n_virt = mo_num - n_up;

/* Quantities involved in the MP2 equation */
   
/*    We assume the TREXIO file contains electron repulsion integrals */
/*    (ERI) in the molecular orbital basis, and orbital energies. */

/*    We first read the orbital energies: */
   

  double*  epsilon = malloc(mo_num * sizeof(double));
  rc = trexio_read_mo_energy(trexio_file, epsilon);
  if (rc != TREXIO_SUCCESS) {
    fprintf(stderr, "Error reading mo_energy");
    exit(1);
  }



/* Then we read the electron repulsion integrals. They are stored in a */
/* sparse data format, so we obtain quartets of indices and values for */
/* non-zero integrals. */


  int64_t  n_integrals;
  rc = trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);
  if (rc != TREXIO_SUCCESS) {
    fprintf(stderr, "Error reading n_integrals");
    exit(1);
  }
  assert(n_integrals > 0);

  int* const index = malloc(4*n_integrals * sizeof(int));
  if (index == NULL) {
    fprintf(stderr, "Malloc failed for index");
    exit(1);
  }

  double* const value = malloc(n_integrals * sizeof(double));
  if (index == NULL) {
    fprintf(stderr, "Malloc failed for value");
    exit(1);
  }

  int64_t count = n_integrals;
  rc = trexio_read_mo_2e_int_eri(trexio_file, 0L, &count, index, value);



/* We transform these arrays of indices and values into an array of */
/* double, where the index of the quartet ~(i,j,a,b)~ is located at */
/* address ~b-n_occ + n_virt*(a-n_occ + n_virt*(j + n_occ*i))~. Only */
/* required integrals are stored. */

/* By symmetry, the integral at ~(i,j,a,b)~ is equal to the integral at */
/* ~(j,i,b,a)~.  As integrals may be stored only once in the integrals */
/* file, to be sure we don't miss integrals we store the value at both */
/* addresses. */


  size_t nmax = n_occ*n_occ*n_virt*n_virt;
  double* integral = malloc(nmax * sizeof(double));
  memset(integral, 0, nmax*sizeof(double));

  for (size_t kk=0; kk<n_integrals ; ++kk) {
    int i = index[4*kk+0];  assert (i >= 0);
    int j = index[4*kk+1];  assert (j >= 0);
    int a = index[4*kk+2];  assert (a >= 0);
    int b = index[4*kk+3];  assert (b >= 0);

    reorder(&i, &j, &a, &b);
    if (i >= n_occ || j >= n_occ || a < n_occ || b < n_occ ) {
      continue;
    } else {

      a -= n_occ;
      b -= n_occ;
      const size_t ijab = b + n_virt*(a + n_virt*(j + n_occ*i));
      const size_t jiba = a + n_virt*(b + n_virt*(i + n_occ*j));
      integral[ijab] = value[kk];
      integral[jiba] = value[kk];
    }
  }

/* MP2 computation */


  double Emp2 = 0.;

  for   (int i=0 ; i<n_occ ; ++i) {
    for (int j=0 ; j<n_occ ; ++j) {

      const size_t shift = n_virt*(j + n_occ*i);

      for   (int a=0 ; a<n_virt ; ++a) {
        for (int b=0 ; b<n_virt ; ++b) {

          const size_t ijab = b + n_virt*(a + shift);
          const size_t ijba = a + n_virt*(b + shift);

          Emp2 += ( integral[ijab]*(2.*integral[ijab]-integral[ijba]) ) /
            (epsilon[i] + epsilon[j] - epsilon[n_occ+a] - epsilon[n_occ+b]);

        }
      }

    }
  }

/* Termination */
  
/*   Print the result: */


  printf("Emp2 = %15.12f\n", Emp2);
}
