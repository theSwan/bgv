#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include "flint/fmpz_vec.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz_poly_mat.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "libbgv.h"

double bgv_get_dvn()
{
	return dvn;
}

void gen_q(fmpz_t q, long len)
{
        mpz_t tmp, hold;
        mpz_init(tmp);
        mpz_init(hold);
        mpz_set_ui(hold, 1);
        mpz_mul_2exp(tmp, hold, len);
        mpz_nextprime(hold, tmp);
        fmpz_set_mpz(q, hold);
}

void hcrypt_random(fmpz_t r)
{
	mpz_t tmp;
	FILE *fp;
	mpz_init(tmp);
	fp = fopen("/dev/urandom", "rb");
        int len = 9;
	if (fp) {
		unsigned char *bytes;
		bytes = (unsigned char *) malloc (len);
              
                if (fread(bytes, 1, len, fp)) {
                        mpz_import(tmp, len, 1, 1, 0, 0, bytes);
                }
                
                else {
                        printf("file read error\n");
                }
                
		fclose(fp);
		free(bytes);
	}
        
        else {
                printf("random number generation error\n");
        }
	

	fmpz_set_mpz(r, tmp);
	mpz_clear(tmp);
}

fmpz *samplez(fmpz *vec, long d)
{
	if ( d == 0 )
		return;
	double tdvn = bgv_get_dvn();
	long a = (long)ceil(-10*tdvn);
	long b = (long)floor(+10*tdvn);
	long x, i;
	double p;
	fmpz_t randseed;
	fmpz_init(randseed);
	hcrypt_random(randseed);
	unsigned long int useed = fmpz_get_ui(randseed);
	srand(useed);
	for( i = 0 ; i < d ; i++) {
		do {
			x = rand()%(b - a) + a;
			p = exp(-pi*x / ( tdvn * tdvn));
		} while ( !( p > 0 && p <= 1) );
                
		vec[i] = x;
	}
	fmpz_clear(randseed);
	return vec;
}

void guassian_poly(fmpz *c, fmpz_poly_t poly, long d)
{
	fmpz *tmp = samplez(c, d);
	long k;
	for( k = 0 ; k < d ; k++ ) {
		fmpz_poly_set_coeff_si(poly, k, tmp[k]);
	}
}

void unif_poly(fmpz_poly_t poly, fmpz_t space, long d)
{
	int i;
	fmpz_t randseed;
	fmpz_init(randseed);
	hcrypt_random(randseed);
	unsigned long int useed = fmpz_get_ui(randseed);
	mpz_t rndnum, rndbd;
	fmpz_t rndfmpz;
	gmp_randstate_t gmpstate;
        
	mpz_init(rndnum);
	mpz_init(rndbd);
	fmpz_get_mpz(rndbd, space);
	fmpz_init(rndfmpz);
	gmp_randinit_default(gmpstate);
	gmp_randseed_ui(gmpstate, useed);
        
	for( i = 0 ; i < d ; i++ ) {
		mpz_urandomm(rndnum, gmpstate, rndbd);
		fmpz_set_mpz(rndfmpz, rndnum);
		fmpz_poly_set_coeff_fmpz(poly, i, rndfmpz);
	}
	fmpz_clear(randseed);
	fmpz_clear(rndfmpz);
	gmp_randclear(gmpstate);
	mpz_clear(rndnum);
	mpz_clear(rndbd);
}

