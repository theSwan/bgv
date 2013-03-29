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

int main(int argc, char *args[]) /*tt, lev, d, {q, n, bign} */
{
        fmpz_t t, tmp;
        fmpz_init(t);
        fmpz_init(tmp);
        fmpz_set_str(t, args[1], 10);
        fmpz_set_str(tmp, args[2], 10);
        
        long lev, d, i, j, row;
        lev = fmpz_get_si(tmp);
        
        fmpz_set_str(tmp, args[3], 10);
        d = fmpz_get_si(tmp);
        
        for( i = 0 ; i <= lev; i++ ) {
                param_node_t *h;
                fmpz_poly_mat_t sk;
                h = param_node_init(h);
                fmpz_set_str(h->q, args[4 + i * 3], 10);
                fmpz_set_str(tmp, args[5 + i * 3], 10);
                h->n = fmpz_get_si(tmp);
                fmpz_set_str(tmp, args[6 + i * 3], 10);
                h->bign = fmpz_get_si(tmp);
                row = 1 + h->n;
                fmpz_poly_mat_init(sk, row, 1);
                e_skeygen(sk, h, d);
                printf("%ld\n1\n", row);
                for( j = 0; j < row; j++ ) {
                        char *ss;
                        ss = fmpz_poly_get_str(fmpz_poly_mat_entry(sk, j, 0));
                        printf("%s\n", ss);
                }
                fmpz_poly_mat_clear(sk);
        }
        fmpz_clear(t);
        fmpz_clear(tmp);
        return 0;
}