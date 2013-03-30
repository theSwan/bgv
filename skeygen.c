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


int main(int argc, char *args[]) /*setup.txt:tt,lev,d, {q, n, bign} */
{
        FILE *fp;
        
        if((fp = fopen(args[1], "r")) == NULL)
        {
                printf("file read error\n");
                exit(0);
        }
        
        char str[100];
        
        fmpz_t t, tmp;
        fmpz_init(t);
        fmpz_init(tmp);
        fgets(str, 100, fp);
        fmpz_set_str(t, str, 10);
        fgets(str, 100, fp);
        fmpz_set_str(tmp, str, 10);
        
        long lev, d, i, j, row;
        lev = fmpz_get_si(tmp);
        fgets(str, 100, fp);
        fmpz_set_str(tmp, str, 10);
        d = fmpz_get_si(tmp);
        
        
        for( i = 0 ; i <= lev; i++ ) {
                param_node_t *h;
                fmpz_poly_mat_t sk;
                h = param_node_init(h);
                fgets(str, 100, fp);
                fmpz_set_str(h->q, str, 10);
                fgets(str, 100, fp);
                fmpz_set_str(tmp, str, 10);
                h->n = fmpz_get_si(tmp);
                fgets(str, 100, fp);
                fmpz_set_str(tmp, str, 10);
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
        fclose(fp);
        fmpz_clear(t);
        fmpz_clear(tmp);
        return 0;
}