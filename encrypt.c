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

char str[100000];

int main(int argc, char *args[])/* setup.txt, ms, pk.txt{row, col, pka, row, col, pkb} */
{
        FILE *fp;
        
        if((fp = fopen(args[1], "r")) == NULL)
        {
                printf("file read error\n");
                exit(0);
        }
        
        fmpz_t t, tmp;
        fmpz_init(t);
        fmpz_init(tmp);
        fgets(str, 100, fp);
        fmpz_set_str(t, str, 10);
        fgets(str, 100, fp);
        fmpz_set_str(tmp, str, 10);
        
        long lev, d, i, j;
        lev = fmpz_get_si(tmp);
        fgets(str, 100, fp);
        fmpz_set_str(tmp, str, 10);
        d = fmpz_get_si(tmp);
        
        fmpz_poly_t fx;
        fmpz_poly_init(fx);
        fmpz_poly_set_coeff_si(fx, 0, 1);
        fmpz_poly_set_coeff_si(fx, d, 1);
        
        param_node_t *param;
        param = param_node_init(param);
        fgets(str, 100, fp);
        fmpz_set_str(param->q, str, 10);
        fgets(str, 100, fp);
        fmpz_set_str(tmp, str, 10);
        param->n = fmpz_get_si(tmp);
        fgets(str, 100, fp);
        fmpz_set_str(tmp, str, 10);
        param->bign = fmpz_get_si(tmp);
        param->next = NULL;

        fclose(fp);
        
        fmpz_poly_t ms;
        fmpz_poly_init(ms);
        fmpz_poly_set_str(ms, args[2]);
                
        if((fp = fopen(args[3], "r")) == NULL)
        {
                printf("file read error\n");
                exit(0);
        }
        
        long row, col;
        
        fgets(str, 30, fp);
        fmpz_set_str(tmp, str,10);
        row = fmpz_get_si(tmp);
        
        fgets(str, 30, fp);
        fmpz_set_str(tmp, str,10);
        col = fmpz_get_si(tmp);
        
        fmpz_poly_mat_t pk;
        fmpz_poly_mat_init(pk, row, col);
        for( i = 0 ; i < row ; i++) {
                for(j = 0; j < col ; j++) {
                        fgets(str, 100000, fp);
                        fmpz_poly_set_str(fmpz_poly_mat_entry(pk, i, j), str);
                }
        }
        
        fclose(fp);
        
        printf("%ld\n", lev);
        
        fmpz_poly_mat_t ct;
        
        row = 1 + param->n;
        col = 1; 
        fmpz_poly_mat_init(ct, row, 1);
        
        e_encrypt(ct, param, pk, ms, t, fx, d);
        printf("%ld\n%ld\n", row, col);
        char *output;
        for( i = 0; i < row; i++) {
                output = fmpz_poly_get_str(fmpz_poly_mat_entry(ct, i, 0));
                printf("%s\n", output);
        }
        
        return 0;
}
/*output: ct->lev, row, col, ct->text */