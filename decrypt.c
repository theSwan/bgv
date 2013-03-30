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

int main(int argc, char *args[])/* setup.txt, ct.txt{ct->lev, ct->row, ct->col, ct}, sk.txt{row, col, poly}*/
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
        
        param_node_t *ph, *pr, *ps, *param;
        ph = param_node_init(ph);
        
        ps = ph;
        for( i = 0 ; i <= lev; i++ ) {
                pr = param_node_init(pr);
                fgets(str, 100, fp);
                fmpz_set_str(pr->q, str, 10);
                fgets(str, 100, fp);
                fmpz_set_str(tmp, str, 10);
                pr->n = fmpz_get_si(tmp);
                fgets(str, 100, fp);
                fmpz_set_str(tmp, str, 10);
                pr->bign = fmpz_get_si(tmp);
                ps->next = pr;
                ps = pr;
        }
        ps->next = NULL;
        fclose(fp);
        
        param = ph->next;

        long row, col, ctlev;
        
        if((fp = fopen(args[2], "r")) == NULL)
        {
                printf("file read error\n");
                exit(0);
        }
        
        fgets(str, 30, fp);
        fmpz_set_str(tmp, str,10);
        ctlev = fmpz_get_si(tmp);
        
        fgets(str, 30, fp);
        fmpz_set_str(tmp, str,10);
        row = fmpz_get_si(tmp);
        
        fgets(str, 30, fp);
        fmpz_set_str(tmp, str,10);
        col = fmpz_get_si(tmp);
        
        fmpz_poly_mat_t ct;
        fmpz_poly_mat_init(ct, row, col);
        for( i = 0 ; i < row ; i++) {
                for(j = 0; j < col ; j++) {
                        fgets(str, 100000, fp);
                        fmpz_poly_set_str(fmpz_poly_mat_entry(ct, i, j), str);
                }
        }
        
        fclose(fp);

        
        if((fp = fopen(args[3], "r")) == NULL)
        {
                printf("file read error\n");
                exit(0);
        }

        long l;
        sk_node_t *sh, *ss, *sr;
        sh = (sk_node_t *)malloc(sizeof(sk_node_t));

        ss = sh;
        for(l = 0; l <= lev ; l++){
                sr = (sk_node_t *)malloc(sizeof(sk_node_t));
                fgets(str, 30, fp);
                fmpz_set_str(tmp, str,10);
                row = fmpz_get_si(tmp);


                fgets(str, 30, fp);
                fmpz_set_str(tmp, str,10);
                col = fmpz_get_si(tmp);

                fmpz_poly_mat_init(sr->sk, row, col);

                for( i = 0 ; i < row ; i++) {
                        for(j = 0; j < col ; j++) {
                                fgets(str, 100000, fp);
                                fmpz_poly_set_str(fmpz_poly_mat_entry(sr->sk, i, j), str);
                        }
                }
                ss->next = sr;
                ss = sr;
        }
        ss->next = NULL;
        
        fclose(fp);

        sh = sh->next;
        while(lev > ctlev){
                sh = sh->next;
                param = param->next;
                lev--;
        }
        
        fmpz_poly_t ms;
        fmpz_poly_init(ms);
        e_decrypt(ms, param, sh->sk, ct, t, fx);
        fmpz_poly_print(ms);
        printf("\n");
        return 0;
        
}