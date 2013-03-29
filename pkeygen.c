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

void print(pk_node_t *node);

int main(int argc, char *args[])/*tt, lev, d,  {q,n, bign},filename{row, col, poly}*/
{
	fmpz_t t, tmp;
        fmpz_init(t);
        fmpz_init(tmp);
        fmpz_set_str(t, args[1], 10);
        fmpz_set_str(tmp, args[2], 10);
        
        long lev, d, i, j;
        lev = fmpz_get_si(tmp);
        
        fmpz_set_str(tmp, args[3], 10);
        d = fmpz_get_si(tmp);
        fmpz_poly_t fx;
        fmpz_poly_init(fx);
        fmpz_poly_set_coeff_si(fx, 0, 1);
        fmpz_poly_set_coeff_si(fx, d, 1);
        param_node_t *ph, *pr, *ps, *param, *pam;
        ph = param_node_init(ph);

        ps = ph;
        for( i = 0 ; i <= lev; i++ ) {
                pr = param_node_init(pr);
                fmpz_set_str(pr->q, args[4 + i * 3], 10);
                fmpz_set_str(tmp, args[5 + i * 3], 10);
                pr->n = fmpz_get_si(tmp);
                fmpz_set_str(tmp, args[6 + i * 3], 10);
                pr->bign = fmpz_get_si(tmp);
                ps->next = pr;
                ps = pr;
        }
        ps->next = NULL;
        
        long ind = 4 + 3 * (lev + 1), row, col;
        FILE *fp;

        if((fp = fopen(args[ind], "r")) == NULL)
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
        pk_node_t *pk;
        
        param = ph->next;
        pam = param->next;
        
        sh = sh->next;
        pk = (pk_node_t *)malloc(sizeof(pk_node_t));
        fmpz_poly_mat_init(pk->pkb, 1, 1);
        fmpz_poly_mat_zero(pk->pkb);
        fmpz_poly_mat_init(pk->pka, param->bign, 1 + (param->n));
        e_pkeygen(pk->pka, param, sh->sk, d, t, fx);
        print(pk);
        free(pk);
        ss = sh;
        sr = ss->next;
        
        
        fmpz_poly_mat_t s1, s2, tensor;
	long row1, row2, len, llog;
        
        for( i = lev ; i > 0 ; i-- ) {
                llog = fmpz_clog(pam->q, t);
                pk = (pk_node_t *)malloc(sizeof(pk_node_t));
		fmpz_poly_mat_init(pk->pka, pam->bign, 1 + (pam->n));
                e_pkeygen(pk->pka, pam, sr->sk, d, t, fx);
                
                row1 = fmpz_poly_mat_nrows(ss->sk);
        	row2 = row1 * row1;
        	fmpz_poly_mat_init(tensor, row2, 1);
                vec_tensor(tensor, ss->sk, param->q, fx);
                
                len = fmpz_clog(param->q, t);
		row2 = row2 * len;
                fmpz_poly_mat_init(s1, row2, 1);
		fmpz_poly_mat_init(s2, row2, 1);
                bitdecomp(s1, tensor, param->q, t, d);
                scale(s2, s1, param->q, pam->q, t);
                
                row1 = fmpz_poly_mat_nrows(s2) * llog;
		row2 = fmpz_poly_mat_nrows(sr->sk);
		fmpz_poly_mat_init(pk->pkb, row1, row2);
                
                switchkeygen(pk->pkb, s2, sr->sk, pam->q, d, t, fx);
                
                fmpz_poly_mat_clear(s1);
                fmpz_poly_mat_clear(s2);
                fmpz_poly_mat_clear(tensor);
                pam = pam->next;
                param = param->next;
                ss = sr;
                sr = sr->next;
                print(pk);
                free(pk);
        }
        
        return 0;
}

void print(pk_node_t *node)
{
        long row, col, i, j;
        row = fmpz_poly_mat_nrows(node->pka);
        col = fmpz_poly_mat_ncols(node->pka);
        printf("%ld\n%ld\n", row, col);
        char *str;
        for( i = 0 ; i < row ; i++ ) {
                for( j = 0 ; j < col ; j++ ) {
                        str = fmpz_poly_get_str(fmpz_poly_mat_entry(node->pka, i, j));
                        printf("%s\n", str);
                }
        }
        
        row = fmpz_poly_mat_nrows(node->pkb);
        col = fmpz_poly_mat_ncols(node->pkb);
        printf("%ld\n%ld\n", row, col);

        for( i = 0 ; i < row ; i++ ) {
                for( j = 0 ; j < col ; j++ ) {
                        str = fmpz_poly_get_str(fmpz_poly_mat_entry(node->pkb, i, j));
                        printf("%s\n", str);
                }
        }
}