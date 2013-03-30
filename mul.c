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

int main(int argc, char *args[])
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
        
        param_node_t *hp, *pr, *ps;
        hp = param_node_init(hp);
        
        ps = hp;
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
        
        hp = hp->next;
        
        long row, col, ctlev1, ctlev2, ctlev;
        
        if((fp = fopen(args[2], "r")) == NULL)
        {
                printf("file read error\n");
                exit(0);
        }
        
        fgets(str, 30, fp);
        fmpz_set_str(tmp, str,10);
        ctlev1 = fmpz_get_si(tmp);
        
        fgets(str, 30, fp);
        fmpz_set_str(tmp, str,10);
        row = fmpz_get_si(tmp);
        
        fgets(str, 30, fp);
        fmpz_set_str(tmp, str,10);
        col = fmpz_get_si(tmp);
        
        fmpz_poly_mat_t ct1;
        fmpz_poly_mat_init(ct1, row, col);
        for( i = 0 ; i < row ; i++) {
                for(j = 0; j < col ; j++) {
                        fgets(str, 100000, fp);
                        fmpz_poly_set_str(fmpz_poly_mat_entry(ct1, i, j), str);
                }
        }
        
        fclose(fp);
        
        
        if((fp = fopen(args[3], "r")) == NULL)
        {
                printf("file read error\n");
                exit(0);
        }
        
        fgets(str, 30, fp);
        fmpz_set_str(tmp, str,10);
        ctlev2 = fmpz_get_si(tmp);
        
        fgets(str, 30, fp);
        fmpz_set_str(tmp, str,10);
        row = fmpz_get_si(tmp);
        
        fgets(str, 30, fp);
        fmpz_set_str(tmp, str,10);
        col = fmpz_get_si(tmp);
        
        fmpz_poly_mat_t ct2;
        fmpz_poly_mat_init(ct2, row, col);
        for( i = 0 ; i < row ; i++) {
                for(j = 0; j < col ; j++) {
                        fgets(str, 100000, fp);
                        fmpz_poly_set_str(fmpz_poly_mat_entry(ct2, i, j), str);
                }
        }
        
        fclose(fp);
        
        
        if((fp = fopen(args[4], "r")) == NULL)
        {
                printf("file read error\n");
                exit(0);
        }
        
        pk_node_t *pktmp, *pks, *pkr;
        long ii, jj;
        pktmp = (pk_node_t *)malloc(sizeof(pk_node_t));
        
        pks = pktmp;
        for(i=0;i<=lev;i++) {
                pkr = (pk_node_t *)malloc(sizeof(pk_node_t));
                fgets(str, 30, fp);
                fmpz_set_str(tmp, str,10);
                row = fmpz_get_si(tmp);
                
                fgets(str, 30, fp);
                fmpz_set_str(tmp, str,10);
                col = fmpz_get_si(tmp);
                
                fmpz_poly_mat_init(pkr->pka, row, col);
                for( ii = 0 ; ii < row ; ii++) {
                        for(jj = 0; jj < col ; jj++) {
                                fgets(str, 100000, fp);
                                fmpz_poly_set_str(fmpz_poly_mat_entry(pkr->pka, ii, jj), str);
                        }
                }
                
                fgets(str, 30, fp);
                fmpz_set_str(tmp, str,10);
                row = fmpz_get_si(tmp);
                
                fgets(str, 30, fp);
                fmpz_set_str(tmp, str,10);
                col = fmpz_get_si(tmp);
                
                fmpz_poly_mat_init(pkr->pkb, row, col);
                for( ii = 0 ; ii < row ; ii++) {
                        for(jj = 0; jj < col ; jj++) {
                                fgets(str, 100000, fp);
                                fmpz_poly_set_str(fmpz_poly_mat_entry(pkr->pkb, ii, jj), str);
                        }
                }
                
                pks->next = pkr;
                pks = pkr;
        }
        pks->next = NULL;
        
        pktmp = pktmp->next;
        
        long high, low, l = lev;
        
        fmpz_poly_mat_t levhigh, levlow;
        if( ctlev1 >= ctlev2 ) {
                fmpz_poly_mat_init_set(levhigh, ct1);
                fmpz_poly_mat_init_set(levlow, ct2);
                high = ctlev1;
                low = ctlev2;
        }
        else {
                fmpz_poly_mat_init_set(levhigh, ct2);
                fmpz_poly_mat_init_set(levlow, ct1);
                high = ctlev2;
                low = ctlev1;
        }
        ctlev = low - 1;
        
        while( l > high ) {
                pktmp = pktmp->next;
                hp = hp->next;
                l--;
        }
        pktmp = pktmp->next;
        l--;
        
        fmpz_poly_mat_t ctmp, tmpp;
        while( l >= low ) {
                fmpz_poly_mat_init(ctmp, 1 + hp->next->n, 1);
                long ltmp = fmpz_poly_mat_nrows(levhigh), k;
                fmpz_poly_mat_init(tmpp, ltmp * ltmp, 1);
                fmpz_poly_mat_zero(tmpp);
                for( k = 0 ; k < ltmp ; k++) {
                        fmpz_poly_set(fmpz_poly_mat_entry(tmpp, k ,0),fmpz_poly_mat_entry(levhigh, k, 0));
                }
                hcrypt_bgv_refresh(ctmp, tmpp, pktmp->pkb, hp->q, hp->next->q, t, fx, d);
                fmpz_poly_mat_swap(ctmp, levhigh);
                fmpz_poly_mat_clear(ctmp);
                fmpz_poly_mat_clear(tmpp);
                pktmp = pktmp->next;
                hp = hp->next;
                l--;
        }
        fmpz_poly_mat_t c3, ct;
        long row1 = fmpz_poly_mat_nrows(levhigh);
        long row2 = fmpz_poly_mat_nrows(levlow);
        row = row1 * row2;
        fmpz_poly_mat_init(c3, row, 1);
        
        for( i = 0 ; i < row1 ; i++ ) {
                for( j = 0 ; j < row2 ; j++ ){
                        fmpz_poly_mul(fmpz_poly_mat_entry(c3, j + i * row1, 0), fmpz_poly_mat_entry(levhigh, i, 0), fmpz_poly_mat_entry(levlow, j, 0));
                        fmpz_poly_rem_basecase(fmpz_poly_mat_entry(c3, j + i * row1, 0), fmpz_poly_mat_entry(c3, j + i * row1, 0), fx);
                        fmpz_poly_scalar_smod_fmpz(fmpz_poly_mat_entry(c3, j + i * row1, 0), fmpz_poly_mat_entry(c3, j + i * row1, 0), hp->q);
                }
        }
        
        row = 1 + hp->next->n;
        fmpz_poly_mat_init(ct, row, 1);
        hcrypt_bgv_refresh(ct, c3, pktmp->pkb, hp->q, hp->next->q, t, fx, d);
        printf("%ld\n%ld\n1\n", ctlev, row);

        char *output;
        for( i = 0; i < row; i++) {
                output = fmpz_poly_get_str(fmpz_poly_mat_entry(ct, i, 0));
                printf("%s\n", output);
        }
        fmpz_poly_mat_clear(levhigh);
        fmpz_poly_mat_clear(levlow);
	return 0;
}