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

ciphertext_t *hcrypt_bgv_mul(ciphertext_t *c, pk_node_t *pbk, ciphertext_t *c1, ciphertext_t *c2)
{
	c = (ciphertext_t *)malloc(sizeof(ciphertext_t));
	pk_node_t *pktmp;
	pktmp = keylist->pubkey;
        param_node_t *hp;
        hp = param;
        
	int l1 = c1->lv, l2 = c2->lv, l, high, low;
        fmpz_poly_mat_t levhigh, levlow;
        if( l1 >= l2 ) {
                fmpz_poly_mat_init_set(levhigh, c1->text);
                fmpz_poly_mat_init_set(levlow, c2->text);
                high = l1;
                low = l2;
        }
        else {
                fmpz_poly_mat_init_set(levhigh, c2->text);
                fmpz_poly_mat_init_set(levlow, c1->text);
                high = l2;
                low = l1;
        }
        l = bgv_get_level();
        c->lv = low - 1;
        
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
                hcrypt_bgv_refresh(ctmp, tmpp, pktmp->pkb, hp->q, hp->next->q, t);
                fmpz_poly_mat_swap(ctmp, levhigh);
                fmpz_poly_mat_clear(ctmp);
                fmpz_poly_mat_clear(tmpp);
                pktmp = pktmp->next;
                hp = hp->next;
                l--;
        }
        fmpz_poly_mat_t c3;
        long row1 = fmpz_poly_mat_nrows(levhigh);
        long row2 = fmpz_poly_mat_nrows(levlow);
        long row = row1 * row2;
        fmpz_poly_mat_init(c3, row, 1);
        long i, j;
        
        for( i = 0 ; i < row1 ; i++ ) {
                for( j = 0 ; j < row2 ; j++ ){
                        fmpz_poly_mul(fmpz_poly_mat_entry(c3, j + i * row1, 0), fmpz_poly_mat_entry(levhigh, i, 0), fmpz_poly_mat_entry(levlow, j, 0));
                        fmpz_poly_rem_basecase(fmpz_poly_mat_entry(c3, j + i * row1, 0), fmpz_poly_mat_entry(c3, j + i * row1, 0), fx);
                        fmpz_poly_scalar_smod_fmpz(fmpz_poly_mat_entry(c3, j + i * row1, 0), fmpz_poly_mat_entry(c3, j + i * row1, 0), hp->q);
                }
        }
        printf("mul\n");
        
        fmpz_poly_mat_init(c->text, hp->next->n + 1, 1);
        hcrypt_bgv_refresh(c->text, c3, pktmp->pkb, hp->q, hp->next->q, t);
        fmpz_poly_mat_clear(levhigh);
        fmpz_poly_mat_clear(levlow);
	return c;
}