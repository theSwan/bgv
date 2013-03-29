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

int main(int argc, char *args[])  /*bb, lam, lev, tt*/
{
        bgv_context_t *ctx;
        ctx = bgv_init(ctx, args[1], args[2], args[3], args[4], 10);
        long bb = fmpz_get_si(ctx->b);
        long lam = fmpz_get_si(ctx->lambda);
        long lev = fmpz_get_si(ctx->level);
        
        fmpz_t mult;
        fmpz_init(mult);
        fmpz_mul(mult, ctx->lambda, ctx->level);
        long miu;
        miu = fmpz_flog(mult, ctx->t);
        param_node_t *param;
        param = param_node_init(param);
        long high;
        high = (lev + 1) * miu;
        hcrypt_random(param->q, high);
        fmpz_t tmp;
        fmpz_init(tmp);
        fmpz_fdiv_q_si(tmp, param->q, bgv_get_bigb());
        long prod, d;
        prod = lam * fmpz_flog(tmp, ctx->t);
        
        if(bb == 0) {
		param->n = prod;
		d = 1;
	}
	else {
		param->n = 1;
		d = prod;
	}
        param->bign = ceil((fmpz_get_si(ctx->t) * param->n + 1) * fmpz_flog(param->q,ctx->t));
        
	param_node_t *r, *pn;
        
        r = param;
        long j;
        for(j = lev - 1 ; j >= 0 ; j--) {
                pn = e_setup((j+1)*miu, lam, bb, ctx->t, pn);
                r->next = pn;
                r = pn;
        }
        r->next = NULL;
        fmpz_print(ctx->t);
        printf(" ");
        fmpz_print(ctx->level);
        printf(" %ld ", d);
        r = param;
        while(r!=NULL) {
                fmpz_print(r->q);
                printf(" %ld %ld ", r->n, r->bign);
                r = r->next;
        }
        
        fmpz_clear(tmp);
        return 0;
}

/* output d {q, n, bign}
 */

