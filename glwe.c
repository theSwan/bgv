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

param_node_t *e_setup(long miu, long lamda, long b, fmpz_t t, param_node_t *param)
{
        param = param_node_init(param);
	gen_q(param->q, miu);
	fmpz_t tmp;
	fmpz_init(tmp);
	fmpz_fdiv_q_si(tmp, param->q, bgv_get_bigb());
	long prod;
	prod = lamda * fmpz_flog(tmp, t);
        
	if(b == 0) {
		param->n = prod;
	}  /* LWE */
	else {
		param->n = 1;
	} /* RLWE */
        
	param->bign = ceil((fmpz_get_si(t) * param->n + 1) * fmpz_flog(param->q,t));
        
	fmpz_clear(tmp);
	return param;
}


void e_skeygen(fmpz_poly_mat_t sk, param_node_t *param, long d)
{
        fmpz *coeffs = _fmpz_vec_init(d);
        fmpz_poly_t poly;
        fmpz_poly_init(poly);
        fmpz_poly_set_coeff_si(poly, 0, 1);
        fmpz_poly_set(fmpz_poly_mat_entry(sk, 0, 0), poly);
        long i;
        for( i = 1 ; i <= param->n ; i++ ) {
                guassian_poly(coeffs, fmpz_poly_mat_entry(sk, i, 0), d);
                fmpz_poly_scalar_smod_fmpz(fmpz_poly_mat_entry(sk, i, 0), fmpz_poly_mat_entry(sk, i, 0), param->q);
        }
        
        _fmpz_vec_clear(coeffs, d);
        fmpz_poly_clear(poly);
}

void e_pkeygen(fmpz_poly_mat_t pk, param_node_t *param, fmpz_poly_mat_t sk, long d, fmpz_t t, fmpz_poly_t fx)
{
        fmpz_poly_mat_t ppk, ee, bb, ss, tmp, tmp1;
        fmpz_poly_mat_init(ppk, param->bign, param->n);
        fmpz_poly_mat_init(ee, param->bign, 1);
        fmpz_poly_mat_init(bb, param->bign, 1);
        fmpz_poly_mat_init(ss, param->n, 1);
        fmpz *coeffs = _fmpz_vec_init(d);
        
        long i, j;
        for( i = 0 ; i < param->n ; i++ ) {
                fmpz_poly_set(fmpz_poly_mat_entry(ss, i, 0), fmpz_poly_mat_entry(sk, i+1, 0 ));
        }
        for( i = 0 ; i < param->bign ; i++ ) {
                guassian_poly(coeffs, fmpz_poly_mat_entry(ee, i, 0),d);
        }
        for( i = 0 ; i < param->bign ; i++ ) {
                for( j = 0; j < param->n; j++ ){
                        unif_poly(fmpz_poly_mat_entry(ppk, i, j), param->q, d);
                }
        }
        fmpz_poly_mat_init(tmp, param->bign, 1);
        fmpz_poly_mat_init(tmp1, param->bign, 1);
        fmpz_poly_mat_mul(tmp, ppk, ss);
        fmpz_poly_mat_scalar_mul_fmpz(tmp1, ee, t);
        fmpz_poly_mat_add(bb, tmp, tmp1);
        for( i = 0 ; i < param->bign ; i++ ) {
                fmpz_poly_set(fmpz_poly_mat_entry(pk, i, 0), fmpz_poly_mat_entry(bb, i, 0));
        }
        for( i = 0 ; i < param->bign ; i++ ) {
                for( j = 1; j <= param->n; j++ ){
                        fmpz_poly_neg(fmpz_poly_mat_entry(pk, i, j), fmpz_poly_mat_entry(ppk, i, j-1));
                }
        }
        for( i = 0 ; i < param->bign ; i++) {
                for( j = 0; j < param->n + 1 ; j++) {
                        fmpz_poly_rem_basecase(fmpz_poly_mat_entry(pk, i, j), fmpz_poly_mat_entry(pk, i, j), fx);
                        fmpz_poly_scalar_smod_fmpz(fmpz_poly_mat_entry(pk, i, j), fmpz_poly_mat_entry(pk, i, j), param->q);
                }
        }
        _fmpz_vec_clear(coeffs, d);
        fmpz_poly_mat_clear(tmp);
        fmpz_poly_mat_clear(tmp1);
        fmpz_poly_mat_clear(ee);
        fmpz_poly_mat_clear(ss);
        fmpz_poly_mat_clear(bb);
        fmpz_poly_mat_clear(ppk);
}

void e_encrypt(fmpz_poly_mat_t ct, param_node_t *param, fmpz_poly_mat_t pk, fmpz_poly_t ms, fmpz_t t, fmpz_poly_t fx, long d)
{
        long i, j;
        fmpz_poly_mat_t mm, rr, tmp, tmp1;
        fmpz_poly_mat_init(mm, 1 + param->n, 1);
        fmpz_poly_mat_init(rr, param->bign, 1);
        fmpz_poly_mat_init(tmp, 1 + param->n, 1);
        fmpz_poly_mat_init(tmp1, 1 + (param->n), param->bign);
        for( i = 0 ; i < 1 + param->n ; i++ ) {
                for( j = 0; j < param->bign; j++ ){
                        fmpz_poly_set(fmpz_poly_mat_entry(tmp1, i, j), fmpz_poly_mat_entry(pk, j, i));
                }
        }
        fmpz_poly_mat_zero(mm);
        fmpz_poly_set(fmpz_poly_mat_entry(mm, 0, 0), ms);
        
        for( i = 0; i < param->bign; i++ ) {
                unif_poly(fmpz_poly_mat_entry(rr, i, 0), t, d);
        }
        fmpz_poly_mat_mul(tmp, tmp1, rr);
        fmpz_poly_mat_add(ct, mm, tmp);
        
        for( i = 0; i < param->n + 1 ; i++) {
                fmpz_poly_rem_basecase(fmpz_poly_mat_entry(ct, i, 0), fmpz_poly_mat_entry(ct, i, 0), fx);
                fmpz_poly_scalar_smod_fmpz(fmpz_poly_mat_entry(ct, i, 0), fmpz_poly_mat_entry(ct, i, 0), param->q);
        }
        
        fmpz_poly_mat_clear(tmp);
        fmpz_poly_mat_clear(tmp1);
        fmpz_poly_mat_clear(mm);
        fmpz_poly_mat_clear(rr);
}

void e_decrypt(fmpz_poly_t ms, param_node_t *param, fmpz_poly_mat_t sk, fmpz_poly_mat_t ct, fmpz_t t, fmpz_poly_t fx)
{
        fmpz_poly_t tmp;
        fmpz_poly_init(tmp);
        fmpz_poly_zero(ms);
        
        long i;
        
        for( i = 0 ; i < param->n + 1 ; i++) {
                fmpz_poly_mul(tmp, fmpz_poly_mat_entry(ct, i, 0), fmpz_poly_mat_entry(sk, i, 0));
                fmpz_poly_add(ms, ms, tmp);
        }
        fmpz_poly_rem_basecase(ms, ms, fx);
        fmpz_poly_scalar_smod_fmpz(ms, ms, param->q);
        fmpz_poly_scalar_smod_fmpz(ms, ms, t);
        
        fmpz_poly_clear(tmp);
}
