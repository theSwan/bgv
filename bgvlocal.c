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


long bgv_get_bigb()
{
        return bigb;
}
param_node_t *param_node_init(param_node_t *pnt)
{
	pnt = (param_node_t *)malloc(sizeof(param_node_t));
	pnt->next = NULL;
	pnt->n = pnt->bign = 0;
	fmpz_init(pnt->q);
	return pnt;
}

void powers(fmpz_poly_mat_t po, fmpz_poly_mat_t x, fmpz_t qq, fmpz_t t)
{
	long xrow = fmpz_poly_mat_nrows(x);
	long len = fmpz_clog(qq, t);
	long qrow = xrow * len;
	long i, j;
	fmpz_poly_mat_init(po, qrow, 1);
	for( i = 0 ; i < xrow ; i++) {
		fmpz_poly_set(fmpz_poly_mat_entry(po, i, 0), fmpz_poly_mat_entry(x, i, 0));
	}
	for( i = 1 ; i < len ; i++) {
		for( j = 0 ; j < xrow ; j++) {
			fmpz_poly_scalar_mul_fmpz(fmpz_poly_mat_entry(po, j + i * xrow, 0), fmpz_poly_mat_entry(po, j + (i-1)*xrow, 0), t);
			fmpz_poly_scalar_smod_fmpz(fmpz_poly_mat_entry(po, j + i * xrow, 0), fmpz_poly_mat_entry(po , j + i * xrow, 0), qq);
		}
	}
}

void bitdecomp(fmpz_poly_mat_t dc, fmpz_poly_mat_t x, fmpz_t qq, fmpz_t t, long d)
{
	long xrow = fmpz_poly_mat_nrows(x);
	long len = fmpz_clog(qq, t);
	long i, j, k;
	fmpz_mat_t bits;
	fmpz_mat_init(bits, d, len);
	fmpz_t hold;
	fmpz_init(hold);
	fmpz_poly_t xtmp;
	for( i = 0 ; i < xrow ; i++ ) {
		fmpz_mat_zero(bits);
		for( j = 0 ; j < d ; j++) {
			fmpz_poly_get_coeff_fmpz(hold, fmpz_poly_mat_entry(x, i, 0), j);
			k = 0;
                        if (fmpz_cmp_si(hold, 0)<0) {
                                fmpz_neg(hold, hold);
                                fmpz_t negval;
                                fmpz_init(negval);
                                while ( !fmpz_is_zero(hold) ) {
                                        fmpz_mod(negval, hold, t);
                                        fmpz_neg(fmpz_mat_entry(bits, j, k), negval);
                                        fmpz_tdiv_q(hold, hold, t);
                                        k++;
                                }
                                fmpz_clear(negval);
                        }
                        else {
                                while ( !fmpz_is_zero(hold) ) {
                                        fmpz_mod(fmpz_mat_entry(bits, j, k), hold, t);
                                        fmpz_tdiv_q(hold, hold, t);
                                        k++;
                                }
                        }
		}
		
		for( j = 0 ; j < len ; j++ ) {
			fmpz_poly_init(xtmp);
			for( k = 0; k < d ; k++ ) {
				fmpz_poly_set_coeff_fmpz(xtmp, k, fmpz_mat_entry(bits, k, j));
			}
                        
			fmpz_poly_set(fmpz_poly_mat_entry(dc, i + j * xrow, 0), xtmp);
			fmpz_poly_clear(xtmp);
		}
	}
	fmpz_clear(hold);
	fmpz_mat_clear(bits);
        
}

void vec_tensor(fmpz_poly_mat_t tensor, fmpz_poly_mat_t x, fmpz_t qq, fmpz_poly_t fx)
{
        long row1 = fmpz_poly_mat_nrows(x);
        
        long i, j;
        for( i = 0 ; i < row1 ; i++ ) {
                for( j = 0 ; j < row1 ; j++ ){
                        fmpz_poly_mul(fmpz_poly_mat_entry(tensor,j+i*row1,0),fmpz_poly_mat_entry(x,i,0),fmpz_poly_mat_entry(x,j,0));
			fmpz_poly_rem_basecase(fmpz_poly_mat_entry(tensor,j+i*row1,0), fmpz_poly_mat_entry(tensor,j+i*row1,0), fx);
        		fmpz_poly_scalar_smod_fmpz(fmpz_poly_mat_entry(tensor,j+i*row1,0), fmpz_poly_mat_entry(tensor,j+i*row1,0), qq);
                }
        }
}


void switchkeygen(fmpz_poly_mat_t mapb, fmpz_poly_mat_t s1, fmpz_poly_mat_t s2, fmpz_t qq, long d, fmpz_t t, fmpz_poly_t fx)
{
	fmpz_poly_mat_t sp1;
	param_node_t *param;
	param = (param_node_t *)malloc(sizeof(param_node_t));
	long n1, n2, i;
	n1 = fmpz_poly_mat_nrows(s1);
	n2 = fmpz_poly_mat_nrows(s2);
	param->n = n2 - 1;
	param->bign = n1 * fmpz_clog(qq, t);
	fmpz_init_set(param->q, qq);
	param->next = NULL;
	e_pkeygen(mapb, param, s2, d, t, fx);
	long s1row = fmpz_poly_mat_nrows(s1);
	long len = fmpz_clog(qq, t);
	long qrow = s1row * len;
	powers(sp1, s1, qq, t);
	for( i = 0 ; i < param->bign ; i++) {
		fmpz_poly_add(fmpz_poly_mat_entry(mapb, i, 0), fmpz_poly_mat_entry(mapb, i, 0), fmpz_poly_mat_entry(sp1, i, 0));
		fmpz_poly_rem_basecase(fmpz_poly_mat_entry(mapb, i, 0), fmpz_poly_mat_entry(mapb, i, 0), fx);
                fmpz_poly_scalar_smod_fmpz(fmpz_poly_mat_entry(mapb, i, 0), fmpz_poly_mat_entry(mapb, i, 0), qq);
        }
	fmpz_poly_mat_clear(sp1);
	free(param);
}

void scale(fmpz_poly_mat_t c2, fmpz_poly_mat_t c1, fmpz_t qq, fmpz_t pp, fmpz_t t)
{
        long row, col, i, j, len, k;
        row = fmpz_poly_mat_nrows(c1);
        col = fmpz_poly_mat_ncols(c1);
        fmpz_poly_t poly;
        fmpz_poly_init(poly);
        fmpz_t coeff, tmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
        fmpz_init(coeff);
        fmpz_init(tmp);
        fmpz_init(tmp1);
        fmpz_init(tmp2);
        fmpz_init(tmp3);
        fmpz_init(tmp4);
        fmpz_init(tmp5);
        fmpz_init(tmp6);
        
        for( i = 0 ; i < row ; i++ ) {
                for( j = 0 ; j < col ; j++ ) {
                        fmpz_poly_set(poly, fmpz_poly_mat_entry(c1, i, j));
                        len = fmpz_poly_length(poly);
                        for( k = 0 ; k < len ; k++ ) {
                                fmpz_poly_get_coeff_fmpz(tmp, poly, k);
                                fmpz_mod(tmp1, tmp, t);   /* tmp1 = base = tmp % r */
                                
                                if(fmpz_cmp_si(tmp, 0) > 0) {
                                        fmpz_mul(coeff, tmp, pp);
                                        fmpz_fdiv_q(tmp2, coeff, qq);
                                        fmpz_fdiv_q(tmp3, tmp2, t);
                                        
                                        fmpz_mul(tmp4, tmp3, t);
                                        fmpz_add(tmp5, tmp4, tmp1);
                                        fmpz_poly_set_coeff_fmpz(poly, k, tmp5);
                                }
                                else if(fmpz_cmp_si(tmp, 0) < 0) {
                                        fmpz_neg(tmp2, tmp);
                                        fmpz_mul(coeff, tmp2, pp);
                                        fmpz_fdiv_q(tmp3, coeff, qq);
                                        fmpz_fdiv_q(tmp4, tmp3, t);
                                        fmpz_mul(tmp5, tmp4, t);
                                        fmpz_sub(tmp6, tmp1, tmp5);
                                        fmpz_poly_set_coeff_fmpz(poly, k, tmp6);
                                }
                                else  {
                                        fmpz_poly_set_coeff_fmpz(poly, k, tmp1);
                                }
                        }
                        fmpz_poly_set(fmpz_poly_mat_entry(c2, i, j), poly);
                }
        }
        fmpz_poly_clear(poly);
        fmpz_clear(coeff);
        fmpz_clear(tmp);
        fmpz_clear(tmp1);
        fmpz_clear(tmp2);
        fmpz_clear(tmp3);
        fmpz_clear(tmp4);
        fmpz_clear(tmp5);
        fmpz_clear(tmp6);
}

void switchkey(fmpz_poly_mat_t c3, fmpz_poly_mat_t mapb, fmpz_poly_mat_t c1, fmpz_t qq, fmpz_t t, fmpz_poly_t fx, long d)
{
	fmpz_poly_mat_t bd, bdt;
	long c1row = fmpz_poly_mat_nrows(c1);
	long len = fmpz_clog(qq, t);
	long qrow = c1row * len;
	fmpz_poly_mat_init(bd, qrow, 1);
	bitdecomp(bd, c1, qq, t, d);
	long bdtrow, bdtcol, i, j;
	bdtrow = fmpz_poly_mat_ncols(bd);
	bdtcol = fmpz_poly_mat_nrows(bd);
	fmpz_poly_mat_init(bdt, bdtrow, bdtcol);
        for( i = 0 ; i < bdtrow ; i++ ) {
		for( j = 0 ; j < bdtcol ; j++ ) {
			fmpz_poly_set(fmpz_poly_mat_entry(bdt, i, j), fmpz_poly_mat_entry(bd, j, i));
		}
	}
	long col = fmpz_poly_mat_ncols(mapb);
	fmpz_poly_mat_mul(c3, bdt, mapb);
	for( i = 0 ; i < bdtrow ; i++ ) {
		for( j = 0 ; j < col ; j++ ) {
			fmpz_poly_rem_basecase(fmpz_poly_mat_entry(c3, i, j), fmpz_poly_mat_entry(c3, i, j), fx);
                	fmpz_poly_scalar_smod_fmpz(fmpz_poly_mat_entry(c3, i, j), fmpz_poly_mat_entry(c3, i, j), qq);
		}
	}
	fmpz_poly_mat_clear(bd);
	fmpz_poly_mat_clear(bdt);
}

void hcrypt_bgv_refresh(fmpz_poly_mat_t c3, fmpz_poly_mat_t c, fmpz_poly_mat_t map, fmpz_t qq, fmpz_t pp, fmpz_t t, fmpz_poly_t fx, long d)
{
        fmpz_poly_mat_t c1;
        powers(c1, c, qq, t);
        fmpz_poly_mat_t c2;
	long row, col, len;
	row = fmpz_poly_mat_nrows(c1);
        col = fmpz_poly_mat_ncols(c1);
        fmpz_poly_mat_init(c2, row, col);
        scale(c2, c1, qq, pp, t);
        switchkey(c3, map, c2, pp, t, fx, d);
        
	fmpz_poly_mat_clear(c1);
	fmpz_poly_mat_clear(c2);
}


bgv_context_t *bgv_init(bgv_context_t *ctx, char *bb, char *lam, char *lev, char *tt, int base)
{
        ctx = (bgv_context_t *)malloc(sizeof(bgv_context_t));
        
        fmpz_init(ctx->b);
        int flag = fmpz_set_str(ctx->b, bb, base);
        if(flag == -1){
                printf("invalid number b\n");
                exit(0);
        }
        long val = fmpz_get_si(ctx->b);
        if(val != 1 && val != 0) {
                printf("the value of first parameter should be 0 or 1\n");
                exit(0);
        }
        
        fmpz_init(ctx->lambda);
        flag = fmpz_set_str(ctx->lambda, lam, base);
        if(flag == -1){
                printf("invalid number lambda\n");
                exit(0);
        }
        
        fmpz_init(ctx->t);
        flag = fmpz_set_str(ctx->t, tt, base);
        if(flag == -1){
                printf("invalid number t\n");
                exit(0);
        }
        
        fmpz_init(ctx->level);
        flag = fmpz_set_str(ctx->level, lev, base);
        if(flag == -1){
                printf("invalid number level\n");
                exit(0);
        }
        
        return ctx;
}