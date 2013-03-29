#ifndef LIBBGV_H
#define LIBBGV_H
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

#ifdef __cplusplus
extern "C" {
#endif
        
typedef struct bgv_context_t {
        fmpz_t lambda;
        fmpz_t level;
        fmpz_t b;
        fmpz_t t;
} bgv_context_t;

typedef struct pk_node_t {
	fmpz_poly_mat_t pka;
	fmpz_poly_mat_t pkb;
	struct pk_node_t *next;
} pk_node_t;

typedef struct sk_node_t {
	fmpz_poly_mat_t sk;
	struct sk_node_t *next;
} sk_node_t;

typedef struct param_node_t {
	fmpz_t q;
	long n;
	long bign;
	struct param_node_t *next;
} param_node_t;

typedef struct ciphertext_t {
        fmpz_poly_mat_t text;
        int lv;
} ciphertext_t;

static const double pi = 3.1415926;
static const double dvn = 1.0;
static const long bigb = 20;
static long chrnd = 0;

double bgv_get_dvn();
long bgv_get_bigb();
bgv_context_t *bgv_init(bgv_context_t *ctx, char *bb, char *lam, char *lev, char *tt, int base);
void hcrypt_random(fmpz_t r, long len);
fmpz *samplez(fmpz *vec, long d);
void guassian_poly(fmpz *c, fmpz_poly_t poly, long d);
void unif_poly(fmpz_poly_t poly, fmpz_t space, long d);
param_node_t *param_node_init(param_node_t *pnt);
void powers(fmpz_poly_mat_t po, fmpz_poly_mat_t x, fmpz_t qq, fmpz_t t);
void bitdecomp(fmpz_poly_mat_t dc, fmpz_poly_mat_t x, fmpz_t qq, fmpz_t t, long d);
void vec_tensor(fmpz_poly_mat_t tensor, fmpz_poly_mat_t x, fmpz_t qq, fmpz_poly_t fx);
void switchkeygen(fmpz_poly_mat_t mapb, fmpz_poly_mat_t s1, fmpz_poly_mat_t s2, fmpz_t qq, long d, fmpz_t t, fmpz_poly_t fx);
void scale(fmpz_poly_mat_t c2, fmpz_poly_mat_t c1, fmpz_t qq, fmpz_t pp, fmpz_t t);
void switchkey(fmpz_poly_mat_t c3, fmpz_poly_mat_t mapb, fmpz_poly_mat_t c1, fmpz_t qq, fmpz_t t, fmpz_poly_t fx, long d);
void hcrypt_bgv_refresh(fmpz_poly_mat_t c3, fmpz_poly_mat_t c, fmpz_poly_mat_t map, fmpz_t qq, fmpz_t pp, fmpz_t t, fmpz_poly_t fx, long d);
param_node_t *e_setup(long miu, long lamda, long b, fmpz_t t, param_node_t *param);
void e_skeygen(fmpz_poly_mat_t sk, param_node_t *param, long d);
void e_pkeygen(fmpz_poly_mat_t pk, param_node_t *param, fmpz_poly_mat_t sk, long d, fmpz_t t, fmpz_poly_t fx);
void e_encrypt(fmpz_poly_mat_t ct, param_node_t *param, fmpz_poly_mat_t pk, fmpz_poly_t ms, fmpz_t t, fmpz_poly_t fx, long d);
void e_decrypt(fmpz_poly_t ms, param_node_t *param, fmpz_poly_mat_t sk, fmpz_poly_mat_t ct, fmpz_t t, fmpz_poly_t fx);
        
#ifdef __cplusplus
}
#endif
#endif
