/*
 * even_order_sample.c - generate an order 2^n difference set
 *
 * Copyright (c) 2023 by Martin Becker, Blaubeuren.
 *
 * This library is free software; you can distribute it and/or modify it
 * under the terms of the Artistic License 2.0 (see the LICENSE file).
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define Long long long
#define L0 0LL
#define L1 1LL

/*
 * numbers of stored conway polynomials and of slots in prime factor lists
 * (higher orders would be doable but need >64 bit arithmetic for exponents)
 */
#define NCP 21
#define NPF 12

typedef Long Fact[NPF];         /* 0-terminated prime factor list */

/* order 2^(i+1) Conway polynomials in binary representation */
static Long CP[NCP] = {
    3,                  /* 1 */
    7,                  /* 2 */
    11,                 /* 3 */
    19,                 /* 4 */
    37,                 /* 5 */
    91,                 /* 6 */
    131,                /* 7 */
    285,                /* 8 */
    529,                /* 9 */
    1135,               /* 10 */
    2053,               /* 11 */
    4331,               /* 12 */
    8219,               /* 13 */
    16553,              /* 14 */
    32821,              /* 15 */
    65581,              /* 16 */
    131081,             /* 17 */
    267267,             /* 18 */
    524327,             /* 19 */
    1050355,            /* 20 */
    2097253,            /* 21 */
};

/* prime factors of order 2^(i+1) moduli without repetition */
static Fact PF[NCP] = {
    {7, 0},
    {3, 7, 0},
    {7, 73, 0},
    {3, 5, 7, 13, 0},
    {7, 31, 151, 0},
    {3, 7, 19, 73, 0},
    {7, 127, 337, 0},
    {3, 5, 7, 13, 17, 241, 0},
    {7, 73, 262657, 0},
    {3, 7, 11, 31, 151, 331, 0},
    {7, 23, 89, 599479, 0},
    {3, 5, 7, 13, 19, 37, 73, 109, 0},
    {7, 79, 8191, 121369, 0},
    {3, 7, 43, 127, 337, 5419, 0},
    {7, 31, 73, 151, 631, 23311, 0},
    {3, 5, 7, 13, 17, 97, 241, 257, 673, 0},
    {7, 103, 2143, 11119, 131071, 0},
    {3, 7, 19, 73, 87211, 262657, 0},
    {7, 32377, 524287, 1212847, 0},
    {3, 5, 7, 11, 13, 31, 41, 61, 151, 331, 1321, 0},
    {7, 73, 127, 337, 92737, 649657, 0},
};

static Long *lg;                        /* logarithm table */
static Long *ex;                        /* inverse logarithm table */
static Long mb;                         /* 2^n - 1 */
static Long mo;                         /* modulus */
static Long m3;                         /* 2^3n - 1 */
static Long p0, p1, p2;                 /* primitive polynomial coefficients */
static Long P0, P1, P2;                 /* coefficient logarithms */
static Long from, until;                /* iteration: from <= i < until */
static Long start0, start1, start2;     /* start polynomial coefficients */

/* print a difference set element */
static void display(Long num)
{
    (void) printf("%Ld\n", num);
    (void) fflush(stdout);
}

/* abort with an error message */
static void bail_out(char *msg)
{
    (void) fprintf(stderr, "%s\n", msg);
    exit(1);
}

/* multiply two 3D vectors (representing polynomials) */
static void multiply(
    Long a, Long b, Long c,
    Long d, Long e, Long f,
    Long *i, Long *j, Long *k)
{
    /*
     * (ax^2 + bx + c) * (dx^2 + ex + f) =
     * ad x^4 + (ae + bd) x^3 + (af + be + cd) x^2 + (bf + ce) x + cf
     */
    Long A, B, C, D, E, F;
    if(a) A = lg[a];
    if(b) B = lg[b];
    if(c) C = lg[c];
    if(d) D = lg[d];
    if(e) E = lg[e];
    if(f) F = lg[f];
    Long gg = (a && d)? ex[A+D]: L0;
    Long hh = (a && e)? ex[A+E]: L0;
    if (b && d) hh ^= ex[B+D];
    Long ii = (a && f)? ex[A+F]: L0;
    if (b && e) ii ^= ex[B+E];
    if (c && d) ii ^= ex[C+D];
    Long jj = (b && f)? ex[B+F]: L0;
    if (c && e) jj ^= ex[C+E];
    Long kk = (c && f)? ex[C+F]: L0;
    /*
     * g x^4 + h x^3 + i x^2 + j x + k ===
     * (h + p2g) x^3 + (i + p1g) x^2 + (j + p0g) x + k === ...
     * (mod x^3 + p2 x^2 + p1 x + p0)
     */
    if (p2) {
        if (gg) {
            Long G = lg[gg];
            hh ^= ex[P2+G];
            ii ^= ex[P1+G];
            jj ^= ex[P0+G];
        }
        if (hh) {
            Long H = lg[hh];
            ii ^= ex[P2+H];
            jj ^= ex[P1+H];
            kk ^= ex[P0+H];
        }
    }
    else if (p1 != L1) {
        if (gg) {
            Long G = lg[gg];
            ii ^= ex[P1+G];
            jj ^= ex[P0+G];
        }
        if (hh) {
            Long H = lg[hh];
            jj ^= ex[P1+H];
            kk ^= ex[P0+H];
        }
    }
    else {
        if (gg) {
            ii ^= gg;
            jj ^= ex[P0+lg[gg]];
        }
        if (hh) {
            jj ^= hh;
            kk ^= ex[P0+lg[hh]];
        }
    }
    *i = ii;
    *j = jj;
    *k = kk;
}

/* calculate the n-th power of the generator x */
static void x_pow_n(
    Long n,
    Long *d, Long *e, Long *f)
{
    Long a = L0;
    Long b = L1;
    Long c = L0;
    Long dd = L0;
    Long ee = L0;
    Long ff = L1;
    while (n) {
        if (n & L1) multiply(dd, ee, ff, a, b, c, &dd, &ee, &ff);
        if (n >>= 1) multiply(a, b, c, a, b, c, &a, &b, &c);
    }
    *d = dd;
    *e = ee;
    *f = ff;
}

/* boolean whether x^n is equal to 1 */
static int x_pow_n_eq_1(Long n)
{
    Long a = L0;
    Long b = L1;
    Long c = L0;
    Long d = L0;
    Long e = L0;
    Long f = L1;
    while (n) {
        if (n & L1) multiply(d, e, f, a, b, c, &d, &e, &f);
        if (n >>= 1) multiply(a, b, c, a, b, c, &a, &b, &c);
    }
    return !d && !e && f == L1;
}

/* boolean whether x^2^n is equal to x */
static int x_pow_2_pow_n_eq_x(int n)
{
    Long a = L0;
    Long b = L1;
    Long c = L0;
    for (int i = 0; i < n; ++i) {
        multiply(a, b, c, a, b, c, &a, &b, &c);
    }
    return !a && !c && b == L1;
}

/* evaluate the polynomial x^3 + p2*x^2 + p1*x + p0 at position x */
static Long evaluate(Long x)
{
    /* precondition: x != 0, p0 != 0 */
    Long X = lg[x];
    Long y;
    if (p2) {
        /* y = ((x + p2) * x + p1) * x + p0 */
        y = x ^ p2;
        y = y? ex[lg[y]+X] ^ p1: p1;
        y = y? ex[lg[y]+X] ^ p0: p0;
    }
    else {
        /* y = ((x * x + p1) * x + p0 */
        y = ex[X+X] ^ p1;
        y = y? ex[lg[y]+X] ^ p0: p0;
    }
    return y;
}

/*
 * find some primitive polynomial of third degree over GF(2^exp),
 * preferably a trinomial x^3 + x + c
 */
static void find_pp(int exp)
{
    Fact es;                    /* exponents for maximal order check */
    Long *fs = PF[exp-1];
    int nf = 0;                 /* number of distinct prime factors */
    while (fs[nf]) ++nf;
    for (int i = 0; i < nf; ++i) {
        es[i] = m3 / fs[nf - 1 - i];
    }
    /* iterate over candidate polynomials */
    for (Long ax = L0; ax <= mb; ++ax) {
        Long a = (ax? mb - ax + L1: L0);
        for (Long b = ax? L0: L1; b <= mb; ++b) {
            for (Long c = mb; c > 0; --c) {
                p2 = a;
                p1 = b;
                p0 = c;
                if (p2) P2 = lg[a];
                if (p1) P1 = lg[b];
                if (p0) P0 = lg[c];
                /* order criteria */
                if (!x_pow_2_pow_n_eq_x(3*exp)) {
                    continue;
                }
                for (int i = 0; i < nf; ++i) {
                    if (x_pow_n_eq_1(es[i])) {
                        goto NEXTPOLY;
                    }
                }
                /* irreducibility (brute-force, Berlekamp would be better) */
                for (Long x = 1; x <= mb; ++x) {
                    if (!evaluate(x)) {
                        goto NEXTPOLY;
                    }
                }
                return;
                NEXTPOLY: ;
            }
        }
    }
    bail_out("found no primitive polynomial");
}

/* initialize global variables */
static void prepare(int exp, int part, int parts)
{
    Long ge = CP[exp-1];
    mb = (L1 << exp) - L1;
    mo = (((L1 << exp) + L1) << exp) + L1;
    m3 = (((mb << exp) | mb) << exp) | mb;
    if (m3 < 0) bail_out("insufficient word size");
    (void) printf("# preparing tables...\n");
    (void) fflush(stdout);
    lg = (Long *) calloc(mb + L1, sizeof (Long));
    if (lg == NULL) bail_out("out of memory");
    ex = (Long *) calloc(mb<<1, sizeof (Long));
    if (ex == NULL) bail_out("out of memory");
    Long y = L1;
    for (Long i = L0; i < mb; ++i) {
        lg[y] = i;
        ex[i] = ex[i+mb] = y;
        y <<= 1;
        if (y > mb) y ^= ge;
    }
    find_pp(exp);
    (void) printf("# primitive polynomial is x^3 + %Ld*x^2 + %Ld*x + %Ld\n",
        p2, p1, p0);
    (void) fflush(stdout);
    if (!part) {
        (void) printf("# dry run: none out of %Ld\n", mo);
    }
    else if (parts == 1) {
        from = L0;
        until = mo;
        start2 = L0;
        start1 = L0;
        start0 = L1;
        (void) printf("# complete run: all out of %Ld\n", mo);
    }
    else {
        from = mo * (part - 1) / parts;
        until = mo * part / parts;
        if (until <= from) {
            (void) printf("# partial run %d/%d none out of %Ld\n",
                part, parts, mo);
            exit(0);
        }
        (void) printf("# partial run %d/%d from %Ld to %Ld out of %Ld\n",
            part, parts, from, until - L1, mo);
        x_pow_n(from, &start2, &start1, &start0);
    }
    (void) printf("# done preparing\n");
    (void) fflush(stdout);
}

/* main loop: generate difference set elements */
static void generate()
{
    Long x2 = start2;
    Long x1 = start1;
    Long x0 = start0;
    (void) printf(
        "# x^%Ld = %Ld %Ld %Ld (mod 1 %Ld %Ld %Ld)\n",
        from, start2, start1, start0, p2, p1, p0);
    if (p2) {
        for (Long i = from; i < until; ++i) {
            if (x2) {
                Long c = lg[x2];
                x2 = x1 ^ ex[c+P2];
                x1 = x0 ^ ex[c+P1];
                x0 = ex[c+P0];
            }
            else {
                display(i);
                x2 = x1;
                x1 = x0;
                x0 = L0;
            }
        }
    }
    else if (p1 != L1) {
        for (Long i = from; i < until; ++i) {
            if (x2) {
                Long c = lg[x2];
                x2 = x1;
                x1 = x0 ^ ex[c+P1];
                x0 = ex[c+P0];
            }
            else {
                display(i);
                x2 = x1;
                x1 = x0;
                x0 = L0;
            }
        }
    }
    else {
        for (Long i = from; i < until; ++i) {
            if (x2) {
                Long x3 = x2;
                x2 = x1;
                x1 = x0 ^ x3;
                x0 = ex[lg[x3]+P0];
            }
            else {
                display(i);
                x2 = x1;
                x1 = x0;
                x0 = L0;
            }
        }
    }
    (void) printf(
        "# x^%Ld = %Ld %Ld %Ld (mod 1 %Ld %Ld %Ld)\n",
        until, x2, x1, x0, p2, p1, p0);
    if (until >= mo) (void) printf("\n");
    (void) fflush(stdout);
}

/* parse command line and run */
int main(int argc, char *argv[])
{
    int exp = 0;
    int part = 1;
    int parts = 1;
    if (argc >= 2) {
        char *eptr;
        exp = strtol(argv[1], &eptr, 0);
        if (*eptr) exp = 0;
    }
    if (argc >= 3) {
        char *mptr, *eptr;
        parts = 0;
        part = strtol(argv[2], &mptr, 0);
        if (*argv[2] && *mptr == '/' && 0 <= part) {
            int num = strtol(++mptr, &eptr, 0);
            if (*mptr && !*eptr && part <= num) parts = num;
        }
    }
    if (exp <= 0 || NCP < exp || !parts || 3 < argc) {
        char *arg0 = argv[0];
        char *sptr = strrchr(arg0, '/');
        if (sptr != NULL && sptr[1]) arg0 = sptr+1;
        (void) fprintf(stderr, "usage: %s exponent [part/parts]\n", arg0);
        exit(1);
    }
    (void) printf("# generating order 2^%d set\n", exp);
    prepare(exp, part, parts);
    if (part) generate();
    return 0;
}
