/*
 * even_order_sample.c - generate an order 2^n difference set
 *
 * Copyright (c) 2023-2024 by Martin Becker, Blaubeuren.
 *
 * This library is free software; you can distribute it and/or modify it
 * under the terms of the Artistic License 2.0 (see the LICENSE file).
 *
 * The licence grants freedom for related software development but does
 * not cover incorporating code or documentation into AI training material.
 * Please contact the copyright holder if you want to use the library whole
 * or in part for other purposes than stated in the licence.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define L0 0LL
#define L1 1LL

/*
 * numbers of stored conway polynomials and of slots in prime factor lists
 */
#define NCP 30
#define NPF 13
#define NEB 60

typedef long long Long;
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
    4202337,            /* 22 */
    8388641,            /* 23 */
    16901801,           /* 24 */
    33554757,           /* 25 */
    67126739,           /* 26 */
    134223533,          /* 27 */
    268443877,          /* 28 */
    536870917,          /* 29 */
    1073948847,         /* 30 */
};

/* exponents, lower 60 bit */
static Fact EL[NCP] = {
    { 1L, 0L},
    { 9L, 21L, 0L},
    { 7L, 73L, 0L},
    { 315L, 585L, 819L, 1365L, 0L},
    { 217L, 1057L, 4681L, 0L},
    { 3591L, 13797L, 37449L, 87381L, 0L},
    { 6223L, 16513L, 299593L, 0L},
    { 69615L, 986895L, 1290555L, 2396745L, 3355443L, 5592405L, 0L},
    { 511L, 1838599L, 19173961L, 0L},
    { 3243933L, 7110873L, 34636833L, 97612893L, 153391689L, 357913941L, 0L},
    { 14329L, 96516119L, 373475417L, 1227133513L, 0L},
    { 630453915L, 941362695L, 1857283155L, 3616814565L, 5286113595L,
      9817068105L, 13743895347L, 22906492245L, 0L},
    { 4529623L, 67117057L, 6958934353L, 78536544841L, 0L},
    { 811597437L, 13050583119L, 34630287489L, 102280151421L,
      628292358729L, 1466015503701L, 0L},
    { 1509346321L, 55759702201L, 233009086681L, 481977699847L,
      1134979744801L, 5026338869833L, 0L},
    { 418239192735L, 1095233372415L, 1167945961455L, 2901803883615L,
      16557351571215L, 21651921285435L, 40210710958665L, 56294995342131L,
      93824992236885L, 0L},
    { 17180000257L, 202518195313L, 1050769861729L, 21862134113449L,
      321685687669321L, 0L},
    { 68585259519L, 206561081853L, 246772582321671L, 948126237341157L,
      2573485501354569L, 6004799503160661L, 0L},
    { 118823881393L, 274878431233L, 4451159405623L, 20587884010836553L, 0L},
    { 872764197279975L, 3483146539597725L, 7635241752363225L,
      18900352534538475L, 28120036697727975L, 37191016277640225L,
      88686269585142075L, 104811045873349725L, 164703072086692425L,
      230584300921369395L, 384307168202282325L, 0L},
    { 14197294936951L, 99457304386111L, 27369056489183311L,
      72624976668147841L, 126347562148695559L, 164703072086692425L, 0L},
    { 123085172783097L, 3537755971368759L, 108033640255985661L,
      829067149380204567L, 1101298153654301589L, 902286394909706329L,
      164703072086692425L, 384307168202282325L, 0L},
    { 58720249L, 3307331370614831L, 1030270280712501553L,
      164703072086692425L, 0L},
    { 121908420447366735L, 529864617590675631L, 1148137597948727279L,
      666367475139737243L, 126347562148695559L, 810161057291297875L,
      667480871088174565L, 1085102592571150095L, 88686269585142075L,
      164703072086692425L, 230584300921369395L, 384307168202282325L, 0L},
    { 3575112450587167L, 374787272576235967L, 224054706614323399L,
      602358323538352663L, 7635241752363225L, 37191016277640225L,
      164703072086692425L, 0L},
    { 13512448149528573L, 184343569761639895L, 4504149450301441L,
      1139412354790875133L, 321066748118362449L, 164703072086692425L,
      384307168202282325L, 0L},
    { 24751301355248209L, 1134907174682624511L, 562543471243439825L,
      892813876301792799L, 126347562148695559L, 164703072086692425L, 0L},
    { 153760103770322799L, 1151219461418647165L, 609943077314748995L,
      27369056489183311L, 72624976668147841L, 1081501588392263535L,
      938424480493945213L, 795118279039204811L, 88686269585142075L,
      164703072086692425L, 230584300921369395L, 384307168202282325L, 0L},
    { 15697568566729L, 652503336579024719L, 864829103742905319L,
      288491691089292625L, 865928168696129703L, 164703072086692425L, 0L},
    { 2005509207195591L, 685145279195171857L, 78566758634064057L,
      3483146539597725L, 7635241752363225L, 126347562148695559L,
      37191016277640225L, 667480871088174565L, 104811045873349725L,
      164703072086692425L, 384307168202282325L, 0L},
};

/* exponents, higher 60 bit */
static Fact EH[NCP] = {
    { 0L, 0L},
    { 0L, 0L, 0L},
    { 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L},
    { 0L, 0L, 0L, 0L, 0L, 1L, 0L},
    { 0L, 0L, 0L, 0L, 0L, 2L, 9L, 21L, 0L},
    { 0L, 0L, 10L, 73L, 0L},
    { 0L, 9L, 16L, 37L, 56L, 110L, 215L, 240L, 315L, 585L, 819L, 1365L, 0L},
    { 0L, 0L, 18L, 54L, 217L, 1057L, 4681L, 0L},
    { 0L, 2L, 32L, 95L, 3318L, 37449L, 87381L, 0L},
    { 0L, 7L, 29L, 808L, 28728L, 299593L, 0L},
    { 1161L, 3095L, 11740L, 49784L, 132104L, 148470L, 390167L, 578524L,
      1290555L, 2396745L, 3355443L, 5592405L, 0L},
    { 0L, 32132L, 64249L, 121684L, 576041L, 19173961L, 0L},
    { 57L, 46061L, 1701651L, 3243933L, 7110873L, 14708792L, 34636833L,
      56512727L, 97612893L, 153391689L, 357913941L, 0L},
};

static Long *lg;                        /* logarithm table */
static Long *ex;                        /* inverse logarithm table */
static Long mb;                         /* 2^n - 1 */
static Long mo;                         /* modulus */
static Long p0, p1, p2;                 /* primitive polynomial coefficients */
static Long P0, P1, P2;                 /* coefficient logarithms */
static Long from, until;                /* iteration: from <= i < until */
static Long start0, start1, start2;     /* start polynomial coefficients */

/* print a difference set element */
static void display(Long num)
{
    (void) printf("%lld\n", num);
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

/* boolean whether x^(nl + 2^60*nh) is equal to 1 */
static int x_pow_n_eq_1(Long nl, Long nh)
{
    Long a = L0;
    Long b = L1;
    Long c = L0;
    Long d = L0;
    Long e = L0;
    Long f = L1;
    for (int r = 0; r < NEB; ++r) {
        if (nl & L1) multiply(d, e, f, a, b, c, &d, &e, &f);
        nl >>= 1;
        if (nl || nh) multiply(a, b, c, a, b, c, &a, &b, &c);
    }
    while (nh) {
        if (nh & L1) multiply(d, e, f, a, b, c, &d, &e, &f);
        if (nh >>= 1) multiply(a, b, c, a, b, c, &a, &b, &c);
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
    Long *el = EL[exp-1];       /* exponents for maximal order check LO */
    Long *eh = EH[exp-1];       /* exponents for maximal order check HI */
    int nf = 0;                 /* number of distinct prime factors */
    while (el[nf] || eh[nf]) ++nf;
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
                    if (x_pow_n_eq_1(el[i], eh[i])) {
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
    if (mo < 0) bail_out("insufficient word size");
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
    (void) printf("# primitive polynomial is x^3 + %lld*x^2 + %lld*x + %lld\n",
        p2, p1, p0);
    (void) fflush(stdout);
    if (!part) {
        (void) printf("# dry run: none out of %lld\n", mo);
    }
    else if (parts == 1) {
        from = L0;
        until = mo;
        start2 = L0;
        start1 = L0;
        start0 = L1;
        (void) printf("# complete run: all out of %lld\n", mo);
    }
    else {
        from = mo * (part - 1) / parts;
        until = mo * part / parts;
        if (until <= from) {
            (void) printf("# partial run %d/%d none out of %lld\n",
                part, parts, mo);
            exit(0);
        }
        (void) printf("# partial run %d/%d from %lld to %lld out of %lld\n",
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
        "# x^%lld = %lld %lld %lld (mod 1 %lld %lld %lld)\n",
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
        "# x^%lld = %lld %lld %lld (mod 1 %lld %lld %lld)\n",
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
