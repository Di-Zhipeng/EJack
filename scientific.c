/* scientific.c
 * --------------------------------
 * The implementation of scientific functions
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numeric.h"
#include "operation.h"
#include "scientific.h"
#include "shift.h"

/* Definitions of some constants */
static Numeric zero;
static Numeric one;
static Numeric two;
static Numeric three;
static Numeric ten;
static Numeric tenth;
Numeric PI;
static Numeric ln10;
static Numeric e;
static Numeric E;
static Udbyte eTable[CALSTACK_FIXED_POINT] = {71828, 18284, 59045, 23536, 2874, 71352, 66249, 77572, 47093, 69995,
                                              95749, 66967,};
static Udbyte ln10Table[CALSTACK_FIXED_POINT] = {30258, 50929, 94045, 68401, 79914, 54684, 36420, 76011, 1488, 62877,
                                                 29760, 33327,};
static Udbyte PItable[CALSTACK_FIXED_POINT] = {14159, 26535, 89793, 23846, 26433, 83279, 50288, 41971, 69399, 37510,
                                               58209, 74945,};
static Numeric ConstTable[50];
static Numeric fact[50];
static Udbyte factable[][50] = {
        {1,  1},
        {1,  1},
        {1,  2},
        {1,  6},
        {1,  24},
        {1,  120},
        {1,  720},
        {1,  5040},
        {1,  40320},
        {2,  62880, 3},
        {2,  28800, 36},
        {2,  16800, 399},
        {2,  1600,  4790},
        {2,  20800, 62270},
        {3,  91200, 71782, 8},
        {3,  68000, 76743, 130},
        {3,  88000, 27898, 2092},
        {3,  96000, 74280, 35568},
        {4,  28000, 37057, 40237, 6},
        {4,  32000, 4088,  64510, 121},
        {4,  40000, 81766, 90200, 2432},
        {4,  40000, 17094, 94217, 51090},
        {5,  80000, 76076, 72777, 24000, 11},
        {5,  40000, 49766, 73888, 52016, 258},
        {5,  60000, 94393, 73323, 48401, 6204},
        {6,  0,     59840, 33098, 10043, 55112, 1},
        {6,  0,     55840, 60563, 61126, 32914, 40},
        {6,  0,     7680,  35216, 50418, 88694, 1088},
        {6,  0,     15040, 86050, 11713, 83446, 30488},
        {7,  0,     36160, 95454, 39701, 19937, 84176, 8},
        {7,  0,     84800, 63630, 91058, 98121, 25285, 265},
        {7,  0,     28800, 72556, 22817, 41779, 83865, 8222},
        {8,  0,     21600, 21801, 30167, 36935, 83693, 63130, 2},
        {8,  0,     12800, 19440, 95518, 18864, 61881, 83317, 86},
        {8,  0,     35200, 60964, 47618, 41408, 3960,  32799, 2952},
        {9,  0,     32000, 33752, 66651, 49296, 38614, 47966, 3331,  1},
        {9,  0,     52000, 15083, 99448, 74679, 90121, 26789, 19933, 37},
        {9,  0,     24000, 58090, 79581, 63159, 34504, 91226, 37530, 1376},
        {9,  0,     12000, 7429,  24100, 72,    11176, 66601, 26174, 52302},
        {10, 0,     68000, 89735, 39902, 2817,  35864, 97443, 20811, 39788, 20},
        {10, 0,     20000, 89427, 96115, 12695, 34561, 97734, 32478, 91528, 815},
        {10, 0,     20000, 66515, 40751, 20534, 17006, 7108,  31638, 52661, 33452},
        {11, 0,     40000, 93638, 11569, 62445, 14260, 98543, 28798, 11775, 5006,  14},
        {11, 0,     20000, 26451, 97507, 85139, 13206, 37355, 38356, 6337,  15263, 604},
        {11, 0,     80000, 63852, 90319, 46158, 81101, 43625, 87680, 78844, 71574, 26582},
        {12, 0,     0,     73376, 64383, 77150, 49565, 63161, 45619, 48019, 20865, 96222, 11},
        {12, 0,     0,     75296, 61651, 48929, 80025, 5428,  98503, 8894,  59812, 26221, 550},
        {12, 0,     0,     38912, 97632, 99691, 61197, 55153, 29643, 18064, 11168, 32415, 25862},
        {13, 0,     0,     67776, 86354, 85214, 37503, 47373, 22890, 67086, 36072, 55925, 41391, 12},
        {13, 0,     0,     21024, 31379, 75528, 37688, 21295, 21633, 87225, 67560, 40342, 28186, 608},
        {13, 0,     0,     51200, 68960, 76415, 84437, 64768, 81660, 61260, 78043, 17133, 9320,  30414},
};

void global_table_init(void) {
    Numeric_new(&zero);
    Numeric_new(&one);
    Numeric_new(&two);
    Numeric_new(&three);
    Numeric_new(&ten);
    Numeric_new(&tenth);
    Numeric_new(&E);
    Numeric_new(&PI);
    Numeric_new(&ln10);
    Numeric_new(&e);
    one.data[0] = 1;
    two.data[0] = 2;
    three.data[0] = 3;
    ten.data[0] = 10;
    tenth.fraction[0] = 10000;
    E.fraction[6] = 1;
    ln10.data[0] = 2;
    PI.data[0] = 3;
    e.data[0] = 2;
    for (int i = 0; i < CALSTACK_FIXED_POINT; i++) {
        PI.fraction[i] = PItable[i];
    }
    for (int i = 0; i < CALSTACK_FIXED_POINT; i++) {
        ln10.fraction[i] = ln10Table[i];
    }
    for (int i = 0; i < CALSTACK_FIXED_POINT; i++) {
        e.fraction[i] = eTable[i];
    }
    for (Udbyte i = 0; i < 50; i++) {
        Numeric_new(&ConstTable[i]);
        ConstTable[i].data[0] = i;
    }
    for (int i = 0; i <= 50; i++) {
        fact[i].nbytes = factable[i][0];
        fact[i].data = &(factable[i][1]);
    }
}


Numeric *Numeric_abs(const Numeric *x) {
    Numeric *result = malloc(sizeof(Numeric));
    Numeric_new(result);
    Numeric_assign(result, x);
    result->sign = POS;
    return result;
}

Numeric *Numeric_sqrt(const Numeric *x) {
    if (Numeric_greater_eq(x, &zero)) {
        /* Estimation */
        Numeric *x0 = malloc(sizeof(Numeric));
        Numeric_new(x0);
        Numeric_assign(x0, x);
        if (x->nbytes) {
            Udbyte n, digit;
            digit = x->data[(x->nbytes) - 1];
            for (n = 0; digit; n++) {
                digit /= 10;
            }
            Numeric *tmp;
            for (int i = 0; i < (n + CALSTACK_UDBYTE_LENGTH * ((x->nbytes) - 1)) / 2; i++) {
                Numeric_assign(x0, tmp = Numeric_mult(x0, &ten));
                Numeric_free_ptr(tmp);
            }
        } else {
            Udbyte n, digit, i, j;
            for (i = 0;; i++) {
                if (x->fraction[i]) {
                    digit = x->fraction[i];
                    break;
                }
            }
            for (n = 0; digit; n++) {
                digit /= 10;
            }
            n = CALSTACK_UDBYTE_LENGTH - n;
            Numeric *tmp;
            for (j = 0; j < (n + i) / 2; j++) {
                Numeric_assign(x0, tmp = Numeric_mult(x0, &ten));
                Numeric_free_ptr(tmp);
            }
        }
        /* Calculation */
        while (Numeric_greater(Numeric_abs(Numeric_sub(Numeric_mult(x0, x0), x)), &E)) {
            Numeric_assign(x0, Numeric_div(Numeric_add(x0, Numeric_div(x, x0)), &two));
        }
        return x0;
    } else {
        fprintf(stderr, "sqrt needs a positive number");
        exit(1);
    }
}

Numeric *Numeric_fact(const Numeric *x) {
    /* Examination */
    if (has_frac(x)) {
        fprintf(stderr, "fact() needs a integer");
        exit(1);
    }
    /* Initialization */
    Numeric *result = malloc(sizeof(Numeric));
    Numeric_new(result);
    Numeric_assign(result, &one);
    /* Calculate the factorial */
    Numeric *i = malloc(sizeof(Numeric));
    Numeric_new(i);
    Numeric *tmp1 = malloc(sizeof(Numeric)), *tmp2;
    Numeric_new(tmp1);
    for (Numeric_assign(i, x); Numeric_greater(i, &one); i = Numeric_sub(tmp1 = i, &one)) {
        result = Numeric_mult(tmp2 = result, i);
        Numeric_free_ptr(tmp1);
        Numeric_free_ptr(tmp2);
    }
    return result;
}

#define SCIENTIFIC_EXP_ITERATE 49

Numeric *Numeric_exp(const Numeric *x) {
    Numeric *x0 = malloc(sizeof(Numeric));
    Numeric_new(x0);
    Numeric_assign(x0, x);
    /* Shrinking */
    if (x->sign == NEG) {
        Numeric_neg(x0);
    }
    Numeric *n = malloc(sizeof(Numeric));
    Numeric_new(n);
    if (x->nbytes) {
        Numeric_assign(n, x0);
        memset(n->fraction, 0, sizeof(Udbyte) * CALSTACK_FIXED_POINT);
        Numeric_assign(x0, Numeric_sub(x0, n));
    }
    /* Taylor */
    Numeric *result = malloc(sizeof(Numeric));
    Numeric_new(result);
    Numeric_assign(result, &one);
    Numeric *power = malloc(sizeof(Numeric));
    Numeric_new(power);
    Numeric_assign(power, x0);
    Numeric *tmp1, *tmp2;
    for (int i = 0; i < SCIENTIFIC_EXP_ITERATE; i++) {
        result = Numeric_add(tmp1 = result, tmp2 = Numeric_div(power, &fact[i + 1]));
        Numeric_free_ptr(tmp1), Numeric_free_ptr(tmp2);
        power = Numeric_mult(tmp1 = power, x0);
        Numeric_free_ptr(tmp1);
    }
    Numeric_free_ptr(power);
    /* Adjustment */
    Numeric_assign(result, Numeric_mult(result, Numeric_pow(&e, n)));
    if (x->sign == NEG) {
        Numeric_assign(result, Numeric_div(&one, result));
        result = Numeric_div(&one, tmp1 = result);
        Numeric_free_ptr(tmp1);
    }
    Numeric_free_ptr(x0);
    return result;
}

#define SCIENTIFIC_LN_ITERATE 50

Numeric *Numeric_ln(const Numeric *x) {
    if (x->sign == POS) {
        /* Scale shrinking */
        Numeric *N = malloc(sizeof(Numeric));
        Numeric_new(N);
        Numeric *x0 = malloc(sizeof(Numeric));
        Numeric_new(x0);
        if (Numeric_less(x, &one)) {
            Numeric_assign(x0, Numeric_div(&one, x));
        } else {
            Numeric_assign(x0, x);
            Udbyte digit, n;
            digit = x->data[(x->nbytes) - 1];
            for (n = 0; digit; n++) {
                digit /= 10;
            }
            Numeric *tmp;
            for (int i = 1; i < n + CALSTACK_UDBYTE_LENGTH * ((x->nbytes) - 1); i++) {
                Numeric_assign(x0, tmp = Numeric_mult(x0, &tenth));
                Numeric_free_ptr(tmp);
                Numeric_assign(N, tmp = Numeric_add(N, &one));
                Numeric_free_ptr(tmp);
            }
        }
        Numeric *tmp;
        for (int i = 0; i < 5; i++) {
            Numeric_assign(x0, tmp = Numeric_sqrt(x0));
            Numeric_free_ptr(tmp);
        }
        /* Taylor */
        Numeric *result = malloc(sizeof(Numeric));
        Numeric_new(result);
        Numeric *x_sub_1 = Numeric_sub(x0, &one);
        Numeric *acc = malloc(sizeof(Numeric));
        Numeric_new(acc);
        Numeric_assign(acc, x_sub_1);
        Numeric *tmp1, *tmp2;
        Numeric *den = malloc(sizeof(Numeric));
        Numeric_new(den);
        den->data[0] = 1;
        for (int64 i = 0; i < SCIENTIFIC_LN_ITERATE; ++i) {
            tmp1 = Numeric_div(acc, den);
            Numeric *(*oper)(const Numeric *, const Numeric *) = i % 2 ? Numeric_sub : Numeric_add;
            result = oper(tmp2 = result, tmp1);
            Numeric_free_ptr(tmp1);
            Numeric_free_ptr(tmp2);
            acc = Numeric_mult(tmp1 = acc, x_sub_1);
            Numeric_free_ptr(tmp1);
            ++den->data[0];
        }
        /* Adjustment */
        for (int i = 0; i < 5; i++) {
            Numeric_assign(result, tmp = Numeric_mult(result, &two));
            Numeric_free_ptr(tmp);
        }
        result = Numeric_add(result, Numeric_mult(N, &ln10));
        if (Numeric_less(x, &one)) {
            Numeric_neg(result);
        }
        return result;
    } else {
        fprintf(stderr, "ln needs a positive number");
        exit(1);
    }
}

Numeric *Numeric_pow(const Numeric *base, const Numeric *power) {
    Numeric *result = malloc(sizeof(Numeric));
    Numeric_new(result);
    /* Change base */
    if (has_frac(power)) {
        Numeric *NewPower = Numeric_ln(base);
        NewPower = Numeric_mult(power, NewPower);
        result = Numeric_exp(NewPower);
        Numeric_free_ptr(NewPower);
    }
        /* Fast Modular Exponentiation */
    else {
        /* Initialization */
        Numeric_assign(result, &one);
        Numeric *base_cpy = malloc(sizeof(Numeric));
        Numeric_new(base_cpy);
        Numeric_assign(base_cpy, base);
        /* Power calculation */
        Numeric *power_cpy = malloc(sizeof(Numeric));
        Numeric_new(power_cpy);
        Numeric_assign(power_cpy, power);
        Numeric *tmp;
        while (!Numeric_is_zero(power_cpy)) {
            if (Numeric_mod2(*power_cpy)) {
                result = Numeric_mult(tmp = result, base_cpy);
                Numeric_free_ptr(tmp);
            }
            power_cpy = Numeric_div_int(tmp = power_cpy, &two);
            Numeric_free_ptr(tmp);
            base_cpy = Numeric_mult(tmp = base_cpy, base_cpy);
            Numeric_free_ptr(tmp);
        }
    }
    return result;
}

Numeric *Numeric_log(const Numeric *base, const Numeric *oper) {
    /* Change base */
    Numeric *tmp1, *tmp2;
    Numeric *result = Numeric_div(tmp1 = Numeric_ln(oper), tmp2 = Numeric_ln(base));
    Numeric_free_ptr(tmp1), Numeric_free_ptr(tmp2);
    return result;
}

Numeric *Numeric_lg(const Numeric *x) {
    /* Change base */
    Numeric *tmp1, *tmp2;
    Numeric *result = Numeric_div(tmp1 = Numeric_ln(x), tmp2 = Numeric_ln(&ten));
    Numeric_free_ptr(tmp1), Numeric_free_ptr(tmp2);
    return result;
}

Numeric *Numeric_log2(const Numeric *x) {
    /* Change base */
    Numeric *tmp1, *tmp2;
    Numeric *result = Numeric_div(tmp1 = Numeric_ln(x), tmp2 = Numeric_ln(&two));
    Numeric_free_ptr(tmp1), Numeric_free_ptr(tmp2);
    return result;
}

Numeric *Numeric_sinh(const Numeric *x) {
    /* Initialization */
    Numeric *result = malloc(sizeof(Numeric));
    Numeric_new(result);
    Numeric *neg = malloc(sizeof(Numeric));
    Numeric_new(neg);
    Numeric_assign(neg, x);
    Numeric_neg(neg);
    /* Calculation */
    Numeric_assign(result, Numeric_div(Numeric_sub(Numeric_exp(x), Numeric_exp(neg)), &two));
    Numeric_free_ptr(neg);
    return result;
}

Numeric *Numeric_cosh(const Numeric *x) {
    /* Initialization */
    Numeric *result = malloc(sizeof(Numeric));
    Numeric_new(result);
    Numeric *neg = malloc(sizeof(Numeric));
    Numeric_new(neg);
    Numeric_assign(neg, x);
    Numeric_neg(neg);
    /* Calculation */
    Numeric_assign(result, Numeric_div(Numeric_add(Numeric_exp(x), Numeric_exp(neg)), &two));
    Numeric_free_ptr(neg);

    return result;
}

Numeric *Numeric_tanh(const Numeric *x) {
    /* Initialization */
    Numeric *result = malloc(sizeof(Numeric));
    Numeric_new(result);
    Numeric *neg = malloc(sizeof(Numeric));
    Numeric_new(neg);
    Numeric_assign(neg, x);
    Numeric_neg(neg);
    /* Calculation */
    Numeric_assign(result, Numeric_div(Numeric_sub(Numeric_exp(x), Numeric_exp(neg)),
                                       Numeric_add(Numeric_exp(x), Numeric_exp(neg))));
    Numeric_free_ptr(neg);
    return result;
}

Numeric *Numeric_sech(const Numeric *x) {
    /* Initialization */
    Numeric *result = malloc(sizeof(Numeric));
    Numeric_new(result);
    Numeric *neg = malloc(sizeof(Numeric));
    Numeric_new(neg);
    Numeric_assign(neg, x);
    Numeric_neg(neg);
    /* Calculation */
    Numeric_assign(result, Numeric_div(&two, Numeric_add(Numeric_exp(x), Numeric_exp(neg))));
    Numeric_free_ptr(neg);
    return result;
}

Numeric *Numeric_csch(const Numeric *x) {
    /* Initialization */
    Numeric *result = malloc(sizeof(Numeric));
    Numeric_new(result);
    Numeric *neg = malloc(sizeof(Numeric));
    Numeric_new(neg);
    Numeric_assign(neg, x);
    Numeric_neg(neg);
    /* Calculation */
    Numeric_assign(result, Numeric_div(&two, Numeric_sub(Numeric_exp(x), Numeric_exp(neg))));
    Numeric_free_ptr(neg);
    return result;
}

Numeric *Numeric_coth(const Numeric *x) {
    /* Initialization */
    Numeric *result = malloc(sizeof(Numeric));
    Numeric_new(result);
    Numeric *neg = malloc(sizeof(Numeric));
    Numeric_new(neg);
    Numeric_assign(neg, x);
    Numeric_neg(neg);
    /* Calculation */
    Numeric_assign(result, Numeric_div(Numeric_add(Numeric_exp(x), Numeric_exp(neg)),
                                       Numeric_sub(Numeric_exp(x), Numeric_exp(neg))));
    Numeric_free_ptr(neg);
    return result;
}

#define SCIENTIFIC_SIN_ITERATE 25

Numeric *Numeric_sin(const Numeric *x) {
    /* Initialization */
    Numeric *x0 = malloc(sizeof(Numeric));
    Numeric_new(x0);
    Numeric_assign(x0, x);
    /* Scale shrinking */
    Numeric *n = malloc(sizeof(Numeric));
    Numeric_new(n);
    Numeric_assign(x0, Numeric_div(x0, &PI));
    Numeric_assign(n, x0);
    memset(n->fraction, 0, sizeof(Udbyte) * CALSTACK_FIXED_POINT);
    Numeric_assign(x0, Numeric_sub(x0, n));
    if (Numeric_greater(Numeric_mult(x0, &two), &one)) {
        Numeric_assign(x0, Numeric_sub(&one, x0));
    }
    Numeric_assign(x0, Numeric_mult(x0, &PI));
    /* Taylor */
    Numeric *tmp1, *tmp2 = malloc(sizeof(Numeric));
    Numeric_new(tmp2);
    Numeric *result = malloc(sizeof(Numeric));
    Numeric_new(result);
    Numeric *pw = malloc(sizeof(Numeric)), *pw2 = malloc(sizeof(Numeric));
    Numeric_new(pw), Numeric_new(pw2);
    Numeric_assign(pw, x0);
    Numeric_move(pw2, tmp1 = Numeric_mult(x0, x0));
    Numeric_free_ptr(tmp1);
    for (int64 i = 0; i < SCIENTIFIC_SIN_ITERATE; i++) {
        Numeric_free_ptr(tmp2);
        tmp1 = Numeric_div(pw, &fact[i * 2 + 1]);
        Numeric *(*oper)(const Numeric *, const Numeric *) = i % 2 ? Numeric_sub : Numeric_add;
        result = oper(tmp2 = result, tmp1);
        Numeric_free_ptr(tmp1);
        pw = Numeric_mult(tmp1 = pw, pw2);
        Numeric_free(tmp1);
    }
    /* Adjustment */
    if (Numeric_mod2(*n)) {
        Numeric_neg(result);
    }
    return result;
}

#define SCIENTIFIC_COS_ITERATE 25

Numeric *Numeric_cos(const Numeric *x) {
    /* Initialization */
    Numeric *x0 = malloc(sizeof(Numeric));
    Numeric_new(x0);
    Numeric_assign(x0, x);
    /* Scale shrinking */
    Numeric *n = malloc(sizeof(Numeric));
    Numeric_new(n);
    Numeric_assign(x0, Numeric_div(x0, &PI));
    Numeric_assign(n, x0);
    memset(n->fraction, 0, sizeof(Udbyte) * CALSTACK_FIXED_POINT);
    Numeric_assign(x0, Numeric_sub(x0, n));
    bool isNeg = 0;
    if (Numeric_greater(Numeric_mult(x0, &two), &one)) {
        Numeric_assign(x0, Numeric_sub(&one, x0));
        isNeg = 1;
    }
    Numeric_assign(x0, Numeric_mult(x0, &PI));
    /* Taylor */
    Numeric *tmp1, *tmp2 = malloc(sizeof(Numeric));
    Numeric_new(tmp2);
    Numeric *result = malloc(sizeof(Numeric));
    Numeric_new(result);
    Numeric *pw = malloc(sizeof(Numeric)), *pw2 = malloc(sizeof(Numeric));
    Numeric_new(pw), Numeric_new(pw2);
    Numeric_assign(pw, &one);
    Numeric_move(pw2, tmp1 = Numeric_mult(x0, x0));
    Numeric_free_ptr(tmp1);
    for (int64 i = 0; i < SCIENTIFIC_COS_ITERATE; i++) {
        Numeric_free_ptr(tmp2);
        tmp1 = Numeric_div(pw, &fact[i * 2]);
        Numeric *(*oper)(const Numeric *, const Numeric *) = i % 2 ? Numeric_sub : Numeric_add;
        result = oper(tmp2 = result, tmp1);
        Numeric_free_ptr(tmp1);
        pw = Numeric_mult(tmp1 = pw, pw2);
        Numeric_free(tmp1);
    }
    /* Adjustment */
    if (Numeric_mod2(*n)) {
        Numeric_neg(result);
    }
    if (isNeg) {
        Numeric_neg(result);
    }
    return result;
}

Numeric *Numeric_tan(const Numeric *x) {
    /* tan(x) = sin(x) / cos(x) */
    Numeric *tmp1, *tmp2;
    Numeric *result = Numeric_div(tmp2 = Numeric_sin(x), tmp1 = Numeric_cos(x));
    Numeric_free_ptr(tmp1), Numeric_free_ptr(tmp2);
    return result;
}

Numeric *Numeric_sec(const Numeric *x) {
    /* sec(x) = 1 / cos(x) */
    Numeric *result = malloc(sizeof(Numeric));
    Numeric_new(result);
    Numeric *tmp;
    Numeric_assign(result, tmp = Numeric_div(&one, Numeric_cos(x)));
    Numeric_free_ptr(tmp);
    return result;
}

Numeric *Numeric_csc(const Numeric *x) {
    /* csc(x) = 1 / sin(x) */
    Numeric *result = malloc(sizeof(Numeric));
    Numeric_new(result);
    Numeric *tmp;
    Numeric_assign(result, tmp = Numeric_div(&one, Numeric_sin(x)));
    Numeric_free_ptr(tmp);
    return result;
}

Numeric *Numeric_cot(const Numeric *x) {
    /* cot(x) = cos(x) / sin(x) */
    Numeric *tmp1, *tmp2;
    Numeric *result = Numeric_div(tmp2 = Numeric_cos(x), tmp1 = Numeric_sin(x));
    Numeric_free_ptr(tmp1), Numeric_free_ptr(tmp2);
    return result;
}

Numeric *Numeric_arcsin(const Numeric *x) {
    if (Numeric_less(Numeric_abs(x), &one)) {
        /* Initialization */

        Numeric *result = malloc(sizeof(Numeric));
        Numeric_new(result);
        Numeric *y = Numeric_abs(x);
        Numeric *angle = malloc(sizeof(Numeric));
        Numeric_new(angle);
        Numeric_assign(angle, &one);
        Numeric *x0 = malloc(sizeof(Numeric));
        Numeric_new(x0);
        Numeric_assign(x0, &one);
        Numeric *y0 = malloc(sizeof(Numeric));
        Numeric_new(y0);
        Numeric *x1 = malloc(sizeof(Numeric));
        Numeric_new(x1);
        Numeric *y1 = malloc(sizeof(Numeric));
        Numeric_new(y1);
        /* Calculation */
        while (Numeric_greater(Numeric_abs(Numeric_sub(y1, y)), &E)) {
            Numeric *tmp1, *tmp2, *tmp3, *tmp4, *tmp5;
            if (Numeric_less(y1, y)) {
                Numeric_assign(x1, tmp5 = Numeric_sub(tmp3 = Numeric_mult(x0, tmp1 = Numeric_cos(angle)),
                                                      tmp4 = Numeric_mult(y0, tmp2 = Numeric_sin(angle))));
                Numeric_free_ptr(tmp1), Numeric_free_ptr(tmp2), Numeric_free_ptr(tmp3), Numeric_free_ptr(
                        tmp4), Numeric_free_ptr(tmp5);
                Numeric_assign(y1, tmp5 = Numeric_add(tmp3 = Numeric_mult(x0, tmp1 = Numeric_sin(angle)),
                                                      tmp4 = Numeric_mult(y0, tmp2 = Numeric_cos(angle))));
                Numeric_free_ptr(tmp1), Numeric_free_ptr(tmp2), Numeric_free_ptr(tmp3), Numeric_free_ptr(
                        tmp4), Numeric_free_ptr(tmp5);
                Numeric_assign(x0, x1);
                Numeric_assign(y0, y1);
                Numeric_assign(result, Numeric_add(result, angle));
            } else {
                Numeric_assign(x1, tmp5 = Numeric_add(tmp3 = Numeric_mult(x0, tmp1 = Numeric_cos(angle)),
                                                      tmp4 = Numeric_mult(y0, tmp2 = Numeric_sin(angle))));
                Numeric_free_ptr(tmp1), Numeric_free_ptr(tmp2), Numeric_free_ptr(tmp3), Numeric_free_ptr(
                        tmp4), Numeric_free_ptr(tmp5);
                Numeric_assign(y1, tmp5 = Numeric_sub(tmp3 = Numeric_mult(y0, tmp1 = Numeric_cos(angle)),
                                                      tmp4 = Numeric_mult(x0, tmp2 = Numeric_sin(angle))));
                Numeric_free_ptr(tmp1), Numeric_free_ptr(tmp2), Numeric_free_ptr(tmp3), Numeric_free_ptr(
                        tmp4), Numeric_free_ptr(tmp5);
                Numeric_assign(x0, x1);
                Numeric_assign(y0, y1);
                Numeric_assign(result, Numeric_sub(result, angle));
            }

            Numeric_assign(angle, Numeric_div(angle, &two));
        }
        /* Adjustment */
        Numeric_free_ptr(y);
        Numeric_free_ptr(x0);
        Numeric_free_ptr(y0);
        Numeric_free_ptr(x1);
        Numeric_free_ptr(y1);
        Numeric_free_ptr(angle);
        if (x->sign == NEG) {
            Numeric_neg(result);
        }
        return result;
    } else if (Numeric_eq(Numeric_abs(x), &one)) {
        Numeric *result = Numeric_div(&PI, &two);
        if (x->sign == NEG) {
            Numeric_neg(result);
        }
        return result;
    } else {
        fprintf(stderr, "arcsin(x) error");
        exit(1);
    }
}

Numeric *Numeric_arccos(const Numeric *x) {
    if (Numeric_less_eq(Numeric_abs(x), &one)) {
        Numeric *result = Numeric_sub(Numeric_div(&PI, &two), Numeric_arcsin(x));
        return result;
    } else {
        fprintf(stderr, "arcsin(x) error");
        exit(1);
    }
}

Numeric *Numeric_arctan(const Numeric *x) {
    /* Initialization */
    Numeric *result = malloc(sizeof(Numeric));
    Numeric_new(result);
    Numeric *n = Numeric_abs(x);
    Numeric *angle = malloc(sizeof(Numeric));
    Numeric_new(angle);
    Numeric_assign(angle, &one);
    Numeric *x0 = Numeric_sqrt(Numeric_div(&one, Numeric_add(Numeric_mult(n, n), &one)));
    Numeric *y0 = Numeric_sqrt(Numeric_sub(&one, Numeric_mult(x0, x0)));
    Numeric *x1 = malloc(sizeof(Numeric));
    Numeric_new(x1);
    Numeric_assign(x1, x0);
    Numeric *y1 = malloc(sizeof(Numeric));
    Numeric_new(y1);
    Numeric_assign(y1, y0);
    /* Calculation */
    while (Numeric_greater(Numeric_abs(y1), &E)) {
        Numeric *tmp1, *tmp2, *tmp3, *tmp4, *tmp5;
        if (Numeric_less(y1, &zero)) {
            Numeric_assign(x1, tmp5 = Numeric_sub(tmp3 = Numeric_mult(x0, tmp1 = Numeric_cos(angle)),
                                                  tmp4 = Numeric_mult(y0, tmp2 = Numeric_sin(angle))));
            Numeric_free_ptr(tmp1), Numeric_free_ptr(tmp2), Numeric_free_ptr(tmp3), Numeric_free_ptr(
                    tmp4), Numeric_free_ptr(tmp5);
            Numeric_assign(y1, tmp5 = Numeric_add(tmp3 = Numeric_mult(x0, tmp1 = Numeric_sin(angle)),
                                                  tmp4 = Numeric_mult(y0, tmp2 = Numeric_cos(angle))));
            Numeric_free_ptr(tmp1), Numeric_free_ptr(tmp2), Numeric_free_ptr(tmp3), Numeric_free_ptr(
                    tmp4), Numeric_free_ptr(tmp5);
            Numeric_assign(x0, x1);
            Numeric_assign(y0, y1);
            Numeric_assign(result, Numeric_sub(result, angle));
        } else {
            Numeric_assign(x1, tmp5 = Numeric_add(tmp3 = Numeric_mult(x0, tmp1 = Numeric_cos(angle)),
                                                  tmp4 = Numeric_mult(y0, tmp2 = Numeric_sin(angle))));
            Numeric_free_ptr(tmp1), Numeric_free_ptr(tmp2), Numeric_free_ptr(tmp3), Numeric_free_ptr(
                    tmp4), Numeric_free_ptr(tmp5);
            Numeric_assign(y1, tmp5 = Numeric_sub(tmp3 = Numeric_mult(y0, tmp1 = Numeric_cos(angle)),
                                                  tmp4 = Numeric_mult(x0, tmp2 = Numeric_sin(angle))));
            Numeric_free_ptr(tmp1), Numeric_free_ptr(tmp2), Numeric_free_ptr(tmp3), Numeric_free_ptr(
                    tmp4), Numeric_free_ptr(tmp5);
            Numeric_assign(x0, x1);
            Numeric_assign(y0, y1);
            Numeric_assign(result, Numeric_add(result, angle));
        }
        Numeric_assign(angle, Numeric_div(angle, &two));
    }
    /* Adjustment */
    Numeric_free_ptr(n);
    Numeric_free_ptr(x0);
    Numeric_free_ptr(y0);
    Numeric_free_ptr(x1);
    Numeric_free_ptr(y1);
    Numeric_free_ptr(angle);
    if (x->sign == NEG) {
        Numeric_neg(result);
    }
    return result;
}

Numeric *Numeric_arcsec(const Numeric *x) {
    Numeric *tmp;
    Numeric *result = Numeric_arccos(tmp = Numeric_div(&one, x));
    Numeric_free_ptr(tmp);
    return result;
}

Numeric *Numeric_arccsc(const Numeric *x) {
    Numeric *tmp;
    Numeric *result = Numeric_arcsin(tmp = Numeric_div(&one, x));
    Numeric_free_ptr(tmp);
    return result;
}

Numeric *Numeric_arccot(const Numeric *x) {
    Numeric *tmp;
    Numeric *result = Numeric_arctan(tmp = Numeric_div(&one, x));
    Numeric_free_ptr(tmp);
    return result;
}

