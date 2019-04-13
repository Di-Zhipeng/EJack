/* mult_fft_helper.c
 * Created by Wang Haoren
 * --------------------------------
 * Implements the Fast Fourier Transformation
 * */

#include <math.h>
#include <string.h>
#include "numeric.h"
#include "mult_fft_helper.h"

#ifndef CMPLX // The macro SHOULD be a part of the C Standard Library but not all compilers implement it
#define CMPLX(x, y) ((long double complex)((long double)(x) + I * (long double)(y)))
#endif


static long double complex cap1[CALSTACK_NBYTES_MAX], cap2[CALSTACK_NBYTES_MAX];
static size_t last_used_length;

size_t align_p2(const Udbyte *seq1, const Udbyte *seq2, size_t len1, size_t len2) {
    size_t roof = 1;
    size_t twolen1 = len1 << 1, twolen2 = len2 << 1;
    while (roof < twolen1 || roof < twolen2) {
        roof <<= 1;
    }
    for (size_t i = 0; i < len1; ++i) {
        cap1[i] = CMPLX(seq1[i], 0);
    }
    for (size_t i = 0; i < len2; ++i) {
        cap2[i] = CMPLX(seq2[i], 0);
    }
    return last_used_length = roof;
}

void bit_rev_shuffle(long double complex *val, size_t len, unsigned int64 *table) {
    for (size_t i = 0; i < len; ++i) {
        if (i < table[i]) {
            long double complex temp = val[i];
            val[i] = val[table[i]];
            val[table[i]] = temp;
        }
    }
}

void bitrevtb(unsigned int64 *table, size_t size) {
    table[0] = 0, table[1] = size >> 1;
    for (int i = 2; i < size; i += 2) {
        table[i] = table[i >> 1] >> 1;
        table[i + 1] = table[i] ^ table[1];
    }
}


void fourier_transform_calc(long double complex *sample, size_t len, bool inv, unsigned int64 *table) {
    bit_rev_shuffle(sample, len, table);
    for (int i = 2; i <= len; i <<= 1) {
        long double conv = ((inv) ? 1 : -1) * 2 * M_PI / i;
        long double complex wn = CMPLX(cosl(conv), sinl(conv));
        for (size_t j = 0; j < len; j += i) {
            long double complex w = CMPLX(1, 0);
            for (size_t k = j; k < j + (i >> 1); ++k) {
                long double complex u = sample[k], t = w * sample[k + (i >> 1)];
                sample[k] = u + t;
                sample[k + (i >> 1)] = u - t;
                w *= wn;
            }
        }
    }
    if (inv) {
        for (int i = 0; i < len; ++i) {
            __real__(sample[i]) /= len;
        }
    }
}

void fourier_transform(size_t len, unsigned int64 *table) {
    fourier_transform_calc(cap1, len, false, table);
    fourier_transform_calc(cap2, len, false, table);
    for (size_t i = 0; i < len; ++i) {
        cap1[i] *= cap2[i];
    }
    fourier_transform_calc(cap1, len, true, table);
}

Uldbyte result_copy(Uldbyte *tar, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        tar[i] = (Uldbyte) llroundl(__real__(cap1[i]));
    }
    Uldbyte carry = 0;
    for (size_t i = 0; i < len; ++i) {
        tar[i] += carry;
        carry = tar[i] / CALSTACK_UDBYTE_MAX;
        tar[i] %= CALSTACK_UDBYTE_MAX;
    }
    memset(cap1, 0, last_used_length * sizeof(long double complex));   // Reinitialize the fft arrays
    memset(cap2, 0, last_used_length * sizeof(long double complex));
    return carry;
}
