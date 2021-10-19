// Adapted from
// http://mumble.net/~campbell/tmp/random_real.c

#include "pcg_variants.h"
#include "stdint.h"

/*
 * random_real: Generate a stream of bits uniformly at random and
 * interpret it as the fractional part of the binary expansion of a
 * number in [0, 1], 0.00001010011111010100...; then round it.
 */
double
random_real01(pcg64_random_t* rng) {
    /*
     * Read zeros into the exponent until we hit a one; the rest
     * will go into the significand.
     */
    int exponent = -64;
    uint64_t significand;
    while ((significand = pcg64_random_r(rng)) == 0) {
        exponent -= 64;
        /*
         * If the exponent falls below -1074 = emin + 1 - p,
         * the exponent of the smallest subnormal, we are
         * guaranteed the result will be rounded to zero.  This
         * case is so unlikely it will happen in realistic
         * terms only if random64 is broken.
         */
        if (exponent < -1074)
            return 0;
    }

    /*
     * There is a 1 somewhere in significand, not necessarily in
     * the most significant position.  If there are leading zeros,
     * shift them into the exponent and refill the less-significant
     * bits of the significand.  Can't predict one way or another
     * whether there are leading zeros: there's a fifty-fifty
     * chance, if random64 is uniformly distributed.
     */
    int shift = __builtin_clzll(significand);
    if (shift != 0) {
        exponent -= shift;
        significand <<= shift;
        significand |= (pcg64_random_r(rng) >> (64 - shift));
    }

    /*
     * Set the sticky bit, since there is almost surely another 1
     * in the bit stream.  Otherwise, we might round what looks
     * like a tie to even when, almost surely, were we to look
     * further in the bit stream, there would be a 1 breaking the
     * tie.
     */
    significand |= 1;

    /*
     * Finally, convert to double (rounding) and scale by
     * 2^exponent.
     */
    uint64_t pow2_exponent_bits = ((uint64_t)(exponent + 1023)) << 52;
    double pow2_exponent = *(double*)&pow2_exponent_bits;
    double result = (double)significand * pow2_exponent;
    return result;
}
