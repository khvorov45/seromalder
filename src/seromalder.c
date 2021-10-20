#include <stdint.h>

#include <immintrin.h>

#include "random_real.h"

typedef void* SmlAllocator(uint64_t size);

typedef struct SmlInputTitre {
    double log2titre;
    double time;
} SmlInputTitre;

typedef enum SmlEventType {
    SmlEvent_Vaccination,
    SmlEvent_Infection,
} SmlEventType;

typedef struct SmlInputEvent {
    SmlEventType type;
    double time;
} SmlInputEvent;

typedef struct SmlInputIndividual {
    uint64_t event_count;
    SmlInputEvent* events;
    uint64_t titre_count;
    SmlInputTitre* titres;
} SmlInputIndividual;

typedef struct SmlInput {
    uint64_t n_individuals;
    SmlInputIndividual* data;
} SmlInput;

typedef struct SmlParameters {
    double vaccination_log2diff;
    double baseline;
    double baseline_sd;
    double wane_rate;
    double residual_sd;
} SmlParameters;

typedef struct SmlOutput {
    uint64_t n_iterations;
    SmlParameters* out;
} SmlOutput;

typedef struct SmlConstants {
    double time_to_peak;
    double lowest_log2titre;
} SmlConstants;

typedef struct SmlMcmcSettings {
    SmlParameters proposal_sds;
} SmlMcmcSettings;

SmlConstants
sml_default_constants() {
    SmlConstants result = {
        .time_to_peak = 14, // NOTE(sen) Days
        .lowest_log2titre = 2.321928, // NOTE(sen) Log2(5)
    };
    return result;
}

SmlMcmcSettings
sml_default_settings() {
    SmlMcmcSettings result;
    result.proposal_sds.vaccination_log2diff = 1;
    result.proposal_sds.baseline = 1;
    result.proposal_sds.baseline_sd = 1;
    result.proposal_sds.wane_rate = 1;
    result.proposal_sds.residual_sd = 1;
    return result;
}

double
sml_pow2(double value) {
    // NOTE(sen) Adapted from
    // https://github.com/pmttavara/pt_math

    double exponent_biased = value + 1023;

    int64_t exponent_biased_floored = (int64_t)(exponent_biased);
    double exponent_biased_frac = exponent_biased - exponent_biased_floored;

    int64_t bits_in_mantissa = 52;
    int64_t closest_pow2_bits = exponent_biased_floored << bits_in_mantissa;
    double closest_pow2 = *(double*)&closest_pow2_bits;

    double two_to_frac = 1 + exponent_biased_frac *
        (0.6931471805599452862268 +
            0.2402265069591006940719 * exponent_biased_frac +
            0.05550410866482157618007 * exponent_biased_frac * exponent_biased_frac);
    double result = two_to_frac * closest_pow2;
    return result;
}

double
sml_log2(double value) {
    // NOTE(sen) Adapted from
    // https://github.com/pmttavara/pt_math

    int64_t value_bits = *(int64_t*)&value;

    int64_t bits_in_mantissa = 52;
    int64_t value_exponent_biased = value_bits >> bits_in_mantissa;

    int64_t value_remainder_bits = (value_bits & 0x000fffffffffffff) | 0x3ff0000000000000;
    double value_remainder = *(double*)&value_remainder_bits;

    double result = (double)(value_exponent_biased)
        -1024 + value_remainder * (value_remainder * -0.33333333333333333333 + 2) - 0.666666666666666666;

    return result;
}

double
sml_log2_std_normal_pdf(double value) {
    double log2_one_over_sqrt_2pi = -1.325748064736159248511;
    double half_log2e = 0.7213475204444816935023;
    double result = log2_one_over_sqrt_2pi - half_log2e * value * value;
    return result;
}

double
sml_sin(double value) {
    // NOTE(sen) Adapted from
    // https://github.com/pmttavara/pt_math

    double one_over_two_pi = 0.1591549430918953456082221009637578390538692474365234375;
    double value_scaled_to_period01 = value * one_over_two_pi;

    double value_scaled_to_period01_rounded = value_scaled_to_period01;
    *(volatile double*)&value_scaled_to_period01_rounded += 6755399441055744.0;
    *(volatile double*)&value_scaled_to_period01_rounded -= 6755399441055744.0;

    double value_neg1_1 = value_scaled_to_period01 - value_scaled_to_period01_rounded;
    double value_neg1_1_abs = value_neg1_1 < 0 ? -value_neg1_1 : value_neg1_1;

    double parabola_div16 = value_neg1_1 * (0.5 - value_neg1_1_abs);
    double parabola_div16_abs = parabola_div16 < 0 ? -parabola_div16 : parabola_div16;

    double result = parabola_div16 * (57.3460872862336 * parabola_div16_abs + 12.4158695446104);

    return result;
}

double
sml_rnorm01(pcg64_random_t* rng) {
    // NOTE(sen) Box-Muller

    double u1 = random_real01(rng);
    double u2 = random_real01(rng);

    double log2_u1 = sml_log2(u1);
    double one_over_log2e = 0.69314718055994528622676398299518041312694549560546875;
    double ln_u1 = log2_u1 * one_over_log2e;
    __m128d neg_ln_u1_w = _mm_set_sd(-ln_u1);
    __m128d sqrt_neg_ln_u1_w = _mm_sqrt_pd(neg_ln_u1_w);
    double sqrt_neg_ln_u1 = *(double*)&sqrt_neg_ln_u1_w;
    double sqrt_2 = 1.4142135623730951454746218587388284504413604736328125;
    double out_a = sqrt_2 * sqrt_neg_ln_u1;

    double two_pi = 6.28318530717958623199592693708837032318115234375;
    double out_b = two_pi * u2;

    double sin_out_b = sml_sin(out_b);

    double result1 = out_a * sin_out_b;

    // NOTE(sen) Can generate a second one here
    //double half_pi = 1.5707963267948965579989817342720925807952880859375;
    //double cos_out_b = sml_sin(out_b + half_pi);
    //double result2 = out_a * cos_out_b;

    return result1;
}

int32_t
sml_rbern(pcg64_random_t* rng, double prop) {
    double rand01 = random_real01(rng);
    double result = rand01 <= prop;
    return result;
}

double
sml_log2_prior_prob(SmlParameters* pars) {
    // TODO(sen) Implement
    return 0;
}

double
sml_log2_likelihood(SmlInput* input, SmlParameters* pars, SmlConstants* consts) {

    double log2_likelihood = 0;

    for (uint64_t individual_index = 0;
        individual_index < input->n_individuals;
        individual_index++) {

        SmlInputIndividual* individual = input->data + individual_index;

        for (uint64_t titre_index = 0;
            titre_index < individual->titre_count;
            titre_index++) {

            SmlInputTitre* titre = individual->titres + titre_index;

            double predicted_titre = consts->lowest_log2titre;

            for (uint64_t event_index = 0;
                event_index < individual->event_count;
                event_index++) {

                SmlInputEvent* event = individual->events + event_index;

                if (titre->time > event->time) {
                    if (event->type == SmlEvent_Vaccination) {
                        double time_since = titre->time - event->time;
                        double up_slope = pars->vaccination_log2diff / consts->time_to_peak;
                        double past_peak = (double)(time_since > consts->time_to_peak);
                        double titre_contribution = up_slope * time_since -
                            past_peak * (up_slope + pars->wane_rate) * (time_since - consts->time_to_peak);
                        predicted_titre += titre_contribution;
                    } else {
                        // TODO(sen) Handle infection
                    }
                }
            } // NOTE(sen) for (event)

            double deviation = titre->log2titre - predicted_titre;
            double titre_prob = sml_log2_std_normal_pdf((deviation - predicted_titre) / pars->residual_sd);
            log2_likelihood += titre_prob;

        } // NOTE(sen) for (titre)

    } // NOTE(sen) for (individual)

    return log2_likelihood;
}

void
sml_mcmc(
    SmlInput* input,
    SmlParameters* pars_init,
    SmlOutput* output,
    SmlConstants* consts,
    SmlMcmcSettings* settings
) {
    SmlParameters pars_cur = *pars_init;
    double log2_prior_prob_cur = sml_log2_prior_prob(&pars_cur);
    double log2_likelihood_cur = sml_log2_likelihood(input, &pars_cur, consts);
    double log2_posterior_cur = log2_prior_prob_cur + log2_likelihood_cur;

    SmlParameters* steps = &settings->proposal_sds;

    pcg64_random_t rng;
    pcg_setseq_128_srandom_r(&rng, 0, 0);

    for (uint64_t iteration = 0; iteration < output->n_iterations; iteration++) {

        SmlParameters pars_next;

        pars_next.vaccination_log2diff =
            sml_rnorm01(&rng) * steps->vaccination_log2diff + pars_cur.vaccination_log2diff;
        pars_next.baseline = sml_rnorm01(&rng) * steps->baseline + pars_cur.baseline;
        pars_next.baseline_sd = sml_rnorm01(&rng) * steps->baseline_sd + pars_cur.baseline_sd;
        pars_next.wane_rate = sml_rnorm01(&rng) * steps->wane_rate + pars_cur.wane_rate;

        double log2_prior_prob_next = sml_log2_prior_prob(&pars_next);
        double log2_likelihood_next = sml_log2_likelihood(input, &pars_next, consts);
        double log2_posterior_next = log2_prior_prob_next + log2_likelihood_next;

        double log2_posterior_diff = log2_posterior_next - log2_posterior_cur;

        if (log2_posterior_diff >= 0) {
            pars_cur = pars_next;
        } else if (log2_posterior_diff > -20) {
            double posterior_ratio = sml_pow2(log2_posterior_diff);
            if (sml_rbern(&rng, posterior_ratio)) {
                pars_cur = pars_next;
            }
        }

        output->out[iteration] = pars_cur;
    } // NOTE(sen) for (iteration)
}
