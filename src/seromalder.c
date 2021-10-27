#include <stdint.h>
#include <stddef.h>

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

typedef union SmlParameters {
    // NOTE(sen) Make sure order is the same as SmlPriors
    struct {
        double baseline;
        double residual_sd;
        double vaccination_log2diff;
        double wane_rate;
        double baseline_sd;
    };
    double par[5];
} SmlParameters;

typedef struct SmlOutput {
    uint64_t n_iterations;
    uint64_t n_burn;
    uint64_t n_accepted_after_burn;
    uint64_t n_accepted_burn;
    SmlParameters* out;
} SmlOutput;

typedef struct SmlConstants {
    double time_to_peak;
    double lowest_log2titre;
} SmlConstants;

typedef enum SmlDistType {
    SmlDist_Normal,
    SmlDist_NormalLeftTrunc,
} SmlDistType;

typedef struct SmlDistNormal {
    double mean;
    double sd;
} SmlDistNormal;

typedef struct SmlDistNormalLeftTrunc {
    double mean;
    double sd;
    double min;
} SmlDistNormalLeftTrunc;

typedef struct SmlDist {
    SmlDistType type;
    union {
        SmlDistNormal normal;
        SmlDistNormalLeftTrunc normal_left_trunc;
    };
} SmlDist;

typedef union SmlPriors {
    // NOTE(sen) Make sure order is the same as SmlParameters
    struct {
        SmlDist baseline;
        SmlDist residual_sd;
        SmlDist vaccination_log2diff;
        SmlDist wane_rate;
        SmlDist baseline_sd;
    };
    SmlDist dist[5];
} SmlPriors;

typedef struct SmlStep {
    uint32_t dim;
    double* mean;
    double* var;
    double* chol;
    double* std_norm;
} SmlStep;

SmlConstants
sml_default_constants() {
    SmlConstants result = {
        .time_to_peak = 14, // NOTE(sen) Days
        .lowest_log2titre = 2.321928, // NOTE(sen) Log2(5)
    };
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
sml_log2_normal_pdf(double value, double mean, double sd) {
    double log2_one_over_sqrt_2pi = -1.325748064736159248511;
    double half_log2e = 0.7213475204444816935023;
    double log2_sd = sml_log2(sd);
    double value_sd = (value - mean) / sd;
    double result = log2_one_over_sqrt_2pi - log2_sd - half_log2e * value_sd * value_sd;
    return result;
}

double
sml_log2_dist_pdf(double value, SmlDist dist) {
    double result;
    switch (dist.type) {
    case SmlDist_Normal: {
        result = sml_log2_normal_pdf(value, dist.normal.sd, dist.normal.mean);
    } break;
    case SmlDist_NormalLeftTrunc: {
        if (value < dist.normal_left_trunc.min) {
            uint64_t bits = (uint64_t)0xFFF << 52;
            result = *(double*)&bits;
        } else {
            result = sml_log2_normal_pdf(
                value, dist.normal_left_trunc.mean, dist.normal_left_trunc.sd
            );
        }
    } break;
    }
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
sml_sqrt(double value) {
    __m128d neg_ln_u1_w = _mm_set_sd(value);
    __m128d sqrt_neg_ln_u1_w = _mm_sqrt_pd(neg_ln_u1_w);
    double result = *(double*)&sqrt_neg_ln_u1_w;
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
    double sqrt_neg_ln_u1 = sml_sqrt(-ln_u1);
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

double
sml_dot(double* vec1, double* vec2, uint32_t dim) {
    double result = 0;
    for (uint32_t index = 0; index < dim; index++) {
        result += vec1[index] * vec2[index];
    }
    return result;
}

double
sml_rnorm(pcg64_random_t* rng, double mean, double sd) {
    double result = sml_rnorm01(rng) * sd + mean;
    return result;
}

double
sml_rnorm_left_trunc(pcg64_random_t* rng, double mean, double sd, double min) {
    for (;;) {
        double result = sml_rnorm(rng, mean, sd);
        if (result > min) {
            return result;
        }
    }
}

void
sml_rnorm_mul(pcg64_random_t* rng, double* dest, double* mean, double* chol, double* std_norm, uint32_t dim) {
    for (uint32_t index = 0; index < dim; index++) {
        std_norm[index] = sml_rnorm01(rng);
    }
    for (uint32_t index = 0; index < dim; index++) {
        double chol_dot_sample = sml_dot(chol + index * dim, std_norm, dim);
        dest[index] = chol_dot_sample + mean[index];
    }
}

int32_t
sml_rbern(pcg64_random_t* rng, double prop) {
    double rand01 = random_real01(rng);
    double result = rand01 <= prop;
    return result;
}

double
sml_log2_prior_prob(SmlParameters* pars, SmlPriors* priors) {
    double result = 0;
    result += sml_log2_dist_pdf(pars->baseline, priors->baseline);
    result += sml_log2_dist_pdf(pars->residual_sd, priors->residual_sd);
    return result;
}

double
sml_log2_likelihood(SmlInput* input, SmlParameters* pars, SmlConstants* consts) {

    double log2_likelihood = 0;

    for (uint64_t individual_index = 0;
        individual_index < input->n_individuals;
        individual_index++) {

        // TODO(sen) Random intercept (baseline_sd)

        SmlInputIndividual* individual = input->data + individual_index;

        for (uint64_t titre_index = 0;
            titre_index < individual->titre_count;
            titre_index++) {

            SmlInputTitre* titre = individual->titres + titre_index;

            double predicted_titre = pars->baseline;

            for (uint64_t event_index = 0;
                event_index < individual->event_count;
                event_index++) {

                SmlInputEvent* event = individual->events + event_index;

                if (titre->time >= event->time) {
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

            double titre_prob = sml_log2_normal_pdf(titre->log2titre, predicted_titre, pars->residual_sd);
            log2_likelihood += titre_prob;

        } // NOTE(sen) for (titre)

    } // NOTE(sen) for (individual)

    return log2_likelihood;
}

double sml_get_mean_of_accepted(SmlParameters* pars, uint64_t n_accept, uint32_t par_index) {
    double first_val = pars[0].par[par_index];
    double sum = first_val;
    for (uint32_t index = 1; index < n_accept; index++) {
        double this_val = pars[index].par[par_index];
        sum += this_val;
    }
    double mean = sum / n_accept;
    return mean;
}

double
sml_get_cov_of_accepted(SmlParameters* pars, uint64_t n_accept, uint32_t par_index1, uint32_t par_index2) {
    double mean1 = sml_get_mean_of_accepted(pars, n_accept, par_index1);
    double mean2 = sml_get_mean_of_accepted(pars, n_accept, par_index2);
    double first_val1 = pars[0].par[par_index1];
    double first_val2 = pars[0].par[par_index2];
    double sum_dev = (first_val1 - mean1) * (first_val2 - mean2);
    for (uint32_t index = 1; index < n_accept; index++) {
        double this_val1 = pars[index].par[par_index1];
        double this_val2 = pars[index].par[par_index2];
        sum_dev += (this_val1 - mean1) * (this_val2 - mean2);
    }
    double cov = sum_dev / (n_accept - 1);
    return cov;
}

double
sml_rdist(pcg64_random_t* rng, SmlDist prior) {
    double result;
    switch (prior.type) {
    case SmlDist_Normal: {
        result = sml_rnorm(rng, prior.normal.mean, prior.normal.sd);
    } break;
    case SmlDist_NormalLeftTrunc: {
        result = sml_rnorm_left_trunc(
            rng,
            prior.normal_left_trunc.mean,
            prior.normal_left_trunc.sd,
            prior.normal_left_trunc.min
        );
    } break;
    }
    return result;
}

double
sml_get_var(SmlDist* dist) {
    double result;
    switch (dist->type) {
    case SmlDist_Normal: {
        result = dist->normal.sd * dist->normal.sd;
    } break;
    case SmlDist_NormalLeftTrunc: {
        result = dist->normal_left_trunc.sd * dist->normal_left_trunc.sd;
    } break;
    }
    return result;
}

double
sml_get_sd(SmlDist* dist) {
    double result;
    switch (dist->type) {
    case SmlDist_Normal: {
        result = dist->normal.sd;
    } break;
    case SmlDist_NormalLeftTrunc: {
        result = dist->normal_left_trunc.sd;
    } break;
    }
    return result;
}

void
sml_cholesky(double* in, double* out, int32_t dim) {
    // Adapted from
    // https://rosettacode.org/wiki/Cholesky_decomposition#C

    // TODO(sen) see if I need to add a small multiple of identity

    for (uint32_t out_index = 0; out_index < dim * dim; out_index++) {
        out[out_index] = 0;
    }

    for (int ii = 0; ii < dim; ii++) {
        for (int jj = 0; jj < (ii + 1); jj++) {
            double sum = 0;
            for (int kk = 0; kk < jj; kk++) {
                sum += out[ii * dim + kk] * out[jj * dim + kk];
            }
            out[ii * dim + jj] = (ii == jj) ?
                sml_sqrt(in[ii * dim + ii] - sum) :
                (1.0 / out[jj * dim + jj] * (in[ii * dim + jj] - sum));
        }
    }
}

void
sml_set_vec(double* dest, double* source, uint32_t dim) {
    for (uint32_t index = 0; index < dim; index++) {
        dest[index] = source[index];
    }
}

void
sml_mcmc(
    SmlInput* input,
    SmlParameters* pars_init,
    SmlOutput* output,
    SmlConstants* consts,
    SmlPriors* priors,
    SmlStep* step
) {
    SmlParameters pars_cur = *pars_init;
    double log2_prior_prob_cur = sml_log2_prior_prob(&pars_cur, priors);
    double log2_likelihood_cur = sml_log2_likelihood(input, &pars_cur, consts);

    pcg64_random_t rng;
    pcg_setseq_128_srandom_r(&rng, 0, 0);

    output->n_accepted_after_burn = 0;
    output->n_accepted_burn = 0;
    output->n_burn = output->n_iterations / 10;

    for (uint32_t step_index = 0; step_index < step->dim * step->dim; step_index++) {
        step->var[step_index] = 0;
    }
    for (uint32_t step_index = 0; step_index < step->dim; step_index++) {
        step->var[step_index * step->dim + step_index] = sml_get_var(priors->dist + step_index);
    }
    sml_cholesky(step->var, step->chol, step->dim);

    for (uint64_t iteration = 1; iteration <= output->n_iterations; iteration++) {

        sml_set_vec(step->mean, pars_cur.par, step->dim);

        SmlParameters pars_next = pars_cur;
        sml_rnorm_mul(&rng, pars_next.par, step->mean, step->chol, step->std_norm, step->dim);

        double log2_prior_prob_next = sml_log2_prior_prob(&pars_next, priors);
        double log2_likelihood_next = sml_log2_likelihood(input, &pars_next, consts);

        double log2_posterior_diff = log2_likelihood_next + log2_prior_prob_next
            - log2_likelihood_cur - log2_prior_prob_cur;

        int32_t accept = 0;
        if (log2_posterior_diff >= 0) {
            accept = 1;
        } else if (log2_posterior_diff > -50) {
            double posterior_ratio = sml_pow2(log2_posterior_diff);
            accept = sml_rbern(&rng, posterior_ratio);
        }

        if (accept) {
            pars_cur = pars_next;
            if (iteration > output->n_burn) {
                ++output->n_accepted_after_burn;
            } else {
                ++output->n_accepted_burn;
            }
            log2_prior_prob_cur = log2_prior_prob_next;
            log2_likelihood_cur = log2_likelihood_next;
        }

        if (iteration > output->n_burn) {
            output->out[iteration - 1 - (output->n_burn - output->n_accepted_burn)] = pars_cur;
        } else if (accept) {
            output->out[output->n_accepted_burn - 1] = pars_cur;
        }

        if (iteration == output->n_burn) {
            if (output->n_accepted_burn < 10) {
                // NOTE(sen) Reduce step variances and continue
                output->n_burn += output->n_iterations / 10;
                if (output->n_burn > output->n_iterations) {
                    output->n_burn = output->n_iterations;
                }
                for (uint32_t step_index = 0; step_index < step->dim; step_index++) {
                    step->var[step_index * step->dim + step_index] *= 0.5;
                }
                sml_cholesky(step->var, step->chol, step->dim);
            } else {
                double var_reduction = 0.0025;
                for (uint32_t index1 = 0; index1 < step->dim; index1++) {
                    for (uint32_t index2 = 0; index2 < step->dim; index2++) {
                        step->var[index1 * step->dim + index2] = var_reduction *
                            sml_get_cov_of_accepted(
                                output->out, output->n_accepted_burn, index1, index2
                            );
                    }
                }
                sml_cholesky(step->var, step->chol, step->dim);

            }
        }
    } // NOTE(sen) for (iteration)
}
