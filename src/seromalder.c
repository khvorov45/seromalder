#include <stdint.h>

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
    double closest_whole = (double)(value_exponent_biased - 1023);

    int64_t value_remainder_bits = (value_bits & 0x000fffffffffffff) | 0x3ff0000000000000;
    double value_remainder = *(double*)&value_remainder_bits;
    double value_remainder_m1 = value_remainder - 1;
    double frac_e = value_remainder_m1 -
        (value_remainder_m1) * (value_remainder_m1) * 0.5 +
        (value_remainder_m1) * (value_remainder_m1) * (value_remainder_m1) * 0.33333333333333333;

    double log2e = 1.442695040888963387005;
    double frac_2 = frac_e * log2e;

    double result = closest_whole + frac_2;

    return result;
}

double
sml_log2_normal_pdf(double value, double mean, double sd) {
    double log2_one_over_sqrt_2pi = -1.325748064736159248511;
    double log2_one_over_sd = -sml_log2(sd);

    double half_log2e = 0.7213475204444816935023;
    double value_standardized = (value - mean) / sd;

    double result = log2_one_over_sqrt_2pi + log2_one_over_sd - half_log2e * value_standardized * value_standardized;

    return result;
}

double
sml_rnorm(double mean, double sd) {
    // TODO(sen) Implement
    return mean;
}

int32_t
sml_rbern(double prop) {
    // TODO(sen) Implement
    return 0;
}

double
sml_log_prior_prob(SmlParameters* pars) {
    // TODO(sen) Implement
    return 0;
}

double
sml_log_likelihood(SmlInput* input, SmlParameters* pars, SmlConstants* consts) {

    double log_likelihood = 0;

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
            double titre_prob = sml_log2_normal_pdf(deviation, predicted_titre, pars->residual_sd);
            log_likelihood += titre_prob;

        } // NOTE(sen) for (titre)

        log_likelihood += 0;
    } // NOTE(sen) for (individual)

    return log_likelihood;
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
    double log_prior_prob_cur = sml_log_prior_prob(&pars_cur);
    double log_likelihood_cur = sml_log_likelihood(input, &pars_cur, consts);
    double log_posterior_cur = log_prior_prob_cur + log_likelihood_cur;

    SmlParameters* steps = &settings->proposal_sds;

    for (uint64_t iteration = 0; iteration < output->n_iterations; iteration++) {

        SmlParameters pars_next;
        pars_next.vaccination_log2diff =
            sml_rnorm(pars_cur.vaccination_log2diff, steps->vaccination_log2diff);
        pars_next.baseline = sml_rnorm(pars_cur.baseline, steps->baseline);
        pars_next.baseline_sd = sml_rnorm(pars_cur.baseline_sd, steps->baseline_sd);
        pars_next.wane_rate = sml_rnorm(pars_cur.wane_rate, steps->wane_rate);

        double log_prior_prob_next = sml_log_prior_prob(&pars_next);
        double log_likelihood_next = sml_log_likelihood(input, &pars_next, consts);
        double log_posterior_next = log_prior_prob_next + log_likelihood_next;

        double log_posterior_diff = log_posterior_next - log_posterior_cur;

        if (log_posterior_diff >= 0) {
            pars_cur = pars_next;
        } else if (log_posterior_diff > -20) {
            double posterior_ratio = sml_pow2(log_posterior_diff);
            if (sml_rbern(posterior_ratio)) {
                pars_cur = pars_next;
            }
        }

        output->out[iteration] = pars_cur;
    } // NOTE(sen) for (iteration)
}
