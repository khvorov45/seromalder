#include <stdint.h>

typedef int32_t i32;
typedef uint64_t u64;
typedef float f32;
typedef double f64;
typedef i32 b32;

typedef void* SmlAllocator(u64 size);

typedef struct SmlInputTitre {
    f64 log2titre;
    f64 time;
} SmlInputTitre;

typedef enum SmlEventType {
    SmlEvent_Vaccination,
    SmlEvent_Infection,
} SmlEventType;

typedef struct SmlInputEvent {
    SmlEventType type;
    f64 time;
} SmlInputEvent;

typedef struct SmlInputIndividual {
    u64 event_count;
    SmlInputEvent* events;
    u64 titre_count;
    SmlInputTitre* titres;
} SmlInputIndividual;

typedef struct SmlInput {
    u64 n_individuals;
    SmlInputIndividual* data;
} SmlInput;

typedef struct SmlParameters {
    f64 vaccination_log2diff;
    f64 baseline;
    f64 baseline_sd;
    f64 wane_rate;
} SmlParameters;

typedef struct SmlOutput {
    u64 n_iterations;
    SmlParameters* out;
} SmlOutput;

typedef struct SmlConstants {
    f64 time_to_peak;
    f64 lowest_log2titre;
} SmlConstants;

SmlConstants
sml_default_constants() {
    SmlConstants result = {
        .time_to_peak = 14, // NOTE(sen) Days
        .lowest_log2titre = 2.321928, // NOTE(sen) Log2(5)
    };
    return result;
}

f64
sml_rnorm(f64 mean, f64 sd) {
    return mean;
}

b32
sml_rbern(f64 prop) {
    return 0;
}

f64
sml_prior_prob(SmlParameters* pars) {
    return 0;
}

f64
sml_log_likelihood(SmlInput* input, SmlParameters* pars, SmlConstants* consts) {

    f64 log_likelihood = 0;

    for (u64 individual_index = 0;
        individual_index < input->n_individuals;
        individual_index++) {

        // SmlInputIndividual* individual = input->data + individual_index;

        log_likelihood += 0;
    } // NOTE(sen) for (individual)

    return log_likelihood;
}

void
sml_mcmc(SmlInput* input, SmlParameters* pars_init, SmlOutput* output, SmlConstants* consts) {

    SmlParameters pars_cur = *pars_init;
    f64 prior_prob_cur = sml_prior_prob(&pars_cur);
    f64 log_likelihood_cur = sml_log_likelihood(input, &pars_cur, consts);
    f64 posterior_cur = prior_prob_cur * log_likelihood_cur;

    for (u64 iteration = 0; iteration < output->n_iterations; iteration++) {

        SmlParameters pars_next;
        pars_next.vaccination_log2diff = sml_rnorm(pars_cur.vaccination_log2diff, 1);
        pars_next.baseline = sml_rnorm(pars_cur.baseline, 1);
        pars_next.baseline_sd = sml_rnorm(pars_cur.baseline_sd, 1);
        pars_next.wane_rate = sml_rnorm(pars_cur.wane_rate, 1);

        f64 prior_prob_next = sml_prior_prob(&pars_next);
        f64 log_likelihood_next = sml_log_likelihood(input, &pars_next, consts);
        f64 posterior_next = prior_prob_next * log_likelihood_next;

        f64 posterior_ratio = posterior_next / posterior_cur;

        if (posterior_ratio > 1 || sml_rbern(posterior_ratio)) {
            pars_cur = pars_next;
        }

        output->out[iteration] = pars_cur;
    } // NOTE(sen) for (iteration)
}
