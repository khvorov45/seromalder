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
} SmlParameters;

typedef struct SmlOutput {
    uint64_t n_iterations;
    SmlParameters* out;
} SmlOutput;

typedef struct SmlConstants {
    double time_to_peak;
    double lowest_log2titre;
} SmlConstants;

SmlConstants
sml_default_constants() {
    SmlConstants result = {
        .time_to_peak = 14, // NOTE(sen) Days
        .lowest_log2titre = 2.321928, // NOTE(sen) Log2(5)
    };
    return result;
}

double
sml_rnorm(double mean, double sd) {
    return mean;
}

int
sml_rbern(double prop) {
    return 0;
}

double
sml_prior_prob(SmlParameters* pars) {
    return 0;
}

double
sml_log_likelihood(SmlInput* input, SmlParameters* pars, SmlConstants* consts) {

    double log_likelihood = 0;

    for (uint64_t individual_index = 0;
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
    double prior_prob_cur = sml_prior_prob(&pars_cur);
    double log_likelihood_cur = sml_log_likelihood(input, &pars_cur, consts);
    double posterior_cur = prior_prob_cur * log_likelihood_cur;

    for (uint64_t iteration = 0; iteration < output->n_iterations; iteration++) {

        SmlParameters pars_next;
        pars_next.vaccination_log2diff = sml_rnorm(pars_cur.vaccination_log2diff, 1);
        pars_next.baseline = sml_rnorm(pars_cur.baseline, 1);
        pars_next.baseline_sd = sml_rnorm(pars_cur.baseline_sd, 1);
        pars_next.wane_rate = sml_rnorm(pars_cur.wane_rate, 1);

        double prior_prob_next = sml_prior_prob(&pars_next);
        double log_likelihood_next = sml_log_likelihood(input, &pars_next, consts);
        double posterior_next = prior_prob_next * log_likelihood_next;

        double posterior_ratio = posterior_next / posterior_cur;

        if (posterior_ratio > 1 || sml_rbern(posterior_ratio)) {
            pars_cur = pars_next;
        }

        output->out[iteration] = pars_cur;
    } // NOTE(sen) for (iteration)
}
