#include "../../src/seromalder.c"

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

int
main() {
    pcg64_random_t rng;
    pcg_setseq_128_srandom_r(&rng, 1, 0);

    SmlConstants constants = sml_default_constants();

    SmlParameters pars_true = {
        .vaccination_log2diff = 2,
        .baseline = constants.lowest_log2titre,
        .baseline_sd = 0,
        .wane_rate = 0,
        .residual_sd = 0.2,
    };

    SmlParameters pars_init = {
        .vaccination_log2diff = 2,
        .baseline = sml_log2(10),
        .baseline_sd = 0,
        .wane_rate = 0,
        .residual_sd = 1,
    };

    SmlInput input;
    input.n_individuals = 1000;
    input.data = malloc(input.n_individuals * sizeof(SmlInputIndividual));

    FILE* input_csv_file = fopen("examples/c1/input.csv", "w");
    char* header = "pid,time,log2titre\n";
    fwrite(header, strlen(header), 1, input_csv_file);

    for (uint64_t individual_index = 0;
        individual_index < input.n_individuals;
        individual_index++) {

        SmlInputIndividual* individual = input.data + individual_index;
        individual->event_count = 1;
        individual->events = malloc(individual->event_count * sizeof(SmlInputEvent));

        for (uint64_t event_index = 0;
            event_index < individual->event_count;
            event_index++) {

            SmlInputEvent* event = individual->events + event_index;
            event->type = SmlEvent_Vaccination;
            event->time = 0;
        }

        individual->titre_count = 1;
        individual->titres = malloc(individual->titre_count * sizeof(SmlInputTitre));
        double titre_times[3] = { 0, 14, 200 };

        for (uint64_t titre_index = 0;
            titre_index < individual->titre_count;
            titre_index++) {

            SmlInputTitre* titre = individual->titres + titre_index;

            titre->time = titre_times[titre_index];

            double time_since = titre->time; // NOTE(sen) Only one event at time 0
            double up_slope = pars_true.vaccination_log2diff / constants.time_to_peak;
            double past_peak = (double)(time_since > constants.time_to_peak);

            double log2titre_expected = pars_true.baseline +
                up_slope * time_since -
                past_peak * (up_slope + pars_true.wane_rate) *
                (time_since - constants.time_to_peak);

            titre->log2titre =
                sml_rnorm01(&rng) * pars_true.residual_sd + log2titre_expected;

            char csv_row[64];
            int csv_row_len = snprintf(csv_row, 64, "%ld,%.0f,%f\n", individual_index, titre_times[titre_index], titre->log2titre);
            fwrite(csv_row, csv_row_len, 1, input_csv_file);
        }
    }

    fclose(input_csv_file);

    SmlOutput output;
    output.n_iterations = 100000;
    output.n_burn = 10000;
    output.out = malloc(output.n_iterations * sizeof(SmlParameters));

    SmlPriors priors;
    priors.baseline.type = SmlDist_Normal;
    priors.baseline.normal.mean = 4.3;
    priors.baseline.normal.sd = 2;
    priors.residual_sd.type = SmlDist_NormalLeftTrunc;
    priors.residual_sd.normal_left_trunc.mean = 2;
    priors.residual_sd.normal_left_trunc.sd = 1;
    priors.residual_sd.normal_left_trunc.min = 0;
    sml_mcmc(&input, &pars_init, &output, &constants, &priors);

    printf(
        "acceptance rate: %f\n",
        (double)output.n_accepted_after_burn / (double)(output.n_iterations - output.n_burn)
    );

    FILE* output_csv = fopen("examples/c1/output.csv", "w");
    char* output_header = "iteration,vaccination_log2diff,baseline,baseline_sd,wane_rate,residual_sd\n";
    fwrite(output_header, strlen(output_header), 1, output_csv);

    for (uint32_t iteration_index = 0; iteration_index < output.n_iterations - (output.n_burn - output.n_accepted_burn); iteration_index++) {
        SmlParameters* pars = output.out + iteration_index;
        char csv_row[64];
        int csv_row_len = snprintf(
            csv_row, 64, "%d,%f,%f,%f,%f,%f\n",
            iteration_index,
            pars->vaccination_log2diff, pars->baseline, pars->baseline_sd, pars->wane_rate, pars->residual_sd
        );
        fwrite(csv_row, csv_row_len, 1, output_csv);
    }

    fclose(output_csv);
}
