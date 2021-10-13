#include "../../src/seromalder.c"

#include "stdlib.h"

int
main() {
    SmlConstants constants = sml_default_constants();

    SmlInput input;
    input.n_individuals = 1;
    input.data = malloc(input.n_individuals * sizeof(SmlInputIndividual));
    for (uint64_t individual_index = 0; individual_index < input.n_individuals; individual_index++) {
        SmlInputIndividual* individual = input.data + individual_index;
        individual->event_count = 1;
        individual->events = malloc(individual->event_count * sizeof(SmlInputEvent));
        for (uint64_t event_index = 0; event_index < individual->event_count; event_index++) {
            SmlInputEvent* event = individual->events + event_index;
            event->type = SmlEvent_Vaccination;
            event->time = 0;
        }
        individual->titre_count = 3;
        individual->titres = malloc(individual->titre_count * sizeof(SmlInputTitre));
        for (uint64_t titre_index = 0; titre_index < individual->titre_count; titre_index++) {
            SmlInputTitre* titre = individual->titres + titre_index;
            titre->log2titre = constants.lowest_log2titre;
        }
    }

    SmlParameters pars_init = {
        .vaccination_log2diff = 2,
        .baseline = constants.lowest_log2titre,
        .baseline_sd = 1,
        .wane_rate = 0,
        .residual_sd = 1,
    };

    SmlOutput output;
    output.n_iterations = 1;
    output.out = malloc(output.n_iterations * sizeof(SmlParameters));

    SmlMcmcSettings settings = sml_default_settings();

    sml_mcmc(&input, &pars_init, &output, &constants, &settings);
}
