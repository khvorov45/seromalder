#include "../../src/seromalder.c"

#include "stdlib.h"

int
main() {
    u64 n_individuals = 1;
    SmlInput input = sml_new_input_alloc(n_individuals, malloc);
    input.data[0].logtitre = 2.321928;
    input.data[0].time_sample = 20;
    input.data[0].time_infection = 30;
    SmlParameters pars_init = {
        .long_term_boost = 0,
        .short_term_boost = 0,
        .time_to_peak = 14,
        .time_to_wane = 50
    };
    i32 n_iterations = 1;
    SmlOutput out = sml_new_output_alloc(n_iterations, malloc);

    sml_mcmc(&input, &pars_init, &out);
}
