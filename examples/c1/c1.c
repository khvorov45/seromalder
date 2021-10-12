#include "../../src/seromalder.c"

#include "stdlib.h"

int
main() {
    u64 n_individuals = 1;
    SmlInput input = sml_new_input_alloc(n_individuals, malloc);
    input.logtitres[0] = 2.321928;
    input.times_sample[0] = 20;
    input.times_infection[0] = 30;
    SmlParameters pars_init = {
        .long_term_boost = 0,
        .short_term_boost = 0,
        .time_to_peak = 14,
        .time_to_wane = 50
    };
    i32 n_iterations = 1;
    SmlOutput out = sml_new_output_alloc(n_iterations, malloc);

    sml_mcmc(input, pars_init, out);
}
