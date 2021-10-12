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
    i32 iterations = 1;
    f64 out_long_term_boost;
    f64 out_short_term_boost;
    f64 out_start_time_to_peak;
    f64 out_time_to_wane;

    sml_mcmc(
        input,
        pars_init,
        iterations,
        &out_long_term_boost,
        &out_short_term_boost,
        &out_start_time_to_peak,
        &out_time_to_wane
    );
}
