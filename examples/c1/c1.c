#include "stdlib.h"
#include "../../src/seromalder.c"

int
main() {
    u64 n_individuals = 1;
    u64 input_bytes = sml_required_input_bytes(n_individuals);
    void* input_memory = malloc(input_bytes);
    SmlInput input = sml_new_input(n_individuals, input_memory);
    input.logtitres[0] = 2.321928;
    input.times_sample[0] = 20;
    input.times_infection[0] = 30;
    f64 start_long_term_boost = 0;
    f64 start_short_term_boost = 0;
    f64 start_time_to_peak = 14;
    f64 start_time_to_wane = 50;
    i32 iterations = 1;
    f64 out_long_term_boost;
    f64 out_short_term_boost;
    f64 out_start_time_to_peak;
    f64 out_time_to_wane;

    sml_mcmc(
        input,
        start_long_term_boost,
        start_short_term_boost,
        start_time_to_peak,
        start_time_to_wane,
        iterations,
        &out_long_term_boost,
        &out_short_term_boost,
        &out_start_time_to_peak,
        &out_time_to_wane
    );
}
