#include "../../src/seromalder.c"

int
main() {
    i32 n_individuals = 1;
    f64 logtitres = 2.321928;
    f64 times_sample = 20;
    f64 times_infection = 30;
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
        n_individuals,
        &logtitres,
        &times_sample,
        &times_infection,
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
