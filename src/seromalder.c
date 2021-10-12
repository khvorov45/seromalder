#include <stdint.h>

typedef int32_t i32;
typedef float f32;
typedef double f64;

void
sml_mcmc(
    i32 n_individuals,
    f64* logtitres, f64* times_sample, f64* times_infection,
    f64 start_long_term_boost, f64 start_short_term_boost,
    f64 start_time_to_peak, f64 start_time_to_wane,
    i32 iterations,
    f64* out_long_term_boost, f64* out_short_term_boost,
    f64* out_start_time_to_peak, f64* out_time_to_wane
) {
    f64 cur_long_term_boost = start_long_term_boost;
    f64 cur_short_term_boost = start_short_term_boost;
    f64 cur_time_to_peak = start_time_to_peak;
    f64 cur_time_to_wane = start_time_to_wane;

    for (i32 iteration = 0; iteration < iterations; iteration++) {

        f64 sum_of_squares = 0;
        for (i32 individual_index = 0;
            individual_index < n_individuals;
            individual_index++) {

            f64 time_sample = times_sample[individual_index];
            f64 time_infection = times_infection[individual_index];
            f64 time_peak = time_infection + cur_time_to_peak;
            f64 time_wane = time_peak + cur_time_to_wane;

            f64 predicted_titre;
            if (time_sample < time_infection) {
                predicted_titre = 2.321928;
            } else if (time_sample < time_peak) {
                predicted_titre = 2.321928 +
                    (cur_long_term_boost + cur_short_term_boost) * (time_sample - time_infection) / cur_time_to_peak;
            } else if (time_sample < time_wane) {
                predicted_titre = 2.321928 + cur_long_term_boost +
                    cur_short_term_boost * (1 - (time_sample - time_peak) / cur_time_to_wane);
            } else {
                predicted_titre = 2.321928 + cur_long_term_boost;
            }

            f64 deviation = predicted_titre - logtitres[individual_index];
            sum_of_squares += deviation * deviation;
        }

        out_long_term_boost[iteration] = sum_of_squares;
        out_short_term_boost[iteration] = iteration;
        out_start_time_to_peak[iteration] = iteration;
        out_time_to_wane[iteration] = iteration;
    }
}

void
sml_r_mcmc(
    i32* n_individuals,
    f64* logtitres, f64* times_sample, f64* times_infection,
    f64* start_long_term_boost, f64* start_short_term_boost,
    f64* start_time_to_peak, f64* start_time_to_wane,
    i32* iterations,
    f64* out_long_term_boost, f64* out_short_term_boost,
    f64* out_start_time_to_peak, f64* out_time_to_wane
) {
    sml_mcmc(
        *n_individuals,
        logtitres, times_sample, times_infection,
        *start_long_term_boost, *start_short_term_boost,
        *start_time_to_peak, *start_time_to_wane,
        *iterations,
        out_long_term_boost, out_short_term_boost,
        out_start_time_to_peak, out_time_to_wane
    );
}
