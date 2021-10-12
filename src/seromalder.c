#include <stdint.h>

typedef int32_t i32;
typedef uint64_t u64;
typedef float f32;
typedef double f64;

typedef void* SmlAllocator(u64 size);

typedef struct SmlInput {
    u64 n_individuals;
    f64* logtitres;
    f64* times_sample;
    f64* times_infection;
} SmlInput;

typedef struct SmlParameters {
    f64 long_term_boost;
    f64 short_term_boost;
    f64 time_to_peak;
    f64 time_to_wane;
} SmlParameters;

u64
sml_required_input_bytes(u64 n_individuals) {
    u64 input_arrays_per_individual = (sizeof(SmlInput) - sizeof(u64)) / sizeof(void*);
    u64 result = n_individuals * input_arrays_per_individual * sizeof(f64);
    return result;
}

SmlInput
sml_new_input(u64 n_individuals, void* memory) {
    u64 one_array_bytes = n_individuals * sizeof(f64);
    SmlInput input;
    input.n_individuals = n_individuals;
    input.logtitres = memory;
    input.times_sample = input.logtitres + one_array_bytes;
    input.times_infection = input.times_sample + one_array_bytes;
    return input;
}

SmlInput
sml_new_input_alloc(u64 n_individuals, SmlAllocator* allocator) {
    u64 input_bytes = sml_required_input_bytes(n_individuals);
    void* input_memory = allocator(input_bytes);
    SmlInput input = sml_new_input(n_individuals, input_memory);
    return input;
}

void
sml_mcmc(
    SmlInput input,
    SmlParameters pars_init,
    i32 iterations,
    f64* out_long_term_boost, f64* out_short_term_boost,
    f64* out_start_time_to_peak, f64* out_time_to_wane
) {
    SmlParameters pars_cur = pars_init;

    for (i32 iteration = 0; iteration < iterations; iteration++) {

        f64 sum_of_squares = 0;
        for (i32 individual_index = 0;
            individual_index < input.n_individuals;
            individual_index++) {

            f64 time_sample = input.times_sample[individual_index];
            f64 time_infection = input.times_infection[individual_index];
            f64 time_peak = time_infection + pars_cur.time_to_peak;
            f64 time_wane = time_peak + pars_cur.time_to_wane;

            f64 predicted_titre;
            if (time_sample < time_infection) {
                predicted_titre = 2.321928;
            } else if (time_sample < time_peak) {
                predicted_titre = 2.321928 +
                    (pars_cur.long_term_boost + pars_cur.short_term_boost) * (time_sample - time_infection) / pars_cur.time_to_peak;
            } else if (time_sample < time_wane) {
                predicted_titre = 2.321928 + pars_cur.long_term_boost +
                    pars_cur.short_term_boost * (1 - (time_sample - time_peak) / pars_cur.time_to_wane);
            } else {
                predicted_titre = 2.321928 + pars_cur.long_term_boost;
            }

            f64 deviation = predicted_titre - input.logtitres[individual_index];
            sum_of_squares += deviation * deviation;
        }

        out_long_term_boost[iteration] = sum_of_squares;
        out_short_term_boost[iteration] = iteration;
        out_start_time_to_peak[iteration] = iteration;
        out_time_to_wane[iteration] = iteration;
    }
}
