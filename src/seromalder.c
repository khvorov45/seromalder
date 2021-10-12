#include <stdint.h>

typedef int32_t i32;
typedef uint64_t u64;
typedef float f32;
typedef double f64;

typedef void* SmlAllocator(u64 size);

typedef struct SmlInputRow {
    f64 logtitre;
    f64 time_sample;
    f64 time_infection;
} SmlInputRow;

typedef struct SmlInput {
    u64 n_individuals;
    SmlInputRow* data;
} SmlInput;

typedef struct SmlParameters {
    f64 long_term_boost;
    f64 short_term_boost;
    f64 time_to_peak;
    f64 time_to_wane;
} SmlParameters;

typedef struct SmlOutput {
    u64 n_iterations;
    SmlParameters* out;
} SmlOutput;

u64
sml_required_input_bytes(u64 n_individuals) {
    u64 result = n_individuals * sizeof(SmlInputRow);
    return result;
}

SmlInput
sml_new_input_alloc(u64 n_individuals, SmlAllocator* allocator) {
    u64 input_bytes = sml_required_input_bytes(n_individuals);
    void* input_memory = allocator(input_bytes);
    SmlInput input = { .n_individuals = n_individuals, .data = input_memory };
    return input;
}

u64
sml_required_output_bytes(u64 n_iterations) {
    u64 result = sizeof(SmlParameters) * n_iterations;
    return result;
}

SmlOutput
sml_new_output_alloc(u64 n_iterations, SmlAllocator* allocator) {
    u64 output_bytes = sml_required_output_bytes(n_iterations);
    void* output_memory = allocator(output_bytes);
    SmlOutput output = { .n_iterations = n_iterations, .out = output_memory };
    return output;
}

void
sml_mcmc(SmlInput* input, SmlParameters* pars_init, SmlOutput* output) {
    SmlParameters pars_cur = *pars_init;

    for (i32 iteration = 0; iteration < output->n_iterations; iteration++) {

        f64 sum_of_squares = 0;
        for (i32 individual_index = 0;
            individual_index < input->n_individuals;
            individual_index++) {

            SmlInputRow* row = input->data + individual_index;

            f64 time_sample = row->time_sample;
            f64 time_infection = row->time_infection;
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

            f64 deviation = predicted_titre - row->logtitre;
            sum_of_squares += deviation * deviation;
        }

        output->out[iteration] = pars_cur;
    }
}
