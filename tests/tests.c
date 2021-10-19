#include "../src/seromalder.c"

#include "math.h"
#include "stdio.h"

void exp_test() {
    int32_t count = 100;
    double min = -20;
    double max = 0;
    double step = (max - min) / (double)count;
    double deviation = 0;
    double deviation_proportion = 0;
    double value = min;
    for (int32_t index = 0; index < count; index++, value += step) {
        double result_sml = sml_pow2(value);
        double result_crt = pow(2.0, value);
        deviation += fabs(result_crt - result_sml);
        deviation_proportion += fabs(result_crt - result_sml) / result_crt;
    }
    double average_deviation = deviation / (double)count;
    double average_deviation_proportion = deviation_proportion / (double)count;
    printf(
        "sml_pow2 average deviation from math.h pow (range %f - %f): %f (%.2f%%)\n",
        min, max, average_deviation, average_deviation_proportion * 100
    );
}

void log_test() {
    int32_t count = 100;
    double min = 0.000001;
    double max = 100000;
    double step = (max - min) / (double)count;
    double deviation = 0;
    double deviation_proportion = 0;
    double value = min;
    for (int32_t index = 0; index < count; index++, value += step) {
        double result_sml = sml_log2(value);
        double result_crt = log2(value);
        deviation += fabs(result_crt - result_sml);
        deviation_proportion += fabs(result_crt - result_sml) / result_crt;
    }
    double average_deviation = deviation / (double)count;
    double average_deviation_proportion = deviation_proportion / (double)count;
    printf(
        "sml_log2 average deviation from math.h log2 (range %f - %f): %f (%.2f%%)\n",
        min, max, average_deviation, average_deviation_proportion * 100
    );
}

void test_log2_normal_pdf() {
    double values[] = { -10, -5, -1, 0, 1, 5, 10 };
    double r_results[] = {
        -73.4605001091843, -19.3594360758482, -2.04709558518064, -1.32574806473616,
        -2.04709558518064, -19.3594360758482, -73.4605001091843
    };
    int32_t count = 7;
    double deviation = 0;
    double deviation_proportion = 0;
    for (int32_t index = 0; index < count; index++) {
        double result_sml = sml_log2_std_normal_pdf(values[index]);
        double result_r = r_results[index];
        deviation += fabs(result_r - result_sml);
        deviation_proportion += fabs(result_r - result_sml) / result_r;
    }
    double average_deviation = deviation / (double)count;
    double average_deviation_proportion = deviation_proportion / (double)count;
    printf(
        "sml_log2_normal_pdf average deviation from r (range -10 - 10): %f (%.2f%%)\n",
        average_deviation, average_deviation_proportion * 100
    );
}

void test_random_real() {
    pcg64_random_t rng;
    pcg_setseq_128_srandom_r(&rng, 0, 0);
    uint32_t n_to_generate = 1000000;
    const uint32_t n_buckets = 100;
    uint32_t bucket_counts[n_buckets];
    for (uint32_t bucket_index = 0; bucket_index < n_buckets; ++bucket_index) {
        bucket_counts[bucket_index] = 0;
    }
    for (uint32_t index = 0; index < n_to_generate; index++) {
        double random01 = random_real01(&rng);
        uint32_t bucket_index = (uint32_t)(random01 * n_buckets);
        bucket_counts[bucket_index] += 1;
    }
    double expected_proportion = 1 / (double)n_buckets;
    double deviation = 0;
    for (uint32_t bucket_index = 0; bucket_index < n_buckets; ++bucket_index) {
        double proportion = (double)bucket_counts[bucket_index] / n_to_generate;
        deviation += fabs(proportion - expected_proportion);
    }
    double average_deviation = deviation / n_buckets;
    printf(
        "average deviation of random_real01: %f (%.2f%%)\n",
        average_deviation, average_deviation / expected_proportion * 100
    );
}

int
main() {
    exp_test();
    log_test();
    test_log2_normal_pdf();
    test_random_real();
}
