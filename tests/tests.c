#include "../src/seromalder.c"

#include "math.h"
#include "stdio.h"

void test_pow2() {
    int32_t count = 100;
    double min = -50;
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

void test_log2() {
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
    double mean = 50;
    double sd = 100;
    double values[] = { -50, -10, 0, 15, 35, 50, 75, 100, 150 };
    double r_results[] = {
        -8.690952, -8.229289, -8.149941, -8.057969,
        -7.985835, -7.969604, -8.014688, -8.149941, -8.690952
    };
    int32_t count = 9;
    double deviation = 0;
    double deviation_proportion = 0;
    for (int32_t index = 0; index < count; index++) {
        double result_sml = sml_log2_normal_pdf(values[index], mean, sd);
        double result_r = r_results[index];
        deviation += fabs(result_r - result_sml);
        deviation_proportion += fabs((result_r - result_sml) / result_r);
    }
    double average_deviation = deviation / (double)count;
    double average_deviation_proportion = deviation_proportion / (double)count;
    printf(
        "sml_log2_normal_pdf average deviation from r: %f (%.2f%%)\n",
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

void
test_sin() {
    double min = -100;
    double max = 100;
    double step = 0.1;
    double deviation = 0;
    double deviation_prop = 0;
    uint32_t count = 0;
    for (double test_value = min; test_value <= max; test_value += step, count++) {
        double sml_result = sml_sin(test_value);
        double crt_result = sin(test_value);
        double abs_deviation = abs(sml_result - crt_result);
        deviation += abs_deviation;
        deviation_prop += abs_deviation / crt_result;
    }
    double av_deviation = deviation / (double)count;
    double av_deviation_prop = deviation_prop / (double)count;
    printf(
        "average deviation of sml_sin from crt sin (range %.0f to %.0f): %f (%.2f%%)\n",
        min, max, av_deviation, av_deviation_prop * 100
    );
}

void
test_rnorm() {
    pcg64_random_t rng;
    pcg_setseq_128_srandom_r(&rng, 0, 0);
    uint32_t n_to_generate = 1000000;
    uint32_t n_buckets = 7;
    double cutoffs[] = {
        -1.959963984540053605343246090342290699481964111328125,
        -0.5244005127080406669648482420598156750202178955078125,
        -0.25334710313579977825071409824886359274387359619140625,
        0,
        0.25334710313579977825071409824886359274387359619140625,
        0.5244005127080406669648482420598156750202178955078125,
        1.959963984540053605343246090342290699481964111328125
    };
    double expected_proportions[] = {
        0.025,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.975
    };
    double bucket_counts[] = {
        0,
        0,
        0,
        0,
        0,
        0,
        0
    };
    for (uint32_t index = 0; index < n_to_generate; index++) {
        double result1 = sml_rnorm01(&rng);
        for (uint32_t bucket_index = 0; bucket_index < n_buckets; ++bucket_index) {
            if (result1 < cutoffs[bucket_index]) {
                bucket_counts[bucket_index] += 1;
            }
        }
    }
    double deviation = 0;
    for (uint32_t bucket_index = 0; bucket_index < n_buckets; ++bucket_index) {
        double proportion = (double)bucket_counts[bucket_index] / n_to_generate;
        double expected_proportion = expected_proportions[bucket_index];
        deviation += fabs(proportion - expected_proportion);
    }
    double average_deviation = deviation / n_buckets;
    printf(
        "average deviation of sml_rnorm: %f\n",
        average_deviation
    );
}

int32_t
test_mat_equal(double* mat1, double* mat2, uint32_t dim) {
    for (uint32_t mat_index = 0; mat_index < dim * dim; mat_index++) {
        if (mat1[mat_index] != mat2[mat_index]) {
            return 0;
        }
    }
    return 1;
}

void
test_mat_print(double* mat, uint32_t dim) {
    for (uint32_t row_index = 0; row_index < dim; row_index++) {
        for (uint32_t col_index = 0; col_index < dim; col_index++) {
            double val = mat[row_index * dim + col_index];
            printf("%f ", val);
        }
        printf("\n");
    }
}

void
test_cholesky() {
    int32_t dim = 2;
    double* in = malloc(dim * dim * sizeof(double));
    double* out = malloc(dim * dim * sizeof(double));
    in[0] = 1;
    in[1] = 0;
    in[2] = 0;
    in[3] = 1;
    sml_cholesky(in, out, dim);
    int32_t result1 = test_mat_equal(in, out, dim);
    printf("cholesky dim2 identity: %d\n", result1);

    //
    //
    //

    double* expected = malloc(dim * dim * sizeof(double));
    in[0] = 25;
    in[1] = 15;
    in[2] = 15;
    in[3] = 18;
    expected[0] = 5;
    expected[1] = 0;
    expected[2] = 3;
    expected[3] = 3;
    sml_cholesky(in, out, dim);
    int32_t result2 = test_mat_equal(expected, out, dim);
    printf("cholesky dim2 test: %d\n", result2);
}

void
test_cov() {
    uint32_t size = 3;
    SmlParameters* pars = malloc(sizeof(SmlParameters) * size);

    pars[0].par[0] = 1;
    pars[1].par[0] = 2;
    pars[2].par[0] = 3;

    pars[0].par[1] = 6;
    pars[1].par[1] = 5;
    pars[2].par[1] = 4;

    SmlParameters means;
    means.par[0] = sml_get_mean(pars, size, 0);
    means.par[1] = sml_get_mean(pars, size, 1);

    int32_t result_mean = means.par[0] == 2 && means.par[1] == 5;
    printf("mean test: %d\n", result_mean);

    double cov00 = sml_get_cov(pars, size, &means, 0, 0);
    double cov01 = sml_get_cov(pars, size, &means, 0, 1);
    double cov10 = sml_get_cov(pars, size, &means, 1, 0);
    double cov11 = sml_get_cov(pars, size, &means, 1, 1);

    int32_t result = cov00 == 1 && cov01 == -1 && cov10 == -1 && cov11 == 1;
    printf("cov test: %d\n", result);
}

int
main() {
    test_pow2();
    test_log2();
    test_log2_normal_pdf();
    test_random_real();
    test_sin();
    test_rnorm();
    test_cholesky();
    test_cov();
}
