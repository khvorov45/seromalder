#include "../src/seromalder.c"

#include "math.h"
#include "stdio.h"

void exp_test() {
    int32_t count = 100;
    double min = -20;
    double max = 0;
    double step = (max - min) / (double)count;
    double deviation = 0;
    double value = min;
    for (int32_t index = 0; index < count; index++, value += step) {
        double result_sml = sml_pow2(value);
        double result_crt = pow(2.0, value);
        deviation += fabs(result_crt - result_sml);
    }
    double average_deviation = deviation / (double)count;
    printf("sml_pow2 average deviation from math.h pow (range %f - %f): %f\n", min, max, average_deviation);
}

void log_test() {
    int32_t count = 100;
    double min = 0.000001;
    double max = 100000;
    double step = (max - min) / (double)count;
    double deviation = 0;
    double value = min;
    for (int32_t index = 0; index < count; index++, value += step) {
        double result_sml = sml_log2(value);
        double result_crt = log2(value);
        deviation += fabs(result_crt - result_sml);
    }
    double average_deviation = deviation / (double)count;
    printf("sml_log2 average deviation from math.h log2 (range %f - %f): %f\n", min, max, average_deviation);
}

void test_log2_normal_pdf() {
    double mean = 20;
    double sd = 5;
    double values[] = { -10, -5, -1, 0, 1, 5, 10 };
    double r_results[] = {
        -29.6161868956249, -21.6813641707356, -16.3722464202642, -15.1892364867352,
        -14.0639343548418, -10.1398038436239, -6.53306624140145
    };
    int32_t count = 7;
    double deviation = 0;
    for (int32_t index = 0; index < count; index++) {
        double result_sml = sml_log2_normal_pdf(values[index], mean, sd);
        double result_r = r_results[index];
        deviation += fabs(result_r - result_sml);
    }
    double average_deviation = deviation / (double)count;
    printf("sml_log2_normal_pdf average deviation from r (range -10 - 10): %f\n", average_deviation);
}

int
main() {
    exp_test();
    log_test();
    test_log2_normal_pdf();
}
