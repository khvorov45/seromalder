#include "../src/seromalder.c"

#include "math.h"
#include "stdio.h"

void exp_test() {
    int32_t count = 100;
    double min = -15;
    double max = 0;
    double step = (max - min) / (double)count;
    double deviation = 0;
    double value = min;
    for (int32_t index = 0; index < count; index++) {
        value += step;
        double result_sml = sml_exp(value);
        double result_crt = exp(value);
        deviation += fabs(result_crt - result_sml);
    }
    double average_deviation = deviation / (double)count;
    printf("sml_exp average deviation from math.h exp (range %f - %f): %f\n", min, max, average_deviation);
}

int
main() {
    exp_test();
}
