gcc -g -c -Wall -Werror -fpic \
    src/seromalder.c -o build/seromalder.o

gcc -shared \
    build/seromalder.o -o build/seromalder.so

gcc -g examples/c1/c1.c -o examples/c1/c1

gcc -g tests/tests.c -o tests/tests -lm

./tests/tests

echo done
