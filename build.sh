gcc -g -c -Wall -Werror -fpic \
    src/seromalder.c -o build/seromalder.o

gcc -shared \
    build/seromalder.o -o build/seromalder.so

echo done
