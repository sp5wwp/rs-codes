all: test.c
	gcc test.c rs.c -Wall -Wextra -O2 -o test -lm

clean:
	rm test
