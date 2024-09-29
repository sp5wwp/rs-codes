#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "rs.h"

rs_t rs;
uint8_t poly[]={1,1,0,0,1}; //dec=19 (ascending exponent order)

uint8_t data[9]={6, 15, 8, 9, 8, 3, 0, 0, 5};
uint8_t pty[6]={0}; //0 12 11 2 0 9
int8_t cword[15]={0};

int main(void)
{
    init_RS(&rs, 15, 9, poly);

    //indexes dump
    /*printf("indx: ");
    for(uint8_t i=0; i<rs.n+1; i++) printf("%d ", rs.index[i]);
    printf("\n");*/

    //alpha dump
    /*printf("alph: ");
    for(uint8_t i=0; i<rs.n+1; i++) printf("%d ", rs.alpha[i]);
    printf("\n");*/

    encode_RS(&rs, pty, data);

    printf("data: ");
    for(uint8_t i=0; i<rs.k; i++) printf("%d ", data[i]);
    printf("\n");

    printf("parity: ");
    for(uint8_t i=0; i<2*rs.t; i++) printf("%d ", pty[i]);
    printf("\n");

    memcpy(&cword[0], data, sizeof(data));
    memcpy(&cword[9], pty, sizeof(pty));

    printf("encd: ");
    for(uint8_t i=0; i<rs.n; i++) printf("%d ", cword[i]);
    printf("\n");

    //error
    cword[3]=0; 
    printf("recd: ");
    for(uint8_t i=0; i<rs.n; i++) printf("%d ", cword[i]);
    printf("\n");

    RS_decode(&rs, cword);

    printf("decd: ");
    for(uint8_t i=0; i<rs.n; i++) printf("%d ", cword[i]);
    printf("\n");

    return 0;
}
