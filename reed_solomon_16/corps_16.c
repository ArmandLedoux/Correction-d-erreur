#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>


#define P 65579
#define alpha 3
#define n 65536 

int tab_alpha[n];
int tab_pow[n];

int deg(int a) {
    int d = -1;
    int pow = 1;
    while(pow <= a) {
        d+=1 ;
        pow = pow<<1;
    }
    return d;
}

int mul_ext (int a, int b) {
    int p = 0;
    for (int i=0; i<=deg(a); i++) {
        for (int j=0; j<=deg(b); j++){
            p = p ^ ((a & (1<<i))*(b & (1<<j)));
        }
    }
    return p;
}

int GF (int a) {
    int dega = deg(a);
    int degb = 16;
    while (degb <= dega) {
        int diff = dega - degb;
        a = a ^ (P<<diff);
        for (int i = dega; (a & (1<<i)) == 0 && i>=0 ; i--) {
            dega -= 1;
        }
    }
    return a;
}

int add(int a, int b) {
    return a^b;
}

void init_tab_alpha(void) {
    // printf("Initialisation de tab_alpha et tab_pow...\n");
    int puiss = 1;
    for (int i=0; i<n; i++) {
        tab_alpha[puiss]=i;
        tab_pow [i] = puiss;
        puiss = GF(mul_ext(puiss, alpha));
    }
    tab_alpha[1] = 0;
    tab_alpha[0] = -1;
    // printf("Initialisation terminée.\n");
}

int mul (int a, int b) {
    assert(a>=0 && b>= 0 && a<n && b<n);
    if (a==0 || b==0) {
        return 0;
    }
    int powa = tab_alpha[a];
    int powb = tab_alpha[b];
    return tab_pow[(powa+powb)%(n-1)];
}

int puiss (int a, int b) {
    assert(a>=0);
    assert(a<n);
    if (a==0) {
        assert(b>=0);
        return (int)(b==0);
    }
    int powa = tab_alpha[a];
    long prod = powa*b;
    return tab_pow[(prod%(n-1) + n-1)%(n-1)]; 
    // on fait cela car (-1 % (n-1) == -1), -1 n'étant pas un indice valide
}

int inv (int a) {
    return (puiss(a, -1));
}

int accede_alpha(int a){
    assert(a>0);
    assert(a<n);
    return tab_alpha[a];
}

int division (int a, int b) {
    return mul(a, inv(b));
}

