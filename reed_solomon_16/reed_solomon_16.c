#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define alpha 3
#define n 65536 


struct poly {
    int* coeff;
    int taille;
};
typedef struct poly poly;

void print_poly (poly p) {
    printf("polynôme de taille %d :", p.taille);
    for (int i=0; i<p.taille;i++) {
        if (i%10 == 0) {
            printf("\n%d : ", i);
        }
        printf("%d, ", p.coeff[i]);
    }
    printf("\n");
}

poly init_poly(int taille){
    poly p;
    p.taille = taille;
    if (taille != 0) {
        p.coeff = (int*)malloc(taille*sizeof(int));
        for (int i=0; i<taille;i++) {
            p.coeff[i] = 0;
        }
    }
    else{
        p.coeff = NULL;
    }
    return p;
}

int max (int a, int b) {
    if (a>b) {
        return a;
    }
    return b;
}

void free_poly(poly p1){
    if (p1.taille != -1) {
        free(p1.coeff);
    }
}

void mul_poly (poly p1, poly p2) {
    // modifie p1 pour qu'il contienne p1*p2
    int* prod = (int*)malloc(p1.taille*sizeof(int));
    for (int i=0; i<p1.taille; i++) {
        prod[i] = 0;
    }
    for (int i=0; i<p1.taille;i++) {
        for (int j=0; j<p2.taille; j++) {
            if(i+j<p1.taille) {
                prod[i+j] = add (prod[i+j] , mul(p1.coeff[i], p2.coeff[j]));
            }
            assert(i+j<p1.taille || mul(p1.coeff[i], p2.coeff[j])== 0);
        }
    }
    for (int i = 0; i<p1.taille; i++) {
        p1.coeff[i] = prod[i];
    }
    free(prod);
}

int evalue(poly p, int x) {
    int somme = 0; 
    for (int i=0; i<p.taille; i++) {
        somme = add(somme, mul(p.coeff[i], puiss(x, i)));
    }
    return somme;
}

poly mul_xi (poly p0, int coefficient, int indice){
    poly p = init_poly(p0.taille+indice);
    for (int i=0; i<p0.taille; i++){
        p.coeff[i+indice] = mul(coefficient, p0.coeff[i]);
    }
    return p;
}

poly add_poly (poly p1, poly p2){
    poly p = init_poly(max(p1.taille, p2.taille));
    for (int i=0; i<p1.taille; i++) {
        p.coeff[i] = p1.coeff[i];
    }
    for (int i=0; i<p2.taille; i++) {
        p.coeff[i] = add(p.coeff[i], p2.coeff[i]);
    }
    return p;
}

poly copy(poly p0) {
    poly p = init_poly(p0.taille);
    for (int i=0; i<p0.taille; i++) {
        p.coeff[i] = p0.coeff[i];
    }
    return p;
}

int div_euc(poly p1, poly p2, poly buffer[2]) {
    // modifie buffer avec le quotient et le reste de la division euclidienne
    int deg2 = p2.taille -1;
    for (int i=deg2; i>=0 && p2.coeff[i]==0; i--) {
        deg2 -=1;
    }
    assert(deg2 != -1);
    poly quotient = init_poly(p1.taille-deg2);
    poly reste = copy(p1);
    int deg1 = p1.taille -1;
    while(deg1 >= deg2) {
        if (reste.coeff[deg1] != 0) {
            int diff = deg1 - deg2 ;
            int facteur = mul(reste.coeff[deg1], inv(p2.coeff[deg2])) ;
            quotient.coeff[diff] = facteur;
            poly temp = mul_xi(p2, facteur, diff);
            poly somme = add_poly(temp, reste);
            free_poly(temp);
            free_poly(reste);
            reste = somme;
        }
        deg1 -=1;
    }
    buffer[0] = quotient;
    buffer[1] = reste;
}

poly poly_prod (poly p1, poly p2) {
    poly p = init_poly(p1.taille+p2.taille-1);
    for (int i = 0; i<p1.taille; i++) {
        for (int j = 0; j<p2.taille; j++) {
            p.coeff[i+j] = add(p.coeff[i+j], mul(p1.coeff[i], p2.coeff[j]));
        }
    }
    return p;
}

void euclide (poly r0, poly r1, int max_deg, poly buffer[2]){
    r0 = copy(r0);
    r1 = copy(r1);
    poly p0 = init_poly(0);
    poly p1 = init_poly(1);
    p1.coeff[0]=1;
    poly div[2];
    poly temp;
    poly p;
    int deg = r1.taille-1;
    for (int i=deg; i>=0 && r1.coeff[i]==0; i--) {
        deg -=1;
    }
    while (deg>=max_deg){
        div_euc(r0, r1, div);
        temp = poly_prod(div[0],p1);
        p = add_poly(temp,p0);
        free_poly(temp);
        free_poly(p0);
        free_poly(div[0]);
        p0 = p1;
        p1 = p;
        free_poly(r0);
        r0 = r1;
        r1 = div[1];
        for (int i=deg; i>=0 && r1.coeff[i]==0; i--) {
            deg -=1;
        }
    }
    buffer[0] = p1;
    buffer[1] = r1;
    free_poly(p0);
    free_poly(r0);
}

poly derive(poly p) {
    // ici on modifie p, puis on renvoie le nouveau p afin de pouvoir
    // utiliser cette fonction dans des calculs
    // (on ne cree pas de polynome pour ne pas avoir a le free)
    if (p.taille == 0){
        return init_poly(0);
    }
    poly prime = init_poly(p.taille -1);
    for (int i=0; i<prime.taille;i++){
        prime.coeff[i] = (int)(i%2==0)*p.coeff[i+1];
    }
    return prime;
}

void* f_thread(void* arg) {
    poly* args = (poly*)arg;
    poly localisation = args[0];
    int i = args[1].taille;
    int * retour = args[1].coeff;
    *retour = evalue(localisation, i);
}

poly construit_g(int t) {
    int taille = 2*t+1;
    poly g = init_poly(taille);
    g.coeff[0] = 1;
    poly prod = init_poly(2);
    prod.coeff[1] = 1;
    for (int i=1; i<taille; i++) {
        prod.coeff[0]=puiss(alpha, i);
        mul_poly(g, prod);
    }
    free_poly(prod);
    return g;
}

poly decodage (poly M, int t) {
    // détermination du syndrôme et du polynôme de localisation des erreurs
    poly r0 = init_poly(2*t+1);
    r0.coeff[2*t] = 1;
    poly syndrome = init_poly(2*t);
    for (int i=0; i<2*t; i++) {
        syndrome.coeff[i] = evalue(M, puiss(alpha,i+1));
    }
    poly retour[2];
    euclide (r0, syndrome, t, retour);
    if (retour[0].coeff[0] == 0){
        free_poly(r0);
        free_poly(syndrome);
        free_poly(retour[0]);
        free_poly(retour[1]);
        return copy(M); // le décodage échoue
    };
    poly localisation = mul_xi(retour[0], inv(retour[0].coeff[0]),0);
    poly sum = mul_xi(retour[1], inv(retour[0].coeff[0]), 0);
    // printf("localisation a été construit\n");


    // détermination des emplacements des erreurs
    int * erreurs = malloc(localisation.taille*sizeof(int));
    for (int i=0; i<localisation.taille; i++) {
        erreurs[i] = 0;
    }
    int nb_erreurs = 0;
    for (int j=0; j<n-1; j++) {
        if (j%10000 == 1) {
            printf("%d\n",j/10000);
        }
        if (evalue(localisation, puiss(alpha, -j))== 0){
            erreurs[nb_erreurs] = j;
            nb_erreurs += 1;
        }
    }
    if (nb_erreurs == 0) {
        free_poly(r0);
        free_poly(syndrome);
        free_poly(retour[0]);
        free_poly(retour[1]);
        free_poly(localisation);
        free_poly(sum);
        return copy(M);
    }
    // int i = 0;
    // // printf("erreurs détectées (%d) : \n", nb_erreurs);
    // for (int j=0; j<n-1; j++) {
    //     if (evalue(localisation, puiss(alpha, -j)) == 0) {
    //         erreurs[i] = j;
    //         // printf("%d ",j);
    //         i += 1;
    //     }
    // }
    // // printf("\n");
    // printf("erreurs ont été trouvees\n");


    // détermination de l'amplitude des erreurs et correction du message
    poly M_corrige = copy(M);
    int epsilon_j;
    poly Uprime = derive(localisation);
    for(int i=0; i<nb_erreurs; i++) {
        int j = erreurs[i];
        int nom = evalue(sum,puiss(alpha, -j));
        int den = evalue(Uprime, puiss(alpha, -j));
        if (den == 0) {
            // printf("div par 0 : %d\n",j);
            epsilon_j = 0;
        }
        else {
            epsilon_j = mul(nom, inv(den));
        } 
        M_corrige.coeff[j] = add(M_corrige.coeff[j], epsilon_j);
    }
    
    free_poly(r0);
    free_poly(syndrome);
    free_poly(retour[0]);
    free_poly(retour[1]);
    free_poly(localisation);
    free_poly(sum);
    free(erreurs);
    free_poly(Uprime);
    return M_corrige;
}



int main() {
    init_tab_alpha();
    int t=1000;
    poly message_original = init_poly(n-2*t-1); 
    for (int i=0; i<message_original.taille;i++) { 
        message_original.coeff[i] = rand()%n;
    }

    // encodage
    poly g0 = construit_g(t); 
    poly message = poly_prod(message_original,g0); 

    // introduction des erreurs 

    // decodage
    poly decode = decodage(message, t);
    poly div[2];
    div_euc(decode, g0, div);
    // le message décodé est div[0]
}
