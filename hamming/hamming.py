from random import randint

def concat (paquets) : 
    texte = [0 for _ in range(len(paquets)*len(paquets[0]))]
    for i in range(len(paquets)) :
        for j in range(len(paquets[0])) :
            texte[i*len(paquets[0])+j] = paquets[i][j]
    return texte


# Encodage

def scinde_paquets_encode(texte, nb_cor) :
    taille_paquet = (1<<nb_cor)
    paquets = []
    i = 0 
    while (i < len(texte)) :
        paquet_actuel = [0]
        puiss = 1
        for j in range(1, taille_paquet) :
            if j==puiss or i>=len(texte) : # On ajoute les bits de parité et on complète par des 0 à la fin
                paquet_actuel.append(0)
                puiss = puiss << 1
            else :
                paquet_actuel.append(texte[i])
                i +=1
        paquets.append(paquet_actuel)
    return paquets


def hamming_code (l,nb_cor) :
    """l de taille 64"""
    q = l.copy()
    x = 0
    tot = 0
    for k in range (1, 1<<(nb_cor)): 
        if l[k] == 1 :
            x = x ^ k               
            tot +=1 
    for i in range (nb_cor) :
        b = 1<<i 
        if x&b != 0 :
            q[b] = 1
            tot += 1                  
    q[0] = tot%2                   
    return q


def encodage (texte, nb_cor) : 
    """scinde en plusieurs paquets puis encode avec hamming"""
    paquets = scinde_paquets_encode(texte,  nb_cor=nb_cor)
    for i in range(len(paquets)) : 
        paquets [i] = hamming_code (paquets[i],nb_cor)
    return concat(paquets)



# Décodage

def scinde_paquets_decode (texte, nb_cor) : 
    """après lecture du qr code, on a créé une liste contenant chacun des qr codes sous forme de liste.
    on les convertit en une liste de paquets contenant nb_mat matrices de correction d'erreur,
    chacune comprenant nb_cor+1 bits de correction d'erreur (soit de taille 1<<nb_cor)"""
    def accede(j) :
        if j < len(texte) : 
            return texte[j]
        else :
            return 0
    taille_paquet = 1<<nb_cor
    paquets = [] 
    i=0
    while i < len(texte) :
        paquets.append([accede(j) for j in range(i, i+taille_paquet)])
        i += taille_paquet
    return paquets 


def retire_bits_de_parite(paquet) :
    """on a corrigé les erreurs, donc on peut retirer les bits de parité"""
    q = []
    puiss2 = 1
    for i in range (1,len(paquet)) :
        if i == puiss2 : # les bits de parité sont des puissances de 2
            puiss2 = puiss2<<1
        else :
            q.append(paquet[i])
    return q


def hamming_decode(l, nb_cor) : 
        """L'erreur, si elle existe et si elle est seule,
        se trouve aux coordonnées fournies par le xor de tous les bits impairs
        (les bits de parités sont là pour s'assurer ce résultat,
        via une méthode dichotomique sur les lignes et les colonnes)"""
        q = l.copy()
        xor = 0
        tot = 0
        for k in range (1, 1<<(nb_cor)):
            if l[k] == 1 :
                xor = xor ^ k
                tot += 1
        if tot % 2 == l[0] and xor != 0:
            # print ("Au moins deux erreurs")
            pass
        else :
            q[xor] = 1-q[xor]
        return q

def decodage (texte, nb_cor) :
    paquets = scinde_paquets_decode (texte, nb_cor)
    for i in range(len(paquets)) : 
        paquets[i] = hamming_decode(paquets[i], nb_cor)
        paquets[i] = retire_bits_de_parite(paquets[i])
    return concat(paquets)


# exemple de fonction de test
def test(nb_cor, nb_erreurs, taille_texte) :
    texte = [randint(0,1) for i in range(taille_texte)]

    encode = encodage(texte, nb_cor)    

    # simulation d'erreur
    for i in range (nb_erreurs) : 
        x = randint(0, len(encode)-1)
        encode [x] = 1 - encode[x]

    result = decodage(encode, nb_cor)


    # comparaison avec le texte d'origine
    difference = 0
    for i in range(min(len(texte), len(result))) :
        if result[i] != texte[i] :
            difference += 1
    return difference
