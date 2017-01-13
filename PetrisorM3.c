#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mersenne.c"
#include <time.h>

typedef struct NodR
{
    int in; // numarul de inlinkuri spre nod
    int out; // numarul de outlinkuri din nod
    int nrviz; // numarul de vizite ale paginii(nodului)
    double PgRkE; // pagerank-ul experimental
    double PgRkA; // pagerank-ul paginii calculat conform algoritmului
} NodR;

void eroare()
{
    printf( "Eroare citire din fisierul sursa.\n" );
    exit(1);
}


double functie(double *p, int k)  //calculeaza functia de repartitie in functie de k, returnat de functia simVarDiscr
{
    int i;
    double s=0;
    for(i=0; i<=k; i++)
        s=s+p[i];
    return s;
}

// functie care simuleaza distributia discreta de probabilitate
int simVarDiscr(int m, double *p)
{
    int k=0;
    double u=genrand_real2();
    while(u>functie(p,k)&&k<m-1)
        k++;
    return k;

}

void LantMarkov( NodR *pag, int m, double *pi0, double **G, int n, int *s )
{
    int i, k,j, aux;
    double *p;
    s=(int*)malloc(n*sizeof(int));
    if(!s)
        eroare();
    pi0=(double*)calloc(m,sizeof(double));
    if(!pi0)
        eroare();
    p=(double*)calloc(n,sizeof(double));
    if(!p)
        eroare();
    for ( j = 0; j < m; j++ )
        pi0[j] = (double)1/m;
    s[0]=simVarDiscr(m,pi0);
    i = s[0];
     for(k=1;k<=n;k++)
    {
        pag[i].nrviz++;
        aux =0;
        p=&G[i][0];
        s[k]=simVarDiscr(m,p);
        i=s[k];


    }
    for ( i = 0; i < m; i++ )
        pag[i].PgRkE = (double)pag[i].nrviz/n;

}

int sortare_descr( int m, NodR *pag, int *v)
{
    int i;
    int j,ok;

    NodR aux;
    int aux_nr;
    for(i=0; i<m; i++)
        v[i]=i;

do{
    ok=1;
    for (i=0; i<m-1; i++)
    if(pag[i].PgRkA<pag[i+1].PgRkA)
       {

           aux=pag[i];
           pag[i]=pag[i+1];
           pag[i+1]=aux;

           aux_nr=v[i];
           v[i]=v[i+1];
           v[i+1]=aux_nr;

           ok=0;
        }

}while(ok==0);
  return v;
}


void afisare_fisier(NodR *pag, int *v, int m)
{
    FILE *po;
    int i;

    po = fopen( "SimulPgrank.out", "w" );
    if ( po == NULL )
    {
        printf( "Eroare citire din fisierul sursa.\n" );
        exit(1);
    }
    printf ( "Pagina \t Nr.inlinkuri \t Nr.outlinkuri \t PageRank Algo \t PageRankExper \n" );
    for(i = 0; i < m; i++ )
        printf("%2d \t\t %d \t\t %d \t\t %.2lf \t\t %.2lf \n",v[i], pag[i].in,pag[i].out, pag[i].PgRkA, pag[i].PgRkE );

    fprintf( po, "Pagina \t Nr.inlinkuri \t Nr.outlinkuri \t PageRank Algo \t PageRankExper \n" );
    for(i = 0; i < m; i++ )
        fprintf( po, "%2d \t\t\t %d \t\t\t\t %d \t\t\t\t %.2lf \t\t\t %.2lf \n", v[i], pag[i].in,pag[i].out, pag[i].PgRkA, pag[i].PgRkE );
    fclose(po);
}

int main()
{
    double **H, **NQ, **Q, **G, *pi0;
    int m, n, i, j, *s, *xx, *v;
    double alfa = 0.85;
    n = 1000;
    NodR *pag = NULL;

    time_t secunde;
    secunde=time(NULL);
    init_genrand(secunde);


//Citirea matricii H din fisier

    FILE *pf;


    pf = fopen( "hyperlink.in", "r" );
    if ( pf == NULL )
    {
        printf( "Eroare deschidere fisier sursa.\n" );
        exit(1);
    }
    do
    {
        fscanf( pf, "%d", &m );
    }
    while( m < 10 || m > 20 );

    H = (double**) malloc( m*sizeof(double*) );
    for ( i = 0; i < m; i++ )
        H[i] = (double*)calloc( m, sizeof(double) );
    if(!H)
        eroare();
    for ( i = 0; i < m; i++ )
        for ( j = 0; j <m; j++ )
            fscanf( pf, "%lf", &H[i][j] );

    fclose(pf);

//Deducerea PageRank-ului experimental
//Prelucrarea nodurilor si construirea matricilor Q, stochastica si Google

    //alocam dinamic memorie pentru pointerul de tip NodR
    pag = (NodR*) malloc( m*sizeof(NodR) );
    if(!pag)
        eroare();
    for ( i = 0; i < m; i++ )
        pag[i].in = pag[i].out = pag[i].nrviz = 0;

    // Calculul numarului de inlink-uri ale paginii cu codul i
    for ( i = 0; i < m; i++ )
        for ( j = 0; j < m; j++ )
            if ( H[j][i] == 1 )
                pag[i].in++;

    // Calculul numarului de outlink-uri ale paginii cu codul i
    for ( i = 0; i < m; i++ )
        for ( j = 0; j < m; j++ )
            if ( H[i][j] == 1 )
                pag[i].out++;

    // Constructia matricii Q
    Q = (double**) malloc( m*sizeof(double*) );
    for ( i = 0; i < m; i++ )
        Q[i] = (double*)calloc( m, sizeof(double) );
    if(!Q)
        eroare();
    for ( i = 0; i < m; i++ )
    {
        for ( j = 0; j < m; j++ )
            if ( pag[i].out != 0 )
                Q[i][j] = H[i][j]/pag[i].out;
            else
                Q[i][j] = 0;
    }

    // Constructia matricii stochastica
    NQ = (double**) malloc( m*sizeof(double*) );
    for ( i = 0; i < m; i++ )
        NQ[i] = (double*)calloc( m, sizeof(double) );
    if(!NQ)
        eroare();
    for ( i = 0; i < m; i++ )
    {
        for ( j = 0; j < m; j++ )
            if ( pag[i].out != 0 )
                NQ[i][j] = H[i][j]/pag[i].out;
            else
                NQ[i][j] = (double)1/m;
    }

    // Constructia matricii Google
    G = (double**) malloc( m*sizeof(double*) );
    for ( i = 0; i < m; i++ )
        G[i] = (double*)calloc( m, sizeof(double) );
    if(!G)
        eroare();
    for ( i = 0; i < m; i++ )
    {
        for ( j = 0; j < m; j++ )
            G[i][j] = alfa * NQ[i][j] + (double)(1 - alfa) / m;
    }



    LantMarkov(pag, m,pi0,G,n,s);


//Calculul PageRank-ului algoritmic

    double *pi, *pip;
    double dist;
    double eps=1e-5;

    int nr = 0; // numarul de iteratii necesare pana la convergenta metodei puterii
    v=(int *) malloc(m* sizeof(int));
    pi = (double*)calloc( m, sizeof(double) );
    if(!pi)
        eroare();
    pip = (double*)calloc( m, sizeof(double) );
    if(!pip)
        eroare();
    xx = (int*)calloc( m, sizeof(int) );
    if(!xx)
        eroare();

    for ( i = 0; i < m; i++ )
        xx[i] = i;

    for ( i = 0; i < m; i++ )
        pi[i] = (double)1/m;

    do
    {
        /* pi' <- pi */
        for( i = 0; i < m; i++ )
        {
            pip[i] = pi[i];
            pi[i] = 0;
        }

        /* pi <- pi'*G */
        for ( j = 0; j < m; j++ )
        {
            for ( i = 0; i < m; i++ )
                pi[j] += pip[i]*G[i][j];
        }
        dist = 0;
        for ( i = 0; i < m; i++ )
            dist += (pi[i] - pip[i])*(pi[i] - pip[i]);
        dist =(double) sqrt(dist);
        nr++;
    }
    while ( dist >= eps );
    printf( "Numarul de iteratii pana la convergenta puterii este: %d\n", nr );
    for( i = 0; i < m; i++ )
        pag[i].PgRkA=pi[i];

//Ordonarea celor m pagini in ordine descrescatoare a PageRank-ului lor algoritmic
    sortare_descr(m,pag,v);
//Afisarea pe ecran sub forma de tabel
    afisare_fisier(pag,v,m);

    free(xx);
    free(pag);
    free(pi0);
    free(H);
    free(Q);
    free(NQ);
    free(G);
    for(i=0;i<m;i++)
        free(H[i]);
    for(i=0;i<m;i++)
        free(Q[i]);
    for(i=0;i<m;i++)
        free(NQ[i]);
    for(i=0;i<m;i++)
        free(G[i]);


    return 0;
}
