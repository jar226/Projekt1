#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <stdio.h>
#include <cmath>


#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL
#define LM 0x7FFFFFFFULL


using namespace std;

static unsigned long long mt[NN];
static int mti=NN+1;

void init_genrand64(unsigned long long seed)
{
    mt[0] = seed;
    for (mti=1; mti<NN; mti++)
        mt[mti] = (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);

}
unsigned long long genrand64_int64(void)
{
    int i;
    unsigned long long x;
    static unsigned long long mag01[2]={0ULL, MATRIX_A};

    if (mti >= NN) {
        if (mti == NN+1)
            init_genrand64(5489ULL);

        for (i=0;i<NN-MM;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        for (;i<NN-1;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        x = (mt[NN-1]&UM)|(mt[0]&LM);
        mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        mti = 0;
    }

    x = mt[mti++];

    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);

    return x;
}
double genrand64_real1(void)
{
    return (genrand64_int64() >> 11) * (1.0/9007199254740991.0);
}
void GenerowanieTab(int N, double *tab)//funkcja generujaca liczy z rok�adu jednorodnego
{
    for(int i=0;i<N;i++)
    {
        tab[i]=genrand64_real1();
    }
}
void GenerowanieTabRozWykladniczy(int N, double *tab)//funkcja generujaca liczy z rozk�adu wyk�adniczego
{
    double lambda=2.0;

    for(int i=0;i<N;i++)
    {
        tab[i]=-(1.0/lambda)*(log(genrand64_real1()));
    }
}
void GenerowanieTabRozNormalny(int N, double *tab)//funkcja generujaca liczy z rozk�adu normalnego
{
    double R,theta;
    double u1,u2;

    for(int i=0;i<N;i++)
    {
        u1=genrand64_real1();
        u2=genrand64_real1();

        R=pow((-2.0*log(u1)),0.5);
        theta=2.0*M_PI*u2;

        tab[i]=R*sin(theta);
    }
}
void GenerowanieTabRozArcusaSinusa(int N, double *tab)//funkcja generujaca liczy a rozk�adu arcusasinusa
{
    double u;
    for(int i=0;i<N;i++)
    {
        u=genrand64_real1();
        tab[i]=pow(sin(M_PI*u/2.0),2);
    }
}

void WyswietlTab(double *tab, int N)
{
    for(int i=0;i<N;i++)
    {
        cout.setf(ios::fixed);
        cout<<tab[i]<<endl;
    }
}
void Swap(double *tab, int a, int b)
{
    double pom;
    pom=tab[a];
    tab[a]=tab[b];
    tab[b]=pom;
}
int median ( double *tab, int left, int right, int &lporownan, int &lzamian)
{
    int mid = ( left + right ) / 2;

    if(tab[right] < tab[left])
	{
        Swap(tab, left, right);
        lzamian++;
    }
    if(tab[mid] < tab[left])
	{
        Swap(tab, mid, left);
        lzamian++;
	}
    if(tab[right] < tab[mid])
	{
        Swap(tab, right, mid);
        lzamian++;
	}
	lporownan+=3;
    return mid;
}
void BubbleSort(double *tab,int p, int N, int &lporownan, int &lzamian)
{
    int zamiana=0;
    int i=0;
    do
    {
        zamiana=0;
        for(int j=p;j<N-1-i;j++)
        {
            if(tab[j]>tab[j+1])
            {
                Swap(tab,j,j+1);
                zamiana++;
                lzamian++;
            }
            lporownan++;
        }
    }while(zamiana>0);
}
void QuickSort(double *tab,int lewy, int prawy, int &lporownan, int &lzamian)
{
    int i,j;
    double piwot;

    int med = median(tab,lewy,prawy,lporownan,lzamian);

    Swap(tab,med,prawy);
    lzamian++;
    piwot=tab[prawy];
    i=lewy;

    for(int j=lewy; j<prawy; j++)
    {
        if(tab[j] < piwot)
        {
            Swap(tab,i,j);
            i++;
            lzamian++;
        }
            lporownan++;
    }

    Swap(tab,prawy,i);
    j=i;
    lzamian++;
    if(lewy < j - 1)
        QuickSort(tab,lewy, j - 1,lporownan,lzamian);
    if(j + 1 < prawy)
        QuickSort(tab,j + 1, prawy,lporownan,lzamian);
}

void HybridSort(double *tab,int lewy, int prawy, int rozmiar_krytyczny, int &lporownan, int &lzamian)
{
    int i,j;
    double piwot;
    int ilosc_el=prawy-lewy+1;

    if(ilosc_el<rozmiar_krytyczny)
    {
        BubbleSort(tab,lewy,lewy+ilosc_el,lporownan,lzamian);
    }
    else
    {
        int med = median(tab,lewy,prawy,lporownan,lzamian);

        Swap(tab,med, prawy);
        lzamian++;
        piwot=tab[prawy];
        i=lewy;

        for(int k=lewy; k<prawy; k++)
        {
            if(tab[k] < piwot)
                {
                    Swap(tab,i,k);
                    i++;
                    lzamian++;
                }
            lporownan++;
        }

        Swap(tab,prawy,i);
        j=i;
        lzamian++;
        if(lewy < j - 1)
        {
            HybridSort(tab,lewy, j - 1, rozmiar_krytyczny,lporownan,lzamian);
            lewy=j+1;
        }
        if(j + 1 < prawy)
        {
            HybridSort(tab,j + 1, prawy, rozmiar_krytyczny,lporownan,lzamian);
            prawy=j-1;
        }
    }
}
void PowtorzQS(int L, int N, double *tab,int lewy, int prawy, int &lporownan, int &lzamian)//funkcja powtarzaj�ca losowanie QS L razy
{
    double sr_porownan=0,sr_zamian=0;
    double kw_porownan=0,kw_zamian=0;
    ofstream plik;
    plik.open("QuickSort.txt",ios::app);
    plik<<N<<" ";

    for(int i=0;i<L;i++)
    {
        lporownan=0, lzamian=0;
        GenerowanieTab(N,tab);
        //GenerowanieTabRozWykladniczy(N,tab);
        //GenerowanieTabRozNormalny(N,tab);
        //GenerowanieTabRozArcusaSinusa(N,tab);

        QuickSort(tab, lewy, prawy, lporownan, lzamian);

        sr_porownan=sr_porownan+lporownan;
        sr_zamian=sr_zamian+lzamian;
        kw_porownan=kw_porownan+(pow(lporownan,2));
        kw_zamian=kw_zamian+(pow(lzamian,2));
    }
    plik.setf(ios::fixed);
    plik.precision(2);
    plik<<sr_porownan/L<<" ";
    plik<<sr_zamian/L<<" ";
    plik<<kw_porownan/L<<" ";
    plik<<kw_zamian/L<<endl;
    plik.close();
}
void PowtorzBS(int L, int N, double *tab, int &lporownan, int &lzamian)//funkcja powtarzaj�ca losowanie BS L razy
{
    ofstream plik;
    plik.open("BubbleSort.txt",ios::app);
    plik<<N<<" ";
    double sr_porownan=0,sr_zamian=0;
    double kw_porownan=0,kw_zamian=0;

    for(int i=0;i<L;i++)
    {
        lporownan=0, lzamian=0;
        GenerowanieTab(N,tab);
        BubbleSort(tab,0,N,lporownan,lzamian);

        sr_porownan=sr_porownan+lporownan;
        sr_zamian=sr_zamian+lzamian;
        kw_porownan=kw_porownan+(pow(lporownan,2));
        kw_zamian=kw_zamian+(pow(lzamian,2));
    }
    plik.precision(2);
    plik.setf(ios::fixed);
    plik<<sr_porownan/L<<" ";
    plik<<sr_zamian/L<<" ";
    plik<<kw_porownan/L<<" ";
    plik<<kw_zamian/L<<endl;
    plik.close();
}

void PowtorzHS(int L, int N, double *tab,int lewy, int prawy, int rozmiar_krytyczny, int &lporownan, int &lzamian)//funkcja powtarzaj�ca losowanie HS L razy
{
    ofstream plik;
    plik.open("HybridSort.txt",ios::app);
    plik<<N<<" ";
    double sr_porownan=0,sr_zamian=0;
    double kw_porownan=0,kw_zamian=0;

    for(int i=0;i<L;i++)
    {
        lporownan=0, lzamian=0;
        GenerowanieTab(N,tab);
        //GenerowanieTabRozWykladniczy(N,tab);
       // GenerowanieTabRozNormalny(N,tab);
       //GenerowanieTabRozArcusaSinusa(N,tab);

        HybridSort(tab, lewy, prawy, rozmiar_krytyczny,lporownan, lzamian);

        sr_porownan=sr_porownan+lporownan;
        sr_zamian=sr_zamian+lzamian;
        kw_porownan=kw_porownan+(pow(lporownan,2));
        kw_zamian=kw_zamian+(pow(lzamian,2));
    }
    plik.precision(2);
    plik.setf(ios::fixed);
    plik<<sr_porownan/L<<" ";
    plik<<sr_zamian/L<<" ";
    plik<<kw_porownan/L<<" ";
    plik<<kw_zamian/L<<" ";
    plik<<rozmiar_krytyczny<<endl;
    plik.close();
}
void fQS(int L, int listaN[], double *tab, int &lporownan, int &lzamian)//funkcja losuj�ca QS dla r�nych wielko�ci tablic z atrybutu listaN
{
    ofstream plik;
    plik.open("QuickSort.txt");
    for(int i=0;i<28;i++)
    {
        PowtorzQS(L,listaN[i],tab,0,listaN[i]-1,lporownan,lzamian);
    }
}
void fBS(int L, int listaN[], double *tab, int &lporownan, int &lzamian)//funkcja losuj�ca BS dla r�nych wielko�ci tablic z atrybutu listaN
{
    ofstream plik;
    plik.open("BubbleSort.txt");
    for(int i=0;i<28;i++)
    {
        PowtorzBS(L,listaN[i],tab,lporownan,lzamian);
    }
}
void fHS(int L, int listaN[], double *tab, int &lporownan, int &lzamian)//funkcja losuj�ca HS dla r�nych wielko�ci tablic z atrybutu listaN
{
    ofstream plik;
    plik.open("HybridSort.txt");
    for(int i=0;i<28;i++)
    {
        PowtorzHS(L,listaN[i],tab,0,listaN[i]-1,7,lporownan,lzamian);
    }
}
void fHS2(int L, int listaN[], double *tab, int &lporownan, int &lzamian)//funkcja powtarzaj�ca losowania HS dla rozmair�w krytycznych od 2 do 50
{
    ofstream plik;
    plik.open("HybridSort.txt");
    for(int i=0;i<28;i++)
    {
        for(int j=2;j<=50;j++)
        {
            PowtorzHS(L,listaN[i],tab,0,listaN[i],j,lporownan,lzamian);
        }
    }
}
void DaneHistogram(int L, int N, double *tab,int lewy, int prawy, int rozmiar_krytyczny, int &lporownan, int &lzamian)//funkcja zabieraj�ca dane z losowa� do stworzenia histogramu
{
    ofstream plik;
    plik.open("DaneHistogramHS10000.txt");

    for(int i=0;i<L;i++)
    {
        lporownan=0, lzamian=0;
        //GenerowanieTab(N,tab);
        //GenerowanieTabRozWykladniczy(N,tab);
        //GenerowanieTabRozNormalny(N,tab);
        GenerowanieTabRozArcusaSinusa(N,tab);

        HybridSort(tab, lewy, prawy-1, rozmiar_krytyczny,lporownan, lzamian);
        //QuickSort(tab, lewy, prawy-1, lporownan, lzamian);

        plik<<N<<" ";
        plik.precision(2);
        plik.setf(ios::fixed);
        plik<<lporownan<<" ";
        plik<<lzamian<<" "<<endl;

    }
plik.close();
}

int main()
{
    int N=10000;
    int L=10000;
    double *tablica=new double[N];
    int lporownan=0, lzamian=0;
    init_genrand64(time(NULL));
    //int listaN[28]={100,200,500,1000,2000,5000,10000};
    int listaN[28]={10,20,30,40,50,60,70,80,90,
                    100,200,300,400,500,600,700,800,900,
                    1000,2000,3000,4000,5000,6000,7000,8000,9000,10000};

    //fQS(L,listaN,tablica,lporownan,lzamian);
    //fHS(L,listaN,tablica,lporownan,lzamian);
    //fHS2(L,listaN,tablica,lporownan,lzamian);
    //fBS(L,listaN,tablica,lporownan,lzamian);
     N=10000;
     DaneHistogram(L,N,tablica,0,N,7,lporownan,lzamian);

    //kolega
    delete[] tablica;
    return 0;
}
