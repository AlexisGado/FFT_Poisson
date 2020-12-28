#include "fft.h"

// Discrete Fourier transform a l'index k du signal f de longueur n.
// s=-1 pour DFT directe, s=+1 pour DFT inverse.
std::complex<float> dft(const std::complex<float> f[], int n, int k, float s) {
    std::complex<float> S = 0;
    for (int j=0;j<n;j++){
        S+=f[j]*std::polar(float(1),float(s*2*j*k*M_PI/n)); //On calcule de manière directe le keme coeff de la DFT de f
    }

    return S/float(std::sqrt(n));
}

// Fast Fourier transform du signal f(deb:pas:fin) (Matlab notation).
// s=-1 pour DFT directe, s=+1 pour DFT inverse.
// Buffer est utilise comme tableau temporaire, sa longueur doit etre au moins
// celle de f.
void fft_main(std::complex<float> f[], int deb, int pas, int fin, float s,
              std::complex<float> buffer[]) {   // Attention ca ne fonctionne que pour des tailles de tableau puissances de 2

    int n = (fin-deb)/pas + 1; //taille du tab

    if (n==1){ // Condition d'arret : si il n y a plus qu'un element dans le tableau considere
        return;
    }

    fft_main(f,deb,2*pas,fin,s,buffer); // On appelle sur les indices pairs
    fft_main(f,deb+pas,2*pas,fin,s,buffer); // et impairs



    std::complex<float> r = std::polar(float(1),float(s*2*M_PI/n)); // r est la raison de la suite (tk)
    std::complex<float> t = 1; //premier terme de la suite

    for (int k=0;k<n;k++){
        buffer[k]=f[deb + k*pas];  // On remplit le buffer
    }

    for (int k=0;k<n/2;k++){  // On recombine les termes selon l'algo de la fft
        f[deb+k*pas] = buffer[2*k] + t*buffer[2*k+1];
        f[deb+(k+n/2)*pas] = buffer[2*k] - t*buffer[2*k+1];

        t*=r; // On passe au terme suivant de la suite
    }
}

// Divise tous les coefficients de f par la racine carree de sa longueur n.
void normalize(std::complex<float> f[], int n, float div) {
    for(int i=0; i<n; i++)
        f[i] /= div;
}

// FFT du signal f de longueur n.
void fft(std::complex<float> f[], int n) {

    std::complex<float>* buffer = new std::complex<float>[n];


    fft_main(f,0,1,n-1,-1,buffer); // On appelle simplement la fft avec les bonnes bornes et un pas de 1

    delete[] buffer;
    normalize(f,n,float(std::sqrt(n))); // on oublie pas de normaliser a la fin

}

// FFT inverse du signal f de longueur n.
void ifft(std::complex<float> f[], int n) {
    std::complex<float>* buffer = new std::complex<float>[n];

    fft_main(f,0,1,n-1,1,buffer); // idem avec s = 1 pour la fft inverse

    delete[] buffer;
    normalize(f,n,float(std::sqrt(n)));
}

// FFT du signal 2D f de dimension wxh.
void fft2(std::complex<float> f[], int w, int h) {

    std::complex<float>* buffer;

    for (int i=0;i<h;i++){ // Les lignes
        buffer = new std::complex<float>[w];
        fft_main(f,i*w,1,(i+1)*w-1,-1,buffer); // On ajuste le début et la fin pour couvrir une ligne
        delete[] buffer;
    }

    for (int j=0;j<w;j++){ // puis les colonnes
        buffer = new std::complex<float>[h];
        fft_main(f,j,w,w*(h-1)+j,-1,buffer); // Cette fois on fait aussi attention au pas !
        delete[] buffer;
    }
    normalize(f,w*h,float(std::sqrt(w*h))); // et on normalise

}

// FFT inverse du signal 2D f de dimentsion wxh.
void ifft2(std::complex<float> f[], int w, int h) {

    std::complex<float>* buffer;

    for (int i=0;i<h;i++){ // Les lignes
        buffer = new std::complex<float>[w];
        fft_main(f,i*w,1,(i+1)*w-1,1,buffer); // On ajuste le début et la fin pour couvrir une ligne
        delete[] buffer;
    }

    for (int j=0;j<w;j++){ // puis les colonnes
        buffer = new std::complex<float>[h];
        fft_main(f,j,w,w*(h-1)+j,1,buffer); // Cette fois on fait aussi attention au pas !
        delete[] buffer;
    }
    normalize(f,w*h,float(std::sqrt(w*h)));
}
