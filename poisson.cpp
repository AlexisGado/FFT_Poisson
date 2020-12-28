#include "poisson.h"
#include "fft.h"
#include <algorithm>

const float POURCENT=.5f;

// Fonction affine amenant les min et max de F (ignorant les POURCENT valeurs
// extremes de chaque cote) a 0 et 255.
void affineContraste(Image<float> F, float& a, float& b) {
    F=F.clone();
    std::sort(F.begin(), F.end());
    int p = (int)(F.totalSize()*POURCENT/100);
    float min = F[p];
    float max = F[F.totalSize()-1-p];
    a = (min<max)? 255.0f/(max-min): 0;
    b = (min<max)? -a*min: 255.0f/2;
}

// Renvoie une image affichable.
Image<byte> affichable(const Image<float>& F) { // Un contraste affine est applique a chaque fois, necessaire ?
    float a,b;
    affineContraste(F, a, b);
    Image<byte> I(F.width(), F.height());
    for(int i=0; i<I.height(); i++)
        for(int j=0; j<I.width(); j++) {
            float f = a*F(j,i)+b + 0.5f; 
            if(f<0) f=0;
            if(f>255) f=255;
            I(j,i) = (byte)f;
        }
    return I;
}

// Affichage d'une image avec intensites reelles.
void affiche(const Image<float>& F) {
    display(affichable(F));
}

// Prend la partie reelle de chaque pixel d'une image complexe.
Image<float> realImage(const Image< std::complex<float> >& F) {
    const int w=F.width(), h=F.height();
    Image<float> I(w,h);
    for(int i=0; i<h; i++)
        for(int j=0; j<w; j++)
            I(j,i) = F(j,i).real();    
    return I;
}

// Puissance de 2 immediatement superieure a i-1.
int puis2(int i) {
    int j=1;
    while(j<i)
        j *= 2;
    return j;
}

// Genere une image plus grande en rajoutant des 0.
Image< std::complex<float> > agrandis(const Image< std::complex<float> >& I,
                                      int w, int h) {
    if(w==I.width() && h==I.height())
        return I.clone();
    Image< std::complex<float> > I2(w,h);
    I2.fill(0.0f);
    for(int i=0; i<I.height(); i++)
        for(int j=0; j<I.width(); j++)
            I2(j,i) = I(j,i);
    return I2;
}

// Gradient de l'image I par differences finies.
void gradient(const Image<float>& I, Image<float>& Vx, Image<float>& Vy) {
    Vx = Image<float>(I.width(), I.height());
    Vy = Image<float>(I.width(), I.height());
    // A completer
    Coords<2> a;
    for (int x=0;x<I.width();x++){
        a[0] = x;
        for (int y=0;y<I.height();y++){
            a[1] = y;
            Vx(x,y) = gradient(I, a)[0];
            Vy(x,y) = gradient(I, a)[1];  // calcul du gradient en chaque point selon les directions x et y
        }
    }
}

void imageToTab(const Image< std::complex<float> > F,std::complex<float> tab[]){  // Fonctions necessaires car fft2 et ifft2 s'appliquent sur des tableaux...
    for (int x=0;x<F.width();x++){
        for (int y=0;y<F.height();y++){
            tab[x+y*F.width()] = F(x,y);
        }
    }
}

void tabToImage(Image< std::complex<float> > F,std::complex<float> tab[]){
    for (int x=0;x<F.width();x++){
        for (int y=0;y<F.height();y++){
            F(x,y) = tab[x+y*F.width()];
        }
    }
}


// Calcul en Fourier de la derivee suivant x.
void Fourier_dx(Image< std::complex<float> >& F) {

    std::complex<float> i(0,1);

    int w = F.width();
    int h = F.height();
    std::complex<float>* tab = new std::complex<float>[w*h];
    imageToTab(F,tab);
    fft2(tab,w,h);
    tabToImage(F,tab);  // conversion dans un sens et dans l'autre necessaire pour appliquer fft2
    delete[] tab;

    // On derive en utilisant la formule de derivation en Fourier
    for (int l=0;l<h;l++){
        for (int k=0;k<w;k++){
            if (k<float(w/2)){
                F(k,l) = F(k,l)*i*float(2*M_PI*k/w);
            }
            else{
                if (k>float(w/2)){
                    F(k,l) = F(k,l)*i*float(2*M_PI*(k-w)/w);
                }
                else{
                    F(k,l) = 0;
                }
            }
        }
    }

}

// Derivee suivant x par DFT.
Image<float> dx(Image< std::complex<float> > F) {

    //F = F.clone();

    int w = F.width();
    int h = F.height();

    int newW = puis2(w);
    int newH = puis2(h);

    Image< std::complex<float> > F2 = agrandis(F,newW,newH);  //On agrandit pour pouvoir appliquer la fft sur toutes les tailles d'images


    Fourier_dx(F2); //On met dans l image nos gradients en Fourier

    std::complex<float>* tab = new std::complex<float>[newW*newH];

    imageToTab(F2,tab);
    ifft2(tab,newW,newH); // transformation inverse
    tabToImage(F2,tab);

    Image< std::complex<float> > F3 = F2.getSubImage(0,0,w,h); // On reduit notre image

    return realImage(F3); // On repasse en reel pour eliminer les restes imaginaires
}

// Calcul en Fourier de la derivee suivant y.
void Fourier_dy(Image< std::complex<float> >& F) {
    std::complex<float> i(0,1);

    int w = F.width();
    int h = F.height();
    std::complex<float>* tab = new std::complex<float>[w*h];
    imageToTab(F,tab);
    fft2(tab,w,h);
    tabToImage(F,tab);
    delete[] tab;

    // On derive
    for (int k=0;k<w;k++){
        for (int l=0;l<h;l++){
            if (l<float(h/2)){
                F(k,l) = F(k,l)*i*float(2*M_PI*l/h);
            }
            else{
                if (l>float(h/2)){
                    F(k,l) = F(k,l)*i*float(2*M_PI*(l-h)/h);
                }
                else{
                    F(k,l) = 0;
                }
            }
        }
    }
}

// Derivee suivant y par DFT.
Image<float> dy(Image< std::complex<float> > F) {

    //F = F.clone();

    int w = F.width();
    int h = F.height();

    int newW = puis2(w);
    int newH = puis2(h);

    Image< std::complex<float> > F2 = agrandis(F,newW,newH);


    Fourier_dy(F2);

    std::complex<float>* tab = new std::complex<float>[newW*newH];

    imageToTab(F2,tab);
    ifft2(tab,newW,newH);
    tabToImage(F2,tab);

    Image< std::complex<float> > F3 = F2.getSubImage(0,0,w,h);

    return realImage(F3);
}

// Resouds l'equation de Poisson.
Image<float> poisson(Image< std::complex<float> > Vx,
                     Image< std::complex<float> > Vy) {

    int w = Vx.width();
    int h = Vx.height();
    std::complex<float> i(0,1);

    int newW = puis2(w);
    int newH = puis2(h);

    Image< std::complex<float> > u(newW,newH); // on cree des le depart une "grande" image

    Image< std::complex<float> > Vx2 = agrandis(Vx,newW,newH); // On agrandit les gradients
    Image< std::complex<float> > Vy2 = agrandis(Vy,newW,newH);


    Fourier_dx(Vx2);  // On derive une nouvelle fois en Fourier
    Fourier_dy(Vy2);
    int kp,lp;

    for (int k=0;k<newW;k++){
        if (k<newW/2) kp = k;
        else{ if (k>newW/2) kp = k-newW;  //On regle la valeur de k selon si l on est a droite, a gauche ou au milieu
              else kp=0;
        }


        for (int l=0;l<newH;l++){
            if (l<newH/2) lp = l;
            else{ if (l>newH/2) lp = l-newH;  //Idem pour l
                  else lp=0;
            }

            if (kp==0 && lp==0){ // On fait attention a la division par 0
                u(0,0) = 0;
            }
            else{  // Derivee 2D en Fourier
                u(k,l) = (Vx2(k,l) +Vy2(k,l))/((i*float(2*M_PI*lp/newH))*(i*float(2*M_PI*lp/newH)) +(i*float(2*M_PI*kp/newW))*(i*float(2*M_PI*kp/newW)));
            }

        }
    }


    std::complex<float>* tab = new std::complex<float>[newW*newH]; // Attentien aux tailles du tableau temporaire

    imageToTab(u,tab);
    ifft2(tab,newW,newH); // retour en reels
    tabToImage(u,tab);

    delete[] tab;

    Image< std::complex<float> > u2 = u.getSubImage(0,0,w,h); //  reduction

    return realImage(u2);
}
