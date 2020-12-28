// Author: Alexis GADONNEIX

#include "poisson.h"


const byte OBSCUR=50; // Seuil d'intensite des pixels dont on change le gradient
const float FACTEUR=3.0f; // Facteur de multiplication du gradient

// Rehausse le gradient des pixels sombres, et montre en rouge ces pixels.
void masque(const Image<float>& I, Image<float>& Vx, Image<float>& Vy) {

    for (int x=0;x<I.width();x++){
        for (int y=0;y<I.height();y++){
            if (I(x,y)<OBSCUR){
                Vx(x,y) = Vx(x,y)*FACTEUR; // On augmente le gradient dans les zones sombres
                Vy(x,y) = Vy(x,y)*FACTEUR;
                drawPoint(x,y,RED); // On dessine ces points sombres
            }
        }
    }

}

int main(int argc, char* argv[]) {
    Image<byte> I;
    if(! load(I, argc>1? argv[1]: srcPath("salon.png"))) {
        std::cout << "Echec de lecture d'image" << std::endl;
        return 1;
    }


    Window w1 = openWindow(I.width(), I.height(),"image originale",400,100);
    setActiveWindow(w1);
    display(I);
    anyClick();
    std::cout << "Contraste simple" << std::endl;
    Window w3 = openWindow(I.width(), I.height(),"Contraste simple",410+I.width(),100);
    setActiveWindow(w3);
    affiche(I);
    anyClick();

    Image<float> Vx,Vy;

    Image<float> Vyf = dy(I);
    Image<float> Vxf = dx(I);

    Window w2 = openWindow(I.width(), I.height(),"Dérivées Fourier",410+I.width(),130+I.height());
    Window w4 = openWindow(I.width(), I.height(),"Différences finies",400,130+I.height());

    gradient(I,Vx,Vy);

    setActiveWindow(w4);
    affiche(Vx);
    setActiveWindow(w2);
    affiche(Vxf);
    anyClick();

    setActiveWindow(w4);
    affiche(Vy);
    setActiveWindow(w2);
    affiche(Vyf);
    anyClick();

    closeWindow(w2);
    closeWindow(w4);
    Window w6 = openWindow(I.width(), I.height(),"Masque",400,130+I.height());
    setActiveWindow(w6);

    display(I);
    masque(I,Vxf,Vyf);

    anyClick();


    Window w5 = openWindow(I.width(), I.height(),"Contraste Poisson",410+I.width(),130+I.height());
    setActiveWindow(w5);

    affiche(poisson(Vxf,Vyf));

    anyClick();



    endGraphics();
    return 0;
}
