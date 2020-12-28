// Author: Alexis GADONNEIX

#include "poisson.h"

// Selection a la souris d'un rectangle. Lui passer la fenetre W1.
bool selectionRect(Window W,
                   int& x1, int& y1, int& x2, int& y2, const Image<byte>& I) {
    setActiveWindow(W);
    Event e;
    do {
        getEvent(0, e);
        if(e.button==3) return false;
        if(e.win != W) continue;
        if(e.type==EVT_BUT_ON) {
            x1 = x2 = e.pix[0];
            y1 = y2 = e.pix[1];
        }
        if(e.type==EVT_MOTION) {
            x2 = e.pix[0];
            y2 = e.pix[1];

            display(I);
            drawRect(std::min(x1,x2), std::min(y1,y2), abs(x1-x2), abs(y1-y2),
                     RED);
        }
    } while(e.win!=W || e.type!=EVT_BUT_OFF || abs(x1-x2)<5 || abs(y1-y2)<5);
    if(x1>x2) std::swap(x1,x2);
    if(y1>y2) std::swap(y1,y2);
    return true;
}

// Place interactivement la petite image fg dans l'image bg. Lui passer W2.
void cloneRect(Window W,
               int& xc, int& yc, const Image<byte>& bg, const Image<byte>& fg) {
    setActiveWindow(W);
    Event e;
    do {
        getEvent(0, e);
        if(e.win != W) continue;
        if(e.type==EVT_BUT_ON || e.type==EVT_MOTION) {
            xc = e.pix[0];
            yc = e.pix[1];
            display(bg);
            display(fg, xc, yc);
        } 
    } while(e.win!=W || e.type!=EVT_BUT_OFF);    
}

// Mets dans Vx2, Vy2 le plus fort gradient.
void maxGradient(const Image<float>& Vx1, const Image<float>& Vy1,
                 Image<float>& Vx2, Image<float>& Vy2,
                 int x1, int y1, int x2, int y2, int xc, int yc) {

    int w = abs(x1-x2);
    int h = abs(y1-y2);

    for (int x=0;x<=w;x++){
        for (int y=0;y<=h;y++){
            if (0<=x1+x<Vx1.width() && 0<=y1+y<Vy1.height() && 0<=xc+x<Vx2.width() && 0<=yc+y<Vy2.height()){ // Ces tests rallongent le temps d'execution mais evite les problemes si le collage deborde
                float nor1 = abs(Vx1(x1+x,y1+y))+abs(Vy1(x1+x,y1+y));
                float nor2 = abs(Vx2(xc+x,yc+y))+abs(Vy2(xc+x,yc+y));
                if (nor1>nor2){
                    Vx2(xc+x,yc+y) = Vx1(x1+x,y1+y);
                    Vy2(xc+x,yc+y) = Vy1(x1+x,y1+y); // On garde le gradient de plus forte norme (norme 1)
                }
            }
        }
    }

}

int main(int argc, char* argv[]) {
    Image<byte> I1, I2, I3;
    Image<float> I4;
    Image<float> Vx1,Vx2,Vy1,Vy2;
    const char* fic1 = srcPath("bateau1.jpg");
    const char* fic2 = srcPath("bateau2.jpg");
    if(argc>2) {
        fic1 = argv[1];
        fic2 = argv[2];
    }
    if(! load(I1, fic1) || ! load(I2, fic2)) {
        std::cout << "Probleme dans le chargement d'images" << std::endl;
        return 1;
    }

    Window W1 = openWindow(I1.width(), I1.height());
    display(I1);
    Window W2 = openWindow(I2.width(), I2.height());
    setActiveWindow(W2);
    I2 = affichable(I2);
    display(I2);


    while(true){

        int x1,y1,x2,y2,xc,yc;

        setActiveWindow(W1);

        if (!selectionRect(W1,x1,y1,x2,y2,I1)){ // Si il ya un clic droit pendant la phase de selection on sort du programme
            break;
        }

        I3 = I1.getSubImage(x1,y1,abs(x2-x1),abs(y2-y1));

        cloneRect(W2,xc,yc,I2,I3); // On clone la portion d image selectionnee

        anyClick();

        Vx1 = dx(I1);
        Vx2 = dx(I2);
        Vy1 = dy(I1);
        Vy2 = dy(I2); // On calcule les gradients avec la FFT

        maxGradient(Vx1,Vy1,Vx2,Vy2,x1,y1,x2,y2,xc,yc); // On garde les plus hauts dans la zone concernee

        I4 = poisson(Vx2,Vy2); // On reconstruit avec Poisson

        setActiveWindow(W2);

        I2 = affichable(I4); // On met a jour l'image 2. Probleme : affichable rehausse le contraste a chaque iteration, notre image s eclaircit au fil des collages
        display(I2);

    }



    endGraphics();
    return 0;
}
