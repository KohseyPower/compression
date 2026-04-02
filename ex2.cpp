#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<math.h>
#include "fichiers.h"
#include "matrix.h" 
#include "pred.h" 
#include "dct.h"

 
int SaveIntImage_pgm(char *nom, int **im, int Height, int Width)
{ 
  int i, j;
  
  unsigned char ** ima_uc = alocamuc(Height, Width);  

  int maxim = 0;

  for(i = 0; i<Height; i++)
    for(j = 0; j<Width; j++)
      {
	if(abs(im[i][j]) > maxim)
	  maxim = abs(im[i][j]); 
      }
  //  fprintf(stderr, "\nnom = %s  maxim = %d\n", nom, maxim); 
  
  for(i = 0; i<Height; i++)
    for(j = 0; j<Width; j++)
      {
	ima_uc[i][j] = (unsigned char) ( (abs(im[i][j]) *255)/ maxim);
      }
  
  ecriture_pgm(nom, ima_uc, Width, Height); 
  
  dalocuc(ima_uc, Height); 
  
  return 1; 
}


int SaveIntImage_pgm_tronc(char *nom, int **im, int Height, int Width)
{
  int i, j;
  
    unsigned char ** ima_uc = alocamuc(Height, Width); 

  
  for(i = 0; i<Height; i++)
    for(j = 0; j<Width; j++)
      {
	if( abs(im[i][j]) >= 255.0 )
	  ima_uc[i][j] = 255; 
	else
	  ima_uc[i][j] = (unsigned char) ( (abs(im[i][j])) );	
      }
  
  ecriture_pgm(nom, ima_uc, Width, Height); 
  
  dalocuc(ima_uc, Height); 
  
  return 1; 
}


/* Q1 b */
void my_codeurDPCM_with_loop(unsigned char **x, int **err, int H, int W, int step)
{
    int i, j;
    int pred, error, q_error, rec;
    
    for(i = 0; i < H; i++) {
        for(j = 0; j < W; j++) {
            // Prédiction DPCM par ligne : le pixel courant est prédit par le précédent
            if(j == 0) {
                // Pas de pixel précédent à gauche pour la première colonne
                pred = 128; 
            } else {
                // Boucle de rétroaction : on utilise le pixel RECONSTRUIT précédent
                pred = rec; 
            }
            
            // Calcul de l'erreur de prédiction
            error = (int)x[i][j] - pred;
            
            // Quantification de l'erreur avec la fonction fournie
            // On suppose que quantiz retourne la valeur quantifiée/déquantifiée
            q_error = quantiz(error, step);
            err[i][j] = q_error;
            
            // Étape de reconstruction locale (ce que fera le décodeur)
            rec = pred + q_error;
            
            // Saturation des valeurs pour rester dans l'espace colorimétrique [0, 255]
            if(rec < 0) rec = 0;
            if(rec > 255) rec = 255;
        }
    }
} 

/* Q1 a */
void my_codeurDPCM_without_loop(unsigned char **x, int **err, int H, int W, int step)
{
    int i, j;
    int pred, error, q_error;
    
    for(i = 0; i < H; i++) {
        for(j = 0; j < W; j++) {
            // Prédiction DPCM par ligne SANS boucle
            if(j == 0) {
                pred = 128; 
            } else {
                // On utilise le VRAI pixel précédent de l'image originale
                pred = (int)x[i][j-1]; 
            }
            
            // Calcul de l'erreur de prédiction
            error = (int)x[i][j] - pred;
            
            // Quantification de l'erreur
            q_error = quantiz(error, step);
            err[i][j] = q_error;
            
            // On n'a plus besoin de calculer "rec" ici puisqu'on ne s'en sert pas !
        }
    }
}

/* Q1 a and b */
void my_decodeurDPCM(int **err, unsigned char **xrec, int H, int W)
{
    int i, j;
    int pred, rec;
    
    for(i = 0; i < H; i++) {
        for(j = 0; j < W; j++) {
            // La logique de prédiction doit être strictement identique à celle du codeur
            if(j == 0) {
                pred = 128; 
            } else {
                pred = rec; 
            }
            
            // Reconstruction du pixel à partir de la prédiction et de l'erreur quantifiée
            rec = pred + err[i][j];
            
            // Saturation
            if(rec < 0) rec = 0;
            if(rec > 255) rec = 255;
            
            // Stockage dans l'image finale
            xrec[i][j] = (unsigned char)rec;
        }
    }
}

void my_codeur_adapt(unsigned char **x, int **err, int H, int W, int step)
{
    int i, j;
    int A, B, C, pred, error, q_error, rec;
    
    // Allocation locale d'une matrice pour stocker l'image reconstruite au fur et à mesure
    // C'est indispensable pour la boucle de rétroaction afin d'avoir accès à la ligne précédente
    unsigned char **xrec_local = alocamuc(H, W);

    for(i = 0; i < H; i++) {
        for(j = 0; j < W; j++) {
            
            // --- GESTION DES BORDS ---
            if (i == 0 && j == 0) {
                // Tout premier pixel (en haut à gauche) : aucun voisin
                pred = 128;
            } else if (i == 0) {
                // Première ligne : B et C n'existent pas, on utilise la prédiction 1D par ligne (A)
                pred = xrec_local[i][j-1];
            } else if (j == 0) {
                // Première colonne : A et B n'existent pas, on utilise le voisin du haut (C)
                pred = xrec_local[i-1][j];
            } else {
                // --- CAS GÉNÉRAL : PRÉDICTION ADAPTATIVE ---
                // On récupère les voisins RECONSTRUITS
                A = xrec_local[i][j-1];
                C = xrec_local[i-1][j];
                B = xrec_local[i-1][j-1];

                // Application de la formule : si |B-C| <= |B-A| alors pred = A, sinon C
                if (abs(B - C) <= abs(B - A)) {
                    pred = A;
                } else {
                    pred = C;
                }
            }

            // Calcul de l'erreur de prédiction
            error = (int)x[i][j] - pred;

            // Quantification de l'erreur
            q_error = quantiz(error, step);
            err[i][j] = q_error;

            // Reconstruction locale pour les prochaines prédictions
            rec = pred + q_error;
            
            // Saturation
            if (rec < 0) rec = 0;
            if (rec > 255) rec = 255;
            
            xrec_local[i][j] = (unsigned char)rec;
        }
    }

    // Libération de la mémoire allouée localement
    dalocuc(xrec_local, H);
}


void my_decodeur_adapt(int **err, unsigned char **xrec, int H, int W)
{
    int i, j;
    int A, B, C, pred, rec;

    // L'ordre de parcours de l'image est de type raster-scan
    for(i = 0; i < H; i++) {
        for(j = 0; j < W; j++) {
            
            // --- La gestion des bords DOIT être identique à celle du codeur ---
            if (i == 0 && j == 0) {
                pred = 128;
            } else if (i == 0) {
                pred = xrec[i][j-1];
            } else if (j == 0) {
                pred = xrec[i-1][j];
            } else {
                // Récupération des voisins (déjà reconstruits car raster-scan)
                A = xrec[i][j-1];
                C = xrec[i-1][j];
                B = xrec[i-1][j-1];

                // Formule adaptative
                if (abs(B - C) <= abs(B - A)) {
                    pred = A;
                } else {
                    pred = C;
                }
            }

            // Reconstruction du pixel
            rec = pred + err[i][j];

            // Saturation
            if (rec < 0) rec = 0;
            if (rec > 255) rec = 255;
            
            xrec[i][j] = (unsigned char)rec;
        }
    }
}

int main (int argc, char *argv[])
{ 
  char nom[200], nom_out[200], nom_err[300]; 
  int W, H; // les dimensions de li'image: H = Height, W = Width
  int i, j; 

  if (argc != 3)
    { 
      fprintf(stderr, "Utilisation: %s <nom_image_pgm> <pas_quantification> \n", argv[0]); exit(0); 
    }

  strcpy(nom, argv[1]); 
  strcpy(nom_out, argv[1]); strcat(nom_out, ".out"); 
  strcpy(nom_err, argv[1]); strcat(nom_err, ".err"); 

  int step = atoi(argv[2]); 
  
  lecture_dim(nom, &W, &H);

  fprintf(stderr, "Width = %d  Height = %d\n", W, H);

  unsigned char **x = alocamuc(H, W);
  unsigned char **xrec = alocamuc(H, W);

  lecture_pgm(nom, x); // x contient l'image initiale

  int **err = alocami(H, W); // allocation de la matric des erreurs de prediction

  for(i=0;i<H;i++)
    for(j=0;j<W;j++)
	err[i][j] = (int)x[i][j]; 
   fprintf(stderr, "\nentropie initiale = %g [bits/pixel]\n", calc_entropie(err, H, W));


  
   /*
  // CODE DCT
  double **tdct = alocamd(H, W); // image transformee
  double **xd = alocamd(H, W);   // image initiale en double
  double **xrecd = alocamd(H, W); // image reconstruite

  

  for(i=0;i<H;i++)
    for(j=0;j<W;j++)
      xd[i][j] = (double)x[i][j]; 


  dct2dim(xd, tdct, H, W);
  fprintf(stderr, "OK DCT directe\n"); 
  
//Quantification des coefficients de la DCT
  for(i=1;i<H;i++)
    for(j=1;j<W;j++)
      {	
	err[i][j] = quantiz(tdct[i][j], step);
	tdct[i][j] = (double)err[i][j];	
      }
   for(j=0;j<W;j++)
     {err[0][j] =  quantiz(tdct[0][j], 1); tdct[0][j] = (double)err[0][j];}
   for(i=0;i<H;i++)
     {err[i][0] =  quantiz(tdct[i][0], 1); tdct[i][0] = (double)err[i][0];}
     
    //FIN Quantification des coefficients de la DCT
 

    double entro = calc_entropie(err, H, W); 
    fprintf(stderr, "\nentropie = %g [bits/pixel]\n", entro);
  
  dct2dim_inv(tdct, xrecd, H, W);
  for(i=0;i<H;i++)
    for(j=0;j<W;j++)
      if(xrecd[i][j] < 0.0)  xrec[i][j] = 0;
      else if(xrecd[i][j] > 255.0)  xrec[i][j] = 255;
      else  xrec[i][j] = (unsigned char)xrecd[i][j]; 
  SaveIntImage_pgm_tronc(nom_err, err, H, W); 
  ecriture_pgm(nom_out, xrec, W, H);
  
  dalocd(xrecd,H);
  dalocd(xd,H);
  dalocd(tdct,H);
  

  //FIN CODE DCT
  */

// DCT PAR BLOCS
/*
double **tdct = alocamd(H, W); // image transformee
double **xd = alocamd(H, W);   // image initiale en double
double **xrecd = alocamd(H, W); // image reconstruite

  for(i=0;i<H;i++)
    for(j=0;j<W;j++)
      xd[i][j] = (double)x[i][j]; 



int Bx=64, By=64; 

dct2dim_bloc(xd, tdct, H, W, Bx, By, step); 

for(i=0;i<H;i++)
  for(j=0;j<W;j++)
    err[i][j] = (int)tdct[i][j];

double entro = calc_entropie(err, H, W); 
fprintf(stderr, "\nentropie = %g [bits/pixel]\n", entro);

dct2dim_bloc_inv(tdct, xrecd, H, W, Bx, By); 

for(i=0;i<H;i++)
    for(j=0;j<W;j++)
      if(xrecd[i][j] < 0.0)  xrec[i][j] = 0;
      else if(xrecd[i][j] > 255.0)  xrec[i][j] = 255;
      else  xrec[i][j] = (unsigned char)xrecd[i][j]; 

  SaveIntImage_pgm_tronc(nom_err, err, H, W); 
  ecriture_pgm(nom_out, xrec, W, H);
   dalocd(xrecd,H);
  dalocd(xd,H);
  dalocd(tdct,H);
*/


// FIN DCT PAR BLOCS


//


//PREDICTION
  
  codeur_adapt(x, err, H, W, step);
  decodeur_adapt(err, xrec, H, W); 
  
   //  codeur_adapt(x, err, H, W, step);  
   //    decodeur_adapt(err, xrec, H, W);  
   
  
   double entro = calc_entropie(err, H, W); 
   fprintf(stderr, "\nentropie = %g [bits/pixel]\n", entro);
  SaveIntImage_pgm(nom_err, err, H, W); 
  

  ecriture_pgm(nom_out, xrec, W, H);
  
//FIN PREDICTION
   
   dalocuc(x,H);
   dalocuc(xrec,H);
  daloci(err,H);

  return 1; 
}


