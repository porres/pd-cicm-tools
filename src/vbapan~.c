/**********************************************************************/
//                                                                    //
// /****************************************************************/ //
// /*                                                              */ //
// /*                         VBAPAN                               */ //
// /*                                                              */ //
// /* Auteur: Rémi MIGNOT                                          */ //
// /*         Elève ingénieur Télécom2 à l'ISPG                    */ //
// /*         (Institut Supérieur Polytechnique de Galilée),       */ //
// /*         Université Paris13.                                  */ //
// /*                                                              */ //
// /* Date de création:   21/07/03                                 */ //
// /* Version: 1.5        16/07/04                                 */ //
// /*                                                              */ //
// /* Réalisé à la MSH Paris Nord (Maison des Sciences de l'Homme) */ //
// /*         en collaboration avec A.Sédes, B.Courribet           */ //
// /*         et J.B.Thiebaut,                                     */ //
// /*         CICM Université Paris8, MSH Paris Nord,              */ //
// /*         ACI Jeunes Chercheurs "Espaces Sonores".             */ //
// /*                                                              */ //
// /****************************************************************/ //
//                                                                    //
/**********************************************************************/


/**
 * Copyright (C) 2003-2004 Rémi Mignot, MSH Paris Nord,
 * 
 * This library is free software; you can redistribute it and/or modify it 
 * under the terms of the GNU Library General Public License as published 
 * by the Free Software Foundation; either version 2 of the License.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public 
 * License for more details.
 * 
 * You should have received a copy of the GNU Library General Public License 
 * along with this library; if not, write to the Free Software Foundation, 
 * Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * 
 * www.mshparisnord.org
 * rmignot@mshparisnord.org
 */


/*Description:  Cet objet permet de spatialiser un son monoral
                à l'aide de N haut-parleurs situés autour
                de l'auditeur. La spatialisation se fait grâce à
                au système VBAP (Vector Base Amplitude Panning).      */ 



/**********************************************************************/
/*                            EN TETE                                 */
/**********************************************************************/
 
#include <stdlib.h> 
#include <math.h>
#include "m_pd.h"

#define   Nmax      16     //Nombre maximum de haut-parleurs,
#define   Ndefaut   4      //Nombre de haut-parleurs par défaut,
#define   Xdefaut   0      //Valeur par défaut de l'abscisse,
#define   Ydefaut   1      //Valeur par défaut de l'ordonnée,
#define   Phidefaut 1.5708F//Valeur par défaut de l'angle,
#define   Rdefaut   1      //Distance par défaut de la source,
#define   Rdisq     .5     //Rayon du disque central par défaut,
#define   T_COS     4096   //Taille du tableau cosinus (puissance de 2),
#define   MASQUE    4095   //Masque pour le modulo = T_COS - 1,
#define   DTIME     10     //Temps d'interpolation par défaut en ms,
#define   Pi        3.1415926535897932384F  //   = Pi
#define   I360      0.0027777777777777778F  //   = 1/360
#define   EPSILON   0.000001 //Critère pour tester si les flottants sont nuls




/**********************************************************/
/*         La classe                                      */
/**********************************************************/

static t_class *vbapan_tilde_class;

typedef struct _vbapan_tilde
{  
  t_object  x_obj;        //Membre indispensable à Max,
	
  int	      base;	        //Type de repére, 1->cartèsienne, 0->polaire, 		
  int 	    controle;	    //1->coordonnées en contrôle, 0->coordonnées en signal,
  int       mute;         //1->entrées signal mutées, 0->non mutées,
  int       Nout;         //Nombre de sorties,
  int       N;            //Nombre de haut-parleurs.
	
  float 	  x,   y;	      //Coordonnées de la source en cartésien,
  float     phi, r;       //Coordonnées polaires de la source,

  int       dtime;        //Temps d'interpolation en échantillons,
  float     r_c;          //Rayon du disque central,
	
  float		  G[Nmax];	    //Tableau contenant les N gains des hp,
  float     teta[Nmax];   //Angles au centre des haut-parleurs,
  float     x_hp[Nmax];   //Abscisse de chaque haut-parleur,
  float     y_hp[Nmax];   //Ordonnée des haut-parleurs,
  float     dst_hp[Nmax]; //distance des haut-parleurs,

  float     dG[Nmax];     //Pas pour l'interpolation,
  float     Gstop[Nmax];  //Cible pour l'interpolation,

  float     G0[Nmax];     //gains des hp pour la position centrale,
  float     L[Nmax][2][2];//Tableaux contenant toute les matrices
                          //inversées des configurations de 2 haut-parleurs,
  int       rev[Nmax];    //Tableau qui permet d'ordonner les hp 
                          //par ordre croissant des angles, 

  float     *cosin;       //Adresse des tableaux pour les cosinus,
	
} t_vbapan_tilde;




//Déclaration de la méthode DSP qui se trouve à la fin du ficher,
void vbapan_tilde_dsp(t_vbapan_tilde *x, t_signal **sp);

//et de la fonction modulo,
float vbapan_tilde_modulo( float anglep );

//et de la fonction qui initialise les matrices inverses.
void vbapan_tilde_init_hp_mat( t_vbapan_tilde *x );




/**********************************************************************/
/*      TRAITEMENT DU SIGNAL avec coordonnées recues en controle.     */
/**********************************************************************/

t_int *vbapan_tilde_perform_controle(t_int *w)
{
  /*Préparation des pointeurs****************************************/
  t_vbapan_tilde *x   = (t_vbapan_tilde *)(w[1]);//adresse de la dataspace
  float     *son_in   = (t_sample *)(w[2]);      //son en entrée     
  int       n         = (int)(w[3+x->Nout]);     //taille des blocs

  /*Déclarations des variables locales*/
  float son;        //échantillon temporaire,
  int   hp;         //hp = n¡ du haut-parleur, 
  int   i;          //i = indice des échantillons,
  int   N = x->N;   //Nombre de haut-parleurs,
  int   retour = 4 + x->Nout;  //pour le pointeur de retour.

  /*Préparation des pointeurs de sorties ****************************/
  float *out[Nmax];
  for( hp = 0; hp < x->Nout; hp++)
    out[hp] = (t_sample *)(w[3+hp]);



  /******************************************************************/
  /*  Traitements  **************************************************/

  for( i=0; i<n; i++)
  {
    //on stocke l'échantillon d'entrée dans une variable temporaire,
    son = son_in[i];


    //Modulation des sorties avec les gains Pn.
    for( hp=N-1; hp >= 0; hp--)
    {
      //Incrémentation des gains pour l'interpolation,
      if( x->G[hp] == x->Gstop[hp] ) /*rien*/  ;
      else if ( fabs(x->Gstop[hp] - x->G[hp]) > fabs(x->dG[hp]) )
        x->G[hp] += x->dG[hp];
      else
        x->G[hp] = x->Gstop[hp];

      /******************************/
      out[hp][i] = son * x->G[hp];/**/
      /******************************/
    }
  }


  /*Initialisation à zéro des sorties inutilisées*/
  if( x->Nout > x->N){
    for( hp = x->Nout-1 ; hp >= N ; hp--){
      for( i=n-1 ; i>=0; i--)
        out[hp][i] = 0;
    }
  }

  return ( w+retour ) ;
}





/**********************************************************************/
/*      TRAITEMENT DU SIGNAL avec coordonnées recues en signal.       */
/**********************************************************************/

t_int *vbapan_tilde_perform_signal(t_int *w)
{
  /*Préparation des pointeurs****************************************/
  t_vbapan_tilde *x = (t_vbapan_tilde *)(w[1]);//adresse de la dataspace
  float     *son_in   = (t_sample *)(w[2]);      //son en entrée
  float     *xp       = (t_sample *)(w[3]);      //réception de x
  float     *yp       = (t_sample *)(w[4]);      //réception de y
  int       n         = (int)(w[5+x->Nout]);     //taille des blocs

  /*Déclarations des variables locales*/
  int   phii;           //angle en entier (de 0 à T_COS),
  float son;            //variable temporaire,
  int   hp;             //indice relatif au haut-parleur,
  int   i;              //indice du n¡ de l'échantillon,
  int   N = x->N;       //nombre de haut-parleurs,
  int   retour = x->Nout + 6; //retour pour pure data,
  float xl, yl, xtemp;  //coordonnéees cartésienne de la source,
  float g[Nmax];        //gains des haut-parleurs (en signal),
  float phi;            //angle de la source,
  float r;              //rayon de la source.
  float Puis, iAmp;     //puissance et amplitude des gains.
  
  
  /*Préparation des pointeurs des sorties ***************************/
  float *out[Nmax];
  for( hp = 0; hp < x->Nout; hp++)
    out[hp] = (t_sample *)(w[hp+5]);



  /******************************************************************/
  /*  Traitements  **************************************************/

  /******************************************************************/
  /*Si entrees signal mutées, on utilise le tableau P de contrôle.***/
  if( x->mute ) 
    for( i=0; i<n; i++)
    {
      //on stocke l'échantillon d'entrée dans une variable temporaire
      //car son_in = out[0].
      son = son_in[i];


      //Modulation des sorties avec les gains Pn.
      for( hp=N-1; hp >= 0; hp--)
      {
        //Incrémentation des gains pour l'interpolation,
        if( x->G[hp] == x->Gstop[hp] ) /*rien*/  ;
        else if ( fabs(x->Gstop[hp] - x->G[hp]) > fabs(x->dG[hp]) )
          x->G[hp] += x->dG[hp];
        else
          x->G[hp] = x->Gstop[hp];

        /******************************/
        out[hp][i] = son * x->G[hp];/**/
        /******************************/
      }
    }

  /******************************************************************/
  /*Si entrées signal non mutées, on utilise les vecteurs de signal**/
  else 
    for( i=n-1; i>=0; i--)
    {
      son = son_in[i];    //On stocke les échantillons d'entrée, dans des
      xl = xp[i];         //variables, car elles vont être écrasées.
      yl = yp[i];
  
      /*Conversion polaires -> cartésiennes,*******/
      if( !x->base )
      {
        phi = vbapan_tilde_modulo(yl)*Pi/180;
        r = xl;

        //ici xl = rayon et yl = angle,
        phii = (int)( yl*T_COS*I360 )&(int)MASQUE; 
        xtemp    = xl*x->cosin[ phii ];
        phii = (int)(.25*T_COS - phii)&(int)MASQUE;
        yl       = xl*x->cosin[ phii ];
        xl = xtemp;
        //maintenant xl = abscisse et yl = ordonnée.

        //Si r est négatif on modifie les coordonnées,
        if( r<0 ){
          r = -r;
          phi = (phi<Pi)?(phi+Pi):(phi-Pi);
        }

      }
      /*Conversion cartésienne -> polaire,**********/
      else
      {
        r = (float)sqrt( xl*xl+yl*yl );
        phi = (float)acos( xl/r );
        if( yl < 0 )
          phi = 2*Pi - phi;
      }


      //remise à zéro des gains:
      for( hp=0; hp<x->N; hp++)
        g[hp] = 0;

      //Recherche des hp encadrant le point:
      for( hp=0; hp<x->N-1; hp++){
        if( phi >= x->teta[x->rev[hp]] && phi < x->teta[x->rev[hp+1]] )
          break;
      }
  
      //Calcul du gain:
      g[x->rev[hp]]          = x->L[hp][0][0]*xl + x->L[hp][0][1]*yl;
      g[x->rev[(hp+1)%x->N]] = x->L[hp][1][0]*xl + x->L[hp][1][1]*yl; 
    
      //Puissance du gain (source et hp sur le cercle) :
      Puis =    ( g[x->rev[hp         ]]*g[x->rev[hp         ]] ) 
              + ( g[x->rev[(hp+1)%x->N]]*g[x->rev[(hp+1)%x->N]] );
      iAmp = (Puis <= EPSILON) ? 0 : (1/(float)sqrt(2*Puis));


      if( r > x->r_c )
      {
        //Normalisation des g, et prise en compte des distances :
        g[x->rev[hp]]            *=  iAmp * x->dst_hp[x->rev[ hp        ]] / r;
        g[x->rev[(hp+1)%x->N]]   *=  iAmp * x->dst_hp[x->rev[(hp+1)%x->N]] / r;
      }
      else
      {
        //Calcul du gain (sur le disque central):
        g[x->rev[hp]]            *=  iAmp * x->dst_hp[x->rev[ hp        ]] * r/(x->r_c*x->r_c);
        g[x->rev[(hp+1)%x->N]]   *=  iAmp * x->dst_hp[x->rev[(hp+1)%x->N]] * r/(x->r_c*x->r_c);
        //ici le gain est proportionnelle à r/r_c dans le disque.
        
        //On mixe linéairement l'effet VBAP avec les gains pour la position centrale:
        for( hp=0; hp<x->N; hp++)
          g[hp] += (1 - r/x->r_c) * x->G0[hp];
      }
      
      for( hp=N-1; hp >= 0 ; hp--)
      {
        /***************************/
        out[hp][i] = son * g[hp];/**/
        /***************************/
      }
    }


  /*Initialisation à zéro des sorties inutilisées*/
  for( hp = x->Nout-1 ; hp >= N ; hp--){
    for( i=n-1 ; i>=0; i--)
      out[hp][i] = 0;
  }

  return ( w+retour );
}



/**********************************************************************/
/*        METHODE RECOIT_x                                            */
/**********************************************************************/
//Cette méthode reçoit x et change tous les gains : 

void vbapan_tilde_recoit_x(t_vbapan_tilde *x, float xp)
{
  float Puis, iAmp;
  int   hp;
  float yp  = x->y;
  float phi = x->phi;
  float r;

  //Conversion polaires -> cartésiennes,
  if( !x->base )
  {
    //ici xp = rayon,
    r  = xp;
    xp = (float)( r * cos( phi ) );
    yp = (float)( r * sin( phi ) );
    //maintenant xp = abscisse.  
  }
  //Conversion cartésienne -> polaire,
  else
  {
    r = (float)sqrt( xp*xp+yp*yp );
    phi = (float)acos( xp/r );
    if( yp < 0 )
      phi = 2*Pi - phi;
  }  
  x->r   = r; 
  x->phi = phi;
  x->y   = yp;
  x->x   = xp;

  //Si r est négatif on modifie les coordonnées mais pas la dataspace,
  if( r<0 )
  {
    r = -r;
    phi = (phi<Pi)?(phi+Pi):(phi-Pi);
  }
  
  //Remise à zéro de tout les gains cibles:
  for(hp=0; hp<x->Nout; hp++)
    x->Gstop[hp] = 0;


  //Recherche des hp encadrant le point:
  for( hp=0; hp<x->N-1; hp++){
    if( phi >= x->teta[x->rev[hp]] && phi < x->teta[x->rev[hp+1]] )
      break;
  }
    

  //Calcul des gains :
  x->Gstop[x->rev[hp]]          = x->L[hp][0][0]*xp + x->L[hp][0][1]*yp;
  x->Gstop[x->rev[(hp+1)%x->N]] = x->L[hp][1][0]*xp + x->L[hp][1][1]*yp;

  //Puissance du gain (source et hp sur le cercle) :
  Puis =    ( x->Gstop[x->rev[hp         ]]*x->Gstop[x->rev[hp         ]] ) 
          + ( x->Gstop[x->rev[(hp+1)%x->N]]*x->Gstop[x->rev[(hp+1)%x->N]] );
  iAmp = (Puis <= EPSILON) ? 0 : (1/(float)sqrt(2*Puis));


  if( r > x->r_c )
  {
    //Normalisation des g, et prise en compte des distances :
    x->Gstop[x->rev[hp]]           *=  iAmp * x->dst_hp[x->rev[ hp        ]] / r;
    x->Gstop[x->rev[(hp+1)%x->N]]  *=  iAmp * x->dst_hp[x->rev[(hp+1)%x->N]] / r;
  }
  else
  {
    //Calcul du gain (sur le disque central):
    x->Gstop[x->rev[hp]]           *=  iAmp * x->dst_hp[x->rev[ hp        ]] * r/(x->r_c*x->r_c);
    x->Gstop[x->rev[(hp+1)%x->N]]  *=  iAmp * x->dst_hp[x->rev[(hp+1)%x->N]] * r/(x->r_c*x->r_c);
    //ici le gain est proportionnelle à r/r_c dans le disque.
    
    //On mixe linéairement l'effet VBAP avec les gains pour la position centrale:
    for( hp=0; hp<x->N; hp++)
      x->Gstop[hp] += (1 - r/x->r_c) * x->G0[hp];
  }
  
  // Interpolation :
  for( hp=0; hp<x->N; hp++)
  {
    x->dG[hp]          = ( x->Gstop[hp]          - x->G[hp]          ) /(float)x->dtime;
    x->dG[(hp+1)%x->N] = ( x->Gstop[(hp+1)%x->N] - x->G[(hp+1)%x->N] ) /(float)x->dtime;
  }

  return;
}



/**********************************************************************/
/*        METHODE RECOIT_y                                            */
/**********************************************************************/
//Cette méthode reçoit y et change tous les paramétres :

void vbapan_tilde_recoit_y(t_vbapan_tilde *x, float yp)
{
  float Puis, iAmp;
  int hp;
  float xp = x->x;
  float r  = x->r;
  float phi;
  
  //Conversion polaire -> cartésienne,
  if( !x->base )
  {
    //ici yp = angle,
    phi = (float)(vbapan_tilde_modulo(yp)*Pi/180.);
    xp  = (float)( r * cos( phi ) );
    yp  = (float)( r * sin( phi ) );
    //maintenant yp = ordonnée,
  }
  //Conversion cartésienne -> polaire,
  else
  {
    r = (float)sqrt( xp*xp+yp*yp );
    phi = (float)acos( xp/r );
    if( yp < 0 )
      phi = 2*Pi - phi;
  }    
  x->y = yp;
  x->x = xp;
  x->r = r;
  x->phi = phi;

  //Si r est négatif on modifie les coordonnées mais pas la dataspace,
  if( r<0 )
  {
    r = -r;
    phi = (phi<Pi)?(phi+Pi):(phi-Pi);
  }
  
  //Remise à zéro de tout les gains cibles:
  for(hp=0; hp<x->Nout; hp++)
    x->Gstop[hp] = 0;
 

  //Recherche des hp encadrant le point:
  for( hp=0; hp<x->N-1; hp++){
    if( phi >= x->teta[x->rev[hp]] && phi < x->teta[x->rev[hp+1]] )
      break;
  }


  //Calcul des gains :
  x->Gstop[x->rev[hp]]          = x->L[hp][0][0]*xp + x->L[hp][0][1]*yp;
  x->Gstop[x->rev[(hp+1)%x->N]] = x->L[hp][1][0]*xp + x->L[hp][1][1]*yp;

  //Puissance du gain (source et hp sur le cercle) :
  Puis =    ( x->Gstop[x->rev[hp         ]]*x->Gstop[x->rev[hp         ]] ) 
          + ( x->Gstop[x->rev[(hp+1)%x->N]]*x->Gstop[x->rev[(hp+1)%x->N]] );
  iAmp = (Puis <= EPSILON) ? 0 : (1/(float)sqrt(2*Puis));


  if( r > x->r_c )
  {
    //Normalisation des g, et prise en compte des distances :
    x->Gstop[x->rev[hp]]           *=  iAmp * x->dst_hp[x->rev[ hp        ]] / r;
    x->Gstop[x->rev[(hp+1)%x->N]]  *=  iAmp * x->dst_hp[x->rev[(hp+1)%x->N]] / r;
  } 
  else
  {
    //Calcul du gain (sur le disque central):
    x->Gstop[x->rev[hp]]           *=  iAmp * x->dst_hp[x->rev[ hp        ]] * r/(x->r_c*x->r_c);
    x->Gstop[x->rev[(hp+1)%x->N]]  *=  iAmp * x->dst_hp[x->rev[(hp+1)%x->N]] * r/(x->r_c*x->r_c);
    //ici le gain est proportionnelle à r/r_c dans le disque.
    
    //On mixe linéairement l'effet VBAP avec les gains pour la position centrale:
    for( hp=0; hp<x->N; hp++)
      x->Gstop[hp] += (1 - r/x->r_c) * x->G0[hp];
  }
  
  // Interpolation :
  for( hp=0; hp<x->N; hp++)
  {
    x->dG[hp]          = ( x->Gstop[hp]          - x->G[hp]          ) /(float)x->dtime;
    x->dG[(hp+1)%x->N] = ( x->Gstop[(hp+1)%x->N] - x->G[(hp+1)%x->N] ) /(float)x->dtime;
  }

  return;
}



/**********************************************************************/
/*        METHODE RECOIT_LISTE                                        */
/**********************************************************************/

void vbapan_tilde_recoit_liste( t_vbapan_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  float xp, yp = 0;
  
  //Récupération du premier paramètre (abscisse ou rayon),
  xp = atom_getfloat( argv );
    
  //Récupération du second paramètre (ordonnée ou angle),
  if( argc >= 2 )
    yp = atom_getfloat( argv+1 );
  else if( x->base )
    yp = x->y;
  else 
    yp = x->phi*180/Pi;
  
  //Initialisation des valeurs,
  if( x->base )    x->y = yp;
  else             x->phi = vbapan_tilde_modulo(yp)*Pi/180;

  vbapan_tilde_recoit_x( x, xp);

  return;
}



/**********************************************************************/
/*      METHODE POSITIONNER HAUT-PARLEUR sur le cercle avec phi       */
/**********************************************************************/

void vbapan_tilde_teta_positionner_hp( t_vbapan_tilde *x, t_symbol *s, int argc, t_atom *argv )
{
  int hp;

  for( hp=0; hp<x->N && hp<argc; hp++)
  {
    //Récupération des angles.
    x->teta[hp] = (float)(vbapan_tilde_modulo(atom_getfloat( argv+hp )) * Pi/180.);

    //Affectations des nouvelles coordonnées:
    x->x_hp[hp] = (float)cos( x->teta[hp] );
    x->y_hp[hp] = (float)sin( x->teta[hp] );

    //La distance sera initialisée dans "init_hp_mat()".
  }
  
  //Modification des matrices relatives aux hp:
  vbapan_tilde_init_hp_mat( x );
  
  /*Affectation des gains des hp*/
  if(x->controle)
    if(x->base)     vbapan_tilde_recoit_x( x, x->x);
    else            vbapan_tilde_recoit_x( x, x->r);

  return;
}



/**********************************************************************/
/*      METHODE POSITIONNER HAUT-PARLEUR  avec phi et le rayon        */
/**********************************************************************/

void vbapan_tilde_dist_teta_positionner_hp( t_vbapan_tilde *x, t_symbol *s, int argc, t_atom *argv )
{
  int hp;
  float rayon;

  for( hp=0; hp<x->N && 2*hp<argc; hp++)
  {
    //Récupération des rayons et des angles.
    rayon = (float)(atom_getfloat( argv+2*hp ));
    if( 2*hp+1<argc )
      x->teta[hp] = vbapan_tilde_modulo(atom_getfloat( argv+2*hp+1 )) * Pi/180;
   
    //Affectations des nouvelles coordonnées:
    x->x_hp[hp] = rayon*(float)cos( x->teta[hp] );
    x->y_hp[hp] = rayon*(float)sin( x->teta[hp] );
  }

  //Modification des matrices relatives aux hp:
  vbapan_tilde_init_hp_mat( x );


  /*Affectation des gains des hp*/
  if(x->base)     vbapan_tilde_recoit_x( x, x->x);
  else            vbapan_tilde_recoit_x( x, x->r);

  return;
}


/**********************************************************************/
/*      METHODE POSITIONNER HAUT-PARLEUR  avec coordonnées x et y     */
/**********************************************************************/

void vbapan_tilde_xy_positionner_hp( t_vbapan_tilde *x, t_symbol *s, int argc, t_atom *argv )
{
  int hp;
  float xp, yp = 0,
        rayon, angle;

  /*Récupération des coordonnées cartésiennes**********/
  for( hp=0; hp<x->N && 2*hp<argc; hp++)
  {
    xp = (float)(atom_getfloat( argv+2*hp ));
    if( 2*hp+1<argc )
      yp = (float)(atom_getfloat( argv+2*hp+1 ));
    else
      yp = x->y_hp[hp];


    /*Conversion en polaires****************************/  
    rayon = (float)sqrt( pow( xp, 2) + pow( yp, 2) );
    if( rayon!=0 )
    {
      angle = (float)acos( (float)xp/rayon );
      if(yp<0)
        angle = 2*Pi - angle;
    }
    else
      angle = 0;

    //Affectation des nouvelles coordonnées:
    x->teta[hp] = angle;
    x->x_hp[hp] = xp;
    x->y_hp[hp] = yp;
  }

  //Modification des matrices relatives aux hp:
  vbapan_tilde_init_hp_mat( x );
    
  /*Affectation des gains des hp*/
  if(x->base)     vbapan_tilde_recoit_x( x, x->x);
  else            vbapan_tilde_recoit_x( x, x->r);
    
  return;
}





/**********************************************************************/
/*              METHODE MODIFIER LE NOMBRE DE HP                      */
/**********************************************************************/

void vbapan_tilde_initialiser_nb_hp( t_vbapan_tilde *x, float f)
{
  int hp;
  int N = (int)f;

  if( x->Nout >= N && 2 < N ) 
  {
    x->N = (int)N;
    
    //Modifications des angles des haut-parleurs et des rayons,
    for( hp=0; hp<x->N; hp++){
      x->teta[hp] = (float)( Pi*( .5 + (1 - 2*hp )/(float)x->N) );
      if( x->teta[hp] < 0)
        x->teta[hp] += 2*Pi;
      x->x_hp[hp] = (float)cos( x->teta[hp] );
      x->y_hp[hp] = (float)sin( x->teta[hp] );
    }

    //Modification des matrices relatives aux hp:
    vbapan_tilde_init_hp_mat( x );
  }
  else
    post("vbapan~: Le nombre de haut-parleurs doit \210tre"
                   "inf\202rieur \205 %d, \n"
                   "  et sup\202rieur \205 3, ou \202gal.", x->Nout);

  /*Affectation des gains */
  if(x->base)     vbapan_tilde_recoit_x( x, x->x);
  else            vbapan_tilde_recoit_x( x, x->r);
   
}



/**********************************************************************/
/*              METHODE POUR MUTER LES ENTREES SIGNAL                 */
/**********************************************************************/

void  vbapan_tilde_muter_entrees_signal(t_vbapan_tilde *x, float mute)
{
  if( mute == 0)   x->mute = 0;
  else             x->mute = 1;
}



/**********************************************************************/
/*                   MODIFICATION DU TYPE DE REPERE                   */
/**********************************************************************/

void vbapan_tilde_changer_type_repere( t_vbapan_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  char   chaine[256];

  //Si pas d'argument, on quitte
  if( argc < 1) 
    return;

  //Récupération de la chaine de caractère.
  atom_string( argv, chaine, 255);

  //Changement de x->base
  if( chaine[0] == 'c')
    x->base = 1;
  else if( chaine[0] == 'p')
    x->base = 0; 
  else 
    post("vbapan~: erreur, ne reconnait le type de rep\212re.");
    
  //Si coordonnées polaires il faut initialiser x->r et x->phi:
  if( x->base == 0 )    //(Inutile en fait).
  {
    //Conversion de x et y en polaires:  
    x->r = (float)sqrt( pow( x->x, 2) + pow( x->y, 2) );
    if( x->r!=0 )
    {
      x->phi = (float)acos( (float)x->x/x->r );
      if(x->y<0)
        x->phi = 2*Pi - x->phi;
    }
    else
      x->phi = 0;
  }  
}




/**********************************************************************/
/*            METHODE CHANGER RAYON DU DISQUE CENTRAL                 */
/**********************************************************************/

void vbapan_tilde_changer_rayon_disq_centrale( t_vbapan_tilde *x, float r_c)
{
  if( r_c <= EPSILON ){
    post("vbapan~: pas de rayon nul ou n\202gatif.");
    return;
  }
  x->r_c = r_c;

  //Modification des matrices relatives aux hp:
  vbapan_tilde_init_hp_mat( x );
    
  /*Affectation des gains des hp*/
  if(x->base)     vbapan_tilde_recoit_x( x, x->x);
  else            vbapan_tilde_recoit_x( x, x->r);

  return;
}



/**********************************************************************/
/*         AFFICHAGE DES INFORMATIONS DANS LA FENETRE DE MAX          */
/**********************************************************************/

void vbapan_tilde_informations( t_vbapan_tilde *x)
{
  int hp;
  
  post("Info Vbapan~: ");
  
  if(x->controle && x->base)
    post("   coordonn\202es cart\202siennes en contr\223le,");
  else if(x->controle && !x->base)
    post("   coordonn\202es polaires en contr\223le,");
  else if(!x->controle && x->base)
    post("   coordonn\202es cart\202siennes en signal,");
  else if(!x->controle && !x->base)
    post("   coordonn\202es polaires en signal,");

  post("   rayon du disque cental = %f,", x->r_c);
    
  if(!x->controle && x->mute )
    post("   avec entr\202es signal inactives,");
  
  if(x->controle)
    post("   temps d'interpolation  = %d ms,", 
         (int)(x->dtime*1000/sys_getsr()));
  post("   nombre d'haut-parleurs = %d," , x->N);
  
  post("   position des haut-parleurs:");
  for( hp=0; hp<x->N-1; hp++)
    post("      hp n\370%d: %f.x + %f.y,", hp+1, x->x_hp[hp], x->y_hp[hp] );
    
  post("      hp n\370%d: %f.x + %f.y.", hp+1, x->x_hp[hp], x->y_hp[hp] );
}
    


/**********************************************************************/
/*                       FONCTION DESTRUCTION                         */
/**********************************************************************/

void vbapan_tilde_dest(t_vbapan_tilde *x)
{
  //Libération de la mémoire allouée pour les cosinus,
  if( !x->controle  )
    free( x->cosin );
 
  return;
}



/**********************************************************************/
/*                       FONCTION CREATION                            */
/**********************************************************************/

void *vbapan_tilde_new( t_symbol *s, int argc, t_atom *argv )
{  
  int    i, hp;
  char   chaine[256];


  /*Allocation de la dataspace *****************************/
  t_vbapan_tilde *x = (t_vbapan_tilde *)pd_new(vbapan_tilde_class);


  /*********************************************************/
  /*Récupération du nombre de sorties et de haut-parleurs. */
  if( argc >= 1 )
    x->N = (int)atom_getint( argv );
  else 
    x->N = Ndefaut;

  if( x->N <= 2 || x->N > Nmax ){
    x->N = Ndefaut;
    post("vbapan~: Il y a un probl\212me dans la d\202claration "
         "du nombre de haut-parleur, ");
    post("  il est de %d par défaut.", Ndefaut);
  }
  x->Nout = x->N;  //Même nombre de sorties et d'haut-parleur.
  

  /*Récupération du type repère ****************************/
  if( argc >=2 ){
    atom_string( argv+1, chaine, 4);

    if( chaine[0] == 'c' )
      x->base=1;
    else if( chaine[0] == 'p' )
      x->base=0;
    else {
      x->base=1;
      post("vbapan~: erreur dans le type des coordonn\202es, ");
      post("  elles sont cart\202siennes par d\202faut.");
    }
  }
  else 
    x->base = 1;

    
  /*Récupération du type des entrées ***********************/
  if( argc >=3 ){    
    atom_string( argv+2, chaine, 4);

    if( chaine[0] == 'c' )
      x->controle=1;
    else if( chaine[0] == 's' )
      x->controle=0;
    else {
      x->controle=1;
      post("vbapan~: erreur dans le type des entr\202es, ");
      post("  elles sont en controle par d\202faut.");
    }
  }
  else 
    x->controle = 1;

  
  /*Récupération du rayon du disque central ****************/
  if( argc >=4 )
    x->r_c = (float)fabs(atom_getfloat( argv+3));
  else
    x->r_c = (float)Rdisq;
  if( x->r_c <= EPSILON ){
    post("vbapan~: Pas de rayon de disque centrale nul s'il vous plait.");
    x->r_c = (float)Rdisq;
  }

  
  /*Récupération du temps d'interpolation ******************/
  if( argc >= 5 )
    x->dtime = (int)(abs(atom_getint( argv+4))*sys_getsr()/1000.);
  else
    x->dtime = (int)(DTIME*sys_getsr()/1000.);
  if( x->dtime == 0 )
    x->dtime = 1;



  /*********************************************************/
  /*Initialisation des données de la classe*****************/
  for( hp = x->N-1; hp >= 0; hp--) //Initialisation des P
    x->G[hp] = 0;

  x->mute= 0;             //entrées signal non mutées
  x->y   = Ydefaut;       //initialisation de y
  x->phi = Phidefaut;     //initialisation de phi
  
  //Initialisation des angles des haut-parleurs et des rayons,
  for( hp=0; hp<x->N; hp++)
  {
    x->teta[hp] = (float)( Pi*( .5 + (1 - 2*hp )/(float)x->N) );
    if( x->teta[hp] < 0)
      x->teta[hp] += 2*Pi;
    x->x_hp[hp] = (float)cos( x->teta[hp]);
    x->y_hp[hp] = (float)sin( x->teta[hp]);
  }

  //Initialisation de tout les tableaux relatifs aux hp:
  vbapan_tilde_init_hp_mat( x );

  //initialisation de x, X, Y, W et Pn par la mŽthode recoit x.
  if( x->base )    vbapan_tilde_recoit_x( x, Xdefaut);
  else             vbapan_tilde_recoit_x( x, Rdefaut);


  /*********************************************************/
  /*CrŽation du tableaux cosinus                           */
  if( !x->controle )
  {
    x->cosin = (float*)malloc(T_COS*sizeof(float));

    //Remplissage du tableau cosinus,
    for( i=0; i<T_COS; i++)
      x->cosin[i] = (float)cos( i*2*Pi/T_COS );
      
    /*Pour avoir besoin de cos( phi ), on cherche:
      cos(phi) = 
      cosin[ (int)( phi*T_COS/360 ))&(((int)T_COS-1) ] */
  }
       

  /*********************************************************/
  /*Création des nouvelles entrées. ************************/
  //si controle on crée des entrées de contrôle,
  if( x->controle )  {
    inlet_new(&x->x_obj, &x->x_obj.ob_pd,     
              gensym("float"), gensym("x"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, 
              gensym("float"), gensym("y"));
  }
  //Sinon on crée des entrées de signal.
  else  {    
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, 
              gensym("signal"), gensym("signal"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, 
              gensym("signal"), gensym("signal"));
  }

  /*Création des sorties************************************/
  for( i=0; i < x->N; i++)
    outlet_new(&x->x_obj, gensym("signal"));

  return (void *)x;
}



/**********************************************************************/
/*                       FONCTION SETUP                               */
/**********************************************************************/

void vbapan_tilde_setup(void) 
{
  /****Définition de la classe**********************************/
  vbapan_tilde_class = class_new( gensym("vbapan~"),
                                (t_newmethod) vbapan_tilde_new,
                                (t_method) vbapan_tilde_dest,
                                sizeof(t_vbapan_tilde),
                                CLASS_DEFAULT, A_GIMME, 0);


  /****Définition des méthodes**********************************/
  class_addmethod(vbapan_tilde_class,
        (t_method)vbapan_tilde_dsp, gensym("dsp"), 0);
  class_addmethod(vbapan_tilde_class, nullfn, gensym("signal"), 0);

  /*Méthodes pour le changement des coordonnées de la source****/
  class_addfloat (vbapan_tilde_class, (t_method)vbapan_tilde_recoit_x);
  class_addmethod(vbapan_tilde_class, 
                 (t_method)vbapan_tilde_recoit_x, gensym("x"), A_DEFFLOAT, 0);
  class_addmethod(vbapan_tilde_class, 
                 (t_method)vbapan_tilde_recoit_y, gensym("y"), A_DEFFLOAT, 0);
  class_addlist(vbapan_tilde_class, 
                 (t_method)vbapan_tilde_recoit_liste);

  /*Méthodes pour le changement des coordonnées des h-p********/
  class_addmethod(vbapan_tilde_class, 
                 (t_method)vbapan_tilde_teta_positionner_hp, 
                  gensym("a_setpos"), A_GIMME, 0);
  class_addmethod(vbapan_tilde_class, 
                 (t_method)vbapan_tilde_dist_teta_positionner_hp, 
                  gensym("ra_setpos"), A_GIMME, 0);
  class_addmethod(vbapan_tilde_class, 
                 (t_method)vbapan_tilde_xy_positionner_hp, 
                  gensym("xy_setpos"), A_GIMME, 0);
  class_addmethod(vbapan_tilde_class, 
                 (t_method)vbapan_tilde_initialiser_nb_hp, 
                  gensym("set_nb_hp"), A_DEFFLOAT, 0);
  class_addmethod(vbapan_tilde_class, 
                 (t_method)vbapan_tilde_muter_entrees_signal, 
                  gensym("mute_sig"), A_DEFFLOAT, 0);
  class_addmethod(vbapan_tilde_class, 
                 (t_method)vbapan_tilde_changer_rayon_disq_centrale, 
                  gensym("set_rc"), A_DEFFLOAT, 0);
  class_addmethod(vbapan_tilde_class, 
                 (t_method)vbapan_tilde_changer_type_repere, 
                  gensym("change_type"), A_GIMME, 0);
  class_addmethod(vbapan_tilde_class, 
                 (t_method)vbapan_tilde_informations , 
                  gensym("get_info"), A_GIMME, 0);
  class_addfloat (vbapan_tilde_class, (t_method)vbapan_tilde_recoit_x);
  

  /****Céation d'un lien pour l'aide****************************/
  class_sethelpsymbol(vbapan_tilde_class, gensym("vbapan~-help"));

}




/**********************************************************************/
/*      FONCTION QUI DONNE LE MODULO 360 D'UN ANGLE EN FLOTTANT       */
/**********************************************************************/

float vbapan_tilde_modulo( float anglep )
{
  float angle;

  if( anglep >= 0 && anglep < 360 )    return anglep;

  angle = anglep;
  while( angle < 0 )     angle += 360;
  while( angle >= 360 )  angle -= 360;
    
  return angle;
}



/**********************************************************************/
/*               FONCTION CHANGER LES MATRICES DES HP                 */
/**********************************************************************/

void vbapan_tilde_init_hp_mat( t_vbapan_tilde *x )
{
  float xhp1, yhp1, xhp2, yhp2; //Coordonnées provisoires des hp,
  int   hp, i, j;               //indices de boucles,
  int   itemp; float ftemp;     //entier et flottant temporaire,
  float g[Nmax];                //gains de la position centrale,
  float Puis, iAmp;             //puissance et amplitude des gains,
  float pteta, pdist;           //angle et distance porovisoires,
  float g1, g2;                 //gains provisoires.
  float px[4]={1,0,-1, 0};      //positions des 4 points pour la position centrale,
  float py[4]={0,1, 0,-1};


  //Initialisation des tableaux:
  for( i=0; i<Nmax; i++)
  {
    x->rev[i] = i;
    g[i] = 0;
  }


  //Arrangement du tableau rev qui donne l'ordre croissant des angles des hp:
  // ordre croissant -> ordre de création.
  for( j=0; j<x->N; j++ ) {
    for( i=0; i<x->N-1; i++ ) {
      if( x->teta[x->rev[i]] > x->teta[x->rev[i+1]] )
      {
        itemp = x->rev[i];
        x->rev[i] = x->rev[i+1];
        x->rev[i+1] = itemp;
      }
    }
  }
  

  //Calcul des distances des haut-parleurs:
  for( hp=0; hp<x->N; hp++)
    x->dst_hp[hp] = (float)sqrt( x->x_hp[hp]*x->x_hp[hp] + x->y_hp[hp]*x->y_hp[hp] );





  //Calcul des matrices inversŽes, dans l'ordre des téta croissants :
  for( hp=0; hp<x->N; hp++)
  {

    // Coodonnées des hp ramenés sur le cercle
    xhp1 = x->x_hp[ x->rev[hp]          ] / x->dst_hp[x->rev[hp]];
    yhp1 = x->y_hp[ x->rev[hp]          ] / x->dst_hp[x->rev[hp]];
    xhp2 = x->x_hp[ x->rev[(hp+1)%x->N] ] / x->dst_hp[x->rev[(hp+1)%x->N]];
    yhp2 = x->y_hp[ x->rev[(hp+1)%x->N] ] / x->dst_hp[x->rev[(hp+1)%x->N]];


    // Calcul du dénominateur
    ftemp = xhp1 * yhp2 - yhp1 * xhp2;

    // Si le dénominateur est nul, c'est que 2 hp consécutifs 
    //sont alignés, c'est pas possible.
    if( fabs( ftemp ) <= EPSILON )    {
      post( "vbapan~: deux haut-parleurs adjacents sont align\202s, \n"
            "  la configuration actuelle n'est pas possible.");
      ftemp = 0;
    }
    else
      ftemp = 1/ftemp;


    // Calcul de la matrice inverse associée aux 2 hp:
    x->L[hp][0][0] =  ftemp*yhp2;
    x->L[hp][0][1] = -ftemp*xhp2;
    x->L[hp][1][0] = -ftemp*yhp1;
    x->L[hp][1][1] =  ftemp*xhp1;
  }



  
  //Boucle sur les quatres points pour calculer les gains de la position centrale.
  for( i = 0; i < 4; i++ )
  {
      pdist = (float)sqrt( px[i]*px[i] + py[i]*py[i] );
  
      pteta = (float)acos( px[i]/pdist );
      if( py[i] < 0 )
        pteta = 2*Pi - pteta;
 
      //Recherche des hp encadrant le point:
      for( hp=0; hp<x->N-1; hp++){
        if( pteta >= x->teta[x->rev[hp]] && pteta < x->teta[x->rev[(hp+1)%x->N]] )
          break;
      }

      //Calcul des gains:
      g1 = x->L[hp][0][0]*px[i] + x->L[hp][0][1]*py[i];
      g2 = x->L[hp][1][0]*px[i] + x->L[hp][1][1]*py[i];

      //Calcul de l'amplitude efficace:
      Puis =  g1*g1 + g2*g2;
      iAmp = (Puis <= EPSILON) ? 0 : (1/(float)sqrt(2*Puis));

      //Normalisation des g, et prise en compte des distances:
      g[x->rev[hp]]           +=  g1 * iAmp * x->dst_hp[x->rev[ hp        ]] /pdist;
      g[x->rev[(hp+1)%x->N]]  +=  g2 * iAmp * x->dst_hp[x->rev[(hp+1)%x->N]] /pdist;
  }

  //Calcul de la puissance global:
  Puis = 0;
  for( hp=0; hp<x->N; hp++)
    Puis += (float)pow( g[hp]/x->dst_hp[x->rev[hp]], 2 );

  //Affectation dans la dataspace avec normalisation :
  for( hp=0; hp<x->N; hp++) 
    x->G0[hp] = (float)(g[hp]/(sqrt(2*Puis)*x->r_c));


  return;
}




/**********************************************************************/
/*                       METHODE DSP                                  */
/**********************************************************************/

void vbapan_tilde_dsp(t_vbapan_tilde *x, t_signal **sp)
{
  if( x->controle )
  {
    switch( x->Nout )
    {
    case 3:
      dsp_add(vbapan_tilde_perform_controle, 6, x,
              sp[0]->s_vec,   //son d'entrŽe,
              sp[1]->s_vec,   //tous les sons de sortie, ...
              sp[2]->s_vec,   // .........
              sp[3]->s_vec,
              sp[0]->s_n);    //taille des blocs.
      break;
    case 4:  
      dsp_add(vbapan_tilde_perform_controle, 7, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[0]->s_n); 
      break;      
    case 5:
      dsp_add(vbapan_tilde_perform_controle, 8, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[0]->s_n); 
      break;      
    case 6:
      dsp_add(vbapan_tilde_perform_controle, 9, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[0]->s_n); 
      break;      
    case 7:
      dsp_add(vbapan_tilde_perform_controle, 10, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[0]->s_n); 
      break;
    case 8:
      dsp_add(vbapan_tilde_perform_controle, 11, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[0]->s_n); 
      break;
    case 9:
      dsp_add(vbapan_tilde_perform_controle, 12, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[9]->s_vec, sp[0]->s_n); 
      break;
    case 10:
      dsp_add(vbapan_tilde_perform_controle, 13, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[9]->s_vec, sp[10]->s_vec, sp[0]->s_n); 
      break;
    case 11:
      dsp_add(vbapan_tilde_perform_controle, 14, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[0]->s_n); 
      break;
    case 12:
      dsp_add(vbapan_tilde_perform_controle, 15, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[0]->s_n); 
      break;
    case 13:
      dsp_add(vbapan_tilde_perform_controle, 16, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[0]->s_n); 
      break;
    case 14:
      dsp_add(vbapan_tilde_perform_controle, 17, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[14]->s_vec, 
              sp[0]->s_n); 
      break;
    case 15:
      dsp_add(vbapan_tilde_perform_controle, 18, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[14]->s_vec, 
              sp[15]->s_vec, sp[0]->s_n); 
      break;
    case 16:
      dsp_add(vbapan_tilde_perform_controle, 19, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[14]->s_vec, 
              sp[15]->s_vec, sp[16]->s_vec, sp[0]->s_n); 
      break;
    }
  }
  else
  {
    switch( x->Nout )
    {
    case 3:
      dsp_add(vbapan_tilde_perform_signal, 8, x,
              sp[0]->s_vec,   //son d'entrŽe,
              sp[1]->s_vec,   //entrŽe de x,
              sp[2]->s_vec,   //entrŽe de y,
              sp[3]->s_vec,   //tous les sons de sortie, ...
              sp[4]->s_vec,   // .........
              sp[5]->s_vec,
              sp[0]->s_n);    //taille des blocs.
      break;
    case 4:
      dsp_add(vbapan_tilde_perform_signal, 9, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[0]->s_n);  
      break;      
    case 5:
      dsp_add(vbapan_tilde_perform_signal, 10, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[0]->s_n);  
      break;      
    case 6:
      dsp_add(vbapan_tilde_perform_signal, 11, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[0]->s_n);
      break;      
    case 7:
      dsp_add(vbapan_tilde_perform_signal, 12, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[0]->s_n);  
      break;      
    case 8:
      dsp_add(vbapan_tilde_perform_signal, 13, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[0]->s_n);  
      break;      
    case 9:
      dsp_add(vbapan_tilde_perform_signal, 14, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[0]->s_n);  
      break;
    case 10:
      dsp_add(vbapan_tilde_perform_signal, 15, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[11]->s_vec, sp[0]->s_n);  
      break;
    case 11:
      dsp_add(vbapan_tilde_perform_signal, 16, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[0]->s_n);  
      break;
    case 12:
      dsp_add(vbapan_tilde_perform_signal, 17, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[14]->s_vec, 
              sp[0]->s_n);  
      break;
    case 13:
      dsp_add(vbapan_tilde_perform_signal, 18, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[14]->s_vec, 
              sp[15]->s_vec, sp[0]->s_n);  
      break;
    case 14:
      dsp_add(vbapan_tilde_perform_signal, 19, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[14]->s_vec, 
              sp[15]->s_vec, sp[16]->s_vec, sp[0]->s_n);  
      break;
    case 15:
      dsp_add(vbapan_tilde_perform_signal, 20, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[14]->s_vec, 
              sp[15]->s_vec, sp[16]->s_vec, sp[17]->s_vec, 
              sp[0]->s_n);  
      break;
    case 16:
      dsp_add(vbapan_tilde_perform_signal, 21, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[14]->s_vec, 
              sp[15]->s_vec, sp[16]->s_vec, sp[17]->s_vec, 
              sp[18]->s_vec, sp[0]->s_n);  
      break;
    }
  }
}
