/**********************************************************************/
//                                                                    //
// /****************************************************************/ //
// /*                                                              */ //
// /*                         AMBIPAN                              */ //
// /*                                                              */ //
// /* Auteur: R�mi MIGNOT                                          */ //
// /*         El�ve ing�nieur T�l�com2 � l'ISPG                    */ //
// /*         (Institut Sup�rieur Polytechnique de Galil�e),       */ //
// /*         Universit� Paris13.                                  */ //
// /*                                                              */ //
// /* Date de cr�ation:   11/07/03                                 */ //
// /* Version: 1.5        14/07/04                                 */ //
// /*                                                              */ //
// /* R�alis� � la MSH Paris Nord (Maison des Sciences de l'Homme) */ //
// /*         en collaboration avec A.S�des, B.Courribet           */ //
// /*         et J.B.Thiebaut,                                     */ //
// /*         CICM Universit� Paris8, MSH Paris Nord,              */ //
// /*         ACI Jeunes Chercheurs "Espaces Sonores".             */ //
// /*                                                              */ //
// /****************************************************************/ //
//                                                                    //
/**********************************************************************/


/**
 * Copyright (C) 2003-2004 R�mi Mignot, MSH Paris Nord,
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
                � l'aide de N haut-parleurs situ�s en cercle autour
                de l'auditeur. La spatialisation se fait gr�ce �
                l'ambisonie de Michael Gerzon.                        */ 



/**********************************************************************/
/*                            EN TETE                                 */
/**********************************************************************/

#include <stdlib.h> 
#include <math.h>
#include "m_pd.h"

#define   Nmax      16     //Nombre maximum de haut-parleurs,
#define   Ndefaut   4      //Nombre de haut-parleurs par d�faut,
#define   Xdefaut   0      //Valeur par d�faut de l'absisse,
#define   Ydefaut   1      //Valeur par d�faut de l'ordonn�e,
#define   Phidefaut 1.5708F//Valeur par d�faut de l'angle (PI/2),
#define   Rdefaut   1      //Distance par d�faut de la source,
#define   OFFSET    .3     //Offset par defaut de l'ambisonie,
#define   T_COS     4096   //Taille du tableau cosinus (puissance de 2),
#define   MASQUE    4095   //Masque pour le modulo = T_COS - 1,
#define   DTIME     10     //Temps d'interpolation par d�faut,
                           //en �chantillons (10ms),
#define   Pi        3.1415926535897932384F  // Pi
#define   I360      0.0027777777777777778F  // 1/360
#define   EPSILON   0.000001 //Crit�re pour tester si les flottants sont nuls



/**********************************************************/
/*         La classe                                      */
/**********************************************************/

static t_class *ambipan_tilde_class;

typedef struct _ambipan_tilde
{
  t_object  x_obj;      //Membre indispensable � pure_data,
	
	int			base;		      //Type de base, 1->cart�sienne, 0->polaire, 		
	int 		controle;	    //1->coordonn�es en contr�le, 0->coordonn�es en signal,
  int     mute;         //1->entr�es signal mut�es, 0->non mut�es,
  int     Nout;         //Nombre de sorties,
  int     N;            //Nombre de haut-parleurs.
	
	float 	x, y;	  	    //Coordonn�es de la source en cart�sien,
  float   phi, r;       //Coordonn�es polaires de la source,

  float   offset;       //Offset de l'ambisonie,
  int     dtime;        //Temps d'interpolation en �chantillons,
	
	float		P[Nmax];	    //Tableau contenant les N coefficients ambisoniques Pn,
  float   teta[Nmax];   //Angle de chaque haut-parleur,
  float   dist[Nmax];   //Distance des haut-parleurs,

  float   dP[Nmax];     //pas pour l'interpolation,
  float   Pstop[Nmax];  //cible pour l'interpolation,   

  float   *cosin;       //Adresse des tableaux pour les cosinus, 
  float   *cos_teta;    //les sinus et cosinus des angles des 
  float   *sin_teta;    //haut-parleurs.
	
} t_ambipan_tilde;


//D�claration de la m�thode DSP qui se trouve � la fin du ficher.
void ambipan_tilde_dsp(t_ambipan_tilde *x, t_signal **sp);




/**********************************************************************/
/*      TRAITEMENT DU SIGNAL avec coordonn�es recues en controle.     */
/**********************************************************************/

t_int *ambipan_tilde_perform_controle(t_int *w)
{
  /*Pr�paration des pointeurs****************************************/
  t_ambipan_tilde *x = (t_ambipan_tilde *)(w[1]);//adresse de la dataspace
  float     *son_in   = (t_sample *)(w[2]);      //son en entr�e     
  int       n         = (int)(w[3+x->Nout]);     //taille des blocs

  /*D�clarations des variables locales*/
  float son;                    //�chantillon temporaire,
  int   hp;                     //hp = n� du haut-parleur, 
  int   i;                      //i = indice des �chantillons,
  int   N = x->N;               //Nombre de haut-parleur,
  int   retour = 4 + x->Nout;   //pour le pointeur de retour.

  /*Pr�paration des pointeurs de sorties ****************************/
  float *out[Nmax];
  for( hp = 0; hp < x->Nout; hp++)
    out[hp] = (t_sample *)(w[3+hp]);



  /******************************************************************/
  /*  Traitements  **************************************************/

  for( i=0; i<n; i++)
  {
    //on stocke l'�chantillon d'entr�e dans une variable temporaire
    //car son_in = out[0].
    son = son_in[i];


    //Modulation des sorties avec les coefficients ambisoniques Pn.
    for( hp = N-1; hp >= 0; hp--)
    {
      //Incr�mentation des P pour l'interpolation,
      if( x->P[hp] == x->Pstop[hp] ) /*rien*/;
      else if ( fabs(x->Pstop[hp] - x->P[hp]) > fabs(x->dP[hp]) )
        x->P[hp] += x->dP[hp];
      else
        x->P[hp] = x->Pstop[hp];

      /******************************/
      out[hp][i] = son * x->P[hp];/**/
      /******************************/
    }
  }


  /*Initialisation � z�ro des sorties inutilis�es*/
  if( x->Nout > x->N){
    for( hp = x->Nout-1 ; hp >= N ; hp--){
      for( i=n-1 ; i>=0; i--)
        out[hp][i] = 0;
    }
  }
 
  return ( w+retour );
}



/**********************************************************************/
/*      TRAITEMENT DU SIGNAL avec coordonn�es recues en signal.       */
/**********************************************************************/

t_int *ambipan_tilde_perform_signal(t_int *w)
{
  /*Pr�paration des pointeurs****************************************/
  t_ambipan_tilde *x = (t_ambipan_tilde *)(w[1]);//adresse de la dataspace
  float     *son_in   = (t_sample *)(w[2]);      //son en entr�e
  float     *xp       = (t_sample *)(w[3]);      //r�ception de x
  float     *yp       = (t_sample *)(w[4]);      //r�ception de y
  int       n         = (int)(w[5+x->Nout]);     //taille des blocs

  /*D�clarations des variables locales*/
  float son;            //variable temporaire,
  int   hp;             //indices relatif au haut-parleur,
  int   i;              //indices du n� de l'�chantillon,
  int   N = x->N;       //Nombre de haut-parleur,_
  int   retour = x->Nout + 6; //retour pour pure data,
  float offset = x->offset;

  float K = (float)(sqrt(1/(float)N)/1.66);   //Facteur correctif.

  //Param�tres ambisoniques:
  float xtemp, xl, yl, ds, dist, X, Y, W, P;
  int   phii;

  /*Pr�paration des pointeurs des sorties ***************************/
  float *out[Nmax];
  for( hp = 0; hp < x->Nout; hp++)
    out[hp] = (t_sample *)(w[hp+5]);



  /******************************************************************/
  /*  Traitements  **************************************************/

  /******************************************************************/
  /*Si entrees signal mut�es, on utilise le tableau P de contr�le.***/
  if( x->mute ) 
    for( i=0; i<n; i++)
    {
      //on stocke l'�chantillon d'entr�e dans une variable temporaire
      //car son_in = out[0].
      son = son_in[i];


      //Modulation des sorties avec les coefficients ambisoniques Pn.
      for( hp = N-1; hp >= 0; hp--)
      {
        //Incr�mentation des P pour l'interpolation,
        if( x->P[hp] == x->Pstop[hp] ) /*rien*/  ;
        else if ( fabs(x->Pstop[hp] - x->P[hp]) > fabs(x->dP[hp]) )
          x->P[hp] += x->dP[hp];
        else
          x->P[hp] = x->Pstop[hp];

        /******************************/
        out[hp][i] = son * x->P[hp];/**/
        /******************************/
      }
    }

  /******************************************************************/
  /*Si entr�es signal non mut�es, on utilise les vecteurs de signal**/
  else 
    for( i=n-1; i>=0; i--)
    {
      son = son_in[i];    //On stocke les �chantillons d'entr�e, dans des
      xl = xp[i];         //constantes, car elles vont �tre �cras�es.
      yl = yp[i];

      //Conversion polaires -> cart�siennes,
      if( !x->base )
      {
        //ici xl = rayon et yl = angle, 
        phii = (int)( yl*T_COS*I360 )&(int)MASQUE;
        xtemp    = xl*x->cosin[ phii ];
        phii = (int)(.25*T_COS - phii)&(int)MASQUE;
        yl       = xl*x->cosin[ phii ];
        xl = xtemp;
        //maintenant xl = abscisse et yl = ordonn�e.  
      }

      //Calcul des distances,
      ds   = xl*xl + yl*yl;
      dist = (float)sqrt(ds);

      //Calcul des param�tres ambisoniques,
      X = (float)( 2*xl / (ds + offset) ); 
      Y = (float)( 2*yl / (ds + offset) ); 
      W = (float)( .707 / (dist + offset) );

      for( hp=N-1; hp >= 0 ; hp--)
      {
        P = K  * ( W + X*x->cos_teta[hp]  
                     + Y*x->sin_teta[hp]  )
               * x->dist[hp];
      
        //Si Pn<0 on les force � 0
        if(P < 0)            P = 0;

        /***********************/
        out[hp][i] = son * P;/**/
        /***********************/
      }
    }


  /*Initialisation � z�ro des sorties inutilis�es*/
  if( x->Nout > x->N){
    for( hp = x->Nout-1 ; hp >= N ; hp--){
      for( i=n-1 ; i>=0; i--)
        out[hp][i] = 0;
    }
  }

  return ( w+retour );
}



/**********************************************************************/
/*        METHODE RECOIT_x                                            */
/**********************************************************************/
//Cette m�thode re�oit x et change tous les param�tres ambisoniques: 
//  X, Y, W et les Pn 
void ambipan_tilde_recoit_x(t_ambipan_tilde *x, t_floatarg xp)
{
  int hp;
  float ds, dist;
  float X, Y, W;
  float yp = x->y;
  float K = (float)(sqrt(1/(float)x->N)/1.66 ); 


  //Conversion polaire -> cart�sienne,
  if( !x->base )
  {
    //ici xp = rayon,
    x->r = xp;
    xp   = (float)( x->r * cos( x->phi ) );
    yp   = (float)( x->r * sin( x->phi ) );
    //maintenant xp = abscisse,
  }
  x->x = xp;
  x->y = yp;

  //Calcul des distances,
  ds = xp*xp + yp*yp;
  dist = (float)sqrt(ds);

  //Calcul des param�tres ambisoniques,
  X = (float)( 2*xp/(ds + x->offset) ); 
  Y = (float)( 2*yp/(ds + x->offset) ); 
  W = (float)( .707/(dist + x->offset) );

  //Calcul des coefficients ambisoniques cibles et des pas,
  for( hp = x->N-1; hp >= 0 ; hp--)
  {
    x->Pstop[hp] = (float)( ( W + X*cos(x->teta[hp]) 
                                + Y*sin(x->teta[hp])
                            ) * K  * x->dist[hp]);
    //Si Pstop_n < 0 on les force � 0.
    if(x->Pstop[hp] < 0)
      x->Pstop[hp] = 0;

    x->dP[hp] = (x->Pstop[hp] - x->P[hp])/(float)x->dtime;
  }

  return;
}



/**********************************************************************/
/*        METHODE RECOIT_y                                            */
/**********************************************************************/
//Cette m�thode re�oit y et change tous les param�tres ambisoniques: 
//  X, Y, W et les Pn 
void ambipan_tilde_recoit_y(t_ambipan_tilde *x, t_floatarg yp)
{
  int hp;
  float ds, dist;
  float X, Y, W;
  float xp = x->x;
  float K = (float)(sqrt(1/(float)x->N)/1.66 ); 


  //Conversion polaires -> cart�sienne,
  if( !x->base )
  {
    //ici yp = angle,
    x->phi = (float)( yp*Pi/180 );
    xp   = (float)( x->r * cos( x->phi ) );
    yp   = (float)( x->r * sin( x->phi ) );
    //maintenant yp = ordonn�e,
  }
  x->y = yp;
  x->x = xp;

  //Calcul des distances,
  ds = xp*xp + yp*yp;
  dist = (float)sqrt(ds);

  //Calcul des param�tres ambisoniques,
  X = (float)( 2*xp/(ds + x->offset) ); 
  Y = (float)( 2*yp/(ds + x->offset) ); 
  W = (float)( .707/(dist + x->offset) );

  //Calcul des coefficients ambisoniques cibles et des pas,
  for( hp=x->N-1; hp >= 0 ; hp--)
  {
    x->Pstop[hp] = (float)( ( W + X*cos(x->teta[hp]) 
                                + Y*sin(x->teta[hp])
                            ) * K  * x->dist[hp]);
    //Si Pstop_n<0 on les force � 0
    if(x->Pstop[hp] < 0)
      x->Pstop[hp] = 0;

    x->dP[hp] = (x->Pstop[hp] - x->P[hp])/(float)x->dtime;    
  }

  return;
}



/**********************************************************************/
/*        METHODE RECOIT_LISTE                                        */
/**********************************************************************/

void ambipan_tilde_recoit_liste( t_ambipan_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  float xp, yp = 0;

  //R�cup�ration du premier param�tre (abscisse ou rayon),
  xp = atom_getfloat( argv );

  //R�cup�ration du second param�tre (ordonn�e ou angle),
  if( argc >= 2 )
    yp = atom_getfloat( argv+1 );
  //si qu'un param�tre, on prend la valeur d�ja pr�sente.
  else if( x->base )     yp = x->y;
  else                   yp = x->phi*180/Pi;

  //Initialisation des valeurs,
  if( x->base )    x->y = yp;
  else             x->phi = yp*Pi/180;

  ambipan_tilde_recoit_x( x, xp);

  return;
}



/**********************************************************************/
/*      METHODE POSITIONNER HAUT-PARLEUR sur le cercle avec phi       */
/**********************************************************************/

void ambipan_tilde_teta_positionner_hp( t_ambipan_tilde *x, t_symbol *s, int argc, t_atom *argv )
{
  int hp;

  for( hp=0; hp<x->N && hp<argc; hp++)
  {
    //R�cup�ration des angles.
    x->teta[hp] = (float)(atom_getfloat( argv+hp )*Pi/180);
   
    x->dist[hp] = 1;

    //Pr�calcul des cos et sin en signal. 
    if( !x->controle )
    {
      x->cos_teta[hp] = (float)cos( x->teta[hp] );
      x->sin_teta[hp] = (float)sin( x->teta[hp] );
    }
  }
  
  /*Affectation des gains de l'ambisonie*/
  if(x->controle)
    if(x->base)     ambipan_tilde_recoit_x( x, x->x);
    else            ambipan_tilde_recoit_x( x, x->r);

  return;
}



/**********************************************************************/
/*      METHODE POSITIONNER HAUT-PARLEUR  avec phi et le rayon        */
/**********************************************************************/

void ambipan_tilde_dist_teta_positionner_hp( t_ambipan_tilde *x, t_symbol *s, int argc, t_atom *argv )
{
  int hp;

  for( hp=0; hp<x->N && 2*hp<argc; hp++)
  {
    //R�cup�ration des rayons et des angles.
    x->dist[hp] = (float)(atom_getfloat( argv+2*hp ));
    if( 2*hp+1<argc )
      x->teta[hp] = (float)(atom_getfloat( argv+2*hp+1 )*Pi/180);
   
    //Pr�calcul des cos et sin en signal.
    if( !x->controle )
    {
      x->cos_teta[hp] = (float)cos( x->teta[hp] );
      x->sin_teta[hp] = (float)sin( x->teta[hp] );
    }
  }
  
  /*Affectation des gains de l'ambisonie*/
  if(x->controle)
    if(x->base)     ambipan_tilde_recoit_x( x, x->x);
    else            ambipan_tilde_recoit_x( x, x->r);

  return;
}



/**********************************************************************/
/*      METHODE POSITIONNER HAUT-PARLEUR  avec coordonn�es x et y     */
/**********************************************************************/

void ambipan_tilde_xy_positionner_hp( t_ambipan_tilde *x, t_symbol *s, int argc, t_atom *argv )
{
  int hp;
  float xp, yp = 0;
  float rayon, angle;

  /*R�cup�ration des coordonn�es cart�siennes**********/
  for( hp=0; hp<x->N && 2*hp<argc; hp++)
  {
    xp = (float)(atom_getfloat( argv+2*hp ));
    if( 2*hp+1<argc )
      yp = (float)(atom_getfloat( argv+2*hp+1 ));


    /*Conversion en polaires***************************/  
    rayon = (float)sqrt( pow( xp, 2) + pow( yp, 2) );
    if( rayon!=0 )
    {
      angle = (float)acos( (float)xp/rayon );
      if(yp<0)
        angle = -angle;
    }
    else
      angle = 0;

    //Affectatoin des nouvelles valeurs
    x->teta[hp] = angle;
    x->dist[hp] = rayon;  
   
    /*Pr�calcul des cos et des sin ********************/
    if( !x->controle )
    {
      x->cos_teta[hp] = (float)(cos( angle ));
      x->sin_teta[hp] = (float)(sin( angle ));
    }
  }
      
  /*Affectation des gains de l'ambisonie*/
  if(x->controle)
    if(x->base)     ambipan_tilde_recoit_x( x, x->x);
    else            ambipan_tilde_recoit_x( x, x->r);

  return;
}



/**********************************************************************/
/*              METHODE MODIFIER LE NOMBRE DE HP                      */
/**********************************************************************/

void ambipan_tilde_initialiser_nb_hp( t_ambipan_tilde *x, float f)
{
  int hp;
  int N = (int)f;

  if( x->Nout >= N && 2 <= N )
  {
    x->N = (int)N;
    
    //Modifications des angles des haut-parleurs et des rayons,
    for( hp=0; hp<x->N; hp++)
    {
      x->teta[hp] = (float)( Pi*( .5 + (1 - 2*hp )/(float)x->N) );
      x->dist[hp] = 1;
    }


    /*****************************************************/
    /*Modification des tableaux cos_teta et sin_teta     */
    /*pour l'audio.                                      */
    if( !x->controle )
    {
      //Pr�calculs des cos et sin des haut-parleurs.
      for( hp=0; hp<x->N; hp++)
      {
        x->cos_teta[hp] = (float)cos( x->teta[hp] );
        x->sin_teta[hp] = (float)sin( x->teta[hp] );
      }
    }
  }  
  else
    post("ambipan~: Le nombre de haut-parleurs doit \210tre"
                    "inf\202rieur \205 %d, \n"
                    "  et sup\202rieur \205 2, ou \202gal.", x->Nout);
    
  /*Affectation des gains de l'ambisonie*/
  if(x->controle)
    if(x->base)     ambipan_tilde_recoit_x( x, x->x);
    else            ambipan_tilde_recoit_x( x, x->r);
   
}



/**********************************************************************/
/*              METHODE POUR MUTER LES ENTREES SIGNAL                 */
/**********************************************************************/

void  ambipan_tilde_muter_entrees_signal(t_ambipan_tilde *x, float mute)
{
  if( mute == 0)   x->mute = 0;
  else             x->mute = 1;
}



/**********************************************************************/
/*                   MODIFICATION DE L'OFFSET                         */
/**********************************************************************/

void ambipan_tilde_changer_offset( t_ambipan_tilde *x, float offset)
{
  if( offset <= EPSILON ){
    post("ambipan~: Pas d'offset n\202gatif ou nul s'il vous plait.");
    return;
  }

  //Changement de l'offset
  x->offset = offset; 
  
  /*Affectation des gains de l'ambisonie*/
  if(x->controle)
    if(x->base)     ambipan_tilde_recoit_x( x, x->x);
    else            ambipan_tilde_recoit_x( x, x->r);
   
}



/**********************************************************************/
/*                   MODIFICATION DU TYPE DE REPERE                   */
/**********************************************************************/

void ambipan_tilde_changer_type_repere( t_ambipan_tilde *x, t_symbol *s, int argc, t_atom *argv )
{
  char   chaine[256];

  //Si pas d'argument, on quitte
  if( argc < 1) 
    return;

  //R�cup�ration de la chaine de caract�re.
  atom_string( argv, chaine, 255);

  //Changement de x->base
  if( chaine[0] == 'c')
    x->base = 1;
  else if( chaine[0] == 'p')
    x->base = 0; 
  else 
    post("ambipan~: erreur, ne reconnait le type de rep\212re.");
    
  //Si coordonn�es polaires il faut initialiser x->r et x->phi:
  if( x->base == 0 )
  {
    //Conversion de x et y en polaires:  
    x->r = (float)sqrt( pow( x->x, 2) + pow( x->y, 2) );
    if( x->r!=0 )
    {
      x->phi = (float)acos( (float)x->x/x->r );
      if(x->y<0)
        x->phi = -x->phi;
    }
    else
      x->phi = 0;
  }  
}



/**********************************************************************/
/*         AFFICHAGE DES INFORMATIONS DANS LA FENETRE DE MAX          */
/**********************************************************************/

void ambipan_tilde_informations( t_ambipan_tilde *x)
{
  int hp;
  
  post("Info Ambipan~: ");
  
  if(x->controle && x->base)
    post("   coordonn\202es cart\202siennes en contr\223le,");
  else if(x->controle && !x->base)
    post("   coordonn\202es polaires en contr\223le,");
  else if(!x->controle && x->base)
    post("   coordonn\202es cart\202siennes en signal,");
  else if(!x->controle && !x->base)
    post("   coordonn\202es polaires en signal,");
    
  if(!x->controle && x->mute )
    post("   avec entr\202es signal inactives,");
  
  post("   offset                 = %f,", x->offset);
  if(x->controle)
    post("   temps d'interpolation  = %d ms,", 
         (int)(x->dtime*1000/sys_getsr()));
  post("   nombre d'haut-parleurs = %d," , x->N);
  
  post("   position des haut-parleurs:");
  for( hp=0; hp<x->N-1; hp++)
    post("      hp n\370%d: %f.x + %f.y,", hp+1, x->dist[hp]*cos(x->teta[hp]), 
                                            x->dist[hp]*sin(x->teta[hp]));
    
  post("      hp n\370%d: %f.x + %f.y.", hp+1, x->dist[hp]*cos(x->teta[hp]), 
                                          x->dist[hp]*sin(x->teta[hp]));
}



/**********************************************************************/
/*                       FONCTION DESTRUCTION                         */
/**********************************************************************/

void ambipan_tilde_dest(t_ambipan_tilde *x)
{
  //Lib�ration de la m�moire allou�e pour les cosinus,
  if( !x->controle  )
    free( x->cos_teta );

  return;
}



/**********************************************************************/
/*                       FONCTION CREATION                            */
/**********************************************************************/

void *ambipan_tilde_new( t_symbol *s, int argc, t_atom *argv )
{  
  int    i, hp;
  char   chaine[256];

 
  /*Allocation de la dataspace******************************/
  t_ambipan_tilde *x = (t_ambipan_tilde *)pd_new(ambipan_tilde_class);


  /*********************************************************/
  /*R�cup�ration du nombre de sorties et de haut-parleurs. */
  if( argc >= 1 )
    x->N = (int)atom_getint( argv );
  else 
    x->N = Ndefaut;

  if( 16 < x->N || 2 > x->N ){
    x->N = Ndefaut;
    post("ambipan~: Il y a un probl\212me dans la d\202claration "
         "du nombre de haut-parleur, \n"
         "  il est de %d par d\202faut.", Ndefaut);
  }
  x->Nout = x->N;  //M�me nombre de sorties et d'haut-parleur.
    
  /*R�cup�ration du type rep�re ****************************/
  if( argc >=2 ){
    atom_string( argv+1, chaine, 4);

    if( chaine[0] == 'c' )
      x->base=1;
    else if( chaine[0] == 'p' )
      x->base=0;
    else {
      x->base=1;
      post("ambipan~: erreur dans le type des coordonn\202es, "
           "\n  elles sont cart\202siennes par d\202faut.");
    }
  }
  else 
    x->base = 1; //Cart�sienne par d�faut.
  
  /*R�cup�ration du type des entr�es ***********************/
  if( argc >=3 ){    
    atom_string( argv+2, chaine, 4);

    if( chaine[0] == 'c' )
      x->controle=1;
    else if( chaine[0] == 's' )
      x->controle=0;
    else {
      x->controle=1;
      post("ambipan~: erreur dans le type des entr\202es, \n"
           "  elles sont en controle par d\202faut.");
    }
  }
  else 
    x->controle = 1;
  
  /*R�cup�ration de l'offset *******************************/
  if( argc >=4 )
    x->offset = (float)fabs(atom_getfloat( argv+3));
  else
    x->offset = (float)OFFSET;
  if( x->offset <= EPSILON ){
    post("ambipan~: Pas d'offset n\202gatif ou nul s'il vous plait.");
    x->offset = (float)OFFSET;;
  }


  /*R�cup�ration du temps d'interpolation ******************/
  if( argc >=5 )
    x->dtime = (int)(abs(atom_getint( argv+4))*sys_getsr()/1000.);
  else
    x->dtime = (int)(DTIME*sys_getsr()/1000.);
  if( x->dtime == 0 ) //Pas de temps d'interpolation nul ou n�gatif.
    x->dtime = 1;



  /*********************************************************/
  /*Initialisation des donn�es de la classe. ***************/
  for( hp = x->N-1; hp >= 0; hp--) //Initialisation des P
    x->P[hp] = 0;

  x->mute= 0;             //entr�es signal non mut�es
  x->y   = Ydefaut;       //initialisation de y
  x->phi = Phidefaut;     //initialisation de phi

  //Initialisation des angles des haut-parleurs et des rayons,
  for( hp=0; hp<x->N; hp++)
  {
    x->teta[hp] = (float)( Pi*( .5 + (1 - 2*hp )/(float)x->N) );
    x->dist[hp] = 1;
  }

  /*********************************************************/
  /* Cr�ation des tableaux cosinus, cos_teta et sin_teta.  */
  /* pour l'audio.                                         */
  if( !x->controle )
  {
    x->cos_teta = malloc( (T_COS+2*x->N)*sizeof(float) );
    x->sin_teta = x->cos_teta+x->N;

    //Pr�calculs des cos et sin des haut-parleurs.
    for( hp=0; hp<x->N; hp++)
    {
      x->cos_teta[hp] = (float)cos( x->teta[hp] );
      x->sin_teta[hp] = (float)sin( x->teta[hp] );
    }

    //Remplissage du tableau cosinus,
    x->cosin = x->sin_teta+x->N;
    /*Pour avoir besoin de cos( phi ), on cherche:
      cos(phi) = 
      cosin[ (int)( phi*T_COS/360 ))&(((int)T_COS-1) ] */
    for( i=0; i<T_COS; i++)
      x->cosin[i] = (float)cos( i*2*Pi/T_COS );
  }
  
  /*********************************************************/
  //initialisation de x, X, Y, W et Pn par la m�thode recoit x.
  if( x->base )    ambipan_tilde_recoit_x( x, Xdefaut);
  else             ambipan_tilde_recoit_x( x, Rdefaut);     




  /*********************************************************/
  /*Cr�ation des nouvelles entr�es. ************************/
  //si controle on cr�e des entr�es de contr�le,
  if( x->controle )  {
    inlet_new(&x->x_obj, &x->x_obj.ob_pd,     
              gensym("float"), gensym("x"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, 
              gensym("float"), gensym("y"));
  }
  //Sinon on cr�e des entr�es de signal.
  else  {    
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, 
              gensym("signal"), gensym("signal"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, 
              gensym("signal"), gensym("signal"));
  }

  /*Cr�ation des sorties************************************/
  for( i=0; i < x->N; i++)
    outlet_new(&x->x_obj, gensym("signal"));

  return (void *)x;
}



/**********************************************************************/
/*                       FONCTION SETUP                               */
/**********************************************************************/

void ambipan_tilde_setup(void) 
{
  /****D�finition de la classe**********************************/
  ambipan_tilde_class = class_new( gensym("ambipan~"),
                                (t_newmethod) ambipan_tilde_new,
                                (t_method) ambipan_tilde_dest,
                                sizeof(t_ambipan_tilde),
                                CLASS_DEFAULT, A_GIMME, 0);

  /****D�finition des m�thodes**********************************/
  class_addmethod(ambipan_tilde_class,
        (t_method)ambipan_tilde_dsp, gensym("dsp"), 0);
  class_addmethod(ambipan_tilde_class, nullfn, gensym("signal"), 0);

  /*M�thodes pour le changement des coordonn�es de la source****/
  class_addmethod(ambipan_tilde_class, 
                 (t_method)ambipan_tilde_recoit_x, gensym("x"), A_DEFFLOAT, 0);
  class_addmethod(ambipan_tilde_class, 
                 (t_method)ambipan_tilde_recoit_y, gensym("y"), A_DEFFLOAT, 0);
  class_addlist(ambipan_tilde_class, 
                 (t_method)ambipan_tilde_recoit_liste);

  /*M�thodes pour le changement des coordonn�es des h-p********/
  class_addfloat (ambipan_tilde_class, (t_method)ambipan_tilde_recoit_x);
  class_addmethod(ambipan_tilde_class, 
                 (t_method)ambipan_tilde_teta_positionner_hp, 
                  gensym("a_setpos"), A_GIMME, 0);
  class_addmethod(ambipan_tilde_class, 
                 (t_method)ambipan_tilde_dist_teta_positionner_hp, 
                  gensym("ra_setpos"), A_GIMME, 0);
  class_addmethod(ambipan_tilde_class, 
                 (t_method)ambipan_tilde_xy_positionner_hp, 
                  gensym("xy_setpos"), A_GIMME, 0);
  class_addmethod(ambipan_tilde_class, 
                 (t_method)ambipan_tilde_initialiser_nb_hp, 
                  gensym("set_nb_hp"), A_DEFFLOAT, 0);
  class_addmethod(ambipan_tilde_class, 
                 (t_method)ambipan_tilde_muter_entrees_signal, 
                  gensym("mute_sig"), A_DEFFLOAT, 0);
  class_addmethod(ambipan_tilde_class, 
                 (t_method)ambipan_tilde_changer_offset, 
                  gensym("set_offset"), A_DEFFLOAT, 0);
  class_addmethod(ambipan_tilde_class, 
                 (t_method)ambipan_tilde_changer_type_repere, 
                  gensym("change_type"), A_GIMME, 0);
  class_addmethod(ambipan_tilde_class, 
                 (t_method)ambipan_tilde_informations , 
                  gensym("get_info"), A_GIMME, 0);
  class_addfloat (ambipan_tilde_class, (t_method)ambipan_tilde_recoit_x);
  

  /****Cr�ation d'un lien pour l'aide****************************/
  class_sethelpsymbol(ambipan_tilde_class, gensym("ambipan~-help"));

  
  /*Publicit� pour le CICM**********************************/
  post("\"AMBIPAN~\"  Auteurs: R.MIGNOT  "
       "et B.COURRIBET, CICM Universit\202 Paris8,");
  post("           MSH Paris Nord, ACI Jeunes "
       "Chercheurs \"Espaces Sonores\".");
  post("           cicm\100univ-paris8.fr.");
}




/**********************************************************************/
/*                       METHODE DSP                                  */
/**********************************************************************/

void ambipan_tilde_dsp(t_ambipan_tilde *x, t_signal **sp)
{
  if( x->controle )
  {
    switch( x->Nout )
    {
    case 2:
      dsp_add(ambipan_tilde_perform_controle, 5, x,
              sp[0]->s_vec,   //son d'entr�e,
              sp[1]->s_vec,   //tous les sons de sortie, ...
              sp[2]->s_vec,   // .........
              sp[0]->s_n);    //taille des blocs.
      break;
    case 3:
      dsp_add(ambipan_tilde_perform_controle, 6, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, 
              sp[3]->s_vec, sp[0]->s_n); 
      break;      
    case 4:
      dsp_add(ambipan_tilde_perform_controle, 7, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[0]->s_n); 
      break;      
    case 5:
      dsp_add(ambipan_tilde_perform_controle, 8, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[0]->s_n); 
      break;      
    case 6:
      dsp_add(ambipan_tilde_perform_controle, 9, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[0]->s_n); 
      break;      
    case 7:
      dsp_add(ambipan_tilde_perform_controle, 10, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[0]->s_n); 
      break;
    case 8:
      dsp_add(ambipan_tilde_perform_controle, 11, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[0]->s_n); 
      break;
    case 9:
      dsp_add(ambipan_tilde_perform_controle, 12, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[9]->s_vec, sp[0]->s_n); 
      break;
    case 10:
      dsp_add(ambipan_tilde_perform_controle, 13, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[9]->s_vec, sp[10]->s_vec, sp[0]->s_n); 
      break;
    case 11:
      dsp_add(ambipan_tilde_perform_controle, 14, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[0]->s_n); 
      break;
    case 12:
      dsp_add(ambipan_tilde_perform_controle, 15, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[0]->s_n); 
      break;
    case 13:
      dsp_add(ambipan_tilde_perform_controle, 16, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[0]->s_n); 
      break;
    case 14:
      dsp_add(ambipan_tilde_perform_controle, 17, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[14]->s_vec, 
              sp[0]->s_n); 
      break;
    case 15:
      dsp_add(ambipan_tilde_perform_controle, 18, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[14]->s_vec, 
              sp[15]->s_vec, sp[0]->s_n); 
      break;
    case 16:
      dsp_add(ambipan_tilde_perform_controle, 19, x,
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
    case 2:
      dsp_add(ambipan_tilde_perform_signal, 7, x,
              sp[0]->s_vec,   //son d'entr�e,
              sp[1]->s_vec,   //entr�e de x,
              sp[2]->s_vec,   //entr�e de y,
              sp[3]->s_vec,   //tous les sons de sortie, ...
              sp[4]->s_vec,   // .........
              sp[0]->s_n);    //taille des blocs.
      break;
    case 3:
      dsp_add(ambipan_tilde_perform_signal, 8, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[0]->s_n);  
      break;      
    case 4:
      dsp_add(ambipan_tilde_perform_signal, 9, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[0]->s_n);  
      break;      
    case 5:
      dsp_add(ambipan_tilde_perform_signal, 10, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[0]->s_n);  
      break;      
    case 6:
      dsp_add(ambipan_tilde_perform_signal, 11, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[0]->s_n);  
      break;      
    case 7:
      dsp_add(ambipan_tilde_perform_signal, 12, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[0]->s_n);  
      break;      
    case 8:
      dsp_add(ambipan_tilde_perform_signal, 13, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[0]->s_n);  
      break;      
    case 9:
      dsp_add(ambipan_tilde_perform_signal, 14, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[0]->s_n);  
      break;
    case 10:
      dsp_add(ambipan_tilde_perform_signal, 15, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[11]->s_vec, sp[0]->s_n);  
      break;
    case 11:
      dsp_add(ambipan_tilde_perform_signal, 16, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[0]->s_n);  
      break;
    case 12:
      dsp_add(ambipan_tilde_perform_signal, 17, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[14]->s_vec, 
              sp[0]->s_n);  
      break;
    case 13:
      dsp_add(ambipan_tilde_perform_signal, 18, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[14]->s_vec, 
              sp[15]->s_vec, sp[0]->s_n);  
      break;
    case 14:
      dsp_add(ambipan_tilde_perform_signal, 19, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[14]->s_vec, 
              sp[15]->s_vec, sp[16]->s_vec, sp[0]->s_n);  
      break;
    case 15:
      dsp_add(ambipan_tilde_perform_signal, 20, x,
              sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
              sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, 
              sp[6]->s_vec, sp[7]->s_vec, sp[8]->s_vec, 
              sp[9]->s_vec, sp[10]->s_vec, sp[11]->s_vec, 
              sp[12]->s_vec, sp[13]->s_vec, sp[14]->s_vec, 
              sp[15]->s_vec, sp[16]->s_vec, sp[17]->s_vec, 
              sp[0]->s_n);  
      break;
    case 16:
      dsp_add(ambipan_tilde_perform_signal, 21, x,
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
