#include "defs.h"

int contlen[NCONTINS] = 
	{NFVALS,NFVALS,NFVALS,NFVALS,NFVALS,NFVALS,NFVALS};

char *phonchar =
	{"pbtdkgsSrlaiu^-"};

char *classchar =
	{"VLDAPBU"}; 
/* Vowel, Liq, Diffuse, Alv-stop, Palat-stop,Bilab-stop,Unv-stop(p/t)*/

char *contname[NCONTINS] = 
{"POW","VOC","DIF","ACU","CNS","VOI","BST"};

double fweight[NCONTINS] = 
  {1., 1., 1., 1., 1., 1., 1. };

double fetval[NPHONS][NCONTINS][NFVALS] = {

/* each phoneme contains specifications for
   the extent to which each value of each feature on
   each continuum tends to excite it;
   continua occur in the order:

   cont:  POWER ,  VOCALIC ,  DIFFUSE ,  ACUTE ,  CONSONANTAL ,  VOICED, BST*/



/*p*/	0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0   , 0 , /* POW */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1.  , 0 , /* VOC */
	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* DIF */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0   , 0 , /* ACU */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* CNS */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1.  , 0 , /* VOI */
	1  ,  .2 ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* BST */


/*b*/	0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0   , 0 , /* POW */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1.  , 0 , /* VOC */
	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* DIF */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0   , 0 , /* ACU */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* CNS */
	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOI */
	.2 ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* BST */

/*t*/	0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0   , 0 , /* POW */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1.  , 0 , /* VOC */
        0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* DIF */
        0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* ACU */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* CNS */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1.  , 0 , /* VOI */
	0  ,  0  ,  1. ,  .2 ,  0  ,  0  ,  0  ,  0   , 0 , /* BST */

/*d*/	0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0   , 0 , /* POW */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1.  , 0 , /* VOC */
	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* DIF */
        0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* ACU */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* CNS */
	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOI */
	0  ,  0  ,  .2 ,  1. ,  0  ,  0  ,  0  ,  0   , 0 , /* BST */

/*k*/	0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0   , 0 , /* POW */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1.  , 0 , /* VOC */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0   , 0 , /* DIF */
	0  ,  0  ,  0  , .1  , .3  ,  1. , .3  , .1   , 0 , /* ACU */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* CNS */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1.  , 0 , /* VOI */
	0  ,  0  ,  0  ,  0  ,  1. ,  .2 ,  0  ,  0   , 0 , /* BST */

/*g*/	0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0   , 0 , /* POW */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1.  , 0 , /* VOC */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0   , 0 , /* DIF */
	0  ,  0  ,  0  , .1  , .3  ,  1. , .3  , .1   , 0 , /* ACU */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* CNS */
	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOI */
	0  ,  0  ,  0  ,  0  ,  .2 ,  1. ,  0  ,  0   , 0 , /* BST */

/*s*/	0  ,  0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* POW */
	0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0   , 0 , /* VOC */
	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* DIF */
	1. , .3  , .1  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* ACU */
	0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0  ,  0   , 0 , /* CNS */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1.  , 0 , /* VOI */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* BST */

/*S*/	0  ,  0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* POW */
	0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0   , 0 , /* VOC */
	0  ,  0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* DIF */
	0  ,  0  , .1  , .3  ,  1. , .3  , .1  ,  0   , 0 , /* ACU */
	0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0  ,  0   , 0 , /* CNS */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1.  , 0 , /* VOI */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* BST */

/*r*/	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* POW */
	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOC */
	0  ,  0  ,  0  ,  0  ,  0  ,  0. ,  .5 ,  1.  , 0 , /* DIF */
        0  ,  0  ,  0  ,  0  ,  .25,  .5 ,  .75,  0   , 0 , /* ACU */
	0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0   , 0 , /* CNS */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOI */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* BST */

/*l*/	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* POW */
	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOC */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  .5  , 0 , /* DIF */
        0  ,  0  ,  0  ,  0  ,  .75,  .5 ,  .25,  0   , 0 , /* ACU */
	0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0   , 0 , /* CNS */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOI */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* BST */

/*a*/	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* POW */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOC */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0   , 0 , /* DIF */
	0  ,  0  ,  0  ,  0  ,  0  , .1  , .3  ,  1   , 0 , /* ACU */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1.  , 0 , /* CNS */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOI */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* BST */

/*i*/	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* POW */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOC */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* DIF */
	1. , .3  , .1  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* ACU */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1.  , 0 , /* CNS */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOI */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* BST */

/*u*/	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* POW */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOC */
	0  ,  0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* DIF */
	0  ,  0  ,  0  ,  0  , .1  , .3  ,  1. , .3   , 0 , /* ACU */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1.  , 0 , /* CNS */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOI */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* BST */

/*^*/	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* POW */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOC */
	0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0  ,  0   , 0 , /* DIF */
	0  ,  0  ,  0  ,  0  ,  0  , .1  , .3  ,  1   , 0 , /* ACU */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1.  , 0 , /* CNS */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOI */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* BST */

/*-*/   0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , 0 ,  1.  , /* POW */
        0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , 0 ,  1.  , /* VOC */
        0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , 0 ,  1.  , /* DIF */
        0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , 0 ,  1.  , /* ACU */
        0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , 0 ,  1.  , /* CNS */
        0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , 0 ,  1.  , /* VOI */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , 0 ,  1.    /* BST */
};

double classval[NCLASSES][NCONTINS][NFVALS] = {
/*vow*/	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0 , /* POW */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0 , /* VOC */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0 , /* DIF */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0 , /* ACU */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0 , /* CNS */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0 , /* VOI */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* BST */

/*L*/	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* POW */
	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOC */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  .5 , .5   , 0 , /* DIF */
	0  ,  0  ,  0  ,  0  ,  .5 ,  .5 ,  .5 ,  0   , 0 , /* ACU */
	0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0   , 0 , /* CNS */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* VOI */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0   , 0 , /* BST */

/*dif*/	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0 , /* POW */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0 , /* VOC */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0 , /* DIF */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0 , /* ACU */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0 , /* CNS */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0 , /* VOI */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0 , /* BST */
	
/*as*/	0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0  ,  0  , /* POW */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0  , /* VOC */
	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , /* DIF */
        0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , /* ACU */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , /* CNS */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , /* VOI */
	0  ,  0  , .5  , .5  ,  0  ,  0  ,  0  ,  0   , 0 , /* BST */

/*vs*/	0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0  ,  0  , /* POW */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0  , /* VOC */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0  , /* DIF */
	0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0  , /* ACU */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , /* CNS */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , /* VOI */
	0  ,  0  ,  0  ,  0  , .5  , .5  ,  0  ,  0   , 0 ,  /* BST */

/*bs*/	0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0  ,  0  , /* POW */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0  , /* VOC */
	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , /* DIF */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0  , /* ACU */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , /* CNS */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , /* VOI */
	.5 ,  .5 ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , /* BST */

/*p/t*/	0  ,  0  ,  0  ,  0  ,  1. ,  0  ,  0  ,  0  ,  0  , /* POW */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0  , /* VOC */
	0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , /* DIF */
	0  ,  0  ,  0  , .5  , .5  ,  0  ,  0  ,  0  ,  0  , /* ACU */
	1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , /* CNS */
	0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1. ,  0  , /* VOI */
        1. ,  0  ,  1. ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  , /* BST */

};

int fetspread[NCONTINS] = {
/* for each dimension, how far it spreads in feature slices */
/* measured in FSLICES 
   POWER   VOCALIC   DIFFUSE   ACUTE   CONSONANTAL   VOICED    BST. */
   3,	    3,		3,	 3,	    3,      3,	   3
};

double wscale[NPHONS] = {
/* measured proporational to the basic phoneme duration */
 1.0, /*p*/
 1.0, /*b*/
 1.0, /*t*/
 1.0, /*d*/
 1.0, /*k*/
 1.0, /*g*/
 1.0, /*s*/
 1.0, /*S*/
 1.0, /*r*/
 1.0, /*l*/
 1.0, /*a*/
 1.0, /*i*/
 1.0, /*u*/
 1.0, /*^*/
 1.0  /*-*/
 };

double phonrsus[NPHONS] = {
 1., /*p*/
 1., /*b*/
 1., /*t*/
 1., /*d*/
 1., /*k*/
 1., /*g*/
 1., /*s*/
 1., /*S*/
 1., /*r*/
 1., /*l*/
 1., /*a*/
 1., /*i*/
 1., /*u*/
 1., /*^*/
 1.  /*-*/
 };

