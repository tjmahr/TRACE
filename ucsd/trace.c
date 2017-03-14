#include "defs.h"

FILE *altiop = NULL;
FILE *infp = stdin;
FILE *outfile = NULL;
char outname[30];
char *re_comp();
int ecomp();

struct word *wordptr[NWORDS + 1];
struct word *wtord[NWORDS + 1];
float *exptr[PSLICES];
int nwptrs;
int printmin, printmax;
float windowcent = 66 + 12;		/* make float for getval routine */

char wordsound[NWORDS][MAXLEN];
int wordfreq[NWORDS];
struct phoneme phoneme[NPHONS];
struct phoneme *phonptr[128]; /* pointers to phoneme structures indexed 
			   by the character corresponding to the phoneme */

struct feature feature[NCONTINS][NFVALS];

int contlen[NCONTINS] = 
	{NFVALS,NFVALS,NFVALS,NFVALS,NFVALS,NFVALS,NFVALS};

char *phonchar =
	{"pbtdkgsSrlaiu^-"};

char *classchar =
	{"VCLNDGBT"}; /*vocalic, consonantal,liquid,noncompact,dental stop,
			velar stop, bilabial stop, bil/alv  */


char *contname[NCONTINS] = 
{"POW","VOC","DIF","ACU","GRD","VOI","BUR"};

float fweight[NCONTINS] = 
  {1., 1., 1., 1., 1., 1., 1. };

int nwords = 0;

char in[BUFSIZ];

float max = { 1.};
float min = {-.3};
float wbase = {.3};
float wchange = {.1};
float wgain = {1.};
float wramp = {1.};
float rate = {.5};

float ymax = .6;
float ydec = .01;
float pthresh[NLEVS]    = {.00,.00,.05};
float linthresh = .2; 
float linthrsq = .04;
float decay[NLEVS]      = {.02,.02,.02}; 
float alpha[NPARAMS] = 
/*   IF   FF    FP  PP  PW  WW  WP  PF  PFC  */
    {1.0,.00,.020,.00,.10,.00,.15,.00,.00};

float ga[NLEVS] = 
/*   FF    PP     WW */
   {.100, .100, .100};

float rest[NLEVS]   = {-.1, -.1, -.01};

float fscale = {.00};
float imax = {3.};

float psum[PSLICES];
float pisum[PSLICES];
int pcount[PSLICES];

float fsum[NCONTINS][FSLICES];
int fcount[NCONTINS][FSLICES];

float wisum[PSLICES];
float wsum[PSLICES];
int wcount[PSLICES];

float ppi[PSLICES];
float wwi[PSLICES];
float ffi[NCONTINS][FSLICES];
int cycleno;

int nreps = 9;
int sigflag;
int echeck;
int outflag;
int sumflag = 1;
int echoflag = 1;
int grace = 0;
int coartflag = 1;
int printfreq = 3;

int stop();

struct phoneme 	*inphone[MAXITEMS];
char 	inspec[MAXITEMS][60];
float 	infet[MAXITEMS][NCONTINS][NFVALS];

float 	instrength[MAXITEMS];
int 	inpeak[MAXITEMS];
int  	indur[MAXITEMS];
float	inrate[MAXITEMS];
float	insus[MAXITEMS];

struct word 	*outword[MAXITEMS];
struct phoneme 	*outphoneme[MAXITEMS];

int ninputs, nwouts, npouts;

float input[FSLICES][NCONTINS][NFVALS];

main(argc,argv)
int argc; char **argv;
{
	int specflag = 0;
	while (--argc > 0 && (*++argv) [0] == '-') {
		switch (argv[0][1]) {
			case 'p':
				load(*++argv);
				argc--;
				break;
			case 'l':
				lexicon(*++argv);
				argc--;
				break;
			case 'i':
				inspecs(*++argv);
				specflag = 1;
				argc--;
				break;
			case 'f':
				infeats(*++argv);
				argc--;
				break;
			case 'a':	/* to make it easy for sdb */
				abort();
				break;
		}
	}
	initialize();
	if (specflag) inputcomp();
	if (signal(2,stop) & 01) signal(2,stop);
	command();
}

cycle() 
{
	int tick,pmin,pmax;
	for (tick = 0; tick < nreps; tick++) {
	    act_features();
	    fpinteract();
	    ppinteract();
	    pfinteract();
	    pwinteract();
	    wpinteract();
	    wwinteract();
	    fupdate();
	    pupdate();
	    wupdate();
	    cycleno++;
	    if ( !(cycleno%printfreq) ) {
			if (sumflag) {
				if (cycleno>windowcent) {
					printmin  = 6;
					printmax += 66+6;
				}
				summarypr(printmin,printmax,FPP,stdout);
				if (outfile) {
				    summarypr(printmin,printmax,FPP,outfile);
				}
			}
	    }
	    if (sigflag) {sigflag = 0;break;}
	}
}

zarrays() 
{
	WORD *wp, **wptr;
	PHONEME *pp;
	register int i,f,c;

	cycleno = 0;
	printmin = 0;
	printmax = 66;
	for (i = 0; i < PSLICES; i++) {
		wwi[i] =0; wsum[i] =0; wisum[i] =0; wcount[i] = 0;
		ppi[i] =0; psum[i] =0; pisum[i] = 0; pcount[i] = 0;
	}

	for (i = 0; i < FSLICES; i++) {
	    for (c = 0; c < NCONTINS; c++) {
		ffi[c][i] =0; fsum[c][i] =0; fcount[c][i] = 0;
	}   }

	for ALLWORDPTRS {
	    wp = *wptr;
		wp->pextot = 0;
	    wp->rest = wp->prime + wp->base;
	    for (i = 0; i < PSLICES; i++) {
		wp->ex[i] = wp->rest;
		wp->nex[i] =0;
		wp->pex[i] = 0;
		wp->wex[i] = 0;
	}   }

	for ALLPHONEMES {
		pp->fextot = 0;
		pp->wextot= 0;
	    for (i = 0; i < PSLICES; i++) {
		pp->ex[i] = rest[P];
		pp->wex[i] = 0; pp->fex[i] = 0; pp->nex[i] = 0;
		pp->pex[i] = 0;
	}   }

	for (c = 0; c < NCONTINS; c++) {
	    for (f=0; f < NFVALS; f++) {
		feature[c][f].pextot = 0;
		for (i = 0; i < FSLICES; i++) {
		    feature[c][f].ex[i] = rest[F];
		    feature[c][f].pex[i] = 0;
		    feature[c][f].fex[i] = 0;
		    feature[c][f].nex[i] = 0;
	}   }   }
}

