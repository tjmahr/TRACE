
#include "defs.h"

/* mostly routines concerned with feature input stuff -- reading,
   writing, compiling, etc included.  also routines concerned
   with changing defs of features and stuff.  these seem useless */

inputcomp()
{
    register int i,j,f,k;
    char *p;
	
    for (i = 0; i < ninputs; i++) {

	for (j = 0; j < NCONTINS; j++)  {
	    for (f = 0; f < NFVALS; f++) {
		infet[i][j][f] = 0;
	    }
	}
	if ( (k = tindex(classchar,inspec[i])) >= 0) {
				/* a class specification */
	    for (j = 0; j < NCONTINS; j++) {
		for (f = 0; f < NFVALS; f++) {
		    infet[i][j][f] = classval[k][j][f];
		}
	    }
	    indur[i] = PWIDTH;
	    insus[i] = 1.0;
	} else {			   /*simply a phoneme */
	    if ( (k = tindex(phonchar,inspec[i])) < 0 ) {
		continue;
	    }
	    inphone[i] = phonptr[phonchar[k]];
	    indur[i] = phonptr[phonchar[k]]->fdur;
	    insus[i] = phonrsus[k];
	    for (j = 0; j < NCONTINS; j++) {
		for (f = 0; f < NFVALS; f++) {
		    infet[i][j][f] = fetval[k][j][f];
		}
	    }
	}
    }
    fillinput();
}

fillinput()
{
	register int inp, c, fet, cyc;
	int midp;
	float t,spread,r_sq,weight;

	for (cyc = 0; cyc < FSLICES; cyc++) {
	  for (c = 0; c < NCONTINS; c++) {
	    for (fet = 0; fet < NFVALS; fet++) {
		input[cyc][c][fet] = 0;
	} } }

	for (inp = 0; inp < ninputs; inp++) {
	    if (inp) {
		inpeak[inp] = inpeak[inp-1] +
		    inrate[inp-1]*insus[inp-1]*indur[inp-1]/2.0
		    + inrate[inp]*insus[inp]*indur[inp]/2.0;
	    }
	    else {
		inpeak[inp] =
		    inrate[inp]*insus[inp]*indur[inp];
	    }
	    for (c = 0; c < NCONTINS; c++) {
		for (fet = 0; fet < NFVALS; fet++) {
		    t = infet[inp][c][fet];
		    spread = inrate[inp]*fetspread[c];
		    for (cyc = nonneg(1 + inpeak[inp] - spread);
			cyc < (inpeak[inp]+spread) && cyc < FSLICES; cyc++) {
			weight = 1 - (abs(inpeak[inp] - cyc)/spread);
			input[cyc][c][fet] += weight*t*instrength[inp];
	}   }   }   }
}


inspecs(str) char *str; 
{
    PHONEME *pp;
    register int i,j;
    int count;
    FILE *iop;
    char phchar[BUFSIZ];

    if (str[0] == '-') iop = infp;
    else if ( (iop = fopen(str,"r")) == NULL) {
	printf("Can't open %s.\n",str);
	BACKTOTTY
	return;
    }

    if (iop == infp) {
	ttyprint("Give one phoneme,strength,rate/line (end with ^D):\n");
    }

    for (ninputs = 0; ninputs < MAXITEMS; ninputs++) {
getphoneme:
	count = fscanf(iop,"%s %f %f",
	    phchar,&instrength[ninputs],&inrate[ninputs]);
	if (count == EOF) break;
	strcpy(inspec[ninputs],phchar);
    }
    clearerr(iop);
    if (iop != infp)
	fclose(iop);
    else
	clearerr(infp);
    inputcomp();
}

sinspec() 
{
	FILE *iop;
	int i;
		getstring("Give filename for spec list: ");
		if ( (iop = fopen(in,"w")) == NULL) {
			printf ("Can't open %s.\n",in);
			BACKTOTTY
			return;
		}
		for (i = 0; i < ninputs; i++) {
			fprintf(iop,"%s %f %f\n",
			  inspec[i],instrength[i],inrate[i]);
		}
		fprintf(iop,"q\n");
		fclose(iop);
}

dinspec() 
{
    int i;
    for (i = 0; i < ninputs; i++) {
	printf("%s %f %f %d %d\n",
	  inspec[i],instrength[i],inrate[i],inpeak[i],indur[i]);
    }
}

sfeatures() 
{
    FILE *iop;
    register int i,f,c;
    if ( (iop = fopen(in,"w")) == NULL) {
	printf ("Can't open %s.\n",in);
	BACKTOTTY
	return;
    }
    for (i = 0; i < FSLICES; i++) {
      for (c = 0; c < NCONTINS; c++) {
	for (f = 0; f < NFVALS; f++) {
	    fprintf(iop,"%8.5f ",input[i][c][f]);
	}
	fprintf(iop,"\n");
      }
    }
    fclose(iop);
}

infeats(str) char *str;
{
	register int i,c,f;
	FILE *iop;

	if ( (iop = fopen(str,"r")) == NULL) {
		printf("Can't open %s.\n",in);
		BACKTOTTY
		return;
	}
	for (i = 0; i < FSLICES; i++) {
	  for (c = 0; c < NCONTINS; c++) {
	    for (f = 0; f < NFVALS; f++) {
		fscanf(iop,"%f ",&input[i][c][f]);
	    }
	    fscanf(iop,"\n");
	  }
	}
	fclose(iop);
}
		
acoustinput(pmin,pmax,pstep,iop)
int pmin, pmax,pstep; 
FILE *iop;
{
	register int i,j,f,cy;
	if (pmax == 0) pmax = cycleno;
	if (pstep == 0) pstep = 1;
	fprintf(iop,"\t\t\tACOUSTIC FEATURES\n");
	fprintf(iop,"Slice->\t");
	for (cy = pmin; cy <= pmax; cy += pstep) {
		fprintf(iop,"%3d",cy);
	}
	fprintf(iop,"\n");
	for (j = 0; j < NCONTINS; j++) {
	    for (f = 0; f < NFVALS; f++) {
		for (cy = pmin; cy <= pmax; cy += pstep) {
		    if (input[cy][j][f] > pthresh[F]) {
			break;
		    }
		}
		if (cy == pmax+pstep) continue;
		fprintf (iop,"%s.%d\t",contname[j],f);
		for (cy = pmin; cy <= pmax; cy += pstep) {
		    if (input[cy][j][f] <= pthresh[F]) 
			fprintf(iop,"   ");
		    else
			fprintf (iop,"%3d",(int) (input[cy][j][f]*100));
		}
		fprintf(iop,"\n");
	    }
	    if (sigflag) {sigflag = 0; return;}
	}
}


/* get new set of feature definitions for the phonemes */
getfeats()
{
    getstring("feature file name: ");
    if (features(in))	{
	initialize();
    }
}	

features(filename)
char *filename;
{
    FILE *iop;
    register int c,f,p;

    if ( (iop = fopen(filename,"r")) == NULL) {
	printf("Can't open %s.\n",filename);
	BACKTOTTY
	return(0);
    }
    
    while (fgets(in,BUFSIZ,iop) != NULL) {
	for (p=0; p<NPHONS; p++) {
	  for (c = 0; c < NCONTINS;c++) {
	    for (f=0; f<NFVALS; f++) {
		sscanf(in,"%f",&fetval[p][c][f]);
	    }	
	    sscanf(in,"\n");
	  }
	}
    }
    fclose(iop);
    return(1);
}

test(str)  char *str; {
	char string[20], word[20];
	strcpy(string,in);
	getstring("words: ");
	strcpy(word,in);
	makein(string);
	zarrays();
	cycle();
}

makein(word) char *word; {
	strcpy(inspec[0],"-");
	instrength[0] = inrate[0] = 1.;
	for (ninputs = 1; ninputs < MAXITEMS-1 && *word; ) {
		if (*word == '.') {
			*word++;
			continue;
		}
		inspec[ninputs][0] = *word++;
		inspec[ninputs][1] = '\0';
		instrength[ninputs] = inrate[ninputs] = 1.0;
		ninputs++;
	}
	strcpy(inspec[ninputs],"-");
	instrength[ninputs] = inrate[ninputs] = 1.;
	ninputs++;
	inputcomp();
}
