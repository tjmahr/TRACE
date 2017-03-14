#include "defs.h"

char *optname[] =	{
	"acoustinput","alpha","cycle","decay","dinspec",
	"expression","echooutput","featacts","fileinput",
	"fweight","freqscale","gamma","grace",
	"help(?)","inspecs","infeatures","imax",
	"lexicon","min","max","mprime","nreps","output",
	"ofile","owtacts","phonacts","pfreq",
	"parameters","pthresh","quit","rc","rest","summarypr","sumopt",
	"spointers","sinspec","sfeatures","test","topdown","windowcent",
	"wordacts","wtacts","wbase","wchange","wgain","wramp",
	"abort","coart", "!(systemcall)", NULL};


command()
{
float tval;
int i,echeck;
int pmin,pmax,pstep;
char tbuf[20];

while (echeck = getstring(": ")) {
	if (echeck == EOF) {
		fileinput(STD);
	}
	else if (startsame(in,"cycle")) {
		cycle();
	}
	else if (startsame(in,"coarticulation")) {
		ttyprint("coarticulation %s: change? ",(coartflag ? "ON":"OFF"));
		fscanf(infp,"%s",in);
		if(in[0] == 'n') continue;
		coartflag = 1 - coartflag;
	}
	else if (startsame(in,"rc")) {
		printf("reset and cycle\n");
		zarrays();
		cycle();
	}
	else if (startsame(in,"wordacts")) {
		printwords(printmin,printmax,3,stdout);
		if (outfile) printwords(printmin,printmax,3,outfile);
	}
	else if (startsame(in,"wtacts")) {
		printwt(printmin,printmax,3,stdout);
		if (outfile) printwt(printmin,printmax,3,outfile);
	}
	else if (startsame(in,"owtacts")) {
		getstring("word: ");
		strcpy(tbuf,in);
		printonewt(tbuf,printmin,printmax,3,stdout);
		if (outfile) printonewt(tbuf,printmin,printmax,3,outfile);
	}
	else if (startsame(in,"phonacts")) {
		printphonemes(printmin,printmax,3,stdout);
		if (outfile) 
		    printphonemes(printmin,printmax,3,outfile);
	}
	else if (startsame(in,"featacts")) {
		printf("feature activations - give min,max,step: ");
		scanf("%d %d %d",&pmin,&pmax,&pstep);
		printfeatures(pmin,pmax,pstep,stdout);
		if (outfile) printfeatures(pmin,pmax,pstep,outfile);
	}
	else if (startsame(in,"expression")) {
		setex();
	}
	else if (startsame(in,"fileinput")) {
		fileinput(NONSTD);
	}
	else if (startsame(in,"inspecs")) {
		ttyprint("File name (- = stdin): ");
		fscanf(infp,"%s",in);
		inspecs(in);
	}
	else if (startsame(in,"infeatures")) {
		ttyprint("File name: ");
		fscanf(infp,"%s",in);
		infeats(in);
	}
	/* NOT PRESENTLY OPERATIVE -- JLM 10-5-82
	else if (startsame(in,"wsubset")) {
		wordsubset();
	}
	*/
	else if (startsame(in,"test")) {
		getstring("test string: ");
		test(in);
	}
	else if (startsame(in,"topdown")) {
		topdown();
	}
	else if (startsame(in,"output")) {
		setout();
	}
	else if (startsame(in,"ofile")) {
		getstring("give filename (or - for none): ");
		setoutfile(in);
	}
	else if (in[0] == '!') {
		system(&in[1]);
	}
	else if (in[0] == '?') {
		help();
	}
	else if (startsame(in,"help")) {
		help();
	}
	else if (startsame(in,"lexicon")) {
		getlex();
	}
	else if (startsame(in,"parameters")) {
		getpars();
	}
	else if (startsame(in,"quit")) {
		quit();
	}
	else if (startsame(in,"decay")) {
		ttyprint ("decay values:\n");
		getfloat(decay,NLEVS,"F","P","W");
	}
	else if (startsame(in,"alpha")) {
		ttyprint ("alpha values:\n");
		getfloat(alpha,NPARAMS,
		 "IF","FF","FP","PP","PW","WW","WP","PF","PFC");
	}
	else if (startsame(in,"gamma")) {
		ttyprint ("gamma values:\n");
		getfloat(ga,NLEVS,"F","P","W");
	}
	else if (startsame(in,"grace")) {
		ttyprint ("grace %s: change? ",(grace ? "ON" : "OFF"));
		fscanf(infp,"%s",in);
		if (in[0] == 'n') continue;
		grace = 1 - grace;
	}
	else if (startsame(in,"rest")) {
		tval = rest[W];
		ttyprint ("rest values:\n");
		getfloat(rest,NLEVS,"F","P","W");
		if (tval != rest[W]) {
			initialize();
		}
	}
	else if (startsame(in,"fweight")) {
		ttyprint ("fweight values:\n");
		getfloat(fweight,NCONTINS,
		"POW","PIT","VOI","ABT","DIF","ACU","CON","LWS",
		"MDS", "HIS","FT1","FT2","FT3","FT4","FT5","FT6");
	}
	else if (startsame(in,"pthresh")) {
		ttyprint ("pthresh values:\n");
		getfloat(pthresh,NLEVS,"F","P","W");
	}
	else if (startsame(in,"nreps")) {
		ttyprint ("%d nreps: change? ",nreps);
		fscanf(infp,"%s",in);
		if (in[0] == 'n') continue;
		sscanf(in,"%d",&nreps);
	}
	else if (startsame(in,"pfreq")) {
		ttyprint ("%d pfreq: change? ",printfreq);
		fscanf(infp,"%s",in);
		if (in[0] == 'n') continue;
		sscanf(in,"%d",&printfreq);
	}
	else if (startsame(in,"summarypr")) {
		ttyprint("give min, max, and step: ");
		fscanf(infp,"%d %d %d",&pmin,&pmax,&pstep);
		summarypr(pmin,pmax,pstep,stdout);
		if (outfile) summarypr(pmin,pmax,pstep,outfile);
	}
	else if (startsame(in,"spointers")) {
		qsort(wordptr,nwptrs,sizeof(struct word *),ecomp);
	}
	else if (startsame(in,"sinspec")) {
		sinspec();
	}
	else if (startsame(in,"sfeatures")) {
		ttyprint("file name: ");
		fscanf(infp,"%s",in);
		sfeatures(in);
	}
	else if (startsame(in,"dinspec")) {
		dinspec();
	}
	else if (startsame(in,"acoustinput")) {
		ttyprint("give min, max, and step: ");
		fscanf(infp,"%d %d %d",&pmin,&pmax,&pstep);
		acoustinput(pmin,pmax,pstep,stdout);
		if (outfile) acoustinput(pmin,pmax,pstep,outfile);
	}
	else if (startsame(in,"sumopt")) {
		ttyprint ("summaryprint %s: change? ",(sumflag ? "ON" : "OFF"));
		fscanf(infp,"%s",in);
		if (in[0] == 'n') continue;
		sumflag = 1 - sumflag;
	}
	else if (startsame(in,"echooutput")) {
		ttyprint ("echooutput %s: change? ",(echoflag ? "ON" : "OFF"));
		fscanf(infp,"%s",in);
		if (in[0] == 'n') continue;
		echoflag = 1 - echoflag;
	}
	else if (startsame(in,"min")) {
		ttyprint ("min = ");
		getval(&min);
	}
	else if (startsame(in,"max")) {
		ttyprint ("max = ");
		getval(&max);
	}
	else if (startsame(in,"windowcent")) {
		ttyprint ("windowcent = ");
		getval(&windowcent);
	}
	else if (startsame(in,"wbase")) {
		ttyprint ("wbase = ");
		getval(&wbase);
	}
	else if (startsame(in,"wchange")) {
		ttyprint ("wchange = ");
		getval(&wchange);
	}
	else if (startsame(in,"wgain")) {
		ttyprint ("wgain = ");
		getval(&wgain);
	}
	else if (startsame(in,"wramp")) {
		ttyprint ("wramp = ");
		getval(&wramp);
	}
	else if (startsame(in,"imax")) {
		ttyprint ("imax = ");
		getval(&imax);
	}
	else if (startsame(in,"freqscale")) {
		ttyprint ("freqscale = ");
		tval = fscale;
		getval(&fscale);
		if (tval != fscale) {
		    initialize();
		}
	}
	else if (startsame(in,"abort")) {
		abort();	/* to get a core dump for sdb */
	}
	else {
		printf("Unrecognized request: For help type ?.\n");
		if (infp != stdin) fileinput(STD);
	}
	wait(0);
}
}

