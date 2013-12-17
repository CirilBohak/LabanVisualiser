var lgetout = function ( allok ) {
/*
    close files and wait

	called by  linter, lcopyfigs, lcopysubs, ldopivot, lchange,
	           lfindstaff, lfindystart, lleggesture, lselectfig,
	           loverlap,
*/
	if ( allok == 0 ) {
		console.log(nudesname + " created OK.");
		if (infile) console.log("Close infile.");
		if (nudesfile) console.log("Close nudesfile.");
		if (figsfile) console.log("Close figsfile.");
	} else {
		console.log("lintel snag, line " + j);
		console.log(lbnline[j]);
	}
};/* lgetout */
/***********************************************************/

var initialise = function () { 
/*
   set up constants to default values

   called by main,
*/
	var a,b;
	var k,m,n;
	
	prdone = FALSE;
	nbar = -1;
	rise = 1;
	prevc = 0;
	previ = 11;
	track = TRUE;
	mspace = false;
	wspace = FALSE;
	
	doub0 = 0;
	doub1 = 1;
	doub2 = 2;
	doub3 = 3;
	doub4 = 4;
	doub10 = 10;
	doub60 = 60;
	doub90 = 90;
	doub150 = 150;
	doub179 = 179;
	doub180 = 180;
	doub181 = 181;
	doub255 = 255;
	doub360 = 360;
	doub500 = 500;
	inv2 = doub1/doub2;
	inv3 = doub1/doub3;
	inv4 = doub1/doub4;
	inv5 = doub1/5.0;
	inv6 = doub1/6.0;
	inv10 = doub1/doub10;
	inv256 = doub1/256.0;
	inv1000 = doub1/1000.0;
	rt3 = Math.sqrt(doub3);
	piby2 = doub2*Math.atan(doub1);
	pi = piby2+piby2 ;
	twopi = pi+pi;
	radten = twopi/3600.0;
	radian = twopi/doub360;
	degree = doub1/radian;
	lg2 = Math.log(doub2);
	freeze = FALSE;
	forward = TRUE;
	single = DONE;
	pause = FALSE;
	shadow = TRUE;
	fnums = TRUE;
	bnums = TRUE;
	hold = NO;
	prevhold = -99;
	prev_time = -1;
	fstart = 0;
	fstop = 0;
	pstart = 0;
	pend = 0;
	fmax = 0;
	vstart = 0;
	vstop = FMAX;
	inmain = TRUE;
	start = -1;
	lline = 0;
	fast = 1;
	slow = 1;
	fslow = 1;
	njts = 0;
	nvars = 0;
	nfiles = 0;
	nvals = 0;
	axlen[0] = 1; axlen[1] = 1; axlen[2] = 1;
	
	for (  j = 0 ; j < EMAX ; ++ j ) {
		if ( j > 2) axlen[j] = -1;
		keylen[j] = 0;
		ellen[j] = -1;
		jntlen[j] = -1;
		fillen[j] = -1;
		figlen[j] = -1;
		sublen[j] = -1;
		varlen[j] = -1;
		knee[j] = 0 ;
		figell[j] = 0;
		ellfig[j] = 0;
		usevar[j] = 0;
		coel[j][0] = -1;
		coel[j][1] = -1;
		subact[j][0] = 0;
		subact[j][1] = 0;
		called[j] = FALSE;
		defined[j] = FALSE;
		val[j] = doub0 ;

		for (  k = 0 ; k < 3 ; ++ k ) {
			cen[j][k] = doub3;
			ax[j][k] = doub2;
			lim[j][k][0] = -doub360;
			lim[j][k][1] =  doub360;
			obs[k][0] = doub0;
			obs[k][1] = doub0;
			obs[k][2] = doub0;
			obs[k][k] = doub1;
		}
		
		col[j][0] = doub255;
		col[j][1] = doub150;
		col[j][2] = doub90;
		quat[j][0] = doub1;
		quat[j][1] = doub0;
		quat[j][2] = doub0;
		quat[j][3] = doub0;
		quat[j][4] = doub1;
		
		axnam[j] = "";
		tname[j] = "";
		fname[j] = "";
		ename[j] = "";
		jname[j] = "";
		vname[j] = "";
		sname[j] = "";
	}
	axnam[0][0] = 'x';
	axnam[1][0] = 'y';
	axnam[2][0] = 'z';
/*
     set all actions by default to stop
*/
	for (  j = 0 ; j < PMAX ; ++ j ) {
		type[j] = stop_keyword_code;
		frstart[j] = 0;
		frstop[j] = 0;
		distrn[j] = 0;
		cline[j] = 0;

		for (  k = 0 ; k < 6 ; ++ k )
			pf[j][k] = 0;
	}
/*
     artificially set up subroutine "nudes",
     file "dummy", figures "every" and "world",
     variable "variable", and ellipsoid "world"-
*/
	nsubs = 1;
	nfigs = 2;
	ne = 1;
	figell[0] = 0;
	figell[1] = 1;
	
	tname[0] = dummy;
	sname[0] = nudes;
	fname[0] = every;
	fname[1] = world;
	ename[0] = world;
	vname[0] = variab;
	
	fillen[0] = 5;
	sublen[0] = 5;
	figlen[0] = 5;
	figlen[1] = 5;
	ellen[0] = 5;
	varlen[0] = 6;
	ax[0][0] = doub1;
	ax[0][1] = doub1;
	ax[0][2] = doub1;
	df = 1;
	f = 0;
	nsph = SSTART;
	anglex = doub0; angley = doub0; anglez = doub0;
	tx = doub0; ty = doub0; tz = doub0;
	scale = doub1;
	alpha = doub3;
	t = 0;
	more = 1;
	ok = 0;

// find bits in double mantissa -	
	b = doub1;
	m = 0;
	for (a = inv2; doub1 + b > doub1 + a; a *= inv2) {
		b = a;
		++ m;
	}
	tolr = b+b;
	j = 2;
	n = 0;

// find bits in integer -
	for (k = 1; k < j; j += j)
	{
		k += k;
		++ n;
	}
	maxint = k;
	console.log("\n tolr " + tolr + " (" + m + " bits), maxint " + maxint + " (" + n + " bits)\n");
}; /* initialise */
/***************************************/

var bell = function ( number, delay) {
	var i, j;
	for ( i = 0; i < number; i++ ) {
		console.log( "\a" );
		for ( j = 0; j < delay; j++ );
	}
}; /* bell */
/***************************************/

var lfindnext = function( c, y1, y2) {
/*
   find next symbol in column 'c' in range 'y1' to 'y2'.

   called by ldostep,
*/
	var k;
	var q;
	var yy;
	
	q = -1;
	yy = y2;
	for (k = sstart; k < ssend; ++k) {
		if ((lbn[k].c == jc) && (lbn[k].y >= y1) && (lbn[k].y <= y2))	{
			if (lbn[k].y < yy)	{
				q = k;
				yy = lbn[k].y;
			}
		}
	}
	return q;
};   /* lfindnext */
/****************************************************/

var lsetframes = function () {
/*
    set the frames over which an action occurs :-
    fstart, fhalf, frange, fend.

    called by laction,
*/
	if (nbar < 1) {
		fstart = 0;
		frange = 1;
		fend = 1;
	} else {
		fstart = Math.floor(inv2 + lbn_fpp * jy-ystart);
		if (fstart < 1) fstart = 1;
		frange = Math.floor(inv2 + lbn_fpp * jh);
		if (frange < 1) frange = 1;
		fend = fstart + frange;
	}
	fhalf = fstart + frange/2;
	if (fend <= fstart) fend = fhalf+1;
	if (fhalf > fend) fend = fhalf+1;
	if (fend > fmax) fmax = fend;
}; /* lsetframes */
/************************************************/

var lcolx = function ( lcentre ) {
/*
    find column number of each symbol
	-5 = L arm
	-3 = L gesture
    -1 = L support
	 1 = R support
	 3 = R gesture
	 5 = R arm
	
    called by linter,
*/
	var k;
	var kc;
	var kwx;
	
	for (k = 0; k < nlabs; ++k) {
	kwx = lbn[k].x + (lbn[k].w/2);
	kc = (kwx - lcentre)/STEP;
	if (kwx < lcentre)
		--kc;
	else
		++kc;
	lbn[k].c = kc;
	}
}; /* lcolx */
/************************************************/


var lbnread = function () {
/*
   read .lbn laban score file

   called by linter,
*/
	var j;
	var i,x,y,s,w,h;
	var d;
	var m0,m1,m2,m3;
	
	j = 0;
	xmax = 0;
	xmin = 10000;
	var done = false;
	getFileFromServer(infile, function(text) {
		if (text === null) {
			console.log("lbnread oops\n");
		} else {
			lbnline = text.split("\n");
			while ( j < LMAX ) {
				buf = lbnline[j];
				sscanf(buf,"%c%c%c%c %d %d %d %d %d %d %c",
					m0,m1,m2,m3,i,x,y,s,w,h,d);
				if (m0 != '#') {
					lbn[j].m = m0;
					if ((m0 == 'P')&&(m1 == 'a'))
						lbn[j].m = Path;
					lbn[j].i = i;
					lbn[j].x = x;
					lbn[j].y = y;
					lbn[j].w = w;
					lbn[j].h = h;
					lbn[j].s = s;
					lbn[j].b = -1;
					lbn[j].l = j;
					lbn[j].a = TODO;
					lbn[j].x2 = x+w;
					lbn[j].y2 = y+h;
					lbn[j].d = BLANK;
					if (d =='M') lbn[j].d = MED;
					if (d =='L') lbn[j].d = LOW;
					if (d =='H') lbn[j].d = HIGH;
					if (x < xmin) xmin = x;
					if (x+w > xmax) xmax = x+w;
					if (j >= LMAX) {
					  console.log("\nBEWARE: score truncated at line "+ j +"\n");
					  console.log("more than "+ LMAX +" laban score items\n");
					}
					++j;
				}
			} /* while reading next line */
			done = true;
		}
	});
	while (done == false) {}
	nlabs = j;
	console.log("\n   lbnread: "+ nlabs +" lbn symbols\n");
};  /* lbnread */
/************************************************/

var lassign = function() {
/*
   assign global variables

   called by laction, lsorty, lbent,
*/
	ja = lbn[j].a;
	jb = lbn[j].b;
	jc = lbn[j].c;
	jd = lbn[j].d;
	jh = lbn[j].h;
	ji = lbn[j].i;
	jl = lbn[j].l;
	jm = lbn[j].m;
	js = lbn[j].s;
	jw = lbn[j].w;
	jx = lbn[j].x;
	jy = lbn[j].y;
	jx2 = lbn[j].x2;
	jy2 = lbn[j].y2;
}; /* lassign */
/**********************************************/

var lsorty = function () {
/*
   sort score symbols into ascending order of 'y'
   (bubble sort)
   find maxy, and fill yj table

   called by linter,
   calls     lassign,
*/
	var k;
	var last;
	var y;
	
	for (j = 0; j < (nlabs-1); ++j)
	{
	  for (k = j; k < nlabs; ++k)
	  {
		 if (lbn[k].y < lbn[j].y)
		 {
			lassign();
			lbn[j].a = lbn[k].a;
			lbn[j].b = lbn[k].b;
			lbn[j].c = lbn[k].c;
			lbn[j].d = lbn[k].d;
			lbn[j].h = lbn[k].h;
			lbn[j].i = lbn[k].i;
			lbn[j].l = lbn[k].l;
			lbn[j].m = lbn[k].m;
			lbn[j].s = lbn[k].s;
			lbn[j].w = lbn[k].w;
			lbn[j].x = lbn[k].x;
			lbn[j].y = lbn[k].y;
			lbn[j].x2 = lbn[k].x2;
			lbn[j].y2 = lbn[k].y2;
			lbn[k].a = ja;
			lbn[k].b = jb;
			lbn[k].c = jc;
			lbn[k].d = jd;
			lbn[k].h = jh;
			lbn[k].i = ji;
			lbn[k].l = jl;
			lbn[k].m = jm;
			lbn[k].s = js;
			lbn[k].w = jw;
			lbn[k].x = jx;
			lbn[k].y = jy;
			lbn[k].x2 = jx2;
			lbn[k].y2 = jy2;
			buf = lbnline[j];
			lbnline[j] = lbnline[k];
			lbnline[k] = buf;
		 }
	  }
	}
	ymax = 0;
	for (j = 0; j < nlabs; ++j)
		if (((lbn[j].y2) > ymax)&&(lbn[j].m != Stav))
			ymax = lbn[j].y2+1;
	for (y = 0; y < ymax; ++y)
		yj[y] = -1;
	for (j = 0; j < nlabs; ++j) {
		y = lbn[j].y;
		if (y < 0) y = 0;
		if (yj[y] < 0) yj[y] = j;
	}
	last = 0;
	for (y = 0; y < ymax; ++y) {
	   if (yj[y] < 0)
			yj[y] = last;
	   else
			last = yj[y];
	}
};   /* lsorty */
/************************************************/

var lsortx = function (stff, nstff) { //TODO !!! ASK about parameters
/*
   sort staff symbols into ascending order of 'x'
   (bubble sort)

   called by lfindstaff,
*/
   var j;
   var k;
   var s0,s1;

   for (j = 0; j < (nstff-1); ++j) {
      for (k = j; k < nstff; ++k) {
         if (stff[k][1] < stff[j][1]) {
            s0 = stff[j][0];
            s1 = stff[j][1];
            stff[j][0] = stff[k][0];
            stff[j][1] = stff[k][1];
            stff[k][0] = s0;
            stff[k][1] = s1;
         }
      }
   }
   return stff;
};   /* lsortx */
/************************************************/

var loverlap = function( p1j, p2j, p1k, p2k) {
/*
   check how much symbols j and k overlap in dimension p

   called by lbent, lleggesture, lhastap, lhasgesture,
             lseeksym, ldopivot,
			 
   calls lgetout,
*/
   var p1max,p2min;
   var lap;

   if ((p1j > p2j)||(p1k > p2k)) {
	   console.log("OOPS: loverlap "+p1j+" "+p2j+" "+p1k+" "+p2k+"\n");
	   lgetout(j);
   }
   lap = FALSE;
   if (p1k < p1j)
      p1max = p1j;
         else
            p1max = p1k; 
   if (p2k < p2j)
      p2min = p2k;
         else
            p2min = p2j;
   lap = p2min - p1max;
   return lap;
} /* loverlap */
/********************************************/

var lfindstaff = function() {
/*
    find the centres of the staves

    called by linter,
	calls     lsortx, lgetout,
*/
	var j,jp,jq;
	var k,kp,kq;
	var staffstart;
	var nstaffstart;
	var nstff;
	var stff = get2DArray(TMAX);  // [TMAX][2];
	
	k = 0;
	staffstart = 0;
	for (j = 0; j < nlabs; ++j) {
		if (lbn[j].m == Stav) {
			stff[k][0] = j;
			stff[k][1] = lbn[j].x;
			if (lbn[j].y > staffstart)
				staffstart = lbn[j].y;
			nstaffstart = j;
			++k;
			lbn[j].a = DONE;
		}
	}
	if (k < 3) {
	  console.log("lfindstaff: only "+k+" staff lines\n");
	  lgetout(1);
	  if (ok == 1)
	  return;
	}
	if (k > TMAX) {
	  console.log("lfindstaff: "+k+" staff lines, max "+TMAX+"\n");
	  lgetout(1);
	  if (ok == 1)
	  return;
	} 
	nstff = k;
	lsortx(stff,nstff);
	k = 0;
	for (j = 1; j < nstff; j += 3)
	{
	  staff[k][0] = stff[j][0];
	  staff[k][1] = stff[j-1][1];
	  staff[k][2] = stff[j][1];
	  staff[k][3] = stff[j+1][1];
	  staff[k][4] = -1;
	  staff[k][5] = TODO;
	  ++k;
	}
	nstaff = k;
	stmiddle = (staff[0][2] + staff[nstaff-1][2])/2;
	npins = 0;
	// seek pins under center stafflines
	for (j = 0; j < nstaffstart; ++j)
	{
	  if (lbn[j].m == Pins)
	  {
		 jp = lbn[j].x;
		 jq = lbn[j].x2;
		 pins[npins][0] = j;
		 pins[npins][1] = -1;
		 for (k = 0; k < nstaff; ++k)
		 {
			kp = staff[k][2] - 1;
			kq = kp+2;
			if (loverlap(jp,jq,kp,kq) > 0)
			{
			   if (lbn[j].d == 0)
			   {
				  staff[k][4] = MAN;
				  pins[npins][1] = k;
				  lbn[j].a = DONE;
			   }
			   else
			   {
				  staff[k][4] = WOMAN;
				  pins[npins][1] = k;
				  lbn[j].a = DONE;
			   } /* empty pin */
			} /* pin under central staff */
		 } /* k : staff lines */
		 ++npins;
	  } /* a pin found */
	} /* j */
	if (nstaff < 1)
	   printf("No staves found\n");
	else
	for (j = 0; j < nstaff; ++j)
	{
		if (j == 0)
			printf("\n");
	  printf("staff %d: ",j+1);
	  if (staff[j][4] == MAN)
		 printf(" man\n");
	  else
	  if (staff[j][4] == WOMAN)
		 printf(" woman\n");
	  else
		 printf(" no gender\n");
	}
	rtrn: ;
}   /* lfindstaff */
/***************************************************/

var lfindystart = function() {
/*
   find y position of first double bar line

   called by linter,
*/
	var j;
	
	ystart = -1;
/*
   seek initial double bar line -
*/
	for (j = 0; ((j < nlabs)&&(ystart < 0)); ++j) {
		if ((lbn[j].m == Bars) && (lbn[j].d == LOW))
			ystart = lbn[j].y + 1;
	}
/*
   if none, seek any bar line -
*/
	if (ystart < 0) {
		for (j = 0; ((j < nlabs)&&(ystart < 0)); ++j)
			if (lbn[j].m == Bars) ystart = lbn[j].y + 1;
	}
/*
   if none, seek any supporting direction symbol -
*/
	if (ystart < 0) {
		for (j = 0; ((j < nlabs)&&(ystart < 0)); ++j)
			if ((lbn[j].m == Dirn) && ((lbn[j].c == 1) || (lbn[j].c == -1)) )
				ystart = lbn[j].y;
		}
		ystart -= 3;
		if (ystart < 0) {
		console.log("linter : findystart finds no direction support symbols\n");
		lgetout(1);
		ystart = 0;
	}
}; /* lfindystart */
/**************************************************/

var output = "";
var lchange = function(d) {
/*
    change bend in ankles,legs, and hips while stepping

    called by ldostep,
*/
   if (d == 'L') {
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set tcls1 tlow1\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set tcls2 tlow2\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set tcls3 tlow3\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set lcls1 llow1\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set lcls2 llow2\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set lcls3 llow3\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set fcls1 flow1\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set fcls2 flow2\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set fcls3 flow3\n";
   } else if (d == 'M') {
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set tcls1 trlx1\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set tcls2 trlx2\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set tcls3 trlx3\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set lcls1 lrlx1\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set lcls2 lrlx2\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set lcls3 lrlx3\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set fcls1 frlx1\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set fcls2 frlx2\n";
      output += 
         "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set fcls3 frlx3\n";
   } else if (d == 'H') {
      if ((ji != 1)&&(ji != 10)) {
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set tcls1 thig1\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set tcls2 thig2\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set tcls3 thig3\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set lcls1 lhig1\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set lcls2 lhig2\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set lcls3 lhig3\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set fcls1 fmed1\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set fcls2 fmed2\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set fcls3 fmed3\n";
      } else if ((ji != 5)&&(ji != 6)) {
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set tcls1 trlx1\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set tcls2 trlx2\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set tcls3 trlx3\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set lcls1 lrlx1\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set lcls2 lrlx2\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set lcls3 lrlx3\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set fcls1 fhig1\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set fcls2 fhig2\n";
         output += 
            "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set fcls3 fhig3\n";
      }
   } else {
      console.log("linter: funny depth parameter, frame "+fstart+"\n");
      lgetout(1);
   }
}; /* lchange */
/*********************************************/

var lsetcoords = function() {
/*
   set the coordinate system for a step
	
	called by ldostep, lleggesture
	*/
	if ( dofig == MAN ) {
		if (mspace == false)
			output += "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set    coords mpelvis\n";
		else
			output += "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set    coords mspace\n";
	} else {
		if (wspace == FALSE)
			output += "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set    coords wpelvis\n";
		else
			output += "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set    coords wspace\n";
	}	/* dofig == WOMAN */
}; /* lsetcoords */
/************************************************/

var ldostep = function() {
/*
   create NUDES calls for steps on score

   called by laction,
   calls     lsetframes, lfindnext, lsetcoords
*/
	var b;
	var havestep;
	var k;
	var n;

	b = 0;
	if ( ( jm == Dirn ) && ( ( jc == -1 ) || ( jc == 1 ) ) ) {
		havestep = TRUE;
		k = lfindnext ( jc, jy + jh, jy + 2 * jh );
		if ( ji > 5 )
			n = ji - 5;
		else
			n = ji + 5;
		output += "*\n";
		if ( ( jc == -1 ) && ( ( ji == 1 ) || ( ji == 5 ) || ( ji == 3 ) ) ) {
			console.log("dostep: funny symbol in left support column, line "+j+", bar "+nbar+"\n");
			console.log(leadingZeros(jm, 3)+" "+leadingZeros(ji, 3)+" "+leadingZeros(jy, 3)+" "+leadingZeros(js, 3)+" "+leadingZeros(jw, 3)+" "+leadingZeros(jh, 3)+" "+jd+"\n");
		} else if ( ( jc == 1 ) && ( ( ji == 10 ) || ( ji == 6 ) || ( ji == 8 ) ) ) {
			console.log("dostep: funny  symbol in right support column, line "+j+", bar "+nbar+"\n");
			console.log(leadingZeros(jm, 3)+" "+leadingZeros(ji, 3)+" "+leadingZeros(jy, 3)+" "+leadingZeros(js, 3)+" "+leadingZeros(jw, 3)+" "+leadingZeros(jh, 3)+" "+jd+"\n");
		} else {
			output += "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set    fend  "+frange+"\n";
			output += "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" call   "+risesub[jd]+"\n";
			lsetcoords();
			if ( jc > 0 )
				output += "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" call   forright\n";
			else
				output += "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" call   forleft\n";
			if ( ( ji == 1 ) || ( ji == 10 ) )
				output += "call      "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" forward\n";
			if ( ji == 3 )
				output += "call      "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" rside\n";
			if ( ji == 8 )
				output += "call      "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" lside\n";
			if ( ( ji == 5 ) || ( ji == 6 ) )
				output += "call      "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" back\n";
			if ( ji == 2 )
				output += "call      "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" rfordiag\n";
			if ( ji == 9 )
				output += "call      "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" lfordiag\n";
			if ( ji == 4 )
				output += "call      "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" rbacdiag\n";
			if ( ji == 7 )
				output += "call      "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" lbacdiag\n";
			if ( ji == 11 )
				output += "call      "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" close\n";
			lbn[j].a = DONE;
			pstart = fstart;
			pend = fend;
		}
		rise = jd;
		prevc = jc;
		previ = ji;
	} /* column OK */
}     /* ldostep */
/****************************************************/

var lhastap = function(j) {
/*
   check if symbol j has overlapping ground contact
   
   called by ldopivot,
*/
   var t;
   var yk;
   var k;
   var kc,ki,ky,ky2;
   var km;

   t = -1;
   yk = jy - jh;
   if (yk < 1) yk = 1;
   for (k = yj[yk]; ((t < 0)&&(k < nlabs)&&(lbn[k].y < jy2)); ++k) {
      km = lbn[k].m;
      if (km == Misc) {
         ki = lbn[k].i;
         if ((ki >= 4) && (ki <= 9)) {
            kc = lbn[k].c;
            if ((kc == -3)||(kc == -2)||(kc == 2)||(kc == 3)) {
               ky = lbn[k].y;
               ky2 = ky + lbn[k].h;
               if (loverlap(jy,jy2,ky,ky2) > 0) {
                  t = k;
               } /* overlap = TRUE*/
            } /* in leg gesture column */
         } /* tap symbol */
      } /* tap menu */
   } /* k */
   return(t);
}; /* lhastap */
/*******************************************************/

var lhasgesture = function(j) {
/*
   check if symbol j has overlapping gesture
   
   called by ldopivot,
   calls     loverlap,
*/
   var kc,ky,ky2;
   var k;
   var g;
   var km;

   g = -1;
   for (k = 0; ((g < 0)&&(k < nlabs)); ++k) {
      km = lbn[k].m;
      if (km == Dirn) {
         kc = lbn[k].c;
         if ((kc == -3)||(kc == 3)) {
            ky = lbn[k].y;
            ky2 = lbn[k].y2;
            if (loverlap(jy,jy2,ky,ky2) > 0)
               g = k;
         }
      }
   }
   return(g);
}; /* lhasgesture */
/*******************************************************/

var lleggesture = function() {
/*
   do gestures of the legs

   called by laction,
   calls     lsetframes, lgetout, lsetcoords,

   Volm   1  RELAX
   Volm   3  BENT
   Volm   2  STRAIGHT
   Volm   4  STRETCH
   Volm   7  hold
*/
   if ((jc == -3)||(jc == 3)) {
         if ((jd < 0) || (jd > 2)) {
            console.log("OOPS: dogesture height problem line "+j+"\n");
            console.log(leadingZeros(jm, 3)+" "+leadingZeros(ji, 3)+" "+leadingZeros(jx, 3)+" "+leadingZeros(jy, 3)+" "+leadingZeros(js, 3)+" "+leadingZeros(jw, 3)+" "+leadingZeros(jh, 3)+" "+jd+"\n");
            lgetout(j);
            if (ok == 1) return;
         } /* level funny */
         output += "*\n";
         if  (ji==11)
            output += "* close without weight";
         else
         if ((ji==1)||(ji==10))
            output += "* forward gesture";
         else
         if ((ji==2)||(ji==9))
            output += "* forward diagonal gesture";
         else
         if ((ji==3)||(ji==8))
            output += "* sideways gesture";
         else
         if ((ji==4)||(ji==7))
            output += "* back diagonal gesture";
         else
         if ((ji==5)||(ji==6))
            output += "* backward gesture";
//
         if (jd == LOW)
            output += " low\n";
         if (jd == MED)
            output += " middle\n";
         if (jd == HIGH)
            output += " high\n";
//
         if (jc < 0) {
            if ((ji <= 1)||(ji == 3)||(ji == 5)||(ji > 11)) {
               console.log("OOPS: dogesture direction problem line "+j+"\n");
               console.log(leadingZeros(jm, 3)+" "+leadingZeros(ji, 3)+" "+leadingZeros(jx, 3)+" "+leadingZeros(jy, 3)+" "+leadingZeros(js, 3)+" "+leadingZeros(jw, 3)+" "+leadingZeros(jh, 3)+" "+jd+"\n");
               lgetout(1);
               if (ok == 1) return;
            } /* i wrong */
            output += "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" call   forleft * left = b\n";
         } /* left side */
         else if (jc > 0) {
            if ((ji < 1)||(ji == 6)||(j == 8)||(ji == 10)||(ji > 11)) {
               console.log("OOPS: dogesture direction problem line "+j+"\n");
               console.log(leadingZeros(jm, 3)+" "+leadingZeros(ji, 3)+" "+leadingZeros(jx, 3)+" "+leadingZeros(jy, 3)+" "+leadingZeros(js, 3)+" "+leadingZeros(jw, 3)+" "+leadingZeros(jh, 3)+" "+jd+"\n");
               lgetout(1);
               if (ok == 1) return;
            } /* i wrong */
            output += "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" call   forright * right = b\n";
         } /* right side */
//
         output += "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" centre afoot  "+xyz[dofig]+"\n";
//
         if (ji == 11) {
            output += "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" call   "+risesub[rise]+"\n";
            output += "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" set    fend  "+frange+"\n",
				lsetcoords();
            output += "call      "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" close\n";
         } /* close without weight */
         else {
            output += "quadratic "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" bendto bthigh bhip   pelvis "+stt[jd][ji][0]+" "+stt[jd][ji][1]+" "+stt[jd][ji][2]+"\n";
            if ((jd == LOW)&&((ji == 1)||(ji == 3)||(ji == 8)||(ji == 10))||(jb == 2)||(jb == 4))
               output += "linear    "+leadingZeros(fhalf, 3)+" "+leadingZeros(fend, 3)+" bendto bleg   bknee  bthigh lhig1 lhig2 lhig3\n";
            else
               output += "linear    "+leadingZeros(fhalf, 3)+" "+leadingZeros(fend, 3)+" bendto bleg   bknee  bthigh lrlx1 lrlx2 lrlx3\n";
            output += "linear    "+leadingZeros(fhalf, 3)+" "+leadingZeros(fend, 3)+" bendto bfoot  bankle bleg   fhig1 fhig2 fhig3\n";
         } /* doing a leg gesture */
         if ((ji != 11)&&(hold == NO)||(st < 1))
            output += "repeat    "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" moveto fig    afoot  "+xyz[dofig]+"\n";
         lbn[j].a = DONE;
      } /* no tap and pivot */ /* c OK */
} /* lleggesture */
/***************************************************/

//PORT ON 2013-12-09 (Errorage)
var ldoarms = function() {
/*
   do movements of the arms

   called by ldolimb,
   calls     lsetframes,
*/
	var absjc;
	
	if (jc < 0) absjc = -jc; else absjc = jc;
	if ((absjc > 3)&&(absjc < 7))
	{
		if (jm == Dirn)
		{
			output += "*\n* arms\n";
			if ((jd < 0) || (jd > 2))
			{
				console.log("ldoarms problem line "+j+" bar "+nbar+"\n");
				console.log(""+ jm+" "+leadingZeros(ji, 3)+" "+leadingZeros(jx, 3)+" "+leadingZeros(jy, 3)+" "+leadingZeros(js, 3)+" "+leadingZeros(jw, 3)+" "+leadingZeros(jh, 3)+" "+jd+"\n");
				lgetout(1);
				if (ok == 1) return;
			}
			if (jc < 0) // left arm
			{
				if (jb == RELAX)
				{
					output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto luarm lshldr shldrs "+leadingZeros(arlx[jd][ji][0], 3)+" "+leadingZeros(arlx[jd][ji][1], 3)+" "+leadingZeros(arlx[jd][ji][2], 3)+"\n";
					output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto llarm lelbow luarm "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+" "+leadingZeros(45, 3)+"\n";
					if (dofig == MAN)
					{  
						output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto lhand lwrist llarm "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+"\n";
					} /* man */
					else
					{
						output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto lhand lwrist llarm "+leadingZeros(270, 3)+" "+leadingZeros(0, 3)+" "+leadingZeros(150, 3)+"\n";
					} /* woman */
				}
				else
				if (jb == BENT)
				{
					if (ji == 11)
					{
						output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto luarm lshldr shldrs "+leadingZeros(abnt[jd][0][0], 3)+" "+leadingZeros(abnt[jd][0][1], 3)+" "+leadingZeros(abnt[jd][0][2], 3)+"\n";
						output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto llarm lelbow luarm "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+" "+leadingZeros(70, 3)+"\n";
					}
					else
					{
						output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto luarm lshldr shldrs "+leadingZeros(abnt[jd][ji][0], 3)+" "+leadingZeros(abnt[jd][ji][1], 3)+" "+leadingZeros(abnt[jd][ji][2], 3)+"\n";
						output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto llarm lelbow luarm "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+" "+leadingZeros(90, 3)+"\n";
					} /* ji != 11 */
				}
				else
				{
					output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto llarm lelbow luarm "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+"\n";
					output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto luarm lshldr shldrs "+leadingZeros(stt[jd][ji][0], 3)+" "+leadingZeros(stt[jd][ji][1], 3)+" "+leadingZeros(stt[jd][ji][2], 3)+"\n";
				}
			}
			else // if (jc > 0) =  right arm
			{
				if (jb == RELAX)
				{
					output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto ruarm rshldr shldrs "+leadingZeros(arlx[jd][ji][0], 3)+" "+leadingZeros(arlx[jd][ji][1], 3)+" "+leadingZeros(arlx[jd][ji][2], 3)+"\n";
					output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto rlarm relbow ruarm "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+" "+leadingZeros(45, 3)+"\n";
					if (dofig == MAN)
					{  
						output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto rhand rwrist rlarm "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+"\n";
					} /* man */	//TODO
					else
					{
						output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto rhand rwrist rlarm "+leadingZeros(270, 3)+" "+leadingZeros(0, 3)+" "+leadingZeros(150, 3)+"\n";
					} /* woman */

				} /* relaxed */
				else if (jb == BENT)
				{		   
					output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto ruarm rshldr shldrs "+leadingZeros(abnt[jd][ji][0], 3)+" "+leadingZeros(abnt[jd][ji][1], 3)+" "+leadingZeros(abnt[jd][ji][2], 3)+"\n";
					if (ji == 11)
						output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto rlarm relbow ruarm "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+" "+leadingZeros(70, 3)+"\n";
					else
						output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto rlarm relbow ruarm "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+" "+leadingZeros(90, 3)+"\n";
				} /* bent */
				else
				{
					output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto rlarm relbow ruarm "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+"\n";
					output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto ruarm rshldr shldrs "+leadingZeros(stt[jd][ji][0], 3)+" "+leadingZeros(stt[jd][ji][1], 3)+" "+leadingZeros(stt[jd][ji][2], 3)+"\n";
					output += "quadratic "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto rhand rwrist rlarm "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+" "+leadingZeros(0, 3)+"\n";
				} /* not bent or relaxed */
			} /* right arm */
		} /* Dirn symbol */
	} /* in arm column */
	lbn[j].a = DONE;
} /* ldoarms */

var lspotturn = function(j, piv, fstart, fend, g)
/*
  maintain straight non-standing foot with ground 
  contact during turn.

  called by ldopivot,
*/
{
   var gc,gi;

   gi = lbn[g].i;
   gc = lbn[g].c;
   fprintf(nudesfile,"*\n* spot turn-\n");
   if (gc < 0)
      output += "repeat "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" call forleft * left = b\n";
   else
      output += "repeat "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" call forright * right = b\n";
   output += "repeat "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" centre afoot "+xyz[dofig]+"\n";
   output += "linear "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" spinby fig afoot pelvis "+piv+" y\n";
   output += "linear "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" bendto bthigh bhip pelvis "+stt[0][ji][0]+" "+stt[0][ji][1]+" "+stt[0][ji][2]+"\n";
   output += "repeat "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" ground fig\n";
   if ((hold == NO)||(st < 1))
      output += "repeat "+leadingZeros( fstart, 3)+" "+leadingZeros(fend, 3)+" moveto fig afoot "+xyz[dofig]+"\n";
   output += "linear "+leadingZeros( fstart, 3)+" "+leadingZeros(fhalf, 3)+" bendto bfoot bankle bleg fhig1 fhig2 fhig3\n";
   output += "repeat "+leadingZeros( fhalf, 3)+" "+leadingZeros(fend, 3)+" drag bfoot bfoot bankle bleg x\n";
   lbn[j].a = DONE;
   lbn[t].a = DONE;
} /* lspotturn */
/******************************************************/


var lgetpin = function()
/*
    seek a pin in a rotation sign
	
	called by ldolimb, ldopivot,
	calls     loverlap,
*/
{
    var k;
	var ki;
    var piv;
    var ymost;
    var xlap,ylap;

    ki = -123;
    ymost = -1;
    for (k = yj[jy-jh]; lbn[k].y < jy2; ++k)
    {
         if (lbn[k].m == Pins)
         {
             xlap = loverlap(jx,jx2,lbn[k].x,lbn[k].x2);
             ylap = loverlap(jy,jy2,lbn[k].y,lbn[k].y2);
             if ((xlap > 0) && (ylap > ymost))
             {
                ki = lbn[k].i;
                ymost = ylap;
             } /* pin overlaps more than previous pins */
         } /* got a pin */
    } /* k loop looking for overlapping pin */
	 piv = 0;
    if ((ki > 0)&&(ki <= 9))
    {
       if (ji == 1) piv = -45*(9-ki);
       if (ji == 2) piv = 45*(ki-1);
       if (ki == 1) piv = 360;
    }
    return(piv);
} /* lgetpin */
/***************************************************/

var ldopivo = function()
/*
   do turns in the support columns

   called by laction,
   calls     lsetframes, lspotturn, lhasgesture, 
             lhastap,    lgetpin,
*/
{
   var g;
   var t;
   var piv;

   if ( (jm == Rotn)&&(nbar > 0)&&
	 ((jc == -2)||(jc == -1)||(jc == 1)||(jc == 2)) )
   {
      piv = lgetpin();
      if (fstart < 1) fstart = 1;
      g = lhasgesture(j);
      t = lhastap(j);
      if ((g > 0)&&(t > 0))
      {
         lspotturn(j,piv,fstart,fend,g);
         pstart = fstart;
         pend = fend;
      }
      else
      {
         output += "*\n* pivot\n";
         if (jc < 0)
            output += "repeat "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" call forleft * b = left\n";
         else
            output += "repeat "+leadingZeros(fstart, 3)+" "+leadingZeros(fend, 3)+" call forright * b = right\n";
         output += "repeat "+leadingZeros( fstart, 3)+" "+leadingZeros(fstart+1, 3)+" centre bfoot "+xyz[dofig]+"\n";
		 //todo
         fprintf(nudesfile,
            "linear    %3d %3d spinby fig    bfoot  pelvis %d y\n",
               fstart,fend,piv);
         if ((hold == NO)||(st < 1))
			 fprintf(nudesfile,
                "repeat    %3d %3d moveto fig    bfoot  %s\n",
                  fstart,fend,xyz[dofig]);
         if (hold == PR)
             hold = NO;
         pstart = fend;
         pend = fend+1;
      } /* spotturn == false */
   }
} /* ldopivot */
/**************************************************/


var lseeksym = function(m, i, x1, x2, y3, y4)
/*
     seek a symbol of menu m, item i,
	 overlapping box x1,x2,y1,y2.

	 called by lbows, lsethold,
	 call      loverlap,
*/
{
   var lap;
   var kstart;
   var k,kx,kx2,ky,ky2;
   var y1,y2;

   lap = -1;
   if (y3 < 0) y1 = 0; else y1 = y3;
   if (y4 < 0) y2 = 0; else y2 = y4;
   kstart = y1 - 2*STEP;
   if (kstart < 1) kstart = 1;
   for (k = yj[kstart]; ((lap < 0)&&(lbn[k].y < y2)); ++k)
   {
      kx = lbn[k].x;
      kx2 = lbn[k].x2;
      ky = lbn[k].y;
      ky2 = lbn[k].y2;
      if ((lbn[k].m == m )&&(lbn[k].i == i))
      {
         if ((loverlap(x1,x2,kx,kx2) > 0)
            && (loverlap(y1,y2,ky,ky2) > 0))
         {
            lap = k;
         }
      } /* m and i true */
   } /* k loop */
   return(lap);
}  /* lseeksym */
/***************************************************/


var lbows = function()
/*
      detect and flag the various contact bows.

	  called by linter,
	  calls     lseeksym,

Relevant symbols:-
     m     i
    Misc   1  bow
    Limb   4  lhand
    Limb   9  rhand
    Area   1  top/front 
    Area   5  back/bottom
    Volm   1  relax
    Volm   2  bent
    Volm   3  straight
    Volm   4  stretch
    Volm   7  hold
*/
{
   var centre;
   var held,front,back;
   var mlhand,mrhand,wlhand,wrhand;

   centre = (staff[0][2] + staff[1][2])/2;
   for (j = 0; j < nlabs; ++j)
   {
      if ((lbn[j].m == Misc)&&(lbn[j].i == 1))
      {
          lassign();
          held = lseeksym(Volm,1,jx,jx2,jy-STEP,jy2);
	  if (held > 0)
          {
              mlhand = lseeksym(Limb,4,jx-STEP/2,jx+STEP/2,jy,jy2+STEP);
              mrhand = lseeksym(Limb,9,jx-STEP/2,jx+STEP/2,jy,jy2+STEP);
              wlhand = lseeksym(Limb,4,jx2-STEP/2,jx2+STEP/2,jy,jy2+STEP);
              wrhand = lseeksym(Limb,9,jx2-STEP/2,jx2+STEP/2,jy,jy2+STEP);
              front  = lseeksym(Area,1,jx-STEP/2,jx+STEP/2,jy,jy2+STEP);
              if (front < 0)
                 front  = lseeksym(Area,1,jx2-STEP/2,jx2+STEP/2,jy,jy2+STEP);
              if (front < 0)
                 front  = lseeksym(Area,2,jx-STEP/2,jx+STEP/2,jy,jy2+STEP);
              if (front < 0)
                 front  = lseeksym(Area,2,jx2-STEP/2,jx2+STEP/2,jy,jy2+STEP);
              back   = lseeksym(Area,5,jx-STEP/2,jx+STEP/2,jy,jy2+STEP);
              if (back < 0)
                 back  = lseeksym(Area,5,jx2-STEP/2,jx2+STEP/2,jy,jy2+STEP);
              jb = 0;
              if (front > 0)
              {
                 jb = FRONT;
                 lbn[front].a = DONE;
              }
              else 
              if (back > 0)
              {
                 jb = BACK;
                 lbn[back].a = DONE;
              }
              if (mlhand > 0) jb += MLHAND;
              if (mrhand > 0) jb += MRHAND;
              if (wlhand > 0) jb += WLHAND;
              if (wrhand > 0) jb += WRHAND;
              if (jb <= 0)
                 fprintf(nudesfile,"* OOPS: lbows: bow %d with no contacts\n",j);
              else
              {
                 if (mlhand > 0) lbn[mlhand].b = jb;
                 if (mrhand > 0) lbn[mrhand].b = jb;
                 if (wlhand > 0) lbn[wlhand].b = jb;
                 if (wrhand > 0) lbn[wrhand].b = jb;
              }
           } /* held */
           output += "* lbowsb "+ j+1+" "+held+" "+front+" "+back+" "+mlhand+" "+mrhand+" "+wlhand+" "+wrhand+" "+jb+"\n";
      } /* contact bow */
   } /* j */
} /* lbows */
/****************************************/


var lstart = function()
/*
   seek pins denoting starting positions.

   called by linter,
*/
{
   var k;
   var p;
   var dx,dz;
   var mx,my;
   var wx,wy;

   mx = -123;
   my = -123;
   wx = -123;
   wy = -123;
   if (nm > 0)
      output += "\n*\ncall 0 1 doman\n";
   else
      output += "\n*\ncall 0 1 dowoman\n";
   output += "call 0 1 forleft\n";
   for (k = 0; k < npins; ++k)
   {
      if (pins[k][1] < 0)
      {
         p = pins[k][0];
         ji = lbn[p].i;
         if (lbn[p].d == LOW)
         {
            mx = lbn[p].x;
            my = lbn[p].y;
            if (nm > 0)
               output += "*\nquadratic 0 1 spinby man mlfoot mpelvis "+ (ji-1)*45+" y\n";
         }
         else
         {
            wx = lbn[p].x;
            wy = lbn[p].y;
            if (nw > 0)
               output += "*\nquadratic 0 1 spinby woman wrfoot wpelvis "+ (ji-1)*45+" y\n";
         }
         if ((wx > 0)&&(mx > 0)&&(wy > 0)&&(my > 0))
         {
            dx = ((wx - mx)*2)/3;
            dz = (wy - my)/2;
            if (nmw > 0)
            {
               output += "*\n";
               output += "repeat 0 1 centre mpelvis kx ky kz\n";
               output += "repeat      0   1 moveto woman    wpelvis kx ky kz\n";
               output += "repeat      0   1 axes   wpelvis  cx cy cz\n";
               output += "linear      0   1 set    dx "+dx+"\n";
               output += "linear      0   1 set    dz "+dz+"\n";
               output += "linear      0   1 mult   wx   dx  cx\n";
               output += "linear      0   1 mult   wz   dz  cx\n";
               output += "repeat      0   1 centre wpelvis  cx cy cz\n";
               output += "repeat      0   1 centre wpelvis  cx cy cz\n";
               output += "repeat      0   1 centre wpelvis  cx cy cz\n";
			   output += "repeat      0   1 set    fpos     1\n";
               output += "repeat      0   1 call   noposn\n";
            }
         }
      }
   }
} /* lstart */
/***********************************************/
   

var lsetrange = function()
/*
   set range of symbols to be interpreted

   called by linter,
*/
{
   var bend;
   var k,kmax;
   var ymax;

   ystart = 0;
   yend = lbn[0].y;
   ymax = yend;
   sstart = 0;
   ssend = nlabs;
   for (k = 0; k < nlabs; ++k)
   {
      if (lbn[k].m == Bars)
      {
         if (lbn[k].i == bstart)
         {
            sstart = k;
            ystart = lbn[k].y;
         }
      }
   }
   bend = bstart + blength;
   for (k = (sstart+1); k < nlabs; ++k)
   {
      if (lbn[k].m == Bars)
      {
         if (lbn[k].i == bend)
            ssend = k;
      }
   }
   for (k = 0; k < nlabs; ++k)
   {
      if (lbn[k].m == Dirn)
      {
         if (lbn[k].y > yend)
            yend = lbn[k].y;
         if ((lbn[k].y+lbn[k].h) > ymax)
		 {
            ymax = lbn[k].y+lbn[k].h;
			kmax = k;
		 }
      }
   }
   fmax = 2 + int(lbn_fpp*double(ymax));
	console.log("\n   lsetrange: pixels "+ymax+", frames "+fmax+"\n");
} /* lsetrange */
/****************************************************/

var lcopyfigs = function()
/*
   finish off

   called by linter,
   calls     lgetout,
*/
{
	
   var figsname = "lintel.n";
   
	var done = false;
	getFileFromServer(figsname, function(text) {
		if (text === null) {
			printf("\n\noops %s not in folder\n",figsname);
			  lgetout(1);
			  if (ok == 1) return;
		} else {
			//TODO REVIEW THIS CODE
			output += text;
			done = true;
		}
	});
	while (done == false) {}
} /* lcopyfigs */
/********************************************/

var lfinish = function()
/*
   finish off

   called by linter,
   calls lgetout,
*/
{
   fmax += 2;
   output += "*\n";
   output += "**************************\n";
   output += "*\n";
   if (nm > 0)
      output += "repeat      0 "+leadingZeros(fmax, 3)+" ground man\n";
   else
      output += "repeat      0   1 moveto man    mlfoot  10000 10000 10000\n";
   if (nw > 0)
      output += "repeat      0 "+leadingZeros(fmax, 3)+" ground woman\n";
   else
      output += "repeat      0   1 moveto woman  wlfoot  10000 10000 10000\n";
   if (nm > 0)
      output += "repeat      0 "+leadingZeros(fmax, 3)+" centre mpelvis fx fy fz\n";
   else
      output += "repeat      0 "+leadingZeros(fmax, 3)+" centre wpelvis fx fy fz\n";
   if (track == TRUE)
   {
	   output += "repeat      0 "+leadingZeros(fmax, 3)+" add     fy -900 fz\n";
       output += "repeat      0 "+leadingZeros(fmax, 3)+" place   fx  500 fy\n";
   }
   output += "repeat      0 "+leadingZeros(fmax, 3)+" observe -9    0  0\n*\n";
   output += "end dance\n****************************\n";
   output += "*\nsubroutine setfmax\n";
   output += "*\nrepeat 0 1 set fmax "+fmax+"\n";
   output += "*\nend setfmax\n";
   output += "****************************\n*\nstop\n";
   if (nbar > 0)
         frperbar = fmax/nbar;
   else
         frperbar = 0;
} /* lfinish */
/********************************************/

var lselectfig = function()
/*
   select figure

   called by linter,
*/
{   
   var k;
   var nf;
   var nogo;
   var st;
   var stv0,stv1,st4;
   var stv = new Array();	//[2]
   var key;

	//TODO REVIEW (Was a lot of GOTOs, now a WHILE loop)
	var again = true;
	while(again == true){
		again = false;
	   for (k = 0; k < nstaff; ++k)
		  staff[k][5] = DONE;
	   nf = 0;
	   nm = 0;
	   nw = 0;
	   nogo = FALSE;
	   if (nstaff < 1)
		  console.log("no staves\n");
	   else
	   if (nstaff == 1)
	   {
		  staff[0][5] = TODO;
		  if (staff[0][4] == MAN) 
			 ++nm;
		  else
			 ++nw;
	   }
	   else
	   if (nstaff > 1)
	   {
		  nmw = 0;
		  if (nstaff > TMAX)
			 console.log("This can only interpret staves from 1 to "+ 		    TMAX+"\n");
		  if (lbn_figures == 2)
		  {
			 stv[0] = 1; stv[1] = 2;
			 track = TRUE;
		  }
		  else // (lbn_figures != 2)
		  {
			 console.log("\nPlease type the number of staves to be interpreted\n");
			 if (gets(buf) == NULL)
			 {
				console.log("OOPS: cannot open standard input\n");
				lgetout(1);
				nogo = TRUE;
				if (ok == 1){
					break;
				}
			 }
			 sscanf(buf,"%d",lbn_figures);
			 if (lbn_figures > 2)
			 {
				console.log("sorry; this program can only interpret 2 staves at a time\n");
				nogo = TRUE;
				again = true;
				continue;
			 }
			 if (lbn_figures == 1)
				console.log("Please enter the staff number to be interpreted\n");
			 else
			 {
				console.log("Please enter staff numbers to be interpreted\n");
				console.log("separated by a space, and followed by the 'enter' key.\n\n");
			 }
			 if (gets(buf) == NULL)
			 {
				console.log("OOPS: cannot read staff numbers\n");
				lgetout(1);
				nogo = TRUE;
				if (ok == 1){
					break;
				}
			 }
			 if (lbn_figures == 1)
			 {
				sscanf(buf,"%d",stv0); 
				stv[0] = stv0; stv[1] = -1;
			 }
			 else
			 {
				sscanf(buf,"%d %d",stv0,stv1); 
				stv[0] = stv0; stv[1] = stv1;
			 }
		  } /* lbn_figures != 2 */
		  for (nf = 0; nf < lbn_figures; ++nf)
		  {
				st = stv[nf]-1;
				if ((st < 0)||(st > nstaff))
				{
					console.log("OOPS: staff number "+(st+1)+" out of range\n");
					again = true;
					continue;
				}
				st4 = staff[st][4];
				if ( ((nm > 0)&&(st4 == MAN))
				   ||((nw > 0)&&(st4 == WOMAN)) )
				{
					console.log("Sorry: can only do one man and/or one woman.");
					console.log("Please select again.\n");
					nogo = TRUE;
				 } /* more than 1 man or woman */
				 else
				 {
					if (st4 == WOMAN) ++nw;
					if (st4 == MAN) ++nm;
					staff[st][5] = TODO;
				 } /* a man or woman */  
				 nmw = nm*nw;
		   } /* nf */
	   } /* nstaff > 1 */
	   if (nogo == TRUE){
			again = true;
			continue;
		}
	}

   if (lbn_figures != 2)
   {
       track = TRUE;
       console.log("Track main figure? Hit 'enter' for Yes, any other key for No\n");
       key = getchar(); 
       if (key != '\n')
          track = FALSE;
   }
   else
       track = TRUE;
   if (track == FALSE)
       console.log("\n   tracking OFF\n");
   else
       console.log("\n   tracking ON\n");
} /* lselectfig */
/***********************************************/
//END PORT ON 2013-12-09 (Errorage)

//PORT ON 2013-12-10 (Errorage)
var ldobar = function()
/*
   write bar number out

   called by laction,
*/
{
   if ((jm == Bars) && (jy < yend))
   {
      ++nbar;
      output += "*\n";
	  output += "***************************\n";
      output += "*\n";
      output += "*   bar "+nbar+"\n";
      output += "*   bar "+nbar+"\n";
   }
} /* ldobar */
/********************************************/

var lbent = function()
/*
   for Volm symbol : flag next 'Dirn' symbol above
   
   called by laction,
   calls lassign,

   Volm   1  RELAX
   Volm   3  BENT
   Volm   2  STRAIGHT
   Volm   4  STRETCH
   Volm   7  hold
*/
{
   var g;
   var k;
   var ki,kx,kx2,ky,ky2;
   var jy2h;
   var km;

   for (j = 0; j < ssend; ++j)
   {
      if ((lbn[j].m == Volm)&&(lbn[j].i <= STRETCH)) 
      {
         lassign();
         jy2h = jy2+jh;
         g = -1;
         for (k = j+1;	((k < nlabs)&&(g < 0)); ++k)
         {
            km = lbn[k].m;
            if ((km == Dirn)&&(lbn[k].a == TODO))
            {
               ky = lbn[k].y;
               if (ky > jy2h)
                  g = 0;
               else
               {
                  ky2 = lbn[k].y2;
                  kx = lbn[k].x;
                  kx2 = lbn[k].x2;
                  if ((loverlap(jx,jx2,kx,kx2) > 0)
                     &&(loverlap(jy2,jy2h,ky,ky2) > 0))
                  {
                      g = k;
                      lbn[j].b = ji;
                      ki = lbn[k].i;
                      lbn[j].m = km;
                      lbn[j].i = ki;
                      lbn[j].y2 = ky2;
                      lbn[j].h = lbn[k].y2 - jy;
                      lbn[j].d = lbn[k].d;
                      lbn[k].a = DONE;
                      if (ji == BENT)
                      {
                         if ((ki == 11)&&(jc < 0))
                            lbn[j].i = 8;
                         else
                         if ((ki == 11)&&(jc > 0))
                            lbn[j].i = 3;
                      } /* ji == BENT */
                  } /* overlapping */
               } /* ky < jy2h */
            } /* km = Dirn */
         } /* k */
      } /* jm = Volm */
   } /* j */
} /* lbent */
/********************************************/

var lrelease = function()
/*
   release the hold when jm = Misc

   called by laction,
   	  Assumes one of the following holds:

	  So far:
   NO  - no hold: arm gestures apply.
   CL  - closed hold: normal ballroom dancing position.
   SS  - semi-shadow hold: both facing same way, bodies touching, 
         man's L hand to lady's L hand,
         man's R hand to front of lady's R hip,
		 lady's R hand free.
   OE  - open extended hold: both facing same way, bodies apart,
         man's R hand to lady's L hand, other hands free.
   CO  - counter open extended hold: both facing same way, bodies apart,
         man's L hand to lady's R hand, other hands free.
   SH  - shadow hold: both facing same way, bodies touching,
         L hand to L hand, R hand to R hand.

      later to do:
   PR  - promenade position: facing partner, bodies touching,
         but both prepared to travel to man's L.
   CP  - counter promenade position: facing partner, bodies touching,
         but both prepared to travel to man's R.
   DB  - double hold: facing partner, bodies apart,
         L hand to R hand, R hand to L hand.
   OP  - open hold: facing partner, bodies apart,
         man's L hand to lady's R hand, other hands free.
   CR  - crossed open hold: facing partner, bodies apart,
         man's R hand to lady's R hand, other hands free.

	Relevant symbols:-
     m     i
    Misc   1  bow
    Misc   2  release1
    Misc   3  release2
    Limb   4  lhand
    Limb   9  rhand
    Area   1  top/front 
    Area   5  back/bottom
    Volm   1  RELAX
    Volm   3  BENT
    Volm   2  STRAIGHT
    Volm   4  STRETCH
    Volm   7  hold

    FRONT   100         // front symbol found
    BACK    200         // back symbol found
    MLHAND    1         // man's left hand symbol found
    MRHAND    2         // man's right hand symbol found
    WLHAND   10         // woman's left hand symbol found
    WRHAND   20         // woman's right hand symbol found
*/
{
   var fdif;
   var fbegin,ffin;

   if ((nmw > 0)&&(ji == 2)) // release
   {
      holdcl = 0;
      holdoe = 0;
      holdco = 0;
      holdpr = 0;
      holdsh = 0;
      holdss = 0;
      fbegin = pend;
      ffin = fend;
      if (ffin <= fbegin) ffin = fbegin + 1;
      fdif = ffin - fbegin;
      if ((st > 0) && (hold != NO))
      {
         output += "repeat    "+leadingZeros(fbegin, 3)+" "+leadingZeros(ffin, 3)+" set fpos "+fdif+"\n";
         output += "call      "+leadingZeros(fbegin, 3)+" "+leadingZeros(ffin, 3)+" noposn\n*\n";
      }
      hold = NO;
      output += "* lreleasea "+             fstart+" "+fend+" "+j+" "+jb+" "+hold+" "+prevhold+"\n";
      keptf = ffin;
   }
} /* lrelease */
/******************************************/

var ldoposn = function()
/*
   set up a couple dance position

   called by lsethold, ldohold
*/
{
	  fbegin = fstart;
	  ffin = fend;
		output += "** ldoposn "+leadingZeros(fbegin,3)+" "+leadingZeros(ffin,3)+", "+leadingZeros(st,3)+" "+leadingZeros(hold,3)+"\n";
		
      if (st > 0)
      {
			flen = ffin - fbegin;
			if (flen < 1) flen = 1;
			if (hold != NO) 
            output += "repeat    "+leadingZeros(                   fbegin, 3)+" "+leadingZeros(ffin, 3)+" set    fpos "+leadingZeros(flen, 3)+"\n";
			if (hold == PR)
            output += "call      "+leadingZeros(fbegin, 3)+" "+leadingZeros(ffin, 3)+" prposn\n*\n";
			else
			if (hold == CO)
            output += "call      "+leadingZeros(fbegin, 3)+" "+leadingZeros(ffin, 3)+" coposn\n*\n";
			else
			if (hold == CL)
            output += "call      "+leadingZeros(fbegin, 3)+" "+leadingZeros(ffin, 3)+" clposn\n*\n";
			else
			if (hold == SS)
            output += "call      "+leadingZeros(fbegin, 3)+" "+leadingZeros(ffin, 3)+" ssposn\n*\n";
			else
			if (hold == OE)
            output += "call      "+leadingZeros(fbegin, 3)+" "+leadingZeros(ffin, 3)+" oeposn\n*\n";
			else
			if (hold == SH)
            output += "call      "+leadingZeros(fbegin, 3)+" "+leadingZeros(ffin, 3)+" shposn\n*\n";
			keptf = ffin;
			prevhold = hold;
      } /* st > 0 */
} /* ldoposn */
/*******************************************/

var ldokeep = function()
/*
   maintain a couple dancing position

   called by dohold,
*/
{
		output += "** ldokeep "+leadingZeros(fbegin, 3)+" "+leadingZeros(ffin, 3)+" "+leadingZeros(hold, 3)+"\n";
         if (hold == PR)
            output += "repeat    "+leadingZeros(fbegin, 3)+" "+leadingZeros(ffin, 3)+" call   prkeep\n*\n";
         else
         if (hold == CL)
		     output += "repeat    "+leadingZeros(fbegin, 3)+" "+leadingZeros(ffin, 3)+" call   clkeep\n*\n";
         else
         if (hold == OE)
            output += "repeat    "+leadingZeros(fbegin, 3)+" "+leadingZeros(ffin, 3)+" call   oekeep\n*\n";
         else
         if (hold == SS)
            output += "repeat    "+leadingZeros(fbegin, 3)+" "+leadingZeros(ffin, 3)+" call   sskeep\n*\n";
		 else
         if (hold == CO)
            output += "repeat    "+leadingZeros(fbegin, 3)+" "+leadingZeros(ffin, 3)+" call   cokeep\n*\n";
         else
         if (hold == SH)
            output += "repeat    "+leadingZeros(fbegin, 3)+" "+leadingZeros(ffin, 3)+" call   shkeep\n*\n";
         keptf = ffin;
}   /* ldokeep */
/******************************************/

var ldohold = function()
/*
    set up and maintain holds
	
	called by laction,
	calls ldokeep, ldoposn,
*/
{
   fbegin = keptf;
   ffin = pend;
		output += "** ldohold "+leadingZeros(fbegin,3)+" "+leadingZeros(ffin ,3)+","+leadingZeros(hold,3)+" "+leadingZeros(prevhold,3)+"\n";
   if (prevhold == hold) 
   {
	   fbegin = keptf;
      if (fbegin < ffin) ldokeep();
   } /* prevhold == hold */
   else
   {
      ldoposn();
   } /* prevhold != hold */
} /* ldohold */
/*************************************************/

var lsethold = function()
/*
   set the hold if jm = Limb or jm = Face

   called by laction,
   calls ldoposn,

      Uses the hand signs to determine the holds if any.
   	  Assumes one of the following holds:

	  So far:
   NO - no hold: arm gestures apply.
   CL - closed hold: normal ballroom dancing position.
   SS - semi-shadow hold: both facing same way, bodies touching, 
        man's L hand to lady's L hand,
        man's R hand to front of lady's R hip,
         ady's R hand free.
   OE - open extended hold: both facing same way, bodies apart,
        man's R hand to lady's L hand, other hands free.
   CO - counter open extended hold: both facing same way, bodies apart,
        man's L hand to lady's R hand, other hands free.
   SH - shadow hold: both facing same way, bodies touching,
        L hand to L hand, R hand to R hand.
   PR - promenade position: diagonally facing partner,
        bodies touching, both travelling to man's L.
   CP - counter promenade position: facing partner, bodies touching,
        but both prepared to travel to man's R.
   DB - double hold: facing partner, bodies apart,
        L hand to R hand, R hand to L hand.
   OP - open hold: facing partner, bodies apart,
        man's L hand to lady's R hand, other hands free.
   CR - crossed open hold: facing partner, bodies apart,
        man's R hand to lady's R hand, other hands free.

#define NO        0        // no hold
#define CL        1        // closed hold
#define PR        2        // promenade position
#define CP        3        // counter promenade position
#define DB        4        // double hold
#define OP        5        // open hold
#define CR        6        // crossed open hold
#define OE        7        // open extended hold
#define CO        8        // counter open extended hold
#define SH        9        // shadow hold
#define SS       10        // semi-shadow hold 

	Relevant symbols:-
     m     i
    Misc   1  bow
    Misc   2  release1
    Misc   3  release2
    Limb   4  lhand
    Limb   9  rhand
    Area   1  top/front 
    Area   5  back/bottom
    Volm   1  RELAX
    Volm   3  BENT
    Volm   2  STRAIGHT
    Volm   4  STRETCH
    Volm   7  hold

#define FRONT   100         // front/top symbol found
#define BACK    200         // back symbol found
#define MLHAND    1         // man's left hand symbol found
#define MRHAND    2         // man's right hand symbol found
#define WLHAND   10         // woman's left hand symbol found
#define WRHAND   20         // woman's right hand symbol found
*/
{
   var i,n;
   var dy,ylap;

   prevhold = hold;
   mface = -1;
   wface = -1;
   facedif = -1;
	if ((jm == Face)&&(oriented == FALSE)&&
		(((dofig == MAN)&&(jc < 0))||(dofig == WOMAN)&&(jc > 0)))
	{
		output += "linear    "+leadingZeros(0,3)+" "+leadingZeros(1,3)+" spinby fig    afoot  pelvis "+leadingZeros(((ji-1)*45),3)+" y\n";
		oriented = TRUE;
	}
   if ((jm == Limb)&&((ji == 4)||(ji == 9)))
   {
      if (jb ==  11) { ++holdss; ++holdsh; }
      if (jb ==  12) ++holdoe;
      if (jb ==  21) { ++holdco; ++holdcl; ++holdpr; }
      if (jb ==  22) ++holdsh;
      if (jb == 110) { ++holdcl; ++holdpr; }
      if (jb == 102) ++holdss;
      if (jb == 120) ++holdss;
      if (jb == 202) { ++holdcl; ++holdpr; }
   } /* jm = a hand */
   else
   if ((jm == Face)&&(jx > stmiddle))
   {
      n = -1;
      ylap = -1;
      wface = ji;
      for (i = 1; i < 9; ++i)
      {
         n = lseeksym(Face,i,xmin,stmiddle,jy,jy2);
         if (n >= 0)
         {
            dy = loverlap(jy,jy2,lbn[n].y,lbn[n].y2);
            if (dy > ylap)
            {
               ylap = dy;
               mface = i;
            }
         } /* found man facing sign */
      }
      if (mface >= 0)
      {
         facedif = mface - wface;
         if (facedif < 0) facedif += 8;
         if (facedif > 7) facedif -= 8;
      }
      else
         facedif = -1;
      if (facedif == 0)
      {
         facecl = 0;
         facepr = 0;
         facesh = 1;
         facess = 1;
      } /* facing same way */
      else
      if (facedif == 2)
      {
         facecl = 0;
         facepr = 1;
         facesh = 0;
         facess = 0;
      } /* facing at right angles */
      else
      if (facedif == 4)
      {
         facecl = 1;
         facepr = 0;
         facesh = 0;
         facess = 0;
      } /* facing opposite ways */
   } /* jm == Face */
   if (holdoe > 1) if (hold != CO) hold = OE;
   if (holdco > 1) if (hold != OE) hold = CO;
   if ((facesh+holdsh) > 4) hold = SH;
   if ((facess+holdss) > 4) hold = SS;
   if ((facepr+holdpr) > 4) hold = PR;
   if ((facecl+holdcl) > 4) hold = CL;
	output += "** lsethold "+hold+" "+prevhold+",  "+facesh+" "+holdsh+",  "+facess+" "+holdss+",  "+facepr+" "+holdpr+",  "+facec1+" "+holdc1+", "+leadingZeros(mface,3)+" "+leadingZeros(wface,3)+" "+leadingZeros(facedif,3)+"\n";
   if (prevhold != hold) ldoposn();
} /* lsethold */
/********************************************/

var ldochest = function(piv)
/*
   rotate the chest and stomach
   
   called by ldolimb,
*/
{
   if (piv == 0)
   {
	   output += "quadratic "+leadingZeros(              fstart, 3)+" "+leadingZeros(fend, 3)+" bendto chest   ribs  stomach 0 0 0\n";
	   output += "quadratic "+leadingZeros(              fstart, 3)+" "+leadingZeros(fend, 3)+" bendto stomach waist pelvis 0 0 0\n";
   } /* piv == 0 */
   else
   {
      if (dofig == MAN)
         output += "quadratic "+leadingZeros(              fstart, 3)+" "+leadingZeros(fend, 3)+" rotate chest   ribs "+leadingZeros(-piv/2, 3)+"\n";
      else
	      output += "quadratic "+leadingZeros(              fstart, 3)+" "+leadingZeros(fend, 3)+" rotate chest   ribs "+leadingZeros(piv/2, 3)+"\n";
      output += "quadratic "+leadingZeros(           fstart, 3)+" "+leadingZeros(fend, 3)+" rotate stomach waist "+leadingZeros(piv/2, 3)+"\n";
   } /* piv != 0 */
} /* ldochest */
/******************************************/

var ldolimb = function()
/*
   do something to some body part
   
   called by laction,
   calls ldoarms, ldochest,

	Volm 7 + Area 9 = chest
*/
{
	var nc;
	var piv;

	nc = jc+8;
	piv = -1;
	if ( (colm[nc] == ARM)&&(jm == Dirn)&&
		((hold == NO)||(hold == OE)||(hold == CO)) )
        ldoarms();
	else
	if (jm == Limb)
		colm[nc] = Limb;
	else
	if ((jm == Volm)&&(ji == 7)
		&&(colm[nc] == Area)&&(jd == BLANK))
	{
		colm[nc] = CHEST;
	   output += "* ldolimba CHEST at column "+nc+"\n";
	}
	else
	if ((jm == Area)&&(ji == 9)
		&&(colm[nc] == Volm)&&(jd == BLANK))
		colm[nc] = CHEST;
	else
	if ((jm == Area)&&(ji == 9))
		colm[nc] = Area;
	else
	if ((jm == Volm)&&(ji == 7))
		colm[nc] = Volm;
	else
	if ((jm == Rotn)&&(colm[nc] == CHEST))
	{
		piv = lgetpin();
		ldochest(piv);
	}
} /* ldolimb */
/*********************************************/

//var lcoords = function(char jm, int ji)
var lcoords = function(jm, ji)
/*
	check for change of coordinates

	called by laction,
	calls lseeksym, lgetpin

	Relevant symbols:-
	m      i
	Volm   5  space hold
	Volm   6  coordinates
	Volm   7  body hold
	Area   9  square
	Pins   1  forward
	Pins   5  backward
	
	 1 Aug 2006 checking piv against maxint
	30 Jul 2006 writing bendtos for mspace and wspace
*/
{
	var k;
	var piv;

	if ((jm == Area)&&(ji == 9))
	{
		piv = lgetpin ( );
		//fprintf(nudesfile,"* lcoordsa %c %d\n",m,piv);
		if (piv != maxint)
		{
			if (piv == 360) piv = 0;
			//coordinates = SPACE;
			if ( dofig == MAN )
			{
				output += "repeat "+ 					fstart+" "+ fend+" bendto mspace jman joist 270 0 "+ piv+"\n";
			   mspace = true;
			}
			else
			{
				output += "repeat "+ 					fstart+" "+ fend+" bendto wspace jwoman joist 270 0 "+ piv+"\n";
				wspace = TRUE;
			}
		} /* space stance found */
	} /* possible space stance found */
	else
	{
		k = lseeksym(Volm,7,jx,jx2,jy,jy2);
		if (k > 0)
		{
			//coordinates = BODY;
			if ( dofig == MAN )
			   mspace = false;
			else
				wspace = FALSE;
		} /* body stance found */
		//fprintf(nudesfile,"* lcoordsb mspace wspace TRUE\n",
			//mspace,wspace,TRUE);
	} /* possible body stance found */
}  /* lcoords */
/*****************************************/

var ldotoetaps = function()
/*

	do toe taps with gestures of the legs
	doing diagonals sideways at present

	Volm   1  RELAX
	Volm   3  BENT
	Volm   2  STRAIGHT
	Volm   4  STRETCH
	Volm   7  hold

	called by laction,
	calls lgetout, lsetframes, bell,

	19 Aug 2006 d076- Don Herbison-Evans
*/
{
	if ( (( jc == -3 )||( jc == 3 )) && ( jd == -1 ) )
	{
			output += "*\n";
			if  ( ji==11 )
				output += "* in place tap\n";
			else if ( ( ji == 1 ) || ( ji == 10 ) )
				output += "* forward tap\n";
			else if ( ( ji == 2 ) || ( ji == 9 ) )
				output += "* forward diagonal tap\n";
			else if ( ( ji == 3 ) || ( ji == 8 ) )
				output += "* sideways tap\n";
			else if ( ( ji == 4 ) || ( ji == 7 ) )
				output += "* back diagonal tap\n";
			else if ( ( ji == 5 ) || ( ji == 6 ) )
				output += "* backward tap\n";
			//
			if ( dofig == MAN )
			{
				if (mspace == false)
					output += "repeat    "+leadingZeros(fstart,3)+" "+leadingZeros(fend,3)+" set    coords mpelvis\n";
				else
					output += "repeat    "+leadingZeros(fstart,3)+" "+leadingZeros(fend,3)+" set    coords mspace\n";
			}
			else
			{
				if (wspace == FALSE)
					output += "repeat    "+leadingZeros(fstart,3)+" "+leadingZeros(fend,3)+" set    coords wpelvis\n";
				else
					output += "repeat    "+leadingZeros(fstart,3)+" "+leadingZeros(fend,3)+" set    coords wspace\n";
			}
			//
			if ( jc < 0 )
			{
				if ( ( ji <= 1 ) || ( ji == 3 ) || ( ji == 5 ) || ( ji > 11 ) )
				{
					console.log("OOPS: ldotoetap left direction problem line "+ j +"\n");
					console.log(""+leadingZeros( jm, 3)+" "+leadingZeros( ji, 3)+" "+leadingZeros( jx, 3)+" "+leadingZeros( jy, 3)+" "+leadingZeros( js, 3)+" "+leadingZeros( jw, 3)+" "+leadingZeros( jh, 3)+" "+leadingZeros( jb, 3)+" "+ jd +"\n");
					lgetout ( 1 );
					if ( ok == 1 ) return;
				} /* i wrong */
				output += "repeat    "+leadingZeros(fstart,3)+" "+leadingZeros(fend,3)+" call   forleft * left = b\n";
			} /* left side */
			else if ( jc > 0 )
			{
				if ( ( ji < 1 ) || ( ji == 6 ) || 
					( ji == 8 ) || ( ji == 10 ) || ( ji > 11 ) )
				{
					console.log("OOPS: ldotoetap right direction problem line "+ j +"\n");
					console.log(""+leadingZeros( jm, 3)+" "+leadingZeros( ji, 3)+" "+leadingZeros( jx, 3)+" "+leadingZeros( jy, 3)+" "+leadingZeros( js, 3)+" "+leadingZeros( jw, 3)+" "+leadingZeros( jh, 3)+" "+leadingZeros( jb, 3)+" "+ jd +"\n");
					lgetout ( 1 );
					if ( ok == 1 ) return;
				} /* i wrong */
				output += "repeat    "+leadingZeros(fstart,3)+" "+leadingZeros(fend,3)+" call   forright * right = b\n";
			} /* right side */
//
			if ( ji == 11 )
			{
			
				output += "repeat    "+leadingZeros( 					fstart, 3)+" "+leadingZeros( fend, 3)+" call   "+ risesub[rise] +"\n";
				output += "repeat    "+leadingZeros( 					fstart, 3)+" "+leadingZeros( fend, 3)+" set    fend  "+ frange +"\n";
				output += "linear    "+leadingZeros( 					fstart, 3)+" "+leadingZeros( fend , 3)+" bendto bleg   bknee  bthigh lrlx1 lrlx2 lrlx3\n";
			} /* close without weight */
			else
				output += "linear    "+leadingZeros( 					fstart, 3)+" "+leadingZeros( fend , 3)+" bendto bleg   bknee  bthigh lhig1 lhig2 lhig3\n";
			output += "linear    "+leadingZeros( 				fstart, 3)+" "+leadingZeros( fend , 3)+" drag   bfoot  bfoot  bankle bleg  x\n";
			lbn[j].a = DONE;
	} /* c OK */
} /* ldotoetaps */
/**************************************/

var laction = function()
/*
   run through and interpret the actions

   called by linter,
   calls     ldobar,   ldosteps, lleggesture, ldolimb,
             ldopivot, lbent,    lassign,  lsetframes,
             lsethold, ldohold,  lrelease, lface,

#define FRONT   100         // front symbol found
#define BACK    200         // back symbol found
#define MLHAND    1         // man's left hand symbol found
#define MRHAND    2         // man's right hand symbol found
#define WLHAND   10         // woman's left hand symbol found
#define WRHAND   20         // woman's right hand symbol found 
	
#define NO        0        // no hold
#define CL        1        // closed hold
#define PR        2        // promenade position
#define CP        3        // counter promenade position
#define DB        4        // double hold
#define OP        5        // open hold
#define CR        6        // crossed open hold
#define OE        7        // open extended hold
#define CO        8        // counter open extended hold
#define SH        9        // shadow hold
#define SS       10        // semi-shadow hold 

Relevant symbols:-
     m     i
    Misc   1  bow
    Misc   2  release1
    Misc   3  release2
    Limb   4  lhand
    Limb   9  rhand
    Area   1  top/front 
    Area   5  back/bottom
	 Area   9  square
    Volm   1  RELAX
    Volm   3  BENT
    Volm   2  STRAIGHT
    Volm   4  STRETCH
	 Volm   6  coordinates
    Volm   7  hold
    Face   n  facing direction

*/
{
	output += "*\n************************************\n";
	oriented = FALSE;
	if ( dofig == MAN )
		output += "*\nrepeat      0 "+leadingZeros( fmax , 3)+" call   doman\n";
	else
		output += "*\nrepeat      0 "+leadingZeros( fmax , 3)+" call   dowoman\n";
	for ( j = 0; j < NCOLM; ++j )
		colm[j] = ARM;
	for ( j = 0; j < ssend; ++j )
	{
		lassign ();
		lsetframes ();
		output += "* "+lbn[j].a+" "+leadingZeros(jc, 3)+" "+lbnline[j]+"";
		if ( lbn[j].a == TODO )
		{
			if ( jm == Bars )
				ldobar ();
			else if ( ( jm == Face ) || ( jm == Limb ) )
				lsethold ();
			else if ( jm == Misc )
			{
				lrelease ();
			}
			else if ( ( jc > -8 ) && ( jc < 8 ) )
			{
				if ( (( jm == Volm )&&( ji == 6 )) 
					||(( jm == Area )&&( ji == 9 )) )
						lcoords(jm, ji);
				if ( ( jm == Rotn ) && ( jc > -4 ) && ( jc < 4 ) )
					ldopivot ();
				else if (( jm == Dirn ) && ( jc > -4 ) && ( jc < 4 ))
				{
					ldostep ();
					lleggesture ();
					ldotoetaps ();
				}
				else
					ldolimb ();
			} /* jc OK */
		} /* ja == TODO */
		if (( (jm == Dirn)||(jm == Rotn) )&&(jc >= -6)&&(jc <= 6)
			&&( nmw > 0 )&&( dofig == WOMAN ) )
			ldohold ();
		pstart = fstart;
		pend = fend;
	} /* j */
} /* laction */
/*************************************************/

var linter = function()
/*
                     linter

      interpret labanotation score into a NUDES file
                version linter50.c

      input : LED Labanotation file:   standard input (led.lbn)
      output: NUDES animation script:  standard output (led.n)

   called by main,
   calls     lbnread, lsorty, lfindstaff, lstart, lhold,
             lfindystart, lcolx, lsetrange, lselectfig,
             lgetout, lcopyfigs, lfinish, lcopysubs,
             lbows,
*/
{
   lbnread();;
   lsorty();
   lfindstaff();
   lsetrange();
   lselectfig();
   lcopyfigs();
   lstart();
   lfindystart();
   lbows(); // flag hand signs
   lbent(); // flag dirn signs
   for (st = 0; st < nstaff; ++st)
   {
      hold = NO;
      holdcl = 0;
      holdco = 0;
      holdoe = 0;
      holdpr = 0;
      holdsh = 0;
      holdss = 0;
      facecl = 0;
      facepr = 0;
      facesh = 0;
      facess = 0;
      prevhold = -9;
      prevc = 0;
      pstart = -1;
      pend = -1;
      keptf = 0;
      gy = -1;
      gh = 0;
      if (staff[st][5] == TODO)
      {
         nbar = -1;
         if (staff[st][4] == MAN)
            dofig = MAN;
         else
            dofig = WOMAN;
         lcolx(staff[st][2]);
         laction();
         staff[st][5] = DONE;
      }
   }
   lfinish();
} /* linter */
/****************************************/

//var shift = function(double x, double y, double z)
var shift = function(x, y, z)
/*
   this adds 'x,y,z' to all centres and joints in lists
   'elist' and 'jlist'.

   called by  action, dogrofig, dogrojt, domovjnt,
              twirl, dodrag,
*/
{
   var e,j,n ;


   for (  n = 0 ; n < ecount ; ++ n )
   {
      e = elist[n] ;
      cen[e][0] += x ;
      cen[e][1] += y ;
      cen[e][2] += z ;
   }
   for (  n = 0 ; n < jcount ; ++ n )
   {
      j = jlist[n] ;
      jnt[j][0] += x ;
      jnt[j][1] += y ;
      jnt[j][2] += z ;
   }
}  /* shift */
/*****************************/

//var rset = function(double r[3][3], double angl, int axis)
var rset = function(r, angl, axis)
/*
   set up the rotation matrix 'r' for a rotation of
   'angl' radians about 'axis'.

   called by  input, setobs, dobalanc, dospinby,
*/
{
	  var v = new Array();	//v[5]
      var i,j,k;

      v[0] = doub0 ;
      v[1] = doub1 ;

/*   fill out values vector with sin and cos- */

	//TODO REVIEW (May be degrees or radians)
      v[2] = Math.cos(angl) ;
      v[3] = Math.sin(angl) ;
      v[4] = -v[3] ;

/*   choose appropriate permutation of values for rotation axis- */

      for (  i = 0 ; i < 3 ; ++ i )
      {
         for (  j = 0 ; j < 3 ; ++ j )
         {
            k = perm[axis][j][i] ;
            r[i][j] = v[k-1] ;
         }
      }
}  /* rset */
/************************************/

//var matmul = function(double a[3][3], double b[3][3], double c[3][3])
var matmul = function(a, b, c)
/*
     this multiplies matrix 'b' by 'a' and puts the product
     in 'ans'.

     called by  dobalanc, matrot, dospinto, dospinby, getwist.
                getaxes, sepn,  getmat,

	  21 Sep 2006  unrolled loops
*/
{
	var ans00,ans01,ans02,ans10,ans11,ans12,ans20,ans21,ans22;
//
	ans00 = a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0];
	ans01 = a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1];
	ans02 = a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2];
	ans10 = a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0];
	ans11 = a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1];
	ans12 = a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2];
	ans20 = a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0];
	ans21 = a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1];
	ans22 = a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2];
//
	c[0][0] = ans00;
	c[0][1] = ans01;
	c[0][2] = ans02;
	c[1][0] = ans10;
	c[1][1] = ans11;
	c[1][2] = ans12;
	c[2][0] = ans20;
	c[2][1] = ans21;
	c[2][2] = ans22;
}  /* matmul */
/**********************************************************/

//var vecmul = function(double v[EMAX][3], double m[3][3], int n)
var vecmul = function(v, m, n)
/*
   multiply the 'n'th vector from array 'v'
   by matrix 'm'.

   called by touch, dogrojnt, domovjnt, domoveby, doabut,
             twirl,
*/
{
      var i,j ;
      var x;
	  var vv = new Array(); //vv[3]

      for (  i = 0 ; i < 3 ; ++ i )
      {
         x = doub0 ;
         for (  j = 0 ; j < 3 ; ++ j )
         {
            x = x+m[i][j]*v[n][j] ;
         }
         vv[i] = x ;
      }
      
      for (  i = 0 ; i < 3 ; ++ i )
      {
         v[n][i] = vv[i] ;
      }
}  /* vecmul */
/**********************************************/

//var rotget = function(double r[3][3], double unr[3][3], int n)
var rotget = function(r, unr, n)
/*
   form a rotation matrix r and its inverse unr
   from the nth entries in quat

   called by  dobalanc, matrot, dospinto, dospinby,
              dogrojnt, domovjnt, doabut, doground,
*/
{
      var i,j ;
      var cp,sp,x,y,z,m,xsp,ysp,zsp,xm,ym,zm,xym,xzm,yzm ;

      x = quat[n][0] ;
      y = quat[n][1] ;
      z = quat[n][2] ;
      sp = quat[n][3] ;
      cp = quat[n][4] ;
      m = doub1-cp ;
      xm = x*m ;
      ym = y*m ;
      zm = z*m ;
      xsp = x*sp ;
      ysp = y*sp ;
      zsp = z*sp ;
      xym = x*ym ;
      xzm = x*zm ;
      yzm = y*zm ;
      r[0][0] = x*xm+cp ;
      r[0][1] = xym+zsp ;
      r[0][2] = xzm-ysp ;
      r[1][0] = xym-zsp ;
      r[1][1] = y*ym+cp ;
      r[1][2] = yzm+xsp ;
      r[2][0] = xzm+ysp ;
      r[2][1] = yzm-xsp ;
      r[2][2] = z*zm+cp ;

      for (  i = 0 ; i < 3 ; ++ i )
      {
         for (  j = 0 ; j < 3 ; ++ j )
         {
            if ((r[j][i] > -tolr) && (r[j][i] < tolr)) r[j][i] = 0;
            unr[i][j] = r[j][i] ;
         }
      }
}  /* rotget */
/**************************************/

//var rotput = function(double r[3][3], int n)
var rotput = function(r, n)
/*
   interpret rotation matrix 'r' as direction cosines of a
   rotation axis, and the sine and cosine of a rotation about
   that axis, and store in array 'quat'.

   uses the fact that any rotation matrix can be written as -

   ( x.x.m+c    x.y.m-z.s  x.z.m+y.s )
   ( x.y.m+z.s  y.y.m+c    y.z.m-x.s )
   ( x.z.m-y.s  y.z.m+x.s  z.z.m+c   )

   where
     x,y,z-components of unit vector along rotation axis
             x=cos(a1)cos(a2)  y=cos(a1)sin(a2)  z=sin(a1)
             a1,a2-azimuth and elevation of axis from x axis
     s,c  -sine and cosine of rotation about that axis
     m     = 1-c

     x,y,z are stored in quat[n,0], quat[n,1], quat[n,2]
     s,c   are stored in quat[n,3], quat[n,4]

   see 'Control of round-off propagation in articulating the
        human figure', D.Herbison-Evans and D.S.Richardson,
        Computer Graphics and Image Processing,
        vol 17, pp. 386-393 (1981)

   called by matrot, dospinto, doangles, dolimb,
             getwist, store3,
*/
{
      var j,k ;
      //double a[3][3],b[3],d[3]
	  var a = get2DArray(3);
	  var b = new Array();
	  var d = new Array();
	  var e,f,g,c,s,trace ;
      var csq;

      b[0] = r[1][2]-r[2][1] ;
      b[1] = r[2][0]-r[0][2] ;
      b[2] = r[0][1]-r[1][0] ;
      e = b[0]*b[0]+b[1]*b[1]+b[2]*b[2] ;
      trace = r[0][0]+r[1][1]+r[2][2] ;
      if (e > doub0) g = sqrt(e); else g = doub0;
      if (e > tolr)
      {
         f = doub1/g ;
         quat[n][0] = f*b[0] ;
         quat[n][1] = f*b[1] ;
         quat[n][2] = f*b[2] ;
/*
     use g=2s, and trace=1+2c to find s and c -
*/
         s = inv2*g;
         csq = doub1-s*s;
         if (csq > doub0) c = sqrt(csq); else c = doub0;
         if (trace < doub1) c = -c;
         quat[n][3] = s ;
         quat[n][4] = c ;
      }
      else
/*
   symmetric matrix (180 or 360 degree rotation) -
*/
      {
         c = inv2*(trace-doub1);
         for (  j = 0 ; j < 3 ; ++ j )
         {
            d[j] = doub0 ;

/*   run across a row- */

            for (  k = 0 ; k < 3 ; ++ k )
            {
               a[j][k] = r[j][k]+r[k][j] ;
               if (j == k) a[j][j] = doub2*(r[j][j]-c) ;
               d[j] = d[j]+a[j][k]*a[j][k] ;
            }
         }

/*   choose most stable row- */

         j = 0 ;
         if (d[1] > d[0]) j = 1 ;
         if (d[2] > d[j]) j = 2 ;
         if (d[j] > doub0) f = doub1/sqrt(d[j]) ;
         else
         {
            f = doub1;
            a[j][0] = doub1;
         }
         quat[n][0] = f*a[j][0] ;
         quat[n][1] = f*a[j][1] ;
         quat[n][2] = f*a[j][2] ;
         quat[n][3] = inv2*g ;
         quat[n][4] = c ;
      }
      for (k = 0; k < 5; ++k)
      {
         if ((quat[n][k] > -tolr) && (quat[n][k] < tolr))
            quat[n][k] = 0;
         if (quat[n][k] >  doub1) quat[n][k] =  doub1;
         if (quat[n][k] < -doub1) quat[n][k] = -doub1;
      }
}  /* rotput */
/********************************************/

//var mkquat = function(int n, double a1, double a2, double a3)
var mkquat = function(n, a1, a2, a3)
/*
   convert angles a1,a2,a3 (in radians) into quat entries

   called by dospinto, inframe, 
*/
{
      var j;
      var s1,c1,s2,c2,s3,c3 ;

      s1 = sin(a1) ;
      c1 = cos(a1) ;
      s2 = sin(a2) ;
      c2 = cos(a2) ;
      s3 = sin(a3) ;
      c3 = cos(a3) ;
      quat[n][0] = c1*c2 ;
      quat[n][1] = s1*c2 ;
      quat[n][2] = s2 ;
      quat[n][3] = s3 ;
      quat[n][4] = c3 ;
      for (j = 0; j < 5; ++j)
         if ((quat[n][j] > -tolr) && (quat[n][j] < tolr)) quat[n][j] = 0;
}  /* mkquat */
/**********************************************************/


//var matrot = function(double r[3][3], int n)
var matrot = function(r, n)
/*
      this rotates the 'n'th ellipsoid by rotation matrix 'r'.

      called by twirl.
      calls     rotget, matmul, rotput,
*/
{
      //double ro[3][3],unro[3][3] ;
	  var ro = get2DArray(3);
	  var unro = get2DArray(3);
	  
      rotget(ro,unro,n) ;
      matmul(r,ro,ro) ;
      rotput(ro,n) ;
}  /* matrot */
/**********************************************/

//var twirl = new function(double x, double y, double z, double r[3][3])
var twirl = new function(x, y, z, r)
/*
   rotates all the rotation matrices 'quat', centres 'cen',
   and joint vectors 'jnt', of ellipsoids and joints in lists
   'elist' and 'jlist' about a point 'x,y,z' using rotation
   matrix 'r'.

   called by  dospinto, dospinby, store3,
   calls      shift, matrot, vecmul, setels,
*/
{
      var e,j,k ;

      shift(-x,-y,-z) ;
      if (ecount >= 0)
      {

/*   rotate the ellipsoids and their centres- */

         for (  e = 0 ; e < ecount ; ++e )
         {
            k = elist[e];
/*  don't rotate world ! :- */
            if (k != 0)
            {
               matrot(r,k) ;
               vecmul(cen,r,k) ;
            }
         }
      }

/*   now for the joints- */

      if (jcount >= 0)
      {
         for (  j = 0 ; j < jcount ; ++j )
         {
            k = jlist[j];
            vecmul(jnt,r,k) ;
         }
      }

/*   put body part back where it came from- */
      shift(x,y,z) ;
}  /* twirl */
/*****************************/

//var dospinto = function(double xx[3], int refell, double ang[3], double pro)
var dospinto = function(xx, refell, ang, pro)
/*
     spins all ellipsoids in 'elist' and joints in 'jlist'
     so that 'ellpsd' is proportion 'pro' of the way to the
     orientation specified as a rotation 'ang' radians
     about axes of the reference ellipsoid 'refell'
     about point 'xx'.

   called by  action, dodrag,
   calls      rotget, rotput, mkquat, matmul, twirl,
*/
{
      var alfa,nualfa;
      //double mv[3][3],unmv[3][3];
      //double rf[3][3],unrf[3][3];
      //double tg[3][3],untg[3][3];
      //double mt[3][3],nu[3][3];
	  
	  var mv = get2DArray(3);
	  var unmv = get2DArray(3);
	  var rf = get2DArray(3);
	  var unrf = get2DArray(3);
	  var tg = get2DArray(3);
	  var untg = get2DArray(3);
	  var mt = get2DArray(3);
	  var nu = get2DArray(3);
/*
   set rotation matrices of moving and reference ellipsoids -
*/
      rotget(mv,unmv,ellpsd);
      rotget(rf,unrf,refell);

/*   find target rotation matrix, and refer to refell- */

      mkquat(EMAX+1,ang[0],ang[1],ang[2]);
      rotget(tg,untg,EMAX+1);
      matmul(rf,tg,tg);

/*   find increment rotation matrix to reach target- */

      matmul(tg,unmv,mt);
      rotput(mt,EMAX+1);
      if (( quat[EMAX+1][3] == doub0 ) 
		  && ( quat[EMAX+1][4] == doub0 ))
      {
         ok = 53;
		 printf("dospinto no sine and cosine");
         alfa = doub0;
      }
      else alfa = atan2(quat[EMAX+1][3],quat[EMAX+1][4]) ;
      nualfa = pro*alfa ;
      if (alfa > pi ) nualfa = pro*(alfa - twopi);
      if (alfa < -pi) nualfa = pro*(alfa + twopi);
      quat[EMAX+1][3] = sin(nualfa);
      quat[EMAX+1][4] = cos(nualfa);
      rotget(nu,mt,EMAX+1);
      twirl(xx[0],xx[1],xx[2],nu);
}  /* dospinto */
/*************************************/

//var dospinby = new function(double xx[3], int refell, double angl, int axis)
var dospinby = function(xx, refell, angl, axis)
/*
   spins all ellipsoids and joints in 'elist' and 'jlist'
   about a point 'x', by an angle 'angl' radians relative to
   an 'axis' of reference ellipsoid 'refell'.

   called by dobalanc, action, dobend, dotouch, fun, dodrag,
   calls     rset, rotget, matmul, twirl,
*/
{
      var j,k = 0;
      //double r[3][3],ro[3][3],unro[3][3];
	  var r = get2DArray(3);
	  var ro = get2DArray(3);
	  var unro = get2DArray(3);
/*
        do transformation on required coordinates
        aligned with axes of the reference ellipsoid-
*/
      rset(r,angl,axis);
      rotget(ro,unro,refell);
      matmul(r,unro,r);
      matmul(ro,r,r);
      for (j = 0; j < 3; ++j)
         for (k = 0; k < 3; ++k)
            if ((r[j][k] > -tolr) && (r[j][k] < tolr)) r[j][k] = 0;
      twirl(xx[0],xx[1],xx[2],r);
}  /* dospinby */
/**********************************************************/

//var mkang = function(int n)
var mkang = function(n)
/*
   get angles in radians from 'n'th entry in 'quat' into 
   array 'ang'.

   called by  doangles, store3, storeang,
*/
{
      var x,y,z,s1,c1,m1 ;
      var j;

      x = quat[n][0] ;
      y = quat[n][1] ;
      z = quat[n][2] ;
      s1 = z ;
      m1 = doub1-z*z ;
      if (m1 > doub0) c1 = sqrt(m1) ;
         else c1 = doub0 ;
      if ((x == doub0 ) && ( y == doub0))
          ang[0] = doub0;
      else
          ang[0] = atan2(y,x) ;
      if ((s1 == doub0 ) && ( c1 == doub0))
      {
          ok = 54;
          printf("mkang: n %d, s1 %f, c1 %f",            
              n,s1,c1);
          ang[1] = doub0;
      }
      else ang[1] = atan2(s1,c1) ;
      if ((quat[n][3] == doub0 ) && ( quat[n][4] == doub0))
      {
          ok = 52;
          printf("mkang: n %d, quat[n][3] %f, quat[n][4] %f",            
              n,quat[n][3],quat[n][4]);
          ang[2] = doub0;
      }
      else 
	  ang[2] = atan2(quat[n][3],quat[n][4]) ;
      for (j = 0; j < 3; ++j)
      {
         if (ang[j] < doub0) ang[j] += twopi;
         if (ang[j] > twopi) ang[j] -= twopi;
      }
}  /* mkang */
/*****************************************/

//var storeang = function(int f, int e, double a1, double a2, double a3)
var storeang = function(f, e, a1, a2, a3)
/*
   convert angles a1,a2,a3 in degrees into quaternions
   and find direction vector of y axis
   for frame f and ellipsoid e

   called by store3,
*/
{
      var s1,c1,s2,c2,s3,c3;

      s1 = sin(a1) ;
      c1 = cos(a1) ;
      s2 = sin(a2) ;
      c2 = cos(a2) ;
      s3 = sin(a3) ;
      c3 = cos(a3) ;

      qu3[f][e][0] = a3*degree ;
      qu3[f][e][1] = c2*c1 ;
      qu3[f][e][2] = c2*s1 ;
      qu3[f][e][3] = -s2 ;

}  /* storeang */
/**********************************************************/

//var doangles = function(int el, int re, double val[EMAX], int var0, int var1, int var2)
var doangles = function(el, re, val, var0, var1, var2)
/*
  store the angles of 'el' relative to 're' in 'val' array.
  in degrees.

  called by action, dodrag,
  calls  matmul, rotget, rotput, mkang,
*/
{
   //double mvro[3][3],mvunro[3][3];
   //double stro[3][3],stunro[3][3];
   //double r[3][3];
   var mvro = get2DArray(3);
   var mvunro = get2DArray(3);
   var stro = get2DArray(3);
   var stunro = get2DArray(3);
   var r = get2DArray(3);

   rotget(stro,stunro,re) ;
   rotget(mvro,mvunro,el) ;
   matmul(stunro,mvro,r) ;
   rotput(r,EMAX) ;
   mkang(EMAX) ;
   val[var0] = ang[0]*degree ;
   val[var1] = ang[1]*degree ;
   val[var2] = ang[2]*degree ;
   if ((val[var0] > doub179)&&(val[var0] < doub181))
   {
	   val[var0] -= doub180;
	   val[var2] = -val[var2];
   }
   if (val[var1] > doub180) val[var1] -= doub360;
}  /* doangles */
/*********************************/

//var dobend = function(double angle, int axis)
var dobend = function(angle, axis)
/*
  implements flex(38), rotate(39), abduct(40).

  called by action,
  calls     dospinby,
*/
{
   var refell ;
   var left ;

   //TODO REVIEW (Was a strangely placed goto)
   refell = ellpsd ;
   if (!(t == rotate_keyword_code)){
	   if (ellpsd == coel[join][0]) refell = coel[join][1] ;
	   if (ellpsd == coel[join][1]) refell = coel[join][0] ;
/*
  assume odd-numbered ellipsoids are on left side of figure-
*/
	}
	if (((ellpsd-figell[fig])%2) == 0)
          left = TRUE; else left = FALSE;
/*
  flex-
*/
   if ((t == flex_keyword_code)&&(knee[join])) angle = -angle ;
/*
  rotate-
*/
   if ((t == rotate_keyword_code)&&( left == FALSE)) angle = -angle ;
/*
  abduct-
*/
   if ((t == abduct_keyword_code)&&(left == TRUE)) angle = -angle ;
   dospinby(xx,refell,angle,axis) ;
}  /* dobend */
/****************************************************/

/*   compl42.h - based on complu

     This translates a NUDES script into a compact
     form for use by 'perfrm'.

   subroutines-
      getout
      llength
      nexts
      match
      value
      addname
      getint
      inells
      injts
      inlims
      inname
      dojoin
      checkin
      valadd
      parset
      inperf
      compl
*/

/***************************************/

//var getout = function(int v)
var getout = function(v)
/*
   exit gracefully

   called by main, inlims, openfile, compl, nexts, doperfrm,
             initsphere, getkeys,
*/
{
   if (v != 0) 
   {
	   printf("lintel problem\nok error number %d\n",ok);
	   printf("line %d, action %d\n%s\n",
		   pp,p,aline[pp]);
   }
   if (infile) fclose(infile);
   ok = 1;
} /* getout */
/********************************************/

var llength = function()
/*
   find length of line

   called by nexts,
*/
{
   var j,sp;

   sp = 0;
   for ( j = 0; line[j] != null; ++j);
   {
      if (line[j] != blank) sp = j;
   }
   return(sp);
} /* llength */
/*******************************************/

//var nexts_a = function( char c ) (returns int)
var nexts_a = function(c)
{
		var astk = '*';
		var tab = '	';
		var code;

		code = 0;
		if ( c == astk ) code = 1;
		if ( c == '\n' ) code = 2; 
		if ( c == blank
				|| c == tab
				|| c == null ) code = 3;
		if ( c == '.'
				|| ( c == '+' ) 
				|| ( c == '-' )
				|| ( c == '_' )
				|| ( c >= '0' ) && ( c <= '9' )
				|| ( c >= 'a' ) && ( c <= 'z' )
				|| ( c >= 'A' ) && ( c <= 'Z' ) ) code = 4;
		return code;
} /* nexts_a */
/*****************************************/

//TODO REVIEW (was a convoluted set of GOTOs, now converted into functions.)
var nexts = function()
/*
     this picks the next continuous string of non-blank integers
     from 'line', starting at 'start'.

     if an asterisk is read, the data is assumed to continue onto
     the next line.

  input -
     line - an image of the current line being scanned.
     start - the start of the scan.

  output -
     line - an image of the current line being scanned.
     start - the start of the next non-blank string in line.
     string - a copy of the first 6 characters of the next
              non-blank string
     length - the length of the string.

     called by inperf, inname, inells, injts, inlims, parset,
     calls     llength, getout,
*/
{
   var j;
   var astk = '*';
   var tab = '	';

   length = 0;

   for (  j = 0 ; j < BMAX ; ++ j )
      string[j] = null ;
/*
     get a new line if required-
*/
   if ((start < lline) && (start > 0)) lab17() ;
	function lab10(){
	   start = 0 ;
	   if (fgets(line,BMAX,infile) == NULL)
	   {
		  console.log("\nOOPS in nexts: unexpected end of file\n");
		  console.log("missing STOP command?\n");
		   ok = 3;
		  getout(ok);
		  return;
	   }
	   ++nline ;
	   lline = llength();
	/*
		 find start of next string-
	*/
	}
	function lab17(){
   for ( j = start ; j < lline ; ++ j )
   {
      if (line[j] == astk) lab10() ;
      if (line[j] == blank) lab1() ;
      if (line[j] == tab) lab1() ;
      if (line[j] == null) lab1() ;
      if (line[j] == '.') lab3();
      if ((line[j] == '+') || (line[j] == '-')) lab3();
      if ((line[j] >= '0') && (line[j] <= '9')) lab3();
      if ((line[j] >= 'a') && (line[j] <= 'z')) lab3();
      if ((line[j] >= 'A') && (line[j] <= 'Z')) lab3();
	  if (line[j] == '_') lab3();
      lab10();
	 }
	}
	function lab1(){
/*
     rest of line empty, so look at next -
*/
   j = 0;
   lab10() ;
   }
/*
     copy up to the characters of the string,
     then move start to next blank-
*/
	function lab3(){
	   for ( start = j ; start < lline ; ++start )
	   {
		  if (line[start] == '\n') return;
		  if (line[start] == blank) return;
		  if (line[start] == tab) return;
		  if (line[start] == null) return;
		  if (line[start] == '.') lab5();
		  if ((line[start] == '+') || (line[start] == '-')) lab5();
		  if ((line[start] >= '0') && (line[start] <= '9')) lab5();
		  if ((line[start] >= 'a') && (line[start] <= 'z')) lab5();
		  if ((line[start] >= 'A') && (line[start] <= 'Z')) lab5();
		  if (line[start] == '_') lab5();
		  return;
		function lab5(){
			  string[length] = line[start] ;
			  ++length ;
		  }
	   }
	   start = -1 ;
   }
} /* nexts */
/***************************************/

//var match = function(int nnames, int lengths[EMAX], char names[EMAX][BMAX]) (returns int)
var match = function(names, lengths, names)
/*
     find which of 'names' fits 'string'.

 input
   string - array of a1 characters to be matched.
   nnames - number of names actually stored in array names.
   names - list of possible strings for which to search.
   lengths - lengths of each of the existing names
   length - length of 'string'

     output
   which of the names that string matches, -1 if none.

     called by inperf, innames, inells, injts, inlims,
               parset, addnam,
*/
{
   var j,k ;
   var found;
   var no;
   var sp = null;

   no = -1 ;
   for (k = 0; k <= nnames; ++k)
   {
      found = TRUE;
      if ((lengths[k] == 0) || (lengths[k] == length))
      {
         for (j = 0; j < length; ++j)
         {
            if (string[j] != names[k][j]) found = FALSE;
         }
         if (found == TRUE)
         {
            no = k ;
            break;
         }
      }
   }
   return(no);
} /* match */
/***************************************/

//TODO REVIEW (GOTOs replaced with continue-s + some logic moved)
//var value = function() (returns double)
var value = function()
/*
     find the value of the number which is encoded as
     'length' characters in array 'string' and put it into 'v'.
     set ok false if string is not a number.

     called by inperf, inells, injts, parset,
*/
{
	var v; //double v;
	var nsign; //double nsign;
	var expon;	//double expon;
	var k; //int k;
	var frac; //int frac;
	var d; //int d;
	var point = '.'; //char point = '.';
	var minus = '-'; //char minus = '-';
	var plus = '+'; //char plus = '+';

	pok = TRUE;
	v = doub0;
	nsign = doub1;
	expon = doub1;
	frac = FALSE;
	if ( ( length < 0 ) || ( length > BMAX ) )
	{
		pok = FALSE;
		return( v );
	}
	for ( k = 0; k < length; ++ k )
	{
		//	if a decimal point encountered, start decimal place counter
		if (string[k] == point )
		{
			frac = TRUE;
			continue;
		}
		if ( string[k] == plus ){
			continue;
		}
		if ( string[k] == minus )
		{
			nsign = -nsign;
			continue;
		}
		for ( d = 0; d < 10; ++ d )
		{
			if ( string[k] == dig[d] ){
				v = v * doub10 + double(d);
				if ( frac == TRUE ) expon = expon / double(10);
				continue;
			}
		}
		pok = FALSE;
		return ( v );
	}
	v = v * expon * nsign;
	return( v );
} /* value */
/***************************************/

//int addnam(int n, char names[EMAX][BMAX], int isvar, int lengths[EMAX])
var addnam = function(n, names, isvar, lengths)
/*
     add a name to a list of names.
	 return the number of names in the list.

     called by inperf, inname, inells, injts, parset,
     calls     match,
*/
{
   var k; //int k;
   var no //int no ;
   var nnames //int nnames;

   nnames = n;
/*
     see if name already exists-
*/
   no = match(nnames,lengths,names);
   if (no >= 0)
   {
      console.log("addnam: string  "+string+" confusable with ");
      for (k = 0; names[no][k] != null; ++k)
         console.log(""+names[no][k]+"");
      console.log("\n");
   }
   if (isvar == FALSE)
/*
     non-variables must first check variable list-
*/
   {
      no = match(nvars,varlen,vname) ;
      if (no > 0)
      {
         console.log("name  "+string+"  confusable with variable ");
         for ( k = 0; vname[no][k] != null ; ++ k)
            console.log(""+vname[no][k]+"");
         console.log("\n");
      }
   }
   else
/*
     variables must check all name lists-
*/
   {
      no = match(nfigs,figlen,fname) ;
      if (no > 0)
      {
         console.log("variable  "+             string+"  confusable with figure ");
         for ( k = 0; fname[no][k] != null; ++ k)
            console.log(""+fname[no][k]+"");
         console.log("\n");
      }
      no = match(ne,ellen,ename) ;
      if (no > 0)
      {
         console.log("variable  "+             string+"  confusable with ellipsoid ");
         for ( k = 0; ename[no][k] != null; ++ k)
            console.log(""+ename[no][k]+"");
         console.log("\n");
      }
      no = match(njts,jntlen,jname) ;
      if (no > 0)
      {
         console.log("variable  "+             string+"  confusable with joint ");
         for ( k = 0; jname[no][k] != null; ++ k)
            console.log(""+jname[no][k]+"");
         console.log("\n");
      }
      no = match(nsubs,sublen,sname) ;
      if (no > 0)
      {
         console.log("variable  "+             string+"  confusable with subroutine ");
         for ( k = 0; sname[no][k] != null; ++ k)
            console.log(""+sname[no][k]+"");
         console.log("\n");
      }
      no = match(nfiles,fillen,tname) ;
      if (no > 0)
      {
         console.log("variable  "+             string+"  confusable with file name ");
         for ( k = 0; tname[no][k] != null; ++ k)
            console.log(""+tname[no][k]+"");
         console.log("\n");
      }
   }
/*
     add name to list-
*/
   if (nnames > EMAX)
   {
      console.log("\nOOPS addnam: "+          string+" makes more than max of "+EMAX+" names\n");
      ok = 84 ;
   }
   else
   {
      for (  k = 0 ; k < length ; ++ k )
         names[nnames][k] = string[k] ;
      for ( k = length ; k < BMAX; ++ k )
         names[nnames][k] = null;
      lengths[nnames] = length;
   }
   ++nnames ;
   return(nnames);
} /* addnam */
/***************************************/

//int getint(void)
var getint = function()
/*
     find value of positive integer which is encoded as 'length'
     characters in array 'string', and put its value into 'k'.
     set 'pok' false if string not a positive integer.

     called by parset, inname,
*/
{
	var j, k, m, ths;
	var plus = '+';

	if ( length <= 0 )
	{
		k = 0;
		pok = FALSE;
	}
	else
	{
		pok = TRUE;
		k = 0;
		for( j = 0; j < length; ++j )
		{
			if ( string[j] != plus )
			{
				ths = -1;
				for( m = 0; m < 10; ++m )
					if ( string[j] == dig[m] ) ths = m;

				if ( ths < 0 )
				{
					pok = FALSE;
					return( k );
				}
				k = 10 * k + ths;
			}
		}
	}
	return( k );
} /* getint */
/***************************************/

//int inells(void)
var inells = function()
/*
     read in next ellipsoid and its axis lengths.

     called by inperf, injts,
     calls     nexts, match, addnam, value,
*/
{
	var el, k;

	nexts();
	el = match ( ne, ellen, ename );
	if ( el < 0 )
	{
		ne = addnam ( ne, ename, 0, ellen );
		el = ne - 1;
	}
	for ( k = 0; k < 3; ++ k )
	{
		nexts ();
		semiax[k] = value ();
		if ( pok == FALSE )
		{
			console.log("\nOOPS inells: ellipsoid snag with  "+ 				string +"\n");
			ok = 83;
			return ( el );
		}
	}
	return ( el );
} /* inells */
/***************************************/

//void injts(void)
var injts = function()
/*
     read in the next joint, the ellipsoids it connects, and the
     position of the joint relative to each ellipsoid centre.

     called by inperf,
     calls     nexts, addnam, inells,
*/
{
   var el,jt,k,e ;//int
   var klet = 'k';//char
   var nlet = 'n';//char
   var elet = 'e';//char

   nexts();
   njts = addnam(njts,jname,0,jntlen);
   jt = njts-1;
   if ( ok > 0 ){
	lab4();
	return;
   }
/*
     is it a knee -
*/
   knee[jt] = FALSE;
   for (  k = 0 ; k < (length-1) ; ++ k )
   {
      if ((string[k] == klet)
       && (string[k+1] == nlet)
       && (string[k+2] == elet)) knee[jt] = TRUE;
   }
/*
  do the two ellipsoids
*/
   for (  e = 0 ; e <= 1 ; ++ e )
   {
      el = inells();
      if ( ok > 0 ){
		lab5();
		return;
	  }

      dcon[jt][e][0] = semiax[0] ;
      dcon[jt][e][1] = semiax[1] ;
      dcon[jt][e][2] = semiax[2] ;
      coel[jt][e] = el ;
   }
   return;
/*
     snags-
*/
	function lab5(){
		console.log("\nOOPS injts with "+string+" \n");
	   njts = njts-1 ;
	   return;
   }

	function lab4(){
		console.log("\nOOPS : injts more joints than max "+EMAX +"\n");
	   ok = 82 ;
   }
} /* injnts */
/***************************************/

//void inlims(void)
var inlims = function()
/*
     read in limits for a joint.

     called by main,
     calls     nexts, match, value,

*/
{
   var k,m,n;

   nexts();
   n = match(njts,jntlen,jname);
   if (n < 0)
   {
      console.log("limits given for nonexistent joint: "+          string+"\n");
      getout(1);
      if (ok == 1) return;
   }
   for (k = 0; k < 3; ++k)
   {
      for (m = 0; m < 2; ++m)
      {
         nexts();
         lim[n][k][m] = value();
      }
   }
} /* inlims */
/***************************************/

//int inname(int n, int isvar, int lengths[EMAX], char names[EMAX][BMAX])
var inname = function(n, isvar, lengths, names)
/*
     read in a number and then that many names.

     called by inperf,
     calls     nexts, getint, match, addnam,
*/
{
   var e;	//int
   var nitems;	//int
   var nnames;	//int
   var no;	//int

   nnames = n;
/*
     get number of names in list
*/
   nexts();
   nitems = getint();
   if ( ok > 0 )
   {
	   console.log("inname problem- number of names not stated on\n");
	   console.log(""+line+"\n");
	   return(nnames);
   }
/*
     get names in list
*/
   if (nitems <= 0){
	   return(nnames);
	}
   for (  e = 0 ; e < nitems ; ++ e )
   {
      nexts();
      if (length < 1){
	   return(nnames);
	}
      no = match(nnames,lengths,names);
      if (no <= 0)
         nnames = addnam(nnames,names,isvar,lengths);
   }
   return(nnames);
} /* inname */
/***************************************/

//TODO REVIEW: GOTO replaced with return and a function, which looks just as ugly as the original goto.
//void dojoin(void)
var dojoin = function()
/*
   this works out the positions of the centres of each ellipsoid
   'cen' and the joints 'jnt', using the data 'dcon'

   called by  compl,
*/
{
      //int e,ecount,newc,old,newel,oldel,j,k;
	  var e,ecount,newc,old,newel,oldel,j,k;
	  var jfound = new Array(); //int jfound[EMAX];	//ERRORAGE
      var efound = new Array() //int efound[EMAX];
      var elist = new Array() //int elist[EMAX];
/*
     clear found and put all ellipsoids at origin -
*/
      for (  e = 0 ; e < ne ; ++ e )
      {
         cen[e][0] = 0; cen[e][1] = 0; cen[e][2] = 0;
         jfound[e] = FALSE ;
         efound[e] = FALSE ;
      }
      if (njts >= 0)
      {
         ecount=0 ;
         elist[ecount]=0 ;
         efound[ecount] = TRUE;
/*
     run through the ellipsoids of current figure -
*/
	lab2();
	function lab2(){
			for (  e = ecount ; e <= ecount ;  e ++ )
			 {

	/*   run through joints, adding to figure's ellipsoids - */

				for (  j = 0 ; j <= njts ; ++ j )
				{
				   if ((jfound[j] == FALSE)
					   && ((coel[j][0] == elist[e])
						 ||(coel[j][1] == elist[e])))
				   {

	/*   found a joint- */

					  oldel = elist[e] ;
					  if (coel[j][1] == oldel) newc = 0 ;
					  if (coel[j][0] == oldel) newc = 1 ;
					  jfound[j] = TRUE ;
					  old = 1-newc ;
					  newel = coel[j][newc] ;

	/*   check for legality- */

					  for (  k = 0 ; k < ecount ; ++ k )
					  {
						 if (newel == elist[k])
						 {
							console.log("cyclic joint structure - perhaps delete doub1 of the joints \n");
							console.log(" "+                         ecount+" "+ne+" "+njts+" "+e+" "+j+" "+newc+" "+old+" "+k+" "+oldel+" "+newel+"\n");
							return;
						 }
					  }
					  ecount = ecount+1 ;
					  elist[ecount] = newel ;
					  efound[newel] = TRUE;

	/*   locate the new joint and ellipsoid- */
	   
					  for (  k = 0 ; k < 3 ; ++ k )
					  {
						  jnt[j][k] = cen[oldel][k]+dcon[j][old][k] ;
						  cen[newel][k] = jnt[j][k]-dcon[j][newc][k] ;
					  }
					}
				}
			}

	/* locate an ellipsoid in some other figure - */

			for (newel = 0; newel < ne; ++ newel)
			{
				if (efound[newel] == FALSE)
				{
				   ++ ecount;
				   elist[ecount] = newel;
				   efound[newel] = TRUE;
				   lab2();
				}
			}
		}
     }
} /* dojoin */
/*******************************************/

//void checkin(void)
var checkin = function()
/*
     check the specifications of the actions
     to be performed.

     called by  compl,
*/
{
   var newa; //int newa;
   var j; //int j;
   var subfirst; //int subfirst;

   newa = TRUE;
   subfirst = TRUE;
/*
     check for snags-
*/
   for (  j = 0 ; j <= ne ; ++ j )
   {
      if (ellfig[j] < 0)
      {
         console.log("\nOOPS checkin: ellipsoid "+       		  j+" "+ename[j]+" defined but not in a figure\n");
	     ok = 79 ;
      }
      if (ax[j][0]*ax[j][1]*ax[j][2] <= doub0)
      {
         console.log("\nOOPS checkin: ellipsoid "+ 		     j+" "+ename[j]+" not dimensioned\n");
	     ok = 80 ;
      }
   }
   if (fstop < fstart)
   {
      console.log("\nOOPS  checkin: view "+           fstart+" "+fstop +" - produces no frames\n");
      ok = 81;
   }
   for (  j = 1 ; j <= nsubs ; ++ j )
   {
      if ((ok == 0) && (called[j] == FALSE))
      {
	     if (subfirst == TRUE)
	     {
	        subfirst = FALSE;
	        newa = FALSE;
	     }
      }
   }
   if ((nvals+nvars) > EMAX)
   {
      console.log("\nOOPS  checkin "+              nvals+"  non-integer values + "+nvars+" variables\n");
      console.log(" give more than max of "+EMAX+" \n");
		ok = 82;
   }
} /* checkin */
/***************************************/

//TODO REVIEW - GOTOs replaced with return and some logic moved.
//int valadd(double v)
var valadd = function(v)
/*
     if 'v' is not in array 'val', then put it at the end.
     wherever it is, put its index into 'j'.

     called by parset,
*/
{
   var j; //int j ;

   for (  j = 1 ; j <= nvals ; ++ j )
      if (val[j] == v){
		return(j);
	  }

   nvals = nvals+1 ;
   if (nvals > EMAX){
		console.log("\nOOPS in valadd: no. of constants "+           nvals+" > max "+EMAX +"\n");
		ok = 90 ;
		return(j);
	}
   j = nvals ;
   val[j] = v ;

/*
     snag-
*/

} /* valadd */
/***************************************/

//TODO REVIEW (A tonne of GOTOs that bounce processing all over the place.) Replaced with functions, which are no prettier than the gotos.
//int parset(int contrl)
var parset = function(contrl)
/*
     decode the parameters of the jth action using:
             0  none expected
             1  numeric value or variable name
             2  ellipsoid name
             3  joint name
             4  figure name
             5  axis name
             6  subroutine name
             7  variable name
             8  anything
             9  image file name

     called by inperf,
     calls  nexts, getint, value, valadd, match, addnam,
*/
{
	var k;//int
	var nax = 2;//int
	var attach = 8;//int
	var detach = 9;//int
	var v;//double

	k = 0;
	if ( contrl == 0 ) return( k );
	nexts();

	//pick an integer constant-
	if ( ( contrl != 1 ) && ( contrl != 8 ) ) return( lab1() );
	k = getint();
	if ( pok == TRUE ) return( k );

	//pick a double constant-
	pok =  TRUE;
	v = value();
	if ( pok == FALSE ) return( lab1() );

	k = valadd(v);
	if ( pok == FALSE ) return( lab1() );

	k = -k;
	return( k );

	function lab1(){
		//pick an axis-
		if ( ( contrl != 5 ) && ( contrl != 8 ) ) return( lab2() );

		k = match ( nax, axlen, axnam );
		if ( k >= 0 ){
			return( k );
		}else{
			return lab2();
		}
	}

	function lab2(){
		//try for a variable-
		k = match ( nvars, varlen, vname );
		if ( k < 0 ) return( lab3() );
		usevar[k] = 1;
		if ( contrl != 7 ) k = k - EMAX + 1;
		return( k );
	}

	function lab3(){
		//pick an ellipsoid-
		if ( ( contrl != 2 ) && ( contrl != 8 ) ) return( lab4() );
		k = match ( ne, ellen, ename );
		if ( k < 0 ) return( lab4() );
		return( k );
	}

	function lab4(){
		//pick a joint-
		if ( ( contrl != 3 ) && ( contrl != 8 ) ) return( lab5() );
		k = match ( njts, jntlen, jname );
		if ( k >= 0 ) return( k );
		if ( ptype != attach ) return( lab5() );

		njts = addnam ( njts, jname, 0, jntlen );
		k = njts - 1;
		return( k );
	}

	function lab5(){
		//pick a figure-
		if ( ( contrl != 4 ) && ( contrl != 8 ) ) return( lab7() );
		k = match ( nfigs, figlen, fname );
		if ( ptype == detach ) return( lab6() );
		if ( k < 0 ) return( lab7() );
		return( k );
	}

	function lab6(){
		//action detach- accept any figure but "all"-
		if ( k == 0 ) return( lab7() );
		if ( k > 0 ) return( k );
		nfigs = addnam ( nfigs, fname, 0, figlen );
		k = nfigs - 1;
		return( k );
	}

	function lab7(){
		//pick a subroutine call-
		if ( ( contrl != 6 ) && ( contrl != 8 ) ) return( lab8() );
		k = match ( nsubs, sublen, sname );
		if ( k <= 0 )
		{
			nsubs = addnam ( nsubs, sname, 0, sublen );
			k = nsubs - 1;
		}
		called[k] = TRUE;
		return( k );
	}

	function lab8(){
		//pick a file name-
		if ( ( contrl != 9 ) && ( contrl != 8 ) ) return( lab9() );
		k = match ( nfiles, fillen, tname );
		if ( k > 0 ) return( k );
		nfiles = addnam ( nfiles, tname, 0, fillen );
		k = nfiles - 1;
		return( k );
	}

	function lab9(){
		//snag-
		printf ( "\nOOPS parset : contrl %d\n", contrl );
		bell ( 1, 1 );
		ok = 89;

		//lab11:
		return( k );
	}
} /* parset */
/***************************************/

//void inperf(void)
var inperf = function()
/*
     this decodes the input text defining the required actions.

  global  variables
     main     -true if in main program.
     nmax     -max number of commands possible.
     ptype    -type of action being read.

     called by compl,
     calls     nexts, match, addnam, inname, inells, injts,
               value, parset, nlims,

	30 Jul 2006 d072 put infinite loop if error
*/
{
	var nells;	//int
	var how;	//int     // type number of current action keyname.
	var j;	//int
	var k;	//int
	var linel;	//int   // length of string 'line'
	var nmax;	//int	  // max number of commands possible
	var s;	//int	     // counter thru subs and along a string
	var thisub;	//int // number of current subroutine
	var v;	//double

	p = 0;
	thisub = 0;
	subact[0][0] = 0;
	nmax = 4 * EMAX + PMAX;

	// run through statements

	for ( comand = 1; ( ( ok == 0 ) && ( comand <= nmax ) ); ++comand )
	{
		start = -1;
		nexts ();
		linel = parseInt(strlen ( line ));	//TODO REVIEW (Was cast to int)
		for ( s = 0; s < linel; ++s )
			aline[comand][s] = line[s];
		how = match ( NKEYS, keylen, keynam );
		ptype = how;
		if ( ( how <= 0 ) || ( how >= NKEYS ) )
		{
			console.log("\nOOPS inperf: how "+ how+" outside range 0 to "+ NKEYS +"\n");
			bell ( 1, 1 );
			ok = 88;
		}
		else if ( how == stop_keyword_code )
		{
			subact[thisub][1] = p - 1;
			break;
		}
		else if ( how == figure_keyword_code  )
		{
			nexts ();
			nfigs = addnam ( nfigs, fname, FALSE, figlen );
			nells = ne;
			figell[nfigs-1] = ne;
			ne = inname ( nells, 0, ellen, ename );
		}
		else if ( how == ellips_keyword_code  )
		{
			j = inells ();
			ax[j][0] = semiax[0];
			ax[j][1] = semiax[1];
			ax[j][2] = semiax[2];
		}
		else if ( how == joint_keyword_code )
		{
			injts ();
		}
		else if ( how  == limits_keyword_code )
		{
			inlims ();
		}
		else if ( how == variable_keyword_code )
		{
			nvars = inname ( nvars, 1, varlen, vname );
			if ( ( nvars + nvals ) > EMAX )
			{
				console.log("\nOOPS inperf nvars "+ nvars+" + nvals "+ nvals+" > EMAX "+ EMAX +"\n");
				bell ( 1, 1 );
				ok = 87;
			}
		}
		else if ( how == speed_keyword_code )
		{
			nexts ();
			v = value ();
			if ( v < doub0 ) slow = int ( -v + inv2 );
			if ( v > doub0 ) fast = int ( v + inv2 );
		}
		else if ( how == view_keyword_code )
		{
			nexts ();
			v = value ();
			if ( pok == TRUE )
			{
				vstart = int ( v ) - 1;
				if ( vstart < 0 ) vstart = 0;
				nexts ();
				v = value ();
				if ( pok == TRUE ) vstop = int ( v );
			}
			if ( pok == FALSE )
			{
				console.log("\nOOPS inperf: view "+ v +"\n");
				bell ( 1, 1 );
				ok = 70;
			}
		}
		else if ( how == subrou_keyword_code  ) // start a subroutine
		{
			inmain = FALSE;
			if ( thisub == 0 ) subact[0][1] = p - 1;
			nexts ();
			thisub = match ( nsubs, sublen, sname );
			if ( thisub <= 0 )
			{
				nsubs = addnam ( nsubs, sname, 0, sublen );
				if ( ok == 0 ) thisub = nsubs - 1;
			}
			defined[thisub] = TRUE;
			subact[thisub][0] = p;
		}
		else if ( how == endsub_keyword_code )  //  end of a subroutine
		{
			nexts ();
			k = match ( nsubs, sublen, sname );
			if ( k == thisub )
			{
				subact[k][1] = p - 1;
			}
			else
			{
				console.log("\nOOPS inperf: k "+ k+" != thisub "+ thisub +"\n");
				bell ( 1, 1 );
				ok = 86;
			}
		}
		else if ( p >= PMAX )//  read an action -
		{
			p = PMAX - 1;
			console.log("beware- more than "+ PMAX +" action specs\n");
			console.log("actions deleted after line  "+ nline+"\n"+ line +"\n");
		}
		else
		{
			distrn[p] = how;
			cline[p] = comand;

			// read frames to which this action refers -

			frstart[p] = parset(1);
			if ( ok == TRUE ) frstop[p] = parset ( 1 );
			if ( inmain == TRUE )
			{
				if ( frstart[p] < fstart )
					fstart = frstart[p];
				if ( frstop[p] > fstop )
					fstop = frstop[p];
			}
			if ( ok == TRUE ) //call of a subroutine
			{
				if ( how == call_keyword_code )
				{
					distrn[p] = call_keyword_code;
					ptype = call_keyword_code;
					type[p] = call_keyword_code;
					nexts ();
					k = match ( nvars, varlen, vname );
					if ( k < 0 )
					{
						k = match ( nsubs, sublen, sname );
						if ( k < 0 )
						{
							nsubs = addnam ( nsubs, sname, 0, sublen );
							k = nsubs - 1;
						}
						called[k] = TRUE;
						pf[p][0] = k;
					}
					else
					{
						pf[p][0] = k - EMAX + 1;
					}
				}
				else // read action done in these frames
				{
					nexts ();
					ptype = match ( NKEYS, keylen, keynam );
					type[p] = ptype;
					if ( ( ptype < 1 ) || ( ptype >= NKEYS ) )
					{
						console.log("\nOOPS inperf: ptype "+ ptype +"\n");
						bell ( 1, 1 );
						ok = 85;
					}
					else // run through parameters of pth action
					{
						for ( j = 0; ( ( ok <= 0 ) && ( j < 6 ) ); ++j )
						{
							pf[p][j] = parset ( par[ptype][j] );
							if ( ok > 0 )
							{
								console.log("\nOOPS in inperf: problem parameter "+ j+" "+ p+" "+ pf[p][j]+" "+ ptype+" "+ par[ptype][j] +"\n");
								bell ( 1, 1 );
							}
						}
					}
				}
			}
		}
		if ( ok > 0 )
		{
			console.log("\nOOPS in inperf: problem near line "+	nline+"\n "+ line +"\n\n");
			bell ( 1, 1 );
			while(true){}//dead: goto dead; WTF?! TODO REVIEW
		}
		if ( distrn[p] > 0 ) ++p;
		npfs = p;
	}
} /* inperf */
/***************************************/

//void compl()
var compl = function()
/*
   calls    inperf, getout, dojoin, checkin,
   called by main,
*/
{ 
   nline = 0;
   inperf();
   if (ok > 0) getout(1);
   if (ok == 1) return;
   dojoin();
   if (ok > 0) getout(1);
   if (ok == 1) return;
   checkin();
   if (ok > 1) getout(1);
} /* compl */
/***************************************/
/*
          actions39.h

      setels      - finds ellipsoids and joints connected to given ellipsoid
      save        - store positions and orientations
      restore     - restore positions and orientations
      store3      - writes data about given frame to arrays
      getvalu     - gets a value from constants or variables
      vecmat      - multiplies a vector by a matrix
      doground    - moves a set of ellipsoids to rest on y = 0
      setjnt      - finds ellipsoids and joints connected to a given joint
      setfrc      - sets proportion of action for current frame
      doscale     - scale a value by some proportion
      findfg      - finds which figure includes a given ellipsoid
      checkpr     - checks parameters for legality
      setper      - decodes the parameters of the current action
      sqr         - square a value
      docolour    - sets colours of an ellipsoid
      doplace     - sets viewing point (array 'pplace')
      setobs      - sets 3*3 matrix for viewing rotation and place
      enquir      - stores values of centres, joints or axis lengths
      doattach    - joins 2 figures into 1
      dodetach    - breaks 1 figure into 2
      domoveby    - moves a set of ellipsoids relative to refell
      dogroell    - scales axes of an ellipsoid
      dogrofig    - scales a set of ellipsoids in size
      dogrojnt    - scales set of ellipsoids keeping a joint fixed
      domovjnt    - moves a joint
      balanc      - balances part of a figure
      dodrag      - keeps an ellipsoid touching ground
      dcen        - find separation of ellipsoid centres
      newton      - solve a polynomial
      getmat      - generate matrix of an ellipsoid
      getaxes     - find axis lengths of an ellipsoid
      surf        - find separation of ellipsoid surfaces
      sepn        - find distance between 2 ellipsoid surfaces
      fun         - used by 'solve' for abut
      solve       - find zero of 'fun'
      angsepn     - find approx angular separation of ell1 and ell2 from x
      dotouch     - bends a figure to make 2 ellipsoids touch
      trying      - 'domoveby' then 'sepn'
      fndmin      - find minimum of function 'trying'.
      doabut      - slide figure to touch another


***************************************/

//TODO REVIEW: Yet another function full of GOTOs that has to be reviewed. Replaced with return-s, continue-s, a big if-statement and a function.
//void setels(int ellpsd, int jthis)
var setels = function(ellpsd, jthis)
/*
     puts into 'elist' and 'jlist' those ellipsoids and joints
     (if any) connected to 'ellpsd'
     (including 'ellpsd' and 'jthis')
     except those connected through joint 'jthis'

     if 'jthis' is negative, puts all joints and ellipsoids
     connected to 'ellpsd' into lists.
     if 'ellpsd' is zero, puts all joints and ellipsoids into lists,
     except ellipsoid zero (world).

     'ecount' is the number of ellipsoids in the list 'elist'.
     'jcount' is the number of joints in the list 'jlist'.

     called by  setper, findfg, dodetach, dogrojnt, dodrag,
                store3, fun,
*/
{
      var change;	//int
      var ell;	//int
      var e,ee,j,jj ;	//int


      if (ellpsd >= ne)
      {
          ok = 79;
          console.log("\nOOPS setels: ellpsd "+ 			  ellpsd+" "+ename[ellpsd]+" >= ne "+ne+"\n");
          return;
      }
      if (ellpsd <= 0){
		lab6();
		return;
	  }
      ecount = 1;
      elist[0] = ellpsd;
      if (njts <= 0) return;
      jcount = 0;
      if (!(jthis < 0)){
		  if (jthis >= njts)
		  {
			  ok = 78;
			  console.log("\nOOPS setels: jthis "+ 			  jthis+"  "+jname[jthis]+" > njts "+njts+"\n");
			  return;
		  }
		  if ((coel[jthis][0] != ellpsd) && (coel[jthis][1] != ellpsd))
		  {
			  ok = 50;
			  console.log("\nOOPS setels: joint "+jthis+" "+name[jthis]+" connected to "+coel[jthis][0]+" "+ename[coel[jthis][0]]+" and "+coel[jthis][1]+" "+ename[coel[jthis][1]]+", not "+ellpsd+" "+ename[ellpsd]+"\n");
			  return;
		  }
		  jcount = 1;
		  jlist[0] = jthis;
	}
	  change == TRUE;
	  while(change == TRUE){
			change = FALSE;
		  for ( e=0; e < ecount; ++e )
		  {
	/*   seek joint not in jlist connected to ellipsoid elist[e]- */

			 for ( j=0; j < njts; ++j )
			 {
				if ((j == jthis) && (jthis > 0)) continue;
				var b = false;//boolean
				for (  jj=0; jj < jcount; ++jj ){
				   if (j == jlist[jj]){
					b = true;
				   }
				}
				if(b == true){
					continue;
				}
	/*   
		j not in list yet-
	*/
				ell = -1;
				if (coel[j][0] == elist[e]) ell = coel[j][1];
				if (coel[j][1] == elist[e]) ell = coel[j][0];
				if (ell < 0) continue;

	/*   store new joint and ellipsoid- */

				jlist[jcount] = j;
				++jcount;
				change = TRUE;
				for (  ee=0; ee < ecount; ++ee )
				   if (ell == elist[ee]) continue;
				elist[ecount] = ell;
				++ecount;
				change = TRUE;
			 } /* j */
		 } /* e */
	 }
	 return;

/*   set all ellipsoids and joints- */

	lab6();
	function lab6(){
		jcount = 0 ;

	/*   all joints with non-null connections- */

		  jcount = 0;
		  for (  j = 0 ; j <= njts ; ++ j )
		  {
			 if (coel[j][0] >= 0)
			 {
				++jcount;
				jlist[jcount-1] = j ;
			 }
		  }

	/*   all ellipsoids except world- */

		  ecount = ne ;
		  for (  e = 1 ; e <= ne ; ++ e )
			 elist[e-1] = e ;
	}
}  /* setels */
/*************************************************/

/*****************************************************/
/*
         shado41.h

     based on shadoq.c

     to add shadows to figures

     26 Apr  2005  adapt to include in drawel
     16 Jan  2003  remove the shadows below ground
     15 Aug  2001  move the shadows below ground
     18 Aug  1993  to accommodate joints
     23 Oct  1992  D.Herbison-Evans  written

***************************************************

   subroutines-
        setcof
        setaxe
        setpro
        setmat
        setnup
        ground
        doshadow

***********************************************/

//void setcof(double coef[7], double el[3][3] )
var setcof = function(coef, el )
/* 
     set up coeffs of outline ellipse of an ellipsoid about 
     its own centre in the form -  
   
     coef(1)*x**2 + coef[2]*z**2 + coef[3]*x*z 
         + coef(4)*x + coef[5]*z + coef[6] = 0 
   
     called by setnup,
*/
{
   var den ;//double

   if (el[1][1] == doub0){
		console.log("setcof "+               el[1][0]+" "+el[1][1]+" "+el[1][2]+"\n");
		ok = 99 ;
		return;
   }
   den = doub1/el[1][1] ;
   coef[1] = el[0][0] - el[0][1]*el[0][1]*den ;
   coef[2] = el[2][2] - el[1][2]*el[1][2]*den ;
   coef[3] = doub2*(el[0][2] - el[0][1]*el[1][2]*den) ;
   coef[4] = doub0 ;
   coef[5] = doub0 ;
   coef[6] =  -doub1 ;
   return;
/* 
     snags -  
*/
} /* setcof */
/************************************************/

//void setaxe(int n, double axe[3], double coef[7])
var setaxe = function(n, axe, coef)
/* 
     find semiminor axis, axe[0], and semimajor 
     axis, axe[2], of ellipse described by coef. 
   
     called by setnup, 
*/
{
   var discrt,lamx,lamz,c12,rtdis ;//double 

   lamx = doub1 ;
   lamz = doub1 ;
   discrt = (coef[1] - coef[2])*(coef[1] - coef[2])+ coef[3]*coef[3];
   if (discrt < doub0){
		lab1();
		return;
   };
   c12 = inv2*(coef[1]+coef[2]) ;
   rtdis = inv2*sqrt(discrt) ;
   lamx = c12 + rtdis ;
   lamz = c12 - rtdis ;
   if (lamx <= doub0){
		lab1();
		return;
   };
   if (lamz <= doub0){
		lab1();
		return;
   };
   axe[0] = doub1/sqrt(lamx) ;
   axe[2] = doub1/sqrt(lamz) ;
   return;
/* 
     snags -  
*/
	function lab1(){
		console.log("setaxe snag "+          lamx+" "+lamz+" "+discrt+"\n");
	   ok = 98 ;
   }
} /* setaxe */
/******************************************/

//double setpro(double coef[7])
var setpro = function(coef)
/*
     for the outline of nth ellipsoid, find 'phi'
     angle between axx axis and scene x axis.
   
     called by setnup, 
*/
{
   var phi ;//double 

   phi = pi-inv2*atan2(coef[3], coef[1]-coef[2]) ;
   if (phi < doub0) phi = phi+twopi ;
   return(phi);
} /* setpro */
/******************************************/

//void setmat ( int n, double el[3][3], double el1[3][3], double unel1[3][3] )
var setmat = function( n, e1, el1, unel1 )
/*

	finds the matrix "el" of the quadratic form of the "n"th
	ellipsoid by taking the diagonal matrix of inverse square
	semiaxes, and doing on it a similarity transform
	for its own rotation.

	called by setnup, cutting,
	calls matmul, rotget,

	12 Aug 2006  returning el1 and unel1
*/
{
	var ii, j;//int
	//double el0[3][3],el2[3][3],el3[3][3];
	//double r[3][3], unr[3][3];
	var e10 = get2DArray(3);
	var e12 = get2DArray(3);
	var e13 = get2DArray(3);
	var r = get2DArray(3);
	var unr = get2DArray(3);

	// initialise diagonal matrix -

	for ( ii = 0; ii < 3; ++ ii )
	{
		for ( j = 0; j < 3; ++ j )
		{
			el0[ii][j] = doub0;
			el3[ii][j] = doub0;
		}
		if ( ax[n][ii] ==  doub0 )
		{
			console.log("setmat  ax["+ n+"]["+ ii +"] = 0\n");
			ok = 97;
			return;
		}
		el0[ii][ii] = doub1 / ax[n][ii];
		el3[ii][ii] = ax[n][ii];
	}
	rotget ( r, unr, n );

	// do similarity transform -

	matmul ( el0, unr, el1 );
	matmul ( r, el0, el2 );
	matmul ( el2, el1, el);
	matmul ( r, el3, unel1);
} /* setmat */
/**********************************************/

//double setnup(int n, double axe[3])
var setnup = function(n, axe)
/* 
     set up parameters of nth ellipsoid relative 
     to own centre. 
   
     called by shadow, 
     calls     setmat, setcof, setaxe, setpro, 
*/
{
   //double el[3][3],el0[3][3],el1[3][3];
   //double con[7];
   var phi;//double
   var e1 = get2DArray(3);
   var e10 = get2DArray(3);
   var e11 = get2DArray(3);
   var con = new Array();

   phi = doub0;
   setmat(n,el,el0,el1) ;
   if ( ok > 0 ){
		lab1();
		return(phi);
	}
   setcof(con,el) ;
   if ( ok > 0 ){
		lab1();
		return(phi);
	}
   setaxe(n,axe,con) ;
   if ( ok > 0 ){
		lab1();
		return(phi);
	}
   phi = setpro(con) ;
   if ( ok > 0 ){
		lab1();
		return(phi);
	}
	return(phi);
/* 
     snag -  
*/
	function lab1(){
		ok = 96;
		console.log("setnup snag in ellipsoid "+n +"\n");
	}
   return(phi);
} /* setnup */
/******************************************/

//double elground(int i)
var elground = function(i)
/*
   find distance of lowest point above the ground
   of the ellipsoid 'i'.

   called by  shadow,
   calls      rotget,
*/
{
   var x,y,z;//double
   var val;//double
   var sumsq;//double
   var sqt;//double
   //double r[3][3],unr[3][3] ;
   var r = get2DArray(3);
   var unr = get2DArray(3);

   val = cen[i][1];

/*   find lowest point- */

   rotget(r,unr,i) ;
   x = unr[0][1]*ax[i][0] ;
   y = unr[1][1]*ax[i][1] ;
   z = unr[2][1]*ax[i][2] ;
   sumsq = x*x+y*y+z*z;
   if (sumsq > doub0)
      sqt = sqrt(sumsq); else sqt = doub0;
   val = cen[i][1] - sqt ;
   return(val);
}  /* elground */
/**********************************************************/

//void doshadow()
var doshadow = function()
/* 
  find the shadow ellipsoids of each ellsoid in the scene

     called by store3,
     calls     setnup, rset, rotput, ground,
*/
{
   var k,n;//int
   var y, phi;//double
   //double r[3][3];
   var r = get2DArray(3);
   //double axe[3];
   var axe = Array();
/* 
     run thru ellipsoids to shadow each in turn -  
*/
    k = ne;
    for (  n = 1 ; n < ne ; ++n )
    {
         phi = setnup(n,axe);
         y = elground(n);
         if (y > doub0)
         {
            cen[k][0] = cen[n][0];
            cen[k][1] = -inv5;
            cen[k][2] = cen[n][2];
            ax[k][0] = axe[0];
            ax[k][1] = inv5;
            ax[k][2] = axe[2];
            rset(r,phi,1);
            rotput(r,k);
            col[k][0] = doub1;
            col[k][1] = doub1;
            col[k][2] = doub1;
			++k;
         } /* y > 0 */
    } /* end n loop */
    ne = k;
} /* doshadow */
/******************************************/

//void save(void)
var save = function()
/*
   save positions and orientations

   called by  store3, doabut, dodrag, dotouch,
*/
{
   var j,n;//int

   nesave = ne;
   for (n = 0; n <= ne; ++n)
   {
      for ( j = 0; j < 3; ++j)
      {
         censav[n][j] = cen[n][j];
         jntsav[n][j] = jnt[n][j];
      }
      for ( j = 0; j < 5; ++j)
         quasav[n][j] = quat[n][j];
   }
} /* save */
/***********************************************/

//void restore(void)
var restore = function()
/*
   restore positions and orientations

   called by  store3, doabut, try, dodrag, fun, dotouch,
*/
{
   var j,n;//int

   ne = nesave;
   for (n = 0; n <= ne; ++n)
   {
      for (j = 0; j < 3; ++j)
      {
         cen[n][j] = censav[n][j];
         jnt[n][j] = jntsav[n][j];
      }
      for ( j = 0; j < 5; ++j)
         quat[n][j] = quasav[n][j];
   }
} /* restore */
/***********************************************/

//void store3(int f)
var store3 = function(f)
/*
    store axes, centres, orientations and colours
    of  nels ellipsoids starting at 1 (avoiding  0 = world),

    called by doframes,
    calls     save, doshadow, setels, shift, twirl, rotput,
	          mkang, storeang, mkquat, rotgrt, restore,
*/
{
   var e,j;//int
   //double invobs[3][3];
   invobs = get2DArray(3);

   save();
   if (shadow == TRUE) doshadow();
   setels(0,-1);
   twirl(pplace[0],pplace[1],pplace[2],obs);
   shift(-pplace[0],-pplace[1],-pplace[2]);
   nels[f] = ne;
   for (e = 0; e < ne; ++e)
   {
      for (j = 0; j < 3; ++j)
      {
         qu3[f][e][j+1] = quat[e][j];
         ax3[f][e][j] = ax[e][j]*inv1000;
         ce3[f][e][j] = cen[e][j]*inv1000;
         co3[f][e][j] = col[e][j]*inv256;
      } /* j*/
      qu3[f][e][3] = -qu3[f][e][3];
      qu3[f][e][0] = degree*atan2(quat[e][3], quat[e][4]);
      ce3[f][e][2] = doub1 - ce3[f][e][2];
      if (col[e][0] < 0)
            sprintf(tn3[f][e],"%s",tname[int(inv2-col[e][0])]);
   } /* e */
   rotput(obs,ne);
   mkang(ne);
   storeang(f,ne,ang[0],ang[1],ang[2]);
   mkquat(ne,ang[0],ang[1],ang[2]);
   rotget(obs,invobs,ne);
   restore();
}  /* store3 */
/***********************************************/

//int getvalu(int p)
var getvalu = function(p)
/*
  get value possibly from array val and put it into v and k.
  if p is negative, get value of variable val(abs(p)),
  if p is positive get p directly.

  called by  doperfrm, setper,
*/
{
   var k;//int
   var ref ;//int

   ref = 0 ;
/*
  is the parameter a variable or direct reference
*/
   if (p < 0)
   {
/*
  parameter is index into array val-
*/
      ref = -p ;
      if ((ref < 0) || (ref >= EMAX))
      {
         ok = 15 ;
         console.log("val index "+ 		     ref+" outside range 0 - "+EMAX+"\n");
      }
      else
      {
         v = val[ref] ;
         k = int(v + inv2) ;
         if (v < doub0) k = int(v -inv2) ;
      }
   }
   else
/*
  parameter is direct reference, use it-
*/
   {
      k = p ;
      v = k ;
   }
   return(k);
}  /* getvalu */
/**************************************/

//void vecmat(double v[3], double m[3][3], double w[3])
var vecmat = function(v, m, w)
/*
   multiply vector 'v' by matrix 'm',
   putting result in 'w'.

   called by  sepn, dobalanc,
*/
{
      var i,j;//int
      //double vv[3]
	  var vv = Arrray();
	  var x;//double

      for (  i = 0 ; i < 3 ; ++ i )
      {
         x = doub0 ;
         for (  j = 0 ; j < 3 ; ++ j )
         {
            x = x+m[i][j]*v[j];
         }
         vv[i] = x ;
      }
      for (  i = 0 ; i < 3 ; ++ i )
      {
         w[i] = vv[i];
      }
}  /* vecmat */
/**********************************************************/

//double elow(int i)
var elow = function(i)
/*
   find height of lowest point of ellipsoid i
   
   called by doground,
   call rotget,
  */
{
   //double r[3][3],unr[3][3];
   var r = get2DArray(3);
   var unr = get2DArray(3);
   var x,y,z;	//double
   var sq,sqt;	//double
   var toty;	//double

   rotget(r,unr,i);
   x = unr[0][1]*ax[i][0];
   y = unr[1][1]*ax[i][1];
   z = unr[2][1]*ax[i][2];
   sq = x*x+y*y+z*z;
   if (sq > doub0) 
      sqt = sqrt(sq); 
   else 
      sqt = doub0;
   toty = cen[i][1] - sqt;
   return(toty);
} /* elow */
/******************************************************/

//double doground(void)
var doground = function()
/*
   find distance of lowest point above the ground
   of the ellipsoids contained in 'elist'.

   called by  action, dodrag, fun, doshadow,
   calls      rotget,
*/
{
   var n ;//int
   var toty;//double
   var val;//double

   if ((ecount < 1) || (ecount > ne))
   {
      ok = 38 ;
      console.log("\nOOPS doground: ecount "+                         ecount+" out of range\n");
   }
   else
   {
      val = cen[elist[0]][1];
/*  run through affected ellipsoids finding lowest point- */
      for (  n = 0 ; n < ecount ; ++ n )
      {
         toty = elow(elist[n]);
         if (toty < val) val = toty;
      } /* n */
   } /* if ecount */
   return(val);
}  /* doground */
/**********************************************************/

//void setjnt(int ellpsd, int jthis)
var setjnt = function(ellpsd, jthis)
/*
     puts into 'elist' and 'jlist' those ellipsoids and joints
     (if any) connected to 'ellpsd'
     (including 'ellpsd' and excluding 'jthis')
     except those connected through joint 'jthis'

     'ecount' is the number of ellipsoids in the list elist.
     'jcount' is the number of joints in the list jlist.

     called by  domovjnt,
*/
{
      var done;//int
      var e,j,i,jt ;//int

      ecount = 1 ;
      elist[0] = ellpsd ;
      jcount = 0 ;
      if ((coel[jthis][0] != ellpsd) && (coel[jthis][1] != ellpsd))
      {
          ok = 64;
          console.log("\nOOPS setjnts: coel "+coel[jthis][0]+"  "+coel[jthis][1]+" out of range "+ellpsd+" "+jthis+"\n");
          return;
      }

	done = FALSE;
	while((done == FALSE))
	{

		for (  e = 0 ; e < ecount ; ++ e )
		{
			done = TRUE;

	/*   seek joint not in jlist but connected to ellipsoid elist[e]- */
		
			for (  jt = 0 ; jt < njts ; ++ jt )
			{
				if (jt == jthis) continue;
				if (jcount > 0)
				{
					for (  j = 0 ; j < jcount ; ++ j )
					{
						if (jt == jlist[j]) continue;
					}
				}
	/*
	   jt not in list yet-
	*/
				i =  -1 ;
				if (coel[jt][0] == elist[e]) i = 1 ;
				if (coel[jt][1] == elist[e]) i = 0 ;
				if (i < 0) continue;

	/*   store new joint and ellipsoid- */

				jlist[jcount] = jt ;
				jcount ++;
				elist[ecount] = coel[jt][i] ;
				ecount ++;
				done = FALSE;
			}
        }
    }
}  /* setjnt */
/*************************************************/
//END PORT ON 2013-12-10 (Errorage)

//PORT ON 2013-12-11 (Errorage)

//void setfrc(int frame, int start, int stp)
var setfrc = function(frame, start, stp)
/*
  set up prop and frac -  proportion of current action time
  to be done for current frame.

  variables -
     a - number of frames done of current action
     aa - number of increments of current action done
     b - number of frames to be done
     bb - number of increments of current action to be done
     even  - 0 if n is even, 1 if n is odd
     nn - total number of increments in current action

  called by  prfrm,
*/
{
   var a,aa,b,bb,h,n,nn,even;//double
   var distrp ;//int

   distrp = distrn[pp];
   n = stp-start ;
   nn = n*(n+doub1)*inv2 ;
   a = frame-start ;
   aa = a*(a+doub1)*inv2 ;
   b = stp-frame+1 ;
   bb = b*(b+doub1)*inv2 ;
   even = parseInt(n)%2 ;	//TODO REVIEW (n was cast to int, parseInt() may round it different)
   h = parseInt(n)/2 ;	//TODO REVIEW (n was cast to int, parseInt() may round it different)
/*
  repeat command-
*/
   if (distrp == repeat_keyword_code)
   {
      prop = doub1 ;
      frac = doub1 ;
   }
   else
/*
  linear command-
*/
   if (distrp == linear_keyword_code)
   {
      if ((b  == doub0) || (n == doub0))
      {
         ok = 10 ;
         console.log("\nOOPS setfrc: linear b "+b+",  n "+n+", start "+start+"\n");
      }
      else
      {
         frac = doub1/n;
         prop = doub1/b;
      }
   }
   else
/*
  acceleration command-
*/
   if (distrp == accele_keyword_code)
   {
      if ((nn == doub0) || ((nn-aa+a) == doub0))
      {
         ok = 11 ;
         console.log("\nOOPS setfrc: accel bb "+bb+",  nn "+nn+", start "+start+"\n");
      }
      else
      {
         frac = a/nn ;
         prop = a/(nn-aa+a) ;
      }
   }
   else
/*
  deceleration command-
*/
   if (distrp == decele_keyword_code)
   {
      if ((nn == doub0) || (bb == doub0))
      {
         ok = 12 ;
         console.log("\nOOPS setfrc: decele bb "+bb+",  nn "+nn+", start "+start+"\n");
      }
      else
      {
         frac = b/nn ;
         prop = b/bb ;
      }
   }
   else
/*
  quadratic command-
*/
   if (distrp == quadra_keyword_code)
   {
      nn = (n*(n+doub2)+even)/doub4 ;
      if (nn == doub0)
      {
         ok = 13 ;
         console.log("\nOOPS setfrc: quadra nn "+nn+", start "+start+"\n");
      }
      else
      {
         frac = a/nn;
         if (a >= b) frac = b/nn ;
      }
      if (bb == doub0)
      {
         ok = 13 ;
         console.log("\nOOPS setfrc: quadra bb "+bb+", start "+start+"\n");
      }
      prop = b/bb ;
      if (a < b)
      {
         if ((nn-aa+a) == doub0)
         {
            ok = 13 ;
         console.log("\nOOPS setfrc: quadra nn-aa+a "+(nn-aa+a,start)+", start %d\n");
         }
         prop = a/(nn-aa+a) ;
      }
   }
   else
/*
   cubic command-
*/
   if (distrp == cubic_keyword_code)
   {
      nn = h*(h+doub1)*(doub2*h+doub1)*inv3 + even*(h+doub1)*(h+doub1);
      aa = a*(a+doub1)*(doub2*a+doub1)*inv6;
      bb = b*(b+doub1)*(doub2*b+doub1)*inv6;
      if ( a < b )
      {
         if (nn > doub0) frac = a*a/nn; else frac = doub1;
         if ((nn-aa+a*a) > doub0) prop = a*a/(nn-aa+a*a); else prop = doub1;
      }
      else
      {
         if (nn > doub0) frac = b*b/nn; else frac = doub1;
         if (bb > doub0) prop = b*b/bb; else prop = doub1;
      }
   }
} /* setfrc */
/************************************************************/

//double doscale(double x)
var doscale = function(x)
/*
   scale x by proportion frac

   called by setper,
*/
{
   var v1,v2,v3,v4,v5;//double

   if (x == doub0)
   {
      ok = 51;
      console.log("\nOOPS scale: scale factor 0\n");
      v5 = doub0;
   }
   else
   {
      if (x > doub0) v1= x; else v1 = -x;
      if ( v1 > doub0) v2 = log(v1); else v2 = doub0;
      v3 = frac*v2;
      v4 = exp(v3);
      if (x > doub0) v5 = v4; else v5 = -v4;
   }
   return(v5);
}  /* doscale */
/***************************************/

//int findfg(int ell)
var findfg = function(ell)
/*
  find the figure (excluding 'every')
  that includes the ellipsoid 'ell'.

  called by setper, doattach, dodetach,
  calls     setels,
*/
{
   var e,f;//int

   setels(ell,-1) ;
   for (  f = 1 ; f <= nfigs ; ++ f )
   {
      for (  e = 0 ; e < ecount ; ++ e )
         if (figell[f] == elist[e]) return(f);
   }
/*
  snag-
*/
   ok = 35 ;
   console.log("\nOOPS findfg: ecount "+ecount+", ell "+ell+" "+ename[ell]+"\n");
   return(-1);
}  /* findfg */
/************************************************/

//void checkpr(void)
var checkpr = function()
/*
  check parameters for legality

  called by  setper,
*/
{
	if ((njts > 0)&&((join < 0) || (join > njts))){
		lab12();
		return
	}
	if ((axis < 0) || (axis > 2)){
		ab13();
		return;
	}
	if ((ellpsd < 0) || (ellpsd > ne)){
		ab14();
		return;
	}
	if ((refell < 0) || (refell > ne)){
		ab15();
		return;
	}
	if ((ell1 < 0) || (ell1 > ne)){
		ab20();
		return;
	}
	if ((ell2 < 0) || (ell2 > ne)){
		ab21();
		return;
	}
	if ((fig < 0) || (fig > nfigs)){
		ab16();
		return;
	}
	if (nvars <= 0) return;
	if ((var0 < 0) || (var0 >= EMAX)){
		ab17();
		return;
	}
	if ((var1 < 0) || (var1 >= EMAX)){
		ab18();
		return;
	}
	if ((var2 < 0) || (var2 >= EMAX)){
		ab19();
		return;
	}
	if (newcol[0] < -nfiles){
		ab45();
		return;
	}
	if (newcol[0] <= 0)
	{
		if (newcol[1] < 0){
			ab48();
			return;
		}
		if (newcol[2] < 0){
			ab49();
			return;
		}
	}
	return;
/*
  data snag-
*/
	function lab12(){
		lab12: ok = 16 ;
		console.log("\nOOPS checkpr: joint "+join+" out of range 0 to  "+njts+"\n");
		return;
	}

	function lab13(){
		lab13: ok = 17 ;
		console.log("\nOOPS checkpr: axis "+axis+" out of range 0 to 2\n");
		return;
	}

	function lab14(){
		lab14: ok = 18 ;
		console.log("\nOOPS checkpr: ellpsd "+ellpsd+" out of range 0 to "+ne+"\n");
		return;
	}
	
	function lab15(){
		lab15: ok = 19 ;
		console.log("\nOOPS checkpr: refell "+refell+" out of range 0 to "+ne+"\n");
		return;
	}

	function lab16(){
		lab16: ok = 20 ;
		console.log("\nOOPS checkpr: fig "+fig+" out of range 0 to "+nfigs+"\n");
		return;
	}

	function lab17(){
		lab17: ok = 21 ;
		console.log("\nOOPS checkpr: var0 "+var0+" out of range 0 to "+nvars+"\n");
		return;
	}
	
	function lab18(){
		lab18: ok = 22 ;
		console.log("\nOOPS checkpr: var1 "+var1+" out of range 0 to "+nvars+"\n");
		return;
	}

	function lab19(){
		ok = 23 ;
		console.log("\nOOPS checkpr: var2 "+var2+" out of range 0 to "+nvars+"\n");
		return;
	}

	function lab20(){
		ok = 36 ;
		console.log("\nOOPS checkpr: ell1 "+ell1+" out of range 0 to "+ne+"\n");
		return;
	}

	function lab21(){
		ok = 37 ;
		console.log("\nOOPS checkpr: ell2 "+ell2+" out of range 0 to "+ne+"\n");
		return;
	}

	function lab45(){
		ok = 45 ;
		console.log("\nOOPS checkpr: newcol[0] "+newcol[0]+" out of range\n");
		return;
	}

	function lab48(){
		ok = 48 ;
		console.log("\nOOPS checkpr: newcol[1] "+                   newcol[1]+" out of range\n");
		return;
	}

	function lab49(){
		ok = 49 ;
		console.log("\nOOPS checkpr: newcol[2] "+  newcol[2]+" out of range\n");
	}


}  /* checkpr */
/***************************************/

//void setper ( int keyword_code )
var setper = function(keyword_code )
/*
  decodes the parameters of the 'pp'th action using -

     0 - none
     1 - x
     2 - y
     3 - z
     4 - ang1
     5 - ang2
     6 - ang3
     7 - x scaling factor
     8 - y scaling factor
     9 - z scaling factor
    10 - value for a variable
    11,12,13 - red, green and blue colour coords
    14 - the debug parameter
    15 - reference to a file name in list 'fname'
    21 - axis
    22 - joint
    23 - reference ellipsoid
    24 - ellpsd (moving or central ellipsoid)
    25 - fig  (figure)
    27,28,29 - var0,var1,var2 (references to variables)
    30 - ell1 (ellipsoid to touch)
    31 - ell2 (ellipsoid to be touched)

  called by doperfrm,
  calls     getvalu, checkpr, findfg, setels,
*/
{
	var c;//int
	var j;//int

	if ( ( keyword_code < 1 ) || ( keyword_code > NKEYS ) )
	{
		ok = 14;
		console.log("\nOOPS setper: keyword_code "+keyword_code+"\n");
		bell ( 1, 1 );
		return;
	}

	for (  j = 0; j < 6; ++ j )//run thru parameters of 'p'th action
	{
		c = code[keyword_code][j];
		if ( c != 0 )
		{
			k = getvalu ( pf[pp][j] );
			if ( ok != 0 ) return;

			//			set real parameters -

			if ( c ==  1 ) xx[0] = v;
			if ( c ==  2 ) xx[1] = v;
			if ( c ==  3 ) xx[2] = v;
			if ( c ==  4 ) ang[0] = v * radian;
			if ( c ==  5 ) ang[1] = v * radian;
			if ( c ==  6 ) ang[2] = v * radian;
			if ( keyword_code != growto_keyword_code )
			{
				if ( c ==  7 ) factor[0] = doscale ( v );
				if ( c ==  8 ) factor[1] = doscale ( v );
				if ( c ==  9 ) factor[2] = doscale ( v );
			}
			if ( keyword_code == growto_keyword_code )
			{
				if ( c ==  7 ) factor[0] = v;
				if ( c ==  8 ) factor[1] = v;
				if ( c ==  9 ) factor[2] = v;
			}
			if ( c == 10 ) varval = v;

			//			set int parameters-

			if ( c == 11 ) newcol[0] = k;
			if ( c == 12 ) newcol[1] = k;
			if ( c == 13 ) newcol[2] = k;
			if ( c == 21 ) axis = k;
			if ( c == 22 ) join = k;
			if ( c == 23 ) refell = k;
			if ( c == 24 ) ellpsd = k;
			if ( c == 25 ) fig = k;
			if ( c == 27 ) var0 = EMAX - k - 1;
			if ( c == 28 ) var1 = EMAX - k - 1;
			if ( c == 29 ) var2 = EMAX - k - 1;
			if ( c == 30 )
			{
				ell1 = k;
				ellpsd = k;
			}
			if ( c == 31 ) ell2 = k;
		} /* c != 0 */
	}	/* j */

	//	check for errors-

	checkpr ( );
	if ( ok != 0 ) return;

	/*
	if appropriate, set up lists of ellipsoids and joints
	in affected figures ( NB figell[every] is -1)  -
	*/
	if ( keyword_code == drag_keyword_code )
	{
		if ( ellpsd == coel[join][0] )
			refell = coel[join][1];
		else
			refell = coel[join][0];
	}

	if ( code[keyword_code][0] == linear_keyword_code )
	{
		setels ( figell[fig], -1 );
		if ( ok != 0 ) return;
	}

	if ( ( keyword_code == rotate_keyword_code )
		|| ( keyword_code == abduct_keyword_code )
		|| ( keyword_code == drag_keyword_code )
		|| ( keyword_code == abut_keyword_code ) )
	{
		fig = findfg ( ellpsd );
		if ( ok != 0 ) return;
	}

	if ( keyword_code == abut_keyword_code )
	{
		setels ( ellpsd, -1 );
		if ( ok != 0 ) return;
	}
	else if ( ( keyword_code == balanc_keyword_code )
		|| ( keyword_code == touch_keyword_code )
		|| ( keyword_code == bendby_keyword_code )
		|| ( keyword_code == bendto_keyword_code )
		|| ( keyword_code == flex_keyword_code )
		|| ( keyword_code == rotate_keyword_code )
		|| ( keyword_code == abduct_keyword_code )
		|| ( keyword_code == drag_keyword_code )
		|| ( keyword_code == linkx_keyword_code ) )
	{
		if ( keyword_code != linkx_keyword_code ) setels ( ellpsd, join );
		if ( ok != 0 ) return;
		xx[0] = jnt[join][0];
		xx[1] = jnt[join][1];
		xx[2] = jnt[join][2];
	}

	if ( ( keyword_code == grofig_keyword_code )
		|| ( keyword_code == spinto_keyword_code )
		|| ( keyword_code == spinby_keyword_code )
		|| ( keyword_code == center_keyword_code )
		|| ( keyword_code == centre_keyword_code ) )
	{
		xx[0] = cen[ellpsd][0];
		xx[1] = cen[ellpsd][1];
		xx[2] = cen[ellpsd][2];
	}
	else if ( keyword_code == moveto_keyword_code )
	{
		xx[0] = xx[0] - cen[ellpsd][0];
		xx[1] = xx[1] - cen[ellpsd][1];
		xx[2] = xx[2] - cen[ellpsd][2];
	}

	if ( ( keyword_code == growto_keyword_code ) )
	{
		xx[0] = cen[ellpsd][0];
		xx[1] = cen[ellpsd][1];
		xx[2] = cen[ellpsd][2];
		return;
	}

	if ( ( keyword_code == opacty_keyword_code ) )
	{
		xx[0] = cen[ellpsd][0];
		xx[1] = cen[ellpsd][1];
		xx[2] = cen[ellpsd][2];
		return;
	}
	if ( ( keyword_code == lghtng_keyword_code ) )
	{
		xx[0] = cen[ellpsd][0];
		xx[1] = cen[ellpsd][1];
		xx[2] = cen[ellpsd][2];
		return;
	}
   for (j = 0; j < ne; ++j)
	{
		minax[j] = ax[j][0];
		if (ax[j][1] < minax[j]) minax[j] = ax[j][1];
		if (ax[j][2] < minax[j]) minax[j] = ax[j][2];
		maxax[j] = ax[j][0];
		if (ax[j][1] > maxax[j]) maxax[j] = ax[j][1];
		if (ax[j][2] > maxax[j]) maxax[j] = ax[j][2];
	}
}  /* setper */
/****************************************************/

//double sqr(double x)
var sqr = function(x)
/*

  called by   surf, angsepn, dotouch,

  16 Sep 2006 called by totouch
*/
{
   return(x*x);
}  /* sqr */
/*******************************************************/

//void docolour(double prop)
var docolour = function(prop)
/*
  sets ellipsoid colour proportionately to the appropriate rgb.

  called by  action,
*/
{
   col[ellpsd][0] += prop*(newcol[0]-col[ellpsd][0]) ;
   col[ellpsd][1] += prop*(newcol[1]-col[ellpsd][1]) ;
   col[ellpsd][2] += prop*(newcol[2]-col[ellpsd][2]) ;
   if (t == 52) col[ellpsd][0] = -newcol[0] ;
}  /* docolour */
/***************************************/

//void doplace(void)
var doplace = function()
/*
  set observers viewing point to new values.

  called by  action,
*/
{
   pplace[0] += prop*(xx[0]-pplace[0]);
   pplace[1] += prop*(xx[1]-pplace[1]);
   pplace[2] += prop*(xx[2]-pplace[2]);
}  /* doplace */
/***************************************/

//void setobs(void)
var setobs = function()
/*
  set up matrix 'obs' for eulerian angles in 'ang'.

  called by action,
  calls     rset, matmul,
*/
{
   //double newang[3];
   //double r1[3][3],r2[3][3],r3[3][3] ;
   
   newang = Array();
   r1 = get2DArray(3);
   r2 = get2DArray(3);
   r3 = get2DArray(3);

   newang[0] = oldang[0] + prop*(ang[0]-oldang[0]);
   newang[1] = oldang[1] + prop*(ang[1]-oldang[1]);
   newang[2] = oldang[2] + prop*(ang[2]-oldang[2]);
   rset(r1,newang[0],0) ;
   rset(r2,newang[1],1) ;
   rset(r3,newang[2],2) ;
   matmul(r1,r2,obs) ;
   matmul(obs,r3,obs) ;
   oldang[0] = newang[0];
   oldang[1] = newang[1];
   oldang[2] = newang[2];
}  /* setobs */
/************************************************************/

//void enquir(int thisp, double array[EMAX][3])
var enquir = function(thisp, array)
/*
  store values from 'array' into variables.

  called by  action,
*/
{
   val[var0] = array[thisp][0] ;
   val[var1] = array[thisp][1] ;
   val[var2] = array[thisp][2] ;
}  /* enquir */
/************************************************************/

//void doattach(void)
var doattach = function()
/*
  create a joint 'join' at point 'x,y,z'
  relative to centre of refell.

  called by action,
  calls     findfg,
*/
{
   var fig1,fig2,low,high ;//int

   if ((coel[join][1] != -1) || (coel[join][0] != -1))
   {
      ok = 42 ;
	  console.log("doattach: join "+join,coel[join][0]+", coel[join][0] %d,  coel[join][1] "+coel[join][1]+"\n");
   }
   else
   {
/*
  find lowest ellipsoid of the figures being connected-
*/
      fig1 = findfg(ellpsd);
      fig2 = findfg(refell);
      if (ok == 0)
      {
         high = fig2 ;
         low = fig1 ;
         if (figell[low] > figell[high])
         {
            low = fig2 ;
            high = fig1 ;
         }
         figell[high] = figell[low] ;
         coel[join][0] = ellpsd ;
         coel[join][1] = refell ;
         jnt[join][0] = xx[0]+cen[ellpsd][0] ;
         jnt[join][1] = xx[1]+cen[ellpsd][1] ;
         jnt[join][2] = xx[2]+cen[ellpsd][2] ;
      }
   }
}  /* doattach */
/*************************************/

//void dodetach(void)
var dodetach = function()
/*
  split 1 figure into 2.

  called by action,
  calls     findfg, setels,
*/
{
	var othere,otherf;//int
	var j,k;//int
	var fgk ;//int
	//int fg[2],el[2] ;
	var fg = Array();
	var el = Array();
/*
  check if the new figure 'fig' is already being used -
*/
	j = figell[fig] ;
	if (j != 0)
	{
		ok = 43;

		console.log("\nOOPS dodetach:  fig "+fig+" "+fname[fig]+",  figell[fig] "+figall[fig]+" "+ename[figell[fig]]+"\n");
		return;
	}
/*
   fig ok, so start detaching-
*/
	othere = 0 ;
	if (coel[join][0] == ellpsd) othere = coel[join][1] ;
	if (coel[join][1] == ellpsd) othere = coel[join][0] ;
	if (othere == 0)
	{
		ok = 44 ;
		console.log("\nOOPS dodetach:  join "+join+" "+jname[join]+",  ellpsd "+ellpsd+" "+ename[ellpsd]+", othere "+othere+" "+ename[othere]+"\n");
		return;
	}
	otherf = findfg(othere);
	if (ok != 0) return;
/*
   move all the joints down one -
*/
	for (j = join; j < njts; ++j)
	{
		coel[j][0] = coel[j+1][0] ;
		coel[j][1] = coel[j+1][1] ;
		jnt[j][0] = jnt[j+1][0];
		jnt[j][1] = jnt[j+1][1];
		jnt[j][2] = jnt[j+1][2];
	}
	--njts;
/*
   find representative ellipsoid of each figure -
*/
	el[0] = ellpsd ;
	el[1] = othere ;
	fg[0] = fig ;
	fg[1] = otherf ;
	for (  k = 0 ; k < 2 ; ++ k )
	{
		  setels(el[k],-1);
		  fgk = fg[k] ;
		  figell[fgk] = EMAX ;
		for (  j = 0 ; j < ecount ; ++ j ){
			if (elist[j] < figell[fgk]) figell[fgk] = elist[j] ;
		}
	}
return;
}  /* dodetach */
/*******************************/

//void domoveby( double x, double y, double z, int refell)
var domoveby = function(x, y, z, refell)
/*
  moves ellipsoids and joints in lists 'elist' and 'jlist'
  by vector 'x,y,z' relative to the axes of 'refell'.

  called by  try, action, abut, fun,
  calls      rotget, vecmul, shift,
*/
{
   //double v[EMAX][3];
   //double r[3][3];
   //double unr[3][3] ;
   var v = get2DArray(EMAX);
   var r = get2DArray(3);
   var unr = get2DArray(3);

   v[0][0] = x ;
   v[0][1] = y ;
   v[0][2] = z ;
   rotget(r,unr,refell) ;
   vecmul(v,r,0) ;
   shift(v[0][0],v[0][1],v[0][2]);
}  /* domoveby */
/************************************************************/

//void dogroell( double f[3], int j, double a[EMAX][3])
var dogroell = function(f, j, a)
/*
  scales jth member of array by factor 'f'.

  called by  action, dogrofig,
*/
{
   a[j][0] *= f[0] ;
   a[j][1] *= f[1] ;
   a[j][2] *= f[2] ;
}  /* dogroell */
/***************************************/

//void dogrofig( double x0, double x1, double x2)
var dogrofig = function(x0, x1, x2)
/*
  scale parts in 'elist' and 'jlist' about the point 'x0,x1,x2',
  multiplying all semi-axes and coords of centres and joints
  by factor[0],factor[1],factor[2] in x,y, and z directions 
  respectively

  called by  action, dogrojnt,
  calls      shift, dogroell,
*/
{
	var e,j,n ;//int

	shift(-x0,-x1,-x2);
	for (  n = 0 ; n < ecount ; ++ n )
	{
		e = elist[n];
		dogroell(factor,e,cen) ;
		dogroell(factor,e,ax) ;
		maxax[e] = ax[e][0];
		minax[e] = ax[e][0];
		for (j = 1; j < 3; ++j)
		{
			if (ax[e][j] > maxax[e]) maxax[e] = ax[e][j];
			if (ax[e][j] < minax[e]) minax[e] = ax[e][j];
		}
	}
	for (  n = 0 ; n < jcount ; ++ n )
		dogroell(factor,jlist[n],jnt) ;
	shift(x0,x1,x2);
}  /* dogrofig */
/********************************************/

//void dogrojnt(void)
var dogrojnt = function()
/*
  scales ellipsoid 'ellpsd' by factors 'f', keeping position
  of joint 'join' fixed and moving all other joints and
  ellipsoids connected to 'ellpsd' appropriately.

  called by  action,
  calls      rotget, twirl, dogrofig, setels, vecmul, shift,
*/
{
   //double d[EMAX][3],r[3][3],unr[3][3] ;
	var d = get2DArray(EMAX);
	var r = get2DArray(3);
	var unr = get2DArray(3);
	var x0,x1,x2;//double
	var jscan,othere,dim ;//int

	if ((njts <= 0) || (njts > EMAX)){
		ok = 40;
		console.log("\nOOPS dogrojnt:  njts "+njts+"\n");
		return;
	}
	if ((coel[join][0] != ellpsd)&&(coel[join][1] != ellpsd)){
		ok = 41;
		console.log("\nOOPS dogrojnt:  ellpsd "+ellpsd+"\n");
		return;
	}

	rotget(r,unr,ellpsd) ;
/*
  scale and shift the growing ellipsoid-
*/
	x0 = jnt[join][0] ;
	x1 = jnt[join][1] ;
	x2 = jnt[join][2] ;
	elist[0] = ellpsd ;
	ecount = 1;
	jcount = 0;
	twirl(x0,x1,x2,unr) ;
	dogrofig(x0,x1,x2) ;
	twirl(x0,x1,x2,r) ;
/*
  now shift everything connected to ellpsd-
*/
	for ( jscan = 0 ; jscan < njts ; ++ jscan )
	{
		if (jscan != join)
		{
			othere = 0 ;
			if (coel[jscan][0] == ellpsd) othere = coel[jscan][1] ;
			if (coel[jscan][1] == ellpsd) othere = coel[jscan][0] ;
			if (othere != 0)
			{
/*
  find parts connected to jscan through othere-
*/
				setels(othere,jscan);
				if (ok != 0) return;
/*
  find out how much to shift things hanging here-
*/
				for ( dim = 0 ; dim < 3 ; ++ dim )
					d[0][dim] = jnt[jscan][dim]-jnt[join][dim] ;
				vecmul(d,unr,0) ;
				for ( dim = 0 ; dim < 3 ; ++ dim )
					d[0][dim] = d[0][dim]*(factor[dim]-doub1) ;
				vecmul(d,r,0) ;
				shift(d[0][0],d[0][1],d[0][2]);
			} /* othere != 0 */
		} /* jscan != join */
	} /* jscan */
	return;
}  /* dogrojnt */
/************************************************************/

//void domovjnt(void)
var domovjnt = function()
/*
  moves joint 'join' by amounts 'x', keeping the position
  of ellipsoid 'ellpsd' fixed and moving all other joints
  and ellipsoids connected to 'join' appropriately.

  called by  action,
  calls      rotget, setjnt, vecmul, shift,
*/
{
   //double d[EMAX][3],r[3][3],unr[3][3] ;
	var d = get2DArray(EMAX);
	var r = get2DArray(3);
	var unr = get2DArray(3);
	var othere,dim ;//int

	if ((njts <= 0) || (njts > EMAX)){
		lab11();
		return;
	}
	othere = 0 ;
	if (coel[join][0] == ellpsd) othere = coel[join][1] ;
	if (coel[join][1] == ellpsd) othere = coel[join][0] ;
	if (othere == 0){
		lab11();
		return;
	}

	rotget(r,unr,ellpsd) ;
/*
  shift the joint -
*/
	for (dim = 0; dim <3; ++dim)
		d[0][dim] = frac*xx[dim];
	setjnt(ellpsd,join);
	vecmul(d,unr,0);
	shift(d[0][0],d[0][1],d[0][2]);
	return;
/*
  snags-
*/
	function lab11(){
		ok = 63 ;
		console.log("\nOOPS domovjnt:  join "+join+",  ellpsd "+ellpsd+"\n");
	}
}  /* domovjnt */
/************************************************************/

//void dobalanc(void)
var dobalanc = function()
/*
  balance part of a figure by bending ellipsoid 'ellpsd' at
  the joint at point 'x' about 'axis' of ellipsoid 'refell'.

  called by  action,
  calls      rotget, rset, matmul, vecmat, dospinby,
*/
{
	var wsum,wmod,uw,vw,alpha,beta,phi,mass;//double
	//double dx[3],u[3],u1[3],u2[3],v[3],w[3],w1[3],ww[3];
	var dx = Array();
	var u = Array();
	var u1 = Array();
	var u2 = Array();
	var v = Array();
	var w = Array();
	var w1 = Array();
	var ww = Array();
	//double ro[3][3],unro[3][3],rph[3][3],ralph[3][3],rb[3][3];
	var ro = get2DArray(3);
	var unro = get2DArray(3);
	var rph = get2DArray(3);
	var ralph = get2DArray(3);
	var rb = get2DArray(3);
	var usq;//double
	var j,k,thise ;//int

	rotget(ro,unro,refell) ;
/*
  form unit vector along rotation axis
*/
	for (  k = 0 ; k < 3 ; ++ k )
	{
		ww[k] = doub0;
		u[k] = ro[k][axis] ;
	}
/*
  run through moving ellipsoids-
*/
	for (  j = 0 ; j < ecount ; ++ j )
	{
		thise = elist[j] ;
		mass = ax[thise][0]*ax[thise][1]*ax[thise][2] ;
/*
  find vector to jth moving ellipsoid centre-
*/
		for (  k = 0 ; k < 3 ; ++ k )
		{
			dx[k] = cen[thise][k]-xx[k];
			ww[k] += mass*dx[k];
		}
	}
/*
  find vector 'w' to centre of mass of moving parts -
*/
	wsum = doub0;
	for (k = 0; k < 3; ++k)
		wsum += ww[k]*ww[k];
	if (wsum > doub0) wmod = sqrt(wsum); else wmod = doub0;
	for (k = 0; k < 3; ++k)
		if (wmod > doub0) w[k] = ww[k]/wmod; else w[k] = doub0;
/*
   find vector 'v' at point on meridien through u
   equal in distance from u as w is from u -
*/
	uw = doub0;
	for (k = 0; k < 3; ++k)
		uw += u[k]*w[k];
	usq = u[0]*u[0] + u[2]*u[2];
	if (usq > doub0) u1[0] = sqrt(usq); else u1[0] = doub0;
	u1[1] = u[1];
	u1[2] = doub0;
	if ((uw < -doub1) || (uw > doub1)) 
		alpha = doub0; 
	else 
		alpha = acos(uw);
	rset(ralph,alpha,2);
	vecmat(u1,ralph,u2);
	phi = -atan2(u[2],u[0]);
	rset(rph,phi,1);
	vecmat(u2,rph,v);
	vw = v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
	if ((vw > -doub1) && (vw < doub1)) 
		beta = acos(vw); 
	else 
		beta = doub0;
	rset(rb,beta,axis);
	matmul(rb,unro,rb);
	matmul(ro,rb,rb);
	vecmat(w,rb,w1);
	if (w[1] > w1[1]) beta = -beta;
	dospinby(xx,refell,beta,axis) ;
}  /* dobalanc */
/************************************************************/

//void dodrag(void)
var dodrag = function()
/*
   bend 'ellpsd' at joint 'join' with coordinates 'x'
   about 'axis' of 'refell' to make 'ell1' touch the ground

   called by  action,
   calls      dospinby, save, restore, shift, doground,
              setels, doangles,
*/
{
	var jb,jc,jd;//int
	var fixde;//int
	var quadrant,qa,qb,qd;//int
	var gap;//double
	var yj;//double
	var proptemp;//double
	var xa,xb,xd;//double
	var dx;//double
	var y,ya,yb,yd;//double

	proptemp = prop;
	prop = doub1;
/*
	set rest of figure section, excluding ellpsd-ell1, 
	to touch the ground -
*/
	fixde = 0;
	if ( coel[join][0] == ellpsd ) fixde = coel[join][1];
	if ( coel[join][1] == ellpsd ) fixde = coel[join][0];
	if ( fixde == 0)
	{
		ok = 65;
		console.log("\nOOPS dodrag: "+ename[ellpsd]+" not connected to joint "+jname[join]+"\n");
		console.log("which joins  "+ename[coel[join][0]]+" and "+ename[coel[join][1]]+"\n");
		bell ( 1, 1 );
		prop = proptemp;
	} /* joint does not connect specified ellipsoid */
	else
	{
		setels(fixde,join);
		gap = doground();
		setels(fixde,-1);
		shift ( doub0, -gap, doub0 );
		yj = jnt[join][1] + minax[ell1];
		if (yj <= minax[ellpsd])
			console.log("dodrag snag : joint too low\n");
		else
		if ( yj > minax[ellpsd] )
		{
/*
	which quadrant is centre of ell1 in wrt join -
*/
			quadrant = 0;
			if ( jnt[join][0] > cen[ell1][0] ) quadrant += 1;
			if ( jnt[join][2] > cen[ell1][2] ) quadrant += 2;
			save ( );
			xx[0] = jnt[join][0];
			xx[1] = jnt[join][1];
			xx[2] = jnt[join][2];
			setels(ellpsd,join);
			y = doground ( );
//printf("dodrag  %d %f %s   %f %f %f\n",
//f,y,ename[ell1],cen[ell1][0],cen[ell1][1],cen[ell1][2]);
//printf("            %s   %f %f %f\n",
//ename[fixde],cen[fixde][0],cen[fixde][1],cen[fixde][2]);
//
			setels(ellpsd,join);
			xa = doub0;
			qa = quadrant;
			ya = y;
//
			xb = doub2;
			qb = -1;
			for (jb = 0; ((jb < 8)&&(qb != quadrant)); ++jb)
			{
				xb *= inv2;
				restore(); 
				dospinby ( xx, refell, xb, axis );
				yb = doground ( );
				qb = 0;
				if ( jnt[join][0] > cen[ell1][0] ) qb += 1;
				if ( jnt[join][2] > cen[ell1][2] ) qb += 2;
//printf("dodragb %d %f %f %d\n",jb,xb,yb,qb);
			}
			xd = doub0;
			if (ya*yb < doub0)
			{
				for (jc = 1; jc < 9; ++jc)
				{
					dx = (xb - xa)*inv2;
					if ((ya*yb < doub0)||(ya*ya < yb*yb))
						xb = xb - dx;
					else 
					{
						xa = xb;
						ya = yb;
						xb = xb + dx;
					}
					restore();
					dospinby ( xx, refell, xb, axis );
					yb = doground ( );
					qb = 0;
					if ( jnt[join][0] > cen[ell1][0] ) qb += 1;
					if ( jnt[join][2] > cen[ell1][2] ) qb += 2;
//printf("dodragc %d %f %f %d\n",jc,xb,yb,qb);
				} /* jc loop */
//printf("dodragd %d %f %s   %f %f %f\n",
//f,xb,ename[ell1],cen[ell1][0],cen[ell1][1],cen[ell1][2]);
				qd = -1;
				xd = xb*doub2;
				for (jd = 0; ((jd < 8)&&(qd != quadrant)); ++jd)
				{
					xd *= inv2;
					restore();
					dospinby ( xx, refell, xd, axis );
					yd = doground ( );
					qd = 0;
					if ( jnt[join][0] > cen[ell1][0] ) qd += 1;
					if ( jnt[join][2] > cen[ell1][2] ) qd += 2;
//printf("dodrage %d %f %f %d\n",jd,xd,yd,qd);
				} /* jd */
			} /* ya*yb < 0 */
			restore();
			prop = proptemp;
//printf("dodragf %d %f %s   %f %f %f\n",
//f,xd,ename[ell1],cen[ell1][0],cen[ell1][1],cen[ell1][2]);
			dospinby ( xx, refell, prop*xd, axis );
//printf("dodragg %d %f %s  %f %f %f\n\n",
//f,prop*xd,ename[ell1],cen[ell1][0],cen[ell1][1],cen[ell1][2]);
//printf("            %s   %f %f %f\n\n",
//ename[fixde],cen[fixde][0],cen[fixde][1],cen[fixde][2]);
			prdone = TRUE;
		} /* joint OK */
	}  /* connections  OK */
} /* dodrag */
/************************************************************/

//double dcen(double d[3])
var dcen = function(d)
/*
  find separation between ellipsoid centres

  called by  doabut, sepn,
*/
{
   var j;//int
   var ans,dsq;//double

   dsq = doub0 ;
   for (j = 0; j < 3; ++ j)
   {
      d[j] = cen[ell1][j]-cen[ell2][j];
      dsq += d[j]*d[j];
   }
   if (dsq > doub0) ans = sqrt(dsq); else ans = doub0;
   return(ans);
}  /* dcen */
/************************************************************/

//double newton(int n, double start, double a[])
var newton = function(n, start, a)
/*
   solve the polynomial of degree (n-1):

           n-1           n-2
   a[n-1]*x    + a[n-2]*x    + ... a[1]*x + a[0] = 0

   using 48 Newton-Raphson iterations starting at 'start'.

   called by  surf,
*/
{
	var x,xold,num,den;//double
	var j,k;//int

	x = start;
	xold = doub2*start+doub1;
	num = doub1;
	for (j = 0;
		((j < 48) && (num != doub0) && (x != xold)); ++j)
	{
		num = doub0;
		den = doub0;
		for (k = n-1; k >= 0; --k)
		{
			num = x*num + a[k];
			if (k > 0) den = x*den + a[k]*double(k);
		}
		if (den == doub0)
		{
			ok = 57;
			//TODO ERRORAGE Correct (mmultiple overloads more than %x marks)
			printf("\nOOPS newton: ell2 %d, den %f\n",            
			join,ell1,coel[join][0],coel[join][1]);
			x = -doub1;
			xold = -doub1;
		}
		else
		{
			xold = x;
			if (den != doub0) x -= num/den;
		}
   }
   return(x);
}  /* newton */
/*******************************************************/

//void getmat(double mat[3][3], double tr[3][3], double untr[3][3],double trtr[3][3], double untrtr[3][3], int ell)
var getmat = function(mat, tr, untr, trtr, untrtr, ell)
/*
   for ellipsoid 'ell', find its matrix 'mat',
   also: its transformation matrix 'tr' and
   its transpose 'trtr', and its inverse 'untr'
   and its inverse transpose 'untrtr'.

   called by  surf,
   calls      rotget, matmul,
*/
{
	//double r[3][3],unr[3][3];
	//double diag[3][3];
	//double undiag[3][3];
	var r = get2DArray(3);
	var unr = get2DArray(3);
	var diag = get2DArray(3);
	var undiag = get2DArray(3);
	var j,k;//int
/*
   make the matrix -
*/
	rotget(r,unr,ell);
	for (j = 0; j < 3; ++j)
	{
		for (k = 0; k < 3; ++k)
		{
			undiag[j][k] = doub0;
			diag[j][k] = doub0;
		}
		if (ax[ell][j] == doub0)
		{
			ok = 60;
			console.log("\nOOPS getmat: ell "+ell+", j "+j,+",  ax[ell][j] "+ax[ell][j]+"n");
			return;
		}
		if (ax[ell][j] > doub0) 
			undiag[j][j] = doub1/ax[ell][j];
		else 
			undiag[j][j] = doub0;
		diag[j][j] = ax[ell][j];
	}
	matmul(undiag,unr,trtr);
	matmul(r,undiag,untrtr);
	matmul(untrtr,trtr,mat);
/*
  make the transformation matrices-
*/
	matmul(diag,unr,untr);
	matmul(r,diag,tr);
} /* getmat */
/***********************************************/

//void getaxes(double m[3][3], double axes[3], double r[3][3])
var getaxes = function(m, axes, r)
/*
   from its matrix 'm', find the 3 semiaxes of an ellipsoid
   and corresponding rotation matrix 'r'.

   called by  surf,
   calls      matmul,
*/
{
	var lambda,mu,nu,nusq,numu;//double
	//double s[3][3],t[3][3];
	var s = get2DArray(3);
	var t = get2DArray(3);
	//double st[3][3],tmp[3][3];
	var st = get2DArray(3);
	var tmp = get2DArray(3);
	var sn,cs,sig;//double
	var a,b,c,abc;//double
	var p,q;
	var j,k,n;

	for (j = 0; j < 3; ++j)
	{
		for (k = 0; k < 3; ++k)
		{
			r[j][k] = doub0;
			t[j][k] = m[j][k];
		}
		r[j][j] = doub1;
	}
/*
   iterate 3 times -
*/
	for (n = 0; n < 3; ++n)
	{
/*
   find largest off-diagonal element -
*/
		a =  t[0][1]; if (a < doub0) a = -a;
		b =  t[0][2]; if (b < doub0) b = -b;
		c =  t[1][2]; if (c < doub0) c = -c;
		abc = a; p = 0; q = 1;
		if (b > abc)
		{
			abc = b;
			q = 2;
		}
		if (c > abc)
		{
			abc = c;
			p = 1;
			q = 2;
		}
		if (abc != doub0)
		{
			lambda = -t[p][q];
			mu = inv2*(t[p][p] - t[q][q]);
			nusq = lambda*lambda + mu*mu;
			if (nusq > doub0) 
				nu = sqrt(nusq); 
			else 
				nu = doub0;
			if (mu > doub0) 
				sig = doub1; 
			else 
				sig = -doub1;
			numu = nu + mu*sig;
			if ((nu > doub0) && (numu > doub0)) 
				cs = sqrt((numu)/(doub2*nu));
			else 
				cs = doub0;
			if ((cs > doub0) && (nu > doub0)) 
				sn = sig*lambda/(cs*doub2*nu);
			else 
				sn = doub0;
			for (j = 0; j < 3; ++j)
			{
				for (k = 0; k < 3; ++k)
				{
					s[j][k] = doub0;
				}
				s[j][j] = doub1;
			}
			s[p][p] = cs;
			s[q][q] = cs;
			s[p][q] = sn;
			s[q][p] = -sn;
			for (j = 0; j < 3; ++j)
			{
				for (k = 0; k < 3; ++k)
				{
					st[j][k] = s[k][j];
				}
			}
			matmul(st,t,tmp);
			matmul(tmp,s,t);
			matmul(st,r,r);
		}
	}
	for (j = 0; j < 3; ++j)
	{
		if (t[j][j] == doub0)
		{
			ok = 60;
			console.log("\nOOPS getaxes: ell1 "+ell1+",  ell2 "+ell2+",  j "+j+",  t[j][j] "+t[j][j]+"\n");
		}
		else
			if (t[j][j] > doub0)
				axes[j] = doub1/sqrt(t[j][j]);
			else
				axes[j] = doub0;
	}
}  /* getaxes */
/******************************************************/

//static void initsphere(void)
var initsphere = function()
/*
    set up a unit sphere centred on the origin

    called by  main, checkeys,
*/
{
	var i,j,k,m,n;//int
	var nlat,nlong;//int
	var longi,lat;//double
	var dlat,dlong;//double
	var xval,yval,zval;//double
	var ssmall;//double

	if (nsph >= SMAX)
	{
		console.log("initsphere oops "+nsph+" > max, reset to "+SMAX+"\n");
		nsph = SMAX-2;
	}
	if (nsph < SMIN)
	{
		console.log("initsphere oops "+nsph+" < min, reset to "+SMIN+"\n");
		nsph = SMIN;
	}
	dlat = pi/nsph;	//TODO REVIEW (nsph was cast to double, might cause integer division)
	dlong = dlat;
	nlong = nsph+nsph;
	nlat = nsph;
	lat = -piby2;
/*
   set up vertices -
*/
	for (i = 0; i <= nlat; i++)
	{
		ssmall = cosf(lat);
		yval = sinf(lat);
		longi = 0;
		for ( j = 0; j <= nlong; ++ j)
		{
			xval = ssmall*cosf(longi);
			point[i][j][0] = xval;
			zval = ssmall*sinf(longi);
			point[i][j][1] = yval;
			point[i][j][2] = zval;
			longi += dlong;
		}
		lat += dlat;
	}
/*
   set up faces and their normals -
*/
	m = 0;
	for (i = 0; i < nlat; ++i)
	{
		for (j = 0; j < nlong; ++j)
		{
			for (k = 0; k < 3; ++k)
			{
				sph[m][0][k] = point[i][j+1][k];
				sph[m][1][k] = point[i+1][j+1][k];
				sph[m][2][k] = point[i+1][j][k];
				sph[m][3][k] = point[i][j][k];
				norm[m][k] = doub0;
				for (n = 0; n < 4; ++n)
					norm[m][k] += inv4*sph[m][n][k];
			}
			++m;
			if (m >= 4*SMAX*SMAX)
			{
				console.log("initsphere faces ("+nsph+")  "+m+" > "+(4*SMAX*SMAX)+"\n");
				getout(1);
				if (ok == 1) return;
			}
		}
	}
	nfaces = m;
}   /* initsphere */
/*************************************************/

//double surf(int ell1, int ell2)
var surf = function(ell1, ell2)
/*
  return 0 if ell1 touches ell2
  return -ve if they intersect
  return +ve if they do not intersect

    5 May 2007 skip if touch found
   15 Sep 2006 return 0 if they touch
   10 Sep 2006 use polyhedral vertices of ell1 
               instead of ellipsoid
   10 Sep 2006 surf() uses 2 parameters
   17 Feb 1993 find distance between surfaces of ell1 and ell2
               using Buckdales algorithm (giving answer in
               transformed coordinates).
    1 Oct 1981 using polyhedral vertices of ell1 and ell2
               written: Don Herbison-Evans

   called by  dotouch, cutting,
   calls      sqr, setmat, vecmul, getout, initsphere
*/
{
	var j;//int
	var d12, dax, dmin;//double
	//double c1[3],c2[3];
	var c1 = Array();
	var c2 = Array();

	dax = minax[ell1] + minax[ell2];
	dmin = doub0;
//
//   check separation of centres -
//
	for (j = 0; j < 3; ++j)
	{
		c1[j] = cen[ell1][j];
		c2[j] = cen[ell2][j];
		d12 = c1[j]-c2[j];
		dmin += d12*d12;
	}
	dmin = sqrt(dmin);
	return( dmin - dax);
} /* surf */
/***********************************/

//void cutting ( void )
var cutting = function()
/*
  find if any ellipsoids intersect that shouldnt.

  called by  doframes,
  calls setmat, vecmul, getout, initsphere, surf

   15 Sep 2006 use surf() = 0 if ellipsoids touch
   10 Sep 2006 moved innards to surf()
   12 Aug 2006 fixed initsphere and setmat
    8 Aug 2006 Don Herbison-Evans
*/
{
	var j,k;//int
	var s,smin,smax;//int
	var ncut,ntemp;//int
	var el1,el2;//int
	var dmin;//double
	//double m[3][3],m1[3][3],unm1[3][3];
	var m = get2DArray(3);
	var m1 = get2DArray(3);
	var unm1 = get2DArray(3);
	//double e1c[EMAX][3];
	var e1c = get2DArray(EMAX);
	//double c1[3];
	var c1 = Array();
	var key;//char

	ncut = 5;
	ntemp = nsph;
	nsph = ncut;
	if (2*nsph*nsph > EMAX)
		nsph = parseInt(sqrt((EMAX-1)*inv2));//TODO REVIEW (EMAX-1 was cast to double, may cause integer division; Entire thing was cast to int, not sure if parseInt() works the same way)
	printf("cuttinga %d %d\n",ntemp,ncut);
	initsphere ( );
	for (el1 = 1; el1 < (ne-1); ++el1)
	{
		setmat ( el1, m, m1, unm1 );
		for (j = 0; j < 3; ++j)
			c1[j] = cen[el1][j];
		s = 0;
		for (j = 0; j <= nsph; ++j)
		{
			for ( k = 0;  k <= (nsph+nsph); ++k)
			{
				e1c[s][0] = point[j][k][0];
				e1c[s][1] = point[j][k][1];
				e1c[s][2] = point[j][k][2];
				vecmul ( e1c, unm1, s);
				e1c[s][0] += c1[0];
				e1c[s][1] += c1[1];
				e1c[s][2] += c1[2];
				++s;
				if (s >= EMAX)
				{
					console.log("snag in 'cutting'\n");
					console.log(" no. vertices "+(2*nsph*nsph)+", max "+EMAX+"\n");
					ok = 66;
					getout(ok);
				}
			} /* k loop */
		} /* j loop */
		smax = s;
		smin = -1;
		for (el2 = el1+1; el2 < ne; ++el2)
		{
			if (forbid[el1][el2] == true)
			{
				dmin = surf(el1, el2);
				if (dmin < doub0)
				{
					console.log("frame "+f+": ellipsoid "+ename[el1]+" cuts ellpsoid "+ename[el2]+"\n");
					console.log("\n carry on regardless? Hit 'enter' for Yes, 'no' then 'enter' for No: ");
					key = getchar();
					if ( key == '\n' )
					{
						for (j = 0; j < ne; ++j)
							for (k = 0; k < ne; ++k)
								forbid[j][k] = false;
					}
					else
					{
						pause = true;
						fstop = f-1;
					}
				} // dmin
			} /* ell2 forbidden */
		} /* el2 */
	} /* el1 */
	nsph = ntemp;
	console.log("cuttingb "+el1+" "+el2+"\n");
	initsphere ( );
} /* cutting */
/******************************************************/

//double sepn(void)
var sepn = function()
/*
   find distance between surfaces of ell1 and ell2.

   called by  dotouch, fun, fndmin, doabut, try,
   calls      dcen, surf,
*/
{
   var dmid;//double
   var ans,minsep;//double
   //double d[3];
   var d = Array();
/*
   find bounds on separation -
*/
	dmid = dcen(d);
	minsep = dmid-minax[ell1]-minax[ell2];
	if (minsep <= doub0)
	{
		ans = -doub1;
	}
	else
	{
		ans = surf(ell1,ell2);
		if (ans < doub0) ans = -doub1;
	}
	return(ans);
} /* sepn */
/**********************************************/

//double fun(double xarg)
var fun = function(xarg)
/*
   called by  solve,
   calls      setels, dospinby, doground, restore, sepn,
              domoveby,
*/
{
   var ans;//double
   //double dx[3];
   var dx = Array();

	ans = doub0;
	if (t == drag_keyword_code)
	{
		setels(ellpsd,join);
		restore();
		dospinby(xx,refell,xarg,axis);
		ecount = 1; elist[0] = ell1;
		ans = doground();
		return(ans);
	}
	else
	if (t == abut_keyword_code)
	{
		dx[0]=doub0;
		dx[1]=doub0;
		dx[2]=doub0;
		dx[axis] = xarg;
		domoveby(dx[0],dx[1],dx[2],refell);
		ans = sepn();
		if (ok != 0) 
			return(doub0);
		restore();
		return(ans);
	}
	else
	{
		ok = 59;
		return(doub0);
	}
} /* fun */
/*******************************************/

//double solve(double a, double b, int n)
var solve = function(a, b, n)
/*
   find a zero of 'fun()' in the range 'a' to 'b'
   to an accuracy of 1.0 on fun, using at most 'n' iterations.

   called by  doabut,
   calls      fun,
*/
{
	var valab,vala,valb;//double
	var angab,anga,angb;//double
	var dval;//double
	var k;//int

	angab = a;
	anga = a; 
	vala = fun(a);
	if ((vala > -doub1) && (vala < doub1)) return(angab);
	if (ok != 0) return(angab);
	angab = b;
	angb = b; 
	valb = fun(b);
	if (ok != 0) return(angab);
	if ((valb > -doub1) && (valb < doub1)) return(angab);
	if (vala*valb > doub0)
	{
		if (vala > doub0)
		{
			if (vala < valb) angab = a;
				else angab = b;
		}
		else
		{
			if (vala < valb) angab = b;
				else angab = a;
		}
		anga = doub0; angb = doub0;
		valab = doub0; vala = doub0; valb = doub0;
		return(angab);
	}
	if (vala > valb)
	{
		valab = vala;
		vala = valb;
		valb = valab;
		angab = anga;
		anga = angb;
		angb = angab;
	}
	for (k = 0; k < n; ++k)
	{
		dval = vala-valb;
		if (dval < doub0) dval = -dval;
		if (dval < doub1) return(angab);
		angab = inv2*(anga+angb);
		valab = fun(angab);
		if (ok != 0) return(angab);
		if (valab < doub0)
		{
			anga = angab;
			vala = valab;
		}
		else
		{
			angb = angab;
			valb = valab;
		}
	}
} /* solve */
/*******************************************/

//double angsepn ( double xx[3], int ell1, int ell2 )
var angsepn = function(xx, ell1, ell2 )
/*
   find approx angular separation in radians at xx
   of ell1 and ell2 using minax 

   called by  dotouch
   calls sqr

	 5 May 2007  zero if overlapping well
   16 Sep 2006  xx,ell1,ell2 parameters instead of global
*/
{
	var dmin,dsep;//double
	var asep;//double
	var dist1, dist2;//double

	dmin = minax[ell1] + minax[ell2];
	dsep = sqrt (sqr ( cen[ell1][0] - cen[ell2][0] )
		+ sqr ( cen[ell1][1] - cen[ell2][1] )
		+ sqr ( cen[ell1][2] - cen[ell2][2] ) );
	if (dsep <= dmin)
		asep = doub0;
	else
	{
		dist1 = sqrt( sqr ( xx[0] - cen[ell1][0] )
			 + sqr ( xx[1] - cen[ell1][1] )
			 + sqr ( xx[2] - cen[ell1][2] ) );

		dist2 = sqrt( sqr ( xx[0] - cen[ell2][0] )
			 + sqr ( xx[1] - cen[ell2][1] )
			 + sqr ( xx[2] - cen[ell2][2] ) );

		asep = (dsep + dsep) / (dist1 + dist2);
	}
	return( asep );
} /* angsepn */
/************************************************/

//void dotouch(void)
var dotouch = function()
/*
  make ellipsoid 'ell1' come as close as possible to 'ell2'
  by bending ellipsoid 'ellpsd' at the joint at point 'x'
  about 'axis' of ellipsoid 'refell'.

  called by  action,
  calls      angsepn, surf, save, restore, dospinby, sqr,

  	5 May 2007  skip if overlapping well already
  16 Sep 2006  find min of square of surf() penetration
  10 Sep 2006  surf() uses 2 parameters
   1 Oct 1981  written Don Herbison-Evans
*/
{
	var angj, angk;//double
	var arange;//double
	var dang;//double
	var gap, g1, gmin;//double
	var pro;//double
	var samples = 10;//int
	var iterations = 5;//int
	var j;//int
	var jend;//int
	var jmin;//int
	var jstart;//int
	var k;//int

	save ();
	pro = prop;
	prop = doub1;
	arange = angsepn ( xx, ell1, ell2 );
	jstart = -samples / 2;
	jend = jstart + samples;
	angk = doub0;
	for ( k = 0; k < iterations; ++k )
	{
		restore ();
		dang = arange / samples;//TODO REVIEW samples was cast to double
		angj = angk + dang * jstart;//TODO REVIEW jstart was cast to double
		dospinby ( xx, refell, angj, axis );
		gap = sqr ( surf ( ell1, ell2 ) );
		g1 = gap;
		gmin = gap;
		jmin = jstart;
		/*
		seek minimum of gap,
		*/
		for ( j = jstart + 1; ( j < jend + 1 ); ++j )
		{
			restore ();
			angj = angk + dang * j;//TODO REVIEW j was cast to double
			dospinby ( xx, refell, angj, axis );
			gap = sqr ( surf ( ell1, ell2 ) );
			if ( gap < gmin )
			{
				jmin = j;
				gmin = gap;
			}
		} /* j */
		if ( ( jmin != jstart ) && ( jmin != jend ) )
			arange = dang + dang;
		angk += dang*jmin;//TODO REVIEW jmin was cast to double
	} /* k */
	prop = pro;
	restore ();
	angj = angk + dang * jmin;//TODO REVIEW jmin was cast to double
	dospinby ( xx, refell, prop*angj, axis );
} /* dotouch */
/************************************************************/

//double trying(double a)
var trying = function(a)
/*
     function to be found a minimum of,
     called from doabut.

     called by  fndmin,
     calls      domoveby, sepn, restore,
*/
{
	//double dx[3]
	var dx = array();
	var s;//double

	if (t == abut_keyword_code)
	{
		dx[0] = doub0; dx[1] = doub0; dx[2] = doub0;
		dx[axis] = a;
		domoveby(dx[0],dx[1],dx[2],refell);
		s = sepn();
		restore();
	}
	else
	{
		ok = 62;
		s = doub0;
		console.log("\nOOPS trying: doabut not calling it!\n");
	}
	return(s);
} /* trying */
/***************************************************/

//double fndmin(double a, double b, int n)
var fndmin = function(a, b, n)
/*
   finds the minimum value of 'try'
   in the range 'a' to 'b' using at most 'n' iterations.

   called by  dotouch,
   calls      try,
*/
{
	var trya,tryb,tryab;//double
	var olda,oldb;//double
	var mina,minb;//double
	var k;//int

	tryab = doub0;
	olda = a; oldb = b;
	for (k = 0; k < n; ++k)
	{
		trya = olda + (oldb-olda)*inv3;
		mina = trying(trya);
		if (ok != 0) return(tryab);
		tryb = oldb - (oldb-olda)*inv3;
		minb = trying(tryb);
		if (ok != 0) return(tryab);
		if (mina < minb) oldb = tryb;
			else olda = trya;
	}
	if (mina < minb) tryab = trya;
		else tryab = tryb;
   
} /* fndmin */
/*******************************************/

//TODO REVIEW (GOTOs replaced with function)
//void doabut(void)
var doabut = function()
/*
   move figure containg ell1 to touch ell1 to ell2
   along direction parallel to given axis of
   reference ellipsoid.

   called by  prfrm,
   calls      save, restore, sepn, dcen,
              rotget, vecmul, domoveby, solve, fndmin,
*/
{
	var j;//int
	var steps;//int
	var min,max;//double
	var mov;//double
	var dold,dnew;//double
	var xold,xnew;//double
	var forward,back;//double
	var shft;//double
	var dist;//double
	var cdist;//double
	//double d[3];
	var d = Array();
	//double dx[3];
	var dx = Array();
	//double v[EMAX][3];
	var v = get2DArray(EMAX);
	//double r[3][3],unr[3][3];
	var r = get2DArray(3);

	save();
	min = minax[ell1];
	if (min < minax[ell2]) min = minax[ell2];
	max = maxax[ell1];
	if (max < maxax[ell2]) max = maxax[ell2];
	if ((max > doub0) && (lg2 != doub0)) 
		steps = int(doub2 + log(max)/lg2);
	else 
		steps = 2;
	dist = sepn();
	if (ok != 0) return;
	for (j = 0; j < 3; ++j)
		dx[j] = doub0;
/*
   do they already just touch -
*/
	if (dist == doub0) 
		return;
	else
	if (dist < doub0)
/*
   they overlap already, so seek shortest way to separate them -
*/
	{
		mov = doub2*max;
		dx[axis] = mov;
		domoveby(dx[0],dx[1],dx[2],refell);
		forward = sepn();
		if (ok != 0) return;
		restore();
		dx[axis] = -mov;
		domoveby(dx[0],dx[1],dx[2],refell);
		back = sepn();
		if (ok != 0) return;
		restore();
		if ((back < doub0) && (forward < doub0))
		{
			ok = 58;
			printf("doabut: ell1 %d,  ell2 %d,  back %f,  forward %f",            
			  ell1,ell2,back,forward);
			return;
		}
		else
		if (back > forward) mov = -mov;
		shft = solve(doub0,mov,steps);
	}
	else
/*
  try to overlap them -
*/
	{
		cdist = dcen(d);
		for (j = 0; j < 3; ++ j)
			v[0][j] = d[j] ;
		rotget(r,unr,refell) ;
		vecmul(v,unr,0) ;
		mov = -v[0][axis];
		dnew = doubmax;
		xold = mov - max - min;
		for (xnew = mov-max; xnew < (mov+max+min); xnew += min*inv2)
		{
			dx[axis] = xnew;
			domoveby(dx[0],dx[1],dx[2],refell);
			dold = dnew;
			dnew = sepn();
			if (ok != 0) return;
			restore();
			if ((mov > doub0) && (dnew < doub0)) gotit();
			if ((mov < doub0) && (dnew > dold)) gotit();
			xold = xnew;
		}
		shft = xnew;
		return;
	}
	
/*
   they wont overlap so just bring them closest together -
*/
	function gotit(){
		if ((dold > doub0) && (dnew > doub0)){
			shft = fndmin(xold-min,xnew,steps);
			restore();
			dx[axis] = prop*shft;
			domoveby(dx[0],dx[1],dx[2],refell);
		}
/*
   they will overlap -
*/
		shft = solve(xold,xnew,steps);
/*
   move proportion to abut -
*/
		restore();
		dx[axis] = prop*shft;
		domoveby(dx[0],dx[1],dx[2],refell);
	}
} /* doabut */
/***************************************************/
/* 
   drawel.cpp version 42

   subroutines -
      getout      - exit keeping window open
      initialise  - initialise variables and constants
      openfile    - read file root and open input file
      action      - selects an action to perform
      doperfrm    - performs each action in turn for a given frame
      doframes    - simulates and stores each frame in turn
      help        - list interactive commands
      initsphere  - set up polyhedral approximation to a sphere
      dopause     - do nothing for a while
      donum       - write numbers on the animation window
      image       - draw the ellipsoids
      animate     - update interactive variables
      checkeys    - respond to keyboard commands
      initgraphics- initialise graphics window
      main       - run the show

/***********************************************/

//TODO REVIEW (File opening over http - I'm FAIRLY sure this isn't needed)
//void openfile(void)
//var openfile = function()
/*
   open the nudes file written by linter routine

   called by main,
*/
/*{
	if ((infile = fopen(nudesname,"r")) == NULL)
	{
		if (infile) fclose(infile);
		printf("\n\n "+nudesname+" oops?\n");
		ok = 2;
		getout(1);
		if (ok == 1) goto rtrn;
	}
	printf("\n   opened %s\n",nudesname);
rtrn: ;
}*/ /* openfile */
/***************************************/

//void dolighting ( double x, double y, double z )
var dolighting = function( x, y, z )
/*
change of lighting
*/
{
	lighting_rgb[0] = x / 255.0;
	lighting_rgb[1] = y / 255.0;
	lighting_rgb[2] = z / 255.0;	
} /* dolighting */
/******************************************/

//void doopacity ( void )
var doopacity = function()
/*
	change of opacity

	calls elground, rotput, rset, setnup
*/
{
	var k,n;//int
	var y, phi;//double
	//double r[3][3];
	var r = get2DArray(3);
	//double axe[3];
	var axe = Array(3);

	// run thru ellipsoids to shadow each in turn -

	k = ne;
	for ( n = 1; n < ne; ++n )
	{
		phi = setnup ( n, axe );
		y = elground ( n );
		if ( y > doub0 )
		{
			cen[k][0] = cen[n][0];
			cen[k][1] = -inv5;
			cen[k][2] = cen[n][2];
			ax[k][0] = axe[0];
			ax[k][1] = inv5;
			ax[k][2] = axe[2];
			rset ( r, phi, 1 );
			rotput ( r, k );
			col[k][0] = doub1;
			col[k][1] = doub1;
			col[k][2] = doub1;
			col[k][3] = 100;  // transparency 0-100 flags 101-255
			
			++k;
		} /* y > 0 */
	}	  /* end n loop */
	ne = k;
} /* doopacity */
/*****************************************/

//void doangles(int el, int re, double v[3])
var doangles = function(el, re, v)
/*
	store the angles of 'el' relative to 're' in 'val' array.
	in degrees.

	called by action, dodrag
	calls  matmul, rotget, rotput, mkang,

	14 Aug 2006  answers put in 3 element array
*/
{
	//double mvro[3][3],mvunro[3][3];
	var mvro = get2DArray(3);
	var mvunro = get2DArray(3);
	//double stro[3][3],stunro[3][3];
	var stro = get2DArray(3);
	var stunro = get2DArray(3);
	//double r[3][3];
	var r = get2DArray(3);

	rotget(stro,stunro,re);
	rotget(mvro,mvunro,el);
	matmul(stunro,mvro,r);
	rotput(r,EMAX);
	mkang(EMAX);
	v[0] = ang[0] * degree;
	v[1] = ang[1] * degree;
	v[2] = ang[2] * degree;
	if ((v[0] > doub179) && (v[0] < doub181))
	{
		v[0] -= doub180;
		v[2] = -v[2];
	}
	if (v[1] > doub180) v[1] -= doub360;
}  /* doangles */
/***********************************************/

//void dogrowto ( double x, double y, double z )
var dogrowto = function(x, y, z )
/*
	find distance between extremal points in each direction 
	and scale to desired output

	growto  (fname)  (referenceellipsoid) (axis) (size)
		where size =
		(variablename)
		(value)
	makes the total extent of the nominated figure in its current articulated
	state equal to the value of size, in the direction parallel to the nominated
	axis of the reference ellipsoid (tangent plane to tangent plane, normal to
	this axis direction).

	called by  action
	calls dogroell, shift, bell,
*/
{
	//double cen_min[3];
	var cen_min = Array();
	//double cen_max[3];
	var cen_max = Array();
	//double cen_dif[3];
	var cen_dif = Array();
	//double cen_scale[3];
	var cen_scale = Array();

	var e;//int
	var j;//int
	var n;//int

	for ( j = 0; j < 3; ++j )
	{
		cen_min[j] = doubmax;
		cen_max[j] = -doubmax;
	}

	shift( -x, -y, -z );
	for ( n = 0; n < ecount; ++ n )
	{
		e = elist[n];
		maxax[e] = ax[e][0];
		minax[e] = ax[e][0];
		for ( j = 1; j < 3; ++j )
		{
			if ( ax[e][j] > maxax[e] ) maxax[e] = ax[e][j];
			if ( ax[e][j] < minax[e] ) minax[e] = ax[e][j];
		}
		for ( j = 0; j < 3; ++j )
		{
			if ( cen[e][j] > cen_min[j] ) cen_min[j] = cen[e][j];
			if ( cen[e][j] < cen_max[j] ) cen_max[j] = cen[e][j];
		}
	}

	for ( j = 0; j < 3; ++j )
	{
		cen_dif[j] = -( cen_max[j] - cen_min[j] );
		cen_scale[j] = 1.0;
		if ( cen_dif[j] > 0.0 ) cen_scale[j] = factor[j] / cen_dif[j];
		factor[j] = cen_scale[j];
	}

	for ( n = 0; n < ecount; ++ n )
	{
		e = elist[n];
		dogroell ( factor, e, cen );
		dogroell ( factor, e, ax );
		maxax[e] = ax[e][0];
		minax[e] = ax[e][0];
		for ( j = 1; j < 3; ++j )
		{
			if ( ax[e][j] > maxax[e] ) maxax[e] = ax[e][j];
			if ( ax[e][j] < minax[e] ) minax[e] = ax[e][j];
		}
	}

	for ( n = 0; n < jcount; ++ n )
		dogroell ( factor, jlist[n], jnt );

	shift ( x, y, z );
	return;
} /* dogrowto */
/******************************************/

//void action ( int keyword_code )
var action = function( keyword_code )
/*
  decode and do an action keyword_code.

  called by  doperfrm,
  calls      doabut, doangles, doattach, dobalanc, dobend,
		docolour, dodetach, dodrag, dogroell, dogrofig, dogrojnt,
		dogrowto, doground, domoveby, domovjnt, dolighting,
		doopacity, doplace, dospinby, dospinto, dotouch,
		enquir, setobs, shift,

	 1 Sep 2006  allow world world means allow any intersection 
	14 Aug 2006  altered parameters of doangles()
	13 Aug 2006  added allow and forbid

*/
{
	var min;//double
	//double v[3];
	v = Array();
	var j,k;//int

	if ( ( keyword_code < 7 ) || ( keyword_code > NKEYS ) )
	{
		/*
		int figure_keyword_code  =  1;
		int ellips_keyword_code  =  2;
		int joint_keyword_code   =  3;

		int accele_keyword_code  =  5;
		int subrou_keyword_code  =  6;
		*/

		ok = 24;
		printf ( "action type %d out of range %d %d\n", keyword_code, 7, NKEYS );
	}

	if ( keyword_code == balanc_keyword_code ) dobalanc();

	if ( keyword_code == attach_keyword_code ) doattach();

	if ( keyword_code == detach_keyword_code ) dodetach();

	if ( keyword_code == grofig_keyword_code ) dogrofig ( xx[0], xx[1], xx[2] );

	if ( keyword_code == spinto_keyword_code ) dospinto ( xx, refell, ang, prop );

	if ( keyword_code == moveby_keyword_code ) domoveby ( frac * xx[0], frac * xx[1], frac * xx[2], refell );

	if ( keyword_code == add_keyword_code ) val[var0] = xx[0] + xx[1];

	if ( keyword_code == touch_keyword_code ) dotouch ();

	if ( keyword_code == spinby_keyword_code ) dospinby ( xx, refell, ang[0] * frac, axis );

	if ( keyword_code == ground_keyword_code )
	{
		min = doground ();
		shift ( doub0, -prop * min, doub0 );
	}
	if ( keyword_code == bendby_keyword_code ) dospinby ( xx, refell, ang[0] * frac, axis );

	if ( keyword_code == set_keyword_code ) val[var0] = varval;

	if ( keyword_code == bendto_keyword_code ) dospinto ( xx, refell, ang, prop );

	if ( keyword_code == repeat_keyword_code )
	{
		ok = 25;
		console.log("non-existent action "+keyword_code+"\n");
	}
	if ( keyword_code == quadra_keyword_code )
	{
		ok = 26;
		console.log("non-existent action "+keyword_code+"\n");
	}
	if ( keyword_code == linear_keyword_code )
	{
		ok = 27;
		console.log("non-existent action "+keyword_code+"\n");
	}
	if ( keyword_code == observ_keyword_code ) setobs();

	if ( keyword_code == moveto_keyword_code ) shift ( prop * xx[0], prop * xx[1], prop * xx[2] );

	if ( keyword_code == invert_keyword_code )
	{
		if ( val[var0] != doub0 )
		{
			val[var0] = doub1 / val[var0];
		}
		else
		{
			val[var0] = doub0;
		}
	}
	if ( keyword_code == groell_keyword_code ) dogroell ( factor, ellpsd, ax );

	if ( keyword_code == grojnt_keyword_code ) dogrojnt();

	if ( keyword_code == growto_keyword_code ) dogrowto( xx[0], xx[1], xx[2] );

	if ( keyword_code == angles_keyword_code )
	{
		doangles ( ellpsd, refell, v );
		val[var0] = v[0]; val[var1] = v[1]; val[var2] = v[2];
	}

	if ( keyword_code == centre_keyword_code || keyword_code == center_keyword_code) enquir ( ellpsd, cen );

	if ( keyword_code == flex_keyword_code ) dobend ( ang[0] * frac, 0 );

	if ( keyword_code == rotate_keyword_code ) dobend ( ang[0] * frac, 1 );

	if ( keyword_code == abduct_keyword_code ) dobend ( ang[0] * frac, 2 );

	if ( keyword_code == negate_keyword_code ) val[var0] = -val[var0];

	if ( keyword_code == subtra_keyword_code ) val[var0] = xx[0] - xx[1];

	if ( keyword_code == divide_keyword_code )
	{
		if ( xx[1] == doub0 )
		{
			ok = 30;
			console.log("action: divide, keyword_code "+keyword_code+",  EMAX-var0-1 "+(EMAX-var0-1)+",  xx[1] "+xx[1]+"");
		}
		else
		{
			val[var0] = xx[0] / xx[1];
		}
	}
	if ( keyword_code == multip_keyword_code ) val[var0] = xx[0] * xx[1];

	if ( keyword_code == 45 )
	{
		ok = 28;
		console.log("action: non-existent action "+keyword_code+"\n");
	}
	if ( keyword_code == cubic_keyword_code )
	{
		ok = 29;
		console.log("action: non-existent action "+keyword_code+"\n");
	}
	if ( keyword_code == place_keyword_code ) doplace();

	if ( keyword_code == axes_keyword_code ) enquir ( ellpsd, ax );

	if ( keyword_code == linkx_keyword_code ) enquir ( join, jnt );

	if ( keyword_code == colour_keyword_code ||
		keyword_code == color_keyword_code ) docolour ( prop );

	if ( keyword_code == print_keyword_code ) {
		console.log("frame "+f+", "+(vname[EMAX-var0-1])+" "+val[var0]+"\n");
	}
	
	if ( keyword_code == textur_keyword_code ) docolour ( doub1 );

	if ( keyword_code == drag_keyword_code ) dodrag ();

	if ( keyword_code == abut_keyword_code ) doabut ();

	if ( keyword_code == movjnt_keyword_code ) domovjnt ();

	if ( keyword_code == opacty_keyword_code ) doopacity ();

	if ( keyword_code == lghtng_keyword_code ) dolighting ( xx[0], xx[1], xx[2] );

	if ( keyword_code == allow_keyword_code ) 
	{
		//printf("actionb allow %s %s\n",ename[ell1],ename[ell2]);
		forbid[ell1][ell2] = false;
		forbid[ell2][ell1] = false;
		if ((ell1+ell2) == 0)
		{
			for (j = 0; j < EMAX; ++j)
				for (k = 0; k < EMAX; ++k)
					forbid[j][k] = false;
		}
	}
	if ( keyword_code == forbid_keyword_code ) 
	{
		forbid[ell1][ell2] = true;
		forbid[ell2][ell1] = true;
		//printf("actionc forbid %s %s\n",ename[ell1],ename[ell2]);
	}
}  /* action */
/*********************************************************/

//void doperfrm(int sub, int fr, int fstart, int fend)
var doperfrm = function(sub, fr, fstart, fend)
/*
   perform actions of subroutine 'sub' for frame 'fr'
   between frames 'fstart' and 'fend'

   called by  doframes,
   calls      getvalu, setfrc, setper, action,
*/
{
	var frame;//int
	var newsub;//int
	var p;//int
	var pstart,pend;v
	var fsstart;//int
	var fstrt,fstp;//int
	var fsubstart ;//int

	pstart = subact[sub][0] ;
	pend = subact[sub][1] ;
/*
  find 'subfrm', the earliest formal frame number in 
  current subroutine ignoring unset variable 
  frame numbers ( == -1) -
*/
	fsubstart = MAXINT ;
	for (p = pstart ; p <= pend ; ++p )
	{
		fsstart = getvalu(frstart[p]) ;
		if (ok != 0){
			snag();
			return;
		}
		if (fsstart >= 0)
		{
			fsstart *= fast ;
			if (fsstart < fsubstart) fsubstart = fsstart ;
		}
	}
	frame = fr + fsubstart - fstart;
	if (fr >= fstop) fstop = fr+1;
/*
  run through actions in the subroutine -
*/
	for (p = pstart; p <= pend; ++p)
	{
		pp = p;
		t = type[p] ;
		if ((t != stop_keyword_code)
			&&(t != subrou_keyword_code)
			&&(t != endsub_keyword_code))
		{
			fstrt = getvalu(frstart[p]) ;
			if (fstrt < 0)
			{
				ok = 46;
				console.log("doperfrm: start "+start+"\n");
			}
			if (ok != 0){
				snag();
				return;
			}
			fstp = getvalu(frstop[p]) ;
			if (fstp < fstrt)
			{
				ok = 47 ;
				console.log("doperfrm: fstrt "+fstrt+",  fstp "+fstp+"\n");
			}
			if (ok != 0){
				snag();
				return;
			}
			if (fstrt == fstp) ++fstp;
			fstrt *= fast ;
			fstp *= fast ;
			if ((fstp > frame) && (fr < fend)) ++more;
			if ((frame > fstrt) && (frame <= fstp))
			{
				if (t == call_keyword_code)
				{
					newsub = getvalu(pf[p][0]);
					if ((newsub <= 0) || (newsub > PMAX))
					{
						ok = 8 ;
						console.log("doperfrm: newsub "+newsub+",  PMAX "+PMAX+"\n");
						snag();
						return;
					}
					if (distrn[p] == call_keyword_code)
						doperfrm(newsub,frame,fstrt,fstp);
					else
						doperfrm(newsub,frame,frame-1,frame);
					if (ok != 0) return; ;
				}
/*
  if not a subroutine call, then do normal action -
*/
				else
				{
					setfrc(frame,fstrt,fstp) ;
					if (ok != 0){
						snag();
						return;
					}
					setper(t);
					if (ok != 0){
						snag();
						return;
					}
					action(t);
					if (ok != 0){
						snag();
						return;
					}
				} /* t not call */
			} /* frame f in range of action p */
		} /* t not start or end of subroutine */
	} /* p */
	return;
/*
  snag-
*/
	function snag(){
		console.log("error in doperfrm, frame "+f+", action "+p+1); 	  getout(ok+"\n");
		getout(ok);
	}
}  /* doperfrm */
/*********************************************/

//void doframes(void)
var doframes = function()
/*
  calls and performs the actions for each frame in turn.

  called by  main,
  calls      doperfrm, store3,
*/
{
	var j;//int

	t = 1;
	axis = 1;
	join = 1;
	var0 = 1;
	var1 = 1;
	var2 = 1;
	fig = 1;
	ellpsd = 1;
	refell = 1;
	ell1 = 1;
	ell2 = 1;
	for (j = 0; j < 3; ++j)
	{
		newcol[j] = 0;
		oldang[j] = doub0;
		ang[j] = doub0;
		xx[j] = doub0;
		factor[j] = doub1;
	}
	varval = doub0;
/*
  simulate and store each frame of movie-
*/
	for (f = 1; more > 0; ++f)
	{
		more = 0;
		doperfrm(0,f,0,fstop);
		if ((fslow == 1) || (f%fslow == 1))
			store3(f);
	}
	if (vstart > fstart) fstart = vstart;
	if ((vstop > 0) && (vstop < fstop)) fstop = vstop;
} /* doframes */
/*********************************************/

//void help(void);
var help = function()
/*
   list the interactive commands
 
   called by  main, checkeys,
*/
{
  console.log("\n\n******* TO ACTIVATE THESE: CLICK IN THE ANIMATION WINDOW FIRST *******\n");
  console.log("\n  Interactive commands :-\n\n");
 
  console.log("\n\n   animation parameters\n");
  console.log("    i - freeze (opp. of 'a')\n"); 
  console.log("    a - continue animating (opp. of 'i')\n");
  console.log("    c - continue through restart at full rate (opp. of 'p')\n");
  console.log("    p - pause on first and last frames (opp. of 'c')\n"); 
  console.log("    b - if frozen, go back 1 frame else run backwards (opp. of 'f')\n");
  console.log("    f - if frozen, go forward 1 frame else run forwards (opp.of 'b')\n");
  console.log("    0 - reset parameters to defaults, and freeze at start\n");
  console.log("    - - slow down the animation \n");
  console.log("    = - speed up the animation \n"); 
 
  console.log("\n\n   scene movement parameters\n");
  console.log("    d - shift down 10 per cent (opp. of 'u')\n");
  console.log("    u - shift up 10 per cent (opp. of 'd')\n");
  console.log("    l - shift scene left 10 per cent (opp. of 'r')\n");
  console.log("    r - shift scene right 10 per cent (opp. of 'l')\n");
  console.log("    t - shift scene away by 10 per cent(opp. of 't')\n");
  console.log("    v - shift away (opp. of 'w')\n");
  console.log("    w - shift nearer (opp. of 'v')\n");
 
  console.log("\n\n   scene rotation parameters\n");
  console.log("    x - rotate 3 degrees about x (left - right) axis (opp. of '1')\n");
  console.log("    y - rotate 3 degrees about y (vertical) axis (opp. of '2')\n");
  console.log("    z - rotate 3 degrees about z (front - back) axis  (opp. of '3')\n");
  console.log("    1 - rotate 10 degrees about x (right - left) axis (opp. of 'x')\n");
  console.log("    2 - rotate 10 degrees about y (vertical) axis (opp. of 'y')\n");
  console.log("    3 - rotate 10 degrees about z (back - front) axis  (opp. of 'z')\n"); 
 
  console.log("\n\n   scene size parameters\n");
  console.log("    g - grow scene by 10 per cent (opp. of 's')\n"); 
  console.log("    s - shrink scene by 10 per cent (opp. of 'g')\n");
 
  console.log("\n\n   polygon parameters\n");
  console.log("    j - increase the number of polygons per sphere by 1 {opp. of 'k')\n"); 
  console.log("    k - decrease the number of polygons per sphere by 1 {opp. of 'j')\n"); 
 
  console.log("\n\n   shading parameters\n");
  console.log("    4 - normal shading/modified shading (toggle)\n"); 
  console.log("    6 - display shadows (toggle) - inoperative at present\n");
  console.log("    8 - display footprints (toggle) - inoperative at present\n");
 
  console.log("\n\n   annotation parameters\n");
  console.log("    n - display of frame numbers (toggle)\n");
  console.log("    o - display of bar numbers (toggle)\n");
  console.log("    5 - display of time (toggle) - inoperative at present\n");
 
  console.log("\n\n   avi output parameters\n");
  console.log("    7 - save display as AVI file - inoperative at present\n");
 
  console.log("\n\n   program usage parameters\n");
  console.log("    h - show these instructions\n");
  console.log("    q - quit\n");
 
  console.log(" \n\n******* TO ACTIVATE THESE: CLICK IN THE ANIMATION WINDOW FIRST *******\n");
} /* help */
/********************************************/

//void dopause(int t)
var dopause = function(t)
/*
   pause for a while

   called by  image, animate,
*/
{
   var a//double;
   var j,k;//int

   a = 1.0;
   for (j = 1; j < t*1000; ++j)
   {
      for (k = 1; k < 1000; ++k)
      {
         a += j/k;	//Both j and k were cast to double
      }
   }
   if (a < 0.0) printf("%d\n",t);
} /* dopause */
/*****************************************/

//void donum(int f)
var donum = function(f)
/*
   draw bar and/or frame numbers 

   called by image,
*/
{
	//char nstr[BMAX];
	var nstr = Array();
	var b,c,nlngth;//int

	if (bnums == TRUE)
	{
		b = int(doub1 + inv2 + (f/frperbar));	//TODO REVIEW (f and frperbar were cast to double)
		nstr = "bar "+b;
		nlngth = strlen(nstr);
		glColor3f(0.9, 0.0, 0.9);
		glRasterPos3f(0.05,0.95,0.999999);
		for (c = 0; c < nlngth; c++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, nstr[c]);
	}
	if (fnums == TRUE)
	{
		nstr = "frame "+f;
		nlngth = strlen(nstr);
		glColor3f(0.9, 0.0, 0.9);
		glRasterPos3f(0.05,0.05,0.999999);
		for (c = 0; c < nlngth; c++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, nstr[c]);

	  //pres_time = clock();
	  //rate = CLOCKS_PER_SEC/(pres_time - prev_time);
	  //prev_time = pres_time;
	  //glRasterPos3f(0.9,0.05,0.999999);
	  //printf("%d fr/sec  %d %d %d\n",
		  //rate,pres_time,prev_time,CLOCKS_PER_SEC);
	  //sprintf(nstr, "%d fr/sec",rate);
	  //nlngth = strlen(nstr);
	  //for (c = 0; c < nlngth; c++)
         //glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, nstr[c]);
	}
} /* donum */
/**********************************************************/

//void ellmat(double r[3][3], int e, int f)
var ellmat = function(r, e, f)
/*
   for ellipsoid 'e' in frame 'f' find its 
   rotation matrix 'r',

   called by  image,
*/
{
	var p,cp,sp;//double
	var x,y,z,m,xsp,ysp,zsp,xm,ym,zm,xym,xzm,yzm ;//double

	p = qu3[f][e][0]*radian ;
	x = qu3[f][e][1] ;
	y = qu3[f][e][2] ;
	z = qu3[f][e][3] ;
	sp = sin(p);
	cp = cos(p);
	m = doub1-cp ;
	xm = x*m ;
	ym = y*m ;
	zm = z*m ;
	xsp = x*sp ;
	ysp = y*sp ;
	zsp = z*sp ;
	xym = x*ym ;
	xzm = x*zm ;
	yzm = y*zm ;
	r[0][0] = x*xm+cp ;
	r[0][1] = xym+zsp ;
	r[0][2] = xzm-ysp ;
	r[1][0] = xym-zsp ;
	r[1][1] = y*ym+cp ;
	r[1][2] = yzm+xsp ;
	r[2][0] = xzm+ysp ;
	r[2][1] = yzm-xsp ;
	r[2][2] = z*zm+cp ;
} /* ellmat */
/***********************************************/

//void image(void) 
var image = function() 
/*
   called by  main,
   calls      donum, dopause, ellmat,
*/
{ 
    var e,j,k;//int
    var amb,shade,illum;//double
    //double r[3][3];
	var r = get2DArray();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    for (e = 1; e < nels[f]; ++e)
    {
		glPushMatrix();

		// interactive parameters -

		glScalef(doub1,doub1,inv10);
		glScalef(scale,scale,scale);
		glTranslatef(tx,ty,tz);
		glTranslatef(inv2,inv2,inv2);
		glRotatef(anglez,0,0,1);
		glRotatef(angley,0,1,0);
		glRotatef(anglex,1,0,0);
         
		// ellipsoid parameters -

		glTranslatef(ce3[f][e][0], ce3[f][e][1], ce3[f][e][2]);
		glRotatef(qu3[f][e][0], qu3[f][e][1], qu3[f][e][2], qu3[f][e][3]);
		glScalef(ax3[f][e][0], ax3[f][e][1], ax3[f][e][2]);
		ellmat(r,e,f);
		glPolygonMode(GL_FRONT, GL_FILL);
		glBegin(GL_QUADS);
		amb = 0.7;
		for (j = 0; j < nfaces; ++j)
		{
            illum = doub0;
            for (k = 0; k < 3; ++k)
                illum += r[k][1]*norm[j][k];
            shade = amb + (doub1 - amb)*illum;
            glColor3f(shade*co3[f][e][0], shade*co3[f][e][1], shade*co3[f][e][2]);
            for (k = 0; k < 4; ++k)
            {
				glVertex3f(sph[j][k][0],sph[j][k][1],sph[j][k][2]);
            } /* k  vertices */
		} /* j faces */
		glEnd();
		glPopMatrix();
    } /* e ellipsoids */
    if ((fnums == TRUE)||(bnums == TRUE)) donum(f);
    glutSwapBuffers();
    glFlush();
    if (slow < 0) slow = 0;
    if (slow > 0) dopause(slow);
} /* image */
/***************************************/

//void animate(void)
var animate = function()
/*
   called by  main,
   calls      dopause,
*/
{
    if (freeze == TRUE)
    {
		if (single == TODO)
		{
            if (forward == TRUE)
				df = 1;
            else
				df = -1;
		}
		else
            df = 0;
    }
    else
    if (forward == TRUE)
		df = 1;
    else
        df = -1;
    f += df;
    if (f < 1) f = fstop-1;
    if (f >= fstop) f = fstart+1;
    if (((f == (fstart+2)) || (f == (fstart+1))) && (pause == TRUE))
        dopause(100);
    anglex += dangx; 
    angley += dangy;
    anglez += dangz;
    if (anglex < -doub180) anglex += doub360;
    if (anglex > doub180) anglex -= doub360;
    if (angley < -doub180) angley += doub360;
    if (angley > doub180) angley -= doub360;
    if (anglez < -doub180) anglez += doub360;
    if (anglez > doub180) anglez -= doub360;
    dangx = doub0;
    dangy = doub0;
    dangz = doub0;
    glutPostRedisplay();
    single = DONE;
} /* animate */
/***************************************/

//void checkeys(unsigned char key, int x, int y) 
var checkeys = function(key, x, y) 
/*
   called by  main,
   calls      initsphere,
*/
{ 
	if ((key == GLUT_KEY_ESCAPE) || (key == 'q'))
	{
		getout(0);
		if (ok == 1) _tmain(0,junk);
	}
	if (key == 'h') help();
	if (key == 'j')
	{
		nsph +=1;
		console.log("'j' more facets "+nsph+" (opp. 'k')\n");
		console.log("checkeysa "+key+"\n");
		initsphere();
		if (ok == 1) _tmain(0,junk);
	}
	if (key == 'k')
	{
		nsph -=1;
		console.log("'k' fewer facets "+nsph+" (opp. 'j')\n");
		console.log("checkeysb "+key+"\n");
		initsphere();
		if (ok == 1)  _tmain(0,junk);
	}
	if (key == 'a') 
	{
		df = 1;
		console.log("'a' animate (opp. 'i')\n");
		if (forward == FALSE) df = -1;
		freeze = FALSE;
	}
	if (key == 'i') 
	{ 
		freeze = TRUE; 
		console.log("'i' freeze (opp. 'a')\n");
	}
	if (key == 'f')
	{
		if (freeze == TRUE) single = TODO;
		console.log("'f' forwards (opp. 'b')\n");
		forward = TRUE;
	}
	if (key == 'b') 
	{
		if (freeze == TRUE) single = TODO;
		console.log("'b' backwards (opp. 'f')\n");
		forward = FALSE; 
	}
	if (key == 'p') { pause = TRUE; printf("'p' pause on first and last frames (opp. 'c')\n"); }
	if (key == 'c') { pause = FALSE; printf("'c' continuous looping (opp. 'p')\n"); }
	if (key == 'g') { scale *= 1.1; printf("'g' grow %f (opp.'s')\n", scale); }
	if (key == 's') { scale /= 1.1; printf("'s' shrink %f (opp. 'g')\n", scale); }
	if (key == 'l') { tx -= 0.1; printf("'l' shift left %f (opp. 'r')\n", tx); }
	if (key == 'r') { tx += 0.1; printf("'r' shift right %f (opp. 'l')\n", tx); }
	if (key == 'd') { ty -= 0.1; printf("'d' shift down %f (opp. 'u')\n", ty); }
	if (key == 'u') { ty += 0.1; printf("'u' shift up %f (opp. 'd')\n", ty); }
	if (key == 'v') { tz -= 0.1; printf("'v' shift in %f (opp. 'w')\n", tz); }
	if (key == 'w') { tz += 0.1; printf("'w' shift away %f (opp. 'v')\n", tz); }
	if (key == 'x') { dangx = alpha; printf("'x' rotate x %f (opp. '1')\n", anglex); }
	if (key == '1') { dangx = -alpha; printf("'1' rotate x %f (opp. 'x')\n", anglex); }
	if (key == 'y') { dangy = alpha; printf("'y' rotate y %f (opp. '2')\n", angley); }
	if (key == '2') { dangy = -alpha; printf("'2' rotate y %f (opp. 'y')\n", angley); }
	if (key == 'z') { dangz = alpha; printf("'z' rotate z %f (opp. '3')\n", anglez); }
	if (key == '3') { dangz = -alpha; printf("'3' rotate z %f (opp. 'z')\n", anglez); }
	if (key == '-') { slow += 2; printf("'-' slower %d (opp. '=')\n", slow); }
	if (key == '=') { slow -= 2; printf("'=' faster %d (opp. '-')\n", slow); }
	if (key == 'n')
	{ 
		if (fnums == TRUE)
		{
			fnums = FALSE;
			console.log("'n' hide frame numbers (toggle)\n");
		}
		else
		{
			fnums = TRUE;
			console.log("'n' show frame numbers (toggle)\n");
		}
	}
	if (key == 'o')
	{ 
		if (bnums == TRUE)
		{
			bnums = FALSE;
			console.log("'n' hide bar numbers (toggle)\n");
		}
		else
		{
			bnums = TRUE;
			console.log("'n' show bar numbers (toggle)\n");
		}
	}
	if (key == '0') 
	{
		console.log("'0' reset parameters and freeze at frame 1\n");
		f = fstart+1;
		freeze = TRUE;
		forward = TRUE;
		df = 1;
		scale = SCALE;
		tx = doub0;
		ty = doub0;
		tz = doub0;
		anglex = doub0;
		angley = doub0;
		anglez = doub0;
		slow = 1;
		nsph = SSTART;
	} /* key = '0' */
} /* checkeys */
/***************************************/

//TODO ERRORAGE FILE *test_File won't work
//void add_id_num ( char name[], char outname[], char ext[] )
var add_id_num = function ( name, outname, ext )
{
	FILE *test_File;
	var j;//int

	for ( j = 0; j <= 999; j++ )
	{
		sprintf ( outname, "%s_%03d%s", name, j, ext );

		if ( ( test_File = fopen ( outname, "r" ) ) != NULL )
		{
			fclose ( test_File );
		}
		else
		{
			return;
		}
	}
	outname = name + "_000" + ext;
} /* add_id_num */
/*******************************************/

//int find_ini_title ( char title[] )
var find_ini_title = function( title )
/*
   called by get_ini_str, get_ini_bool, 
	          get_ini_char,
*/
{
	var value = -1;//int
	var ini_no;//int
	var j;//int
	var plen;//int
	var iplen;//int
	
	if ( number_ini <= 0 ) return( NULL );
	for ( ini_no = 0; ini_no < number_ini; ini_no++ )
	{
		plen = parseInt(strlen ( title ));//TODO REVIEW (was cast to int)
		iplen = 0;
		for ( j = 0; j < plen; j++ )
		{
			if ( title[j] == ini_title[ini_no][j] )
			{
				iplen = iplen + 1;
			}
		}
		if ( iplen == plen )
		{
			return( ini_no );
		}
	}
	return( value );
} /* find_ini_title */
/************************************************/

//void get_ini_dump ( void )
var get_ini_dump = function()
{
	var ini_no;//int
	
	console.log("  number ini "+number_ini+"\n");
	if ( number_ini <= 0 ) return;
	for ( ini_no = 0; ini_no < number_ini; ini_no++ )
	{
		console.log(" ini_no "+leadingZeros(ini_no, 2)+" title "+ini_title[ini_no][0]+" value "+ini_value[ini_no][0]+"\n");	//ini_title and ini_value had & symbols next to them.
	}
} /* get_ini_dump */
/************************************************/

//bool get_if_ini ( void )
var get_if_ini = function()
{
	if ( number_ini > 0 ) return( true );
	return( false );
} /* get_if_ini */
/************************************************/

//bool get_ini_bool ( char title[] )
var get_ini_bool = function(title)
{
	var value;//bool
	var ini_no;//int
	value = -1;NULL;
	if ( number_ini <= 0 ) return( NULL );
	ini_no = find_ini_title ( title );
	if ( ini_no < 0 ) return( NULL );
	if ( toupper( ini_value[ini_no][0] ) == 'T' )
		return( true );
	if ( toupper( ini_value[ini_no][0] ) == 'F' )
		return( false );
	return( NULL );

} /* get_ini_bool */
/************************************************/

//char* get_ini_char ( char title[] )
var get_ini_char = function( title )
/*
   calls find_ini_title, ini_value
*/
{
	var value;//char*
	var ini_no;//int
	value = NULL;
	if ( number_ini <= 0 ) return( NULL );
	ini_no = find_ini_title ( title );
	if ( ini_no < 0 ) return( NULL );
	return( ini_value[ini_no][0] );//TODO REVIEW (ini_value had &)

} /* get_ini_char */
/************************************************/

//int get_ini_int ( char title[] )
var get_ini_int = function( title )
/*
   calls find_ini_title, ini_value
*/
{
	var value = 0;//int
	var ini_no;//int
	if ( number_ini <= 0 ) return( NULL );
	ini_no = find_ini_title ( title );
	if ( ini_no < 0 ) return( NULL );
	value = atoi ( ini_value[ini_no][0] );//TODO REVIEW (ini_value had &)
	//printf ( " ini_no %d value %d\n", ini_no, value );
	return( value );

} /* get_ini_int */
/************************************************/

//float get_ini_float ( char title[] )
var get_ini_float = function( title )
/*
   calls find_ini_title, ini_value
*/
{
	var value = 0.0;//float	//TODO REVIEW (was 0.0f)
	var ini_no;//int
	if ( number_ini <= 0 ) return( NULL );
	ini_no = find_ini_title ( title );
	if ( ini_no < 0 ) return( NULL );
	value = atof ( ini_value[ini_no][0] );	//TODO REVIEW (ini_value has &)
	return( value );
} /* get_ini_float */
/************************************************/

//double get_ini_double ( char title[] )
var get_ini_double = function( title )
/*
   calls find_ini_title, ini_value
*/
{
	var value = 0.0;//double
	var ini_no;//int
	if ( number_ini <= 0 ) return( NULL );
	ini_no = find_ini_title ( title );
	if ( ini_no < 0 ) return( NULL );
	value = strtod ( ini_value[ini_no][0], NULL );//TODO REVIEW (ini_value had &)
	return( value );
} /* get_ini_double */
/************************************************/

//bool get_ini_str ( char title[], char value[] )
var get_ini_str = function(title, value)
/*
   calls find_ini_title, ini_value
*/
{
	//char* value;
	var ini_no;//int
	var i;//int
	var len;//int

	value[0] = NULL;
	if ( number_ini <= 0 ) return( false );
	ini_no = find_ini_title ( title );
	if ( ini_no < 0 ) return( false );
	len = parseInt(strlen( ini_value[ini_no][0] ));//TODO REVIEW (ini_value had &; strlen was cast to int)
	if ( len <= 0 )  return( false );
	i = -1;
	do
	{
		i = i + 1;
		value[i] = ini_value[ini_no][i];		
	} while ( ini_value[ini_no][i] != NULL );

	return( true );

} /* get_ini_str */
/************************************************/

//TODO ERRORAGE (will not work with the FILE* pointer)
//void get_ini ( int dump )
var get_ini = function( dump )
/*
   open  and decode .ini file

   called by main,
*/
{
	FILE *ini_file_unit;
	var ini_no;//int
	var j;//int
	var k;//int
	var len;//int
	var loc_comma;//int
	var loc_semi;//int
	var asterisk = '*';	//chat
	var blank = ' ';	//chat
	var comma = ',';	//chat
	var semi = ';';	//chat
	var getout;//bool

	for ( ini_no = 0; ini_no < max_ini; ini_no++ )
	{
		ini_title[ini_no][0] = NULL;
		ini_value[ini_no][0] = NULL;
	}

	number_ini = -1;
	ini_diag = 0;
	ini_file_unit = NULL;

	// open ini file - check if it exists

	if ( ( ini_file_unit = fopen ( "lintel.ini", "r" ) ) == NULL )
	{
		if ( ini_file_unit ) fclose ( ini_file_unit );
		console.log( "\n\n    lintel.ini not available - will continue\n\n");
	}
	else
	{
		ini_no = 0;
		number_ini = ini_no;
		len = -1;
		do
		{
			ini_title[ini_no][0] = NULL;
			ini_value[ini_no][0] = NULL;

			if ( fgets ( buf, BMAX, ini_file_unit ) != NULL )
			{
				if ( ini_diag >= 1 )
					console.log(" ini_no "+leadingZeros(ini_no, 2)+" buf "+buf+"");
				if ( buf[0] != asterisk )
				{
					if ( ini_diag >= 1 )
						console.log(" ini_no "+leadingZeros(ini_no, 2)+" buf "+buf+"");
					loc_comma = -1;
					loc_semi = -1;
					getout = false;
					len = parseInt(strlen( buf ));	//TODO REVIEW (Was cast to int)
					if ( ini_diag >= 1 ) printf ( " len %d\n", len );
					for ( j = 0; j < len; j++ )
					{
						if ( buf[j] == comma && loc_semi == -1 ) loc_comma = j;
						if ( buf[j] == semi  )
						{
							loc_semi = j;
							getout = true;
						}
						if ( getout == true ) break;
					}
					if ( ini_diag >= 1 )
						console.log(" loc_comma "+loc_comma+" loc_semi "+loc_semi+"\n");

					// get parameter title

					k = 0;
					for ( j = 0; j < loc_comma; j++ )
					{
						if ( buf[j] != blank )
						{
							ini_title[ini_no][k] = buf[j];
							k = k + 1;
						}

					}
					ini_title[ini_no][k] = NULL;

					// get parameter value

					k = 0;
					for ( j = loc_comma + 1; j < loc_semi; j++ )
					{
						if ( buf[j] != blank )
						{
							ini_value[ini_no][k] = buf[j];
							if ( ini_diag > 1 )
							{
								console.log(" j "+j+" k "+k+" buf[j] "+buf[j]+" ini "+ini_value[ini_no][k]+"\n");
							}
							k = k + 1;
						}
					}
					ini_value[ini_no][k] = NULL;
					ini_no = ini_no + 1;
				}
				else
				{
					if ( buf[1] == 'd' && buf[2] == 'u' 
						&& buf[3] == 'm' && buf[4] == 'p' )
						dump = 1;
				}
			}
		}
		while ( !feof( ini_file_unit ) && len != 0 );
		number_ini = ini_no;
	}
	if ( dump == 1 ) get_ini_dump ();
} /* get_ini */
/************************************************/

//bool strcmpend ( char str1[], char str2[] )
var strcmpend = function( str1, str2 )
/*
	compare strings to see if str2 is included at end of str1
*/
{
	var len1, len2;//int
	var i1, i2;//int
	var cnt;//int

	len1 = parseInt(strlen( str1 ));//TODO REVIEW (was cast to int)
	len2 = parseInt(strlen( str2 ));//TODO REVIEW (was cast to int)

	cnt = 0;
	i2 = len2 - 1;
	for ( i1 = len1 - 1; i1 >= len1 - len2; i1-- )
	{
		if ( str1[i1] == str2[i2] ) cnt = cnt + 1;
		i2 = i2 - 1;
	}
	if ( cnt == len2 )		return( true );

	return( false );
}/* strcmpend */

//void get_filesa ( bool lbn_type, int error )
var get_filesa = function( lbn_type, error )
{
	console.log( "\n" );
	if ( error == 0 )
	{
		console.log( "    Please type input filename followed by pressing the 'enter' key\n\n" );
	}
	else
	{
		console.log( "\n" );
		console.log( "    Please type a correct input filename\n\n" );
	}

	if ( lbn_type == true )
	{
		console.log( "      NUDES file (xxx.nud or xxx.n)\n" );
		console.log( "        - full filename (xxx.nud or xxx.n)\n" );
		console.log( "      LBN file (yyy.lbn)\n" );
		console.log( "        - root portion (yyy) of filename\n" );
		console.log( "           (interprets staves 1 and 2 with figure tracking)\n" );
		console.log( "        - full filename (yyy.lbn)\n" );
		console.log( "           (choice of staves, choice of tracking)\n\n" );
		console.log( "    Filename:  " );
	}


	if ( lbn_type == false )
	{
		console.log( "    Filename:  " );
	}
} /* get_filesa */
/************************************************/

//void get_files ( char file[] )
var get_files = function( file )
/*
   called by main
	calls get_filesa, strcmpend, bell, add_id_num,
*/
{
	var c;//int
	var i;//int
	var len;//int
	var last;//int
	var err_count;//int
	var error;//int
	var loc_dot;//int
	var from_ini;//int
	var key;//char
	var get_out;//bool
	var ini_ok;//bool
	var file_ok;//bool
	var dir_ok;//bool
	var lbn_type;//bool
	//char dir[BMAX];
	dir = Array();

	from_ini = 0;
	err_count = 0;
	error = 0;
	get_out = false;
	ini_ok = false;
	file_ok = false;
	dir_ok = false;
	lbn_type = true;

	var again = true;
	while(again){
		again = false;
		err_count = err_count + 1;
		if ( err_count >= 25 ) 
		{
			printf( " Limit: tried %d times for input file %s\n",
				err_count,name );
			ok = -1;
			return;
		}
		input_file_type = -1;
		for ( c = 0; c < BMAX; ++c )
		{
			name[c] = NULL;
			finname[c] = NULL;
			nudesname[c] = NULL;
		}

		if ( file == NULL )
		{
			file_ok = false;
			if ( from_ini == 0 )
			{
				if ( get_if_ini () == true )
				{
					ini_ok = get_ini_bool ( "input_file_default" );
					if ( ini_ok == true ) 
					{
						file_ok = get_ini_str ( "input_file_name", name );
						dir_ok = get_ini_str ( "input_file_dir", dir );
						if ( dir[0] == NULL ) dir_ok = false;
						len = parseInt(strlen( dir ));	//TODO REVIEW (Was cast to int)
						if ( dir_ok == true && dir[len - 1] != '\\' )
							dir[len - 1] = '\\';
						lbn_type = get_ini_bool ( "lbn_file_encoded" );
						from_ini = 1;
					}
				}
			}

			if ( file_ok == false )
			{
				name[0] = NULL;
				get_filesa ( lbn_type, error );
				if ( gets ( name ) != NULL )
				{
					len = parseInt(strlen( name ));	//TODO REVIEW (Was cast to int)
					if ( len == 0 )
					{
						get_out = true;
						error = 1;
						again = true;
						continue;
					}
				}
				else
				{
					get_out = true;
					error = 1;
					again = true;
					continue;
				}
			}
		}
		else
		{
			strcpy ( name, file );
		}

		len = parseInt(strlen( name )); //TODO REVIEW (was cast to int)
		last = len - 1;

		loc_dot = -1;
		i = -1;
		do
		{
			i = i + 1;
			key = name[i];
			if ( key == '.' ) loc_dot = i;
		} while ( key != NULL );

		if ( loc_dot >= 0 ) loc_dot = last - loc_dot;

		input_file_type = -1;
		haslbn = FALSE;
		get_out = false;

		if ( lbn_type == true ) // use filename to decide lbn type
		{
			switch ( loc_dot )
			{
			case 3:
				// .nud extention
				if ( strcmpend ( name, ".nud" ) ) 
				{
					input_file_type = 0;
					haslbn = FALSE;
				}
				// .lbn extention
				if ( strcmpend ( name, ".lbn" ) )
				{
					input_file_type = 1;
					haslbn = TRUE;
				}
				if ( input_file_type < 0 ) get_out = true;
				break;
			case 2:
				// problem
				get_out = true;
				break;
			case 1:
				// .n extention
				if ( strcmpend ( name, ".n" ) )
				{
					input_file_type = 0;
					haslbn = FALSE;
				}
				else
				{
					get_out = true;
				}
				break;
			case 0:
				// . extention
				if ( strcmpend ( name, "." ) )
				{
					input_file_type = 2;
					strcat( name, "lbn" );
					haslbn = TRUE;
				}
				else
				{
					get_out = true;
				}
				break;
			case -1:
				// no extention
				if ( len > 0 && !strcmpend ( name, "." ) )
				{
					input_file_type = 2;
					strcat( name, ".lbn" );
					haslbn = TRUE;
				}
				else
				{
					get_out = true;
				}
				break;
			default:
				get_out = true;
				break;
			}
		}
		if ( get_out == true )
		{
			if ( from_ini == 1 ) 
			{
				console.log("\n\nFile: "+name+" from lintel.ini is not available.\n");
				from_ini = -1;
			}
			error = 1;
			name[0] = NULL;
			again = true;
			continue;
		}

		// add directory to filename
		if ( dir_ok == true )
		{
			strcat ( dir, name );
			strcpy ( name, dir );
		}
		console.log("\n    ");
		if ( input_file_type == 0 )
		{
			sprintf ( nudesname, "%s", name );
			if ( ( infile = fopen( nudesname, "r" ) ) == NULL )
			{
				if ( infile ) fclose ( infile );
				console.log("\n\n "+nudesname+" OOPS?\n");
				bell ( 1, 1 );
				if ( from_ini == 1 ) 
				{
					from_ini = -1;
				}
				again = true;
				continue;
			}
			console.log("  Opened "+nudesname+"\n");
		}
		else if ( input_file_type > 0 )
		{
			strcpy( finname, name );

			if ( ( infile = fopen ( finname, "r" ) ) == NULL )
			{
				if ( infile ) fclose ( infile );
				console.log("\n   "+finname+" ?  OOPS - file not found.\n");
				bell ( 1, 1 );
				if ( from_ini == 1 ) 
				{
					console.log("\n\n    File: "+name+" from lintel.ini is not available.\n");
					from_ini = -1;
				}
				again = true;
				continue;
			}

			console.log("\n   opened input file "+finname+"\n");

			add_id_num ( name, nudesname, ".n" );
			if ( ( nudesfile = fopen ( nudesname, "w" ) ) == NULL )
			{
				if ( nudesfile ) fclose ( nudesfile );
				console.log("\n\n "+nudesname+" OOPS?\n");
				bell ( 1, 1 );
				again = true;
				continue;
			}
			console.log("\n   created nudes file "+nudesname+"\n");
		}
		if ( ( infile = fopen(nudesname, "r" ) ) == NULL )
		{
			if ( infile ) fclose ( infile );
			console.log("\n\n "+nudesname+" OOPS?\n");
			bell ( 1, 1 );
			again = true;
			continue;
		}
	}
} /* get_files */
/************************************************/

//bool led_opena ( int min_fps, int max_fps, int min_beats, int max_beats )
var led_opena = function( min_fps, max_fps, min_beats, max_beats )
{
	var get_out;//bool
	get_out = true;
	if ( lbn_fps < min_fps || lbn_fps > max_fps )
	{
		console.log("\n   Oops: fps value is "+lbn_fps+" but must be between "+min_fps+" and "+max_fps+"\n");
		get_out = false;
	}
	if ( lbn_bpm < min_beats || lbn_bpm > max_beats )
	{
		if ( lbn_bpm < 0 )
		{
			console.log("\n   Oops: bpm value missing\n");
		}
		else
		{
			console.log("\n   Oops: bpm value is "+lbn_bpm+" but must be between "+min_beats+" and "+max_beats+"\n");
		}
		get_out = false;
	}
	return( get_out );
} /* led_opena */
/********************************************/

//void led_param ( void )
var led_param = function()
/*
   set up parameters of .lbn interpretation from .ini file

   called by main
*/
{
	var get_out;//bool
	var min_fps;//int
	var max_fps;//int
	var min_beats;//int
	var max_beats;//int
	var lbn_figures_in;//int
	var lbn_default;//bool
	var lbn_fps_in;//int
	var lbn_bpm_in;//int

	lbn_fps = -1;
	lbn_bpm = -1;
	lbn_ppb = 23;
	min_fps = 1;
	max_fps = 250;
	min_beats = 25;
	max_beats = 250;	
	lbn_default = false;
	lbn_figures = 1;
	if ( get_if_ini () == true )
	{
		lbn_figures_in = get_ini_int ( "lbn_figures" );
		lbn_default = get_ini_bool ( "lbn_default" );
		lbn_fps_in = get_ini_int ( "lbn_fps" );
		lbn_bpm_in = get_ini_int ( "lbn_bpm" );

		if ( lbn_fps_in < min_fps || lbn_fps_in > max_fps 
			|| lbn_bpm_in < min_beats || lbn_bpm_in > max_beats )
				lbn_default = false;
		if ( lbn_default == true )
		{
			lbn_fps = lbn_fps_in;
			lbn_bpm = lbn_bpm_in;
			lbn_figures = lbn_figures_in;
		}
	}

	if ( lbn_default == false )
	{
		get_out = false;
		do
		{
			console.log("\n   Please enter frames/second ("+min_fps+"-"+max_fps+")");
			console.log("\n            and beats/minute ("+min_beats+"-"+max_beats+")");
			console.log("\n            separated by a space\n   :");
			if ( gets ( buf ) != NULL && buf[0] != 0 )
			{
				buf = ""+lbn_fps+" "+lbn_bpm+"";//TODO REVIEW (lbn_fps and lbn_bpm both had &)
				get_out = led_opena ( min_fps, max_fps, min_beats, max_beats);
			}
			else
			{
				lbn_fps = 25;
				lbn_bpm = 120;
				console.log("\n   Oops: cannot read fps and bpm");
				console.log("\n   values set to "+lbn_fps+" and "+lbn_bpm+" respectively\n");
				get_out = true;
			}
		} while ( get_out == false );
	}
	lbn_fpp = lbn_fps*doub60 / lbn_bpm*lbn_ppb;//TODO REVIEW (All except doub60 were cast to double)
	
	console.log("\n   frames/pixel "+lbn_fpp+", fps "+lbn_fps+", bpm "+lbn_bpm+", ppb "+lbn_ppb+"\n");
	
	console.log("   number of figures "+lbn_figures+"\n");
}/* led_param */
/************************************************/

//void initgraphics(void) 
var initgraphics = function() 
/*
   called by  main,
*/
{ 
   //char title[BMAX];
   var title = Array();

   title = ptitle + "  -  " + finname;
   glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB); 
   glutInitWindowSize(width, height); 
   glutInitWindowPosition(xw,yw); 
   glutCreateWindow(title);
/* run in full screen if WINDOW_MODE macro undefined */  
//#ifndef WINDOW_MODE 
//   glutFullScreen(); 	//TODO REVIEW (These three lines commented out due to error)
//#endif 
/* background color */ 
   glClearColor(1.0, 1.0, 1.0, 0.5); 

/* init viewing matrix */ 
   glMatrixMode(GL_PROJECTION); 
   glLoadIdentity(); 
   glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0); 
   glEnable(GL_DEPTH_TEST);
} /* initgraphics */
/***************************************/ 

//int _tmain(int argc, _TCHAR* argv[])
var _tmain = function(argc, argv)
/*
   calls initialise, lgetfiles, linter, openfile, compl,
         doframes, initsphere, initgraphics,
         checkeys, image, animate, getout,
         help, gluInit, glutKeyboardFunc, glutDisplayFunc,
         glutIdleFunc, glutMainLoop,
			get_ini, get_files, led_param
*/
{
	ptitle = "lintel084";
	console.log("\n   "+ptitle+" running\n");

	var more = true;
	while(more){
		more = false;
		initialise();
		get_ini ( 0 );
		led_param();
		get_files ( NULL );
		if ( ok != 0 ){
			more = true;
			continue;
		}
		if (haslbn == TRUE)
		{
			 output += "*\n* created "+nudesname+" from "+name+" using "+ptitle+"\n*\n";
			 linter();
		}
		fstart = 0;
		if (ok == 0) openfile(); 
		   else 
			  if (ok != 1) getout(1);
				 else 
				 {
					more = true;
					continue;
				}
		compl();
		if (ok == 0) doframes();
		   else 
			  if (ok != 1) getout(1);
				 else{
						more = true;
						continue;
					}
		if (ok == 0)
		  initsphere();
		   else 
			  if (ok != 1) getout(1);
				 else{
						more = true;
						continue;
					}
		glutInit(argc, argv); //TODO REVIEW (argc had &)
		initgraphics(); 
		console.log("For interactive command list:\n");
		console.log("    click in animation window, then press 'h' key\n");
		glutKeyboardFunc(checkeys); // register Keyboard handler 
		glutDisplayFunc(image);     // register Display handler  
		glutIdleFunc(animate);
		glutMainLoop();

		more = true;
		continue;
	}
}
//END PORT ON 2013-12-11 (Errorage)


















































