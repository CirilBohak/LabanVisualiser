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
