// ini file variables-
var ini_title = get2DArray(256); // [256][32]
var ini_value = get2DArray(256); // [256][128]
var max_ini = 256;
var max_ini_len = 32;
var number_ini;
var ini_diag = 0;
var input_file_type;
var lbn_fps = -1;
var lbn_bpm = -1;
var lbn_ppb = 23;
var lbn_figures = 1;     // number of staves to be varerpreted
var time;
var nextcount = 0;
var output_file_id; // boolean

// symbol menus-
var Area ='A';
var Bars ='B';
var Dirn ='D';
var Face ='F';
var Path ='H';
var Keys ='K';
var Limb ='L';
var Misc ='M';
var Pins ='P';
var Rotn ='R';
var Stav ='S';
var Volm ='V';

var forbid = get2DArray(EMAX); // [EMAX][EMAX] // booleans

// doubles
var  doub0;
var  doub1;
var  doub2;
var  doub3;
var  doub4;
var  doub10;
var  doub60;
var  doub90;
var  doub150;
var  doub179;
var  doub180;
var  doub181;
var  doub255;
var  doub360;
var  doub500;
var  doubmax;
var  inv2;
var  inv3;
var  inv4;
var  inv5;
var  inv6;
var  inv10;
var  inv256;
var  inv1000;
var  lg2;                   //logarithm of 2 
var  rt3;
var  tolr;
var  twopi;
var  pi;                    // 3.142...etc /
var  piby2;
var  degree;                // number of degrees in a radian 
var  radian;                // number of radians in a degree
var  rad;                   // conversion factor from degrees to radians 
var  radten;                // conversion factor from tenths of a degree to radians 
var  alpha;                 // basic interactive angle increment
var  anglex,angley,anglez;  // interactive rotation angles 
var  dangx,dangy,dangz;     // interactive rotation angle increments 
//var  fac;                   // lbn conversion factor from y to frames 
var  scale;                 // interactive scaling factor 
var  SCALE = 1.0;           // default scaling to fit window
var  tx,ty,tz;              // interactive translations
//var  x1a,x1b,x2a,x2b;
//var  y1a,y1b,y2a,y2b;
//var  x1s,x2s,y1s,y2s;
var  frac;                  // fraction of action to be done in this frame
var  prop;                  // proportion of action yet to be done 
var  step1,step2;
var  v;
var  varval;                // varval - value of variable
var  ang = [];              // the 3 eulerian angles [3]
var  oldang = [];			// [3]
var  obsang = [];			// [3]
var  factor = [];           // factors in x,y, and z directions [3]
var  lighting_rgb = [];		// [3]
var  pplace = [];           // position of centre of observers attention [3]
var  semiax = [];           // coordinate triple read from input [3]
var  xx = [];               // x,y and z values (x&y used for commands add, subtract,multiply,divide) [3]
var  val = [];              // val[i] - if i <= nvals then value of ith constant
                            //           otherwise (s - i+1)th variable [VMAX]
var  maxax = [];            // maxax[j] - largest semiaxis of jth ellipsoid [EMAX]
var  minax = [];
var  pres = [];;
var  obs = get2DArray(3);              // rotation matrix applied to all ellipsoids to obtain observers view  [3][3]
var  col = get2DArray(EMAX);              // col[i][j] - colour of ell i; j = 0 red, 1 green, 2 blue [EMAX][3]
var  quasav = get2DArray(EMAX+2);			// [EMAX+2][5]
var  ax = get2DArray(EMAX);               // ax3[i][j] - length of jth semiaxis of ith ellipsoid [EMAX][3]
var  cen = get2DArray(EMAX);              // cen[i][j] - jth coordinate of centre of ith ellipsoid [EMAX][3]
var  censav = get2DArray(EMAX);			// [EMAX][3]
var  dcon = get3DArray(EMAX, 2);  		    // distances of joint from ellipsoid centres [EMAX][2][3]
var  jnt = get2DArray(EMAX);              // coordinates of joints between ellipsoids [EMAX][3]
var  jntsav = get2DArray(EMAX);			// [EMAX][3]
var  ob3 = get2DArray(FMAX);       		// observation angles of each frame [FMAX][3]
var  pl3 = get2DArray(FMAX);     		    // centre of view of each frame [FMAX][3]
var  norm = get2DArray(4*SMAX*SMAX);  			// normals at sphere faces [4*SMAX*SMAX][3]
var  sph = get3DArray(4*SMAX*SMAX, 4);				// vertices of facets of  sphere [4*SMAX*SMAX][4][3]
var  lim = get3DArray(EMAX, 3);				// [EMAX][3][2]
var  co3 = get3DArray(FMAX, EMAX);	    		// colours of ellipsoids [FMAX][EMAX][3]
var  ce3 = get3DArray(FMAX, EMAX);    			// coordinates of ellipsoid centres [FMAX][EMAX][3]
var  ax3 = get3DArray(FMAX, EMAX);    			// ellipsoid semiaxis lengths [FMAX][EMAX][3]
var  qu3 = get3DArray(FMAX, EMAX);    			// quaternions of ellipsoids [FMAX][EMAX][4]
var  quat = get2DArray(EMAX+2);       		// quat[i][j] - angle of ith ellipsoid [EMAX+2][5]
						    // j = 0,1,2 - components of direction of rotation axis  
						    // j = 3,4   - sine and cosine of rotation about axis 
var point = get3DArray(SMAX, 2*SMAX+1);				// [SMAX][2*SMAX+1][3]

// Structs
function Symbol() {
	var a;       // TRUE = 0 if already done
	var b;       // bent indicator
	var c;       // column position relative to right support column
	var h;       // height
	var i;       // item in menu
	var l;       // lbn file line number
	var m;     	 // menu (char)
	var s;       // drawing step size
	var w;       // width
	var x;       // horizontal position of left side
	var y;       // vertical position of bottom
	var x2;      // horizontal position of right side
	var y2;      // vertical position of top
	var d;       // height indicator
   };

var lbn = [];	 // [LMAX] of Symbol
var ja;                  // TRUE = 0 if already done
var jb;                  // bendedness of current symbol
var jc;                  // current symbol column
var jh;                  // current symbol height
var ji;                  // current symbol item in menu
var jl;                  // line of current symbol
var jm;                 // current symbol menu
var js;                  // current symbol step size
var jw;                  // current symbol width
var jx;                  // current symbol x bottom
var jy;                  // current symbol y bottom 
var jx2;                 // current symbol x top
var jy2;                 // current symbol y top
var jd;                  // current symbol shading

/*************************************************************/

// linter variables -
var mspace,wspace;

var lbn_fpp;      // frames per pixel

var blength;         // number of bars to interpret
var bpm;             // beats per minute
var bstart;          // bar to start at
var lcentre;         // x position of centre staff line
var complete;        // true if Gloria and Frank to be used
var dofig;           // required gender of current staff
var dostaff;         // index in staff[] of current staff
var facedif;         // difference between facing directions of man and lwoman
var facecl;          // facing score of closed position
var facepr;          // facing score of promenade position
var facesh;          // facing score of shadow position
var facess;          // facing score of semishadow position
var fbegin,ffin,flen;// start,end, and length of a position
var fend;            // frame number of end of current movement
var fhalf;           // frame halfway through a movement
var fmax;            // maximum frame number
var fps;             // frames/second
var frange;          // number of frames in an action
var frperbar;        // frames per bar;
var fstart;          // first frame number of action
var gy,gh;           // arm gesture range disabled by contact bow
var haslbn;          // TRUE if input is lbn file, FALSE for .n file
var hold;            // one of the defined holds NO,CL,PR,CP,DB,OP,CR,OE,CO,SH,SS
var holdcl;          // closed hold counter 
var holdco;          // counter open extended hold counter
var holdoe;          // open extended hold counter
var holdpr;          // promenade hold counter
var holdsh;          // shadow hold counter
var holdss;          // semishadow hold counter
var j;               // counter through symbols
var keptf;           // last frame when last position kept
var mface,wface;     // facing directions of man and woman
var nbar;            // number of current bar
var nlabs;           // number of laban score entries
var npins;           // number of pins below first barline
var nm;              // number of men
var nw;              // number of women
var nmw;             // nm * nw
var nstaff;          // number of staves
var oriented;        // true after orientation
var pend;            // last frame of previous action
var pstart;          // first fame of previous action
var ppb;             // pixels per beat (= 23 );
var prev_time;       // clock reading of previous frame
var pres_time;       // clock reading of current frame
var previ;           // item of previous support symbol
var prevc;           // column of previous support symbol
var prevhold;        // previous hold
var rise;            // height of previous step;
var ssend;           // ending score symbol
var sstart;          // starting score symbol
var st;              // current staff number
var stmiddle;        // halfway across between L and R staves
var track;           // TRUE when tracking viewpoint on main figure
var xmin,xmax;       // width range of score symbols
var ymax;            // top of score
var yend;            // y position of last movement
var ystart;          // y position of start of movement
var yj = [];      	 // symbols starting at given y positions [5*FMAX]
var pins = get2DArray(TMAX);   	 // index and use of initial pins [TMAX][2]
var staff = get2DArray(TMAX);  	 // index, x positions, gender, and use of staves [TMAX][6]
var colm = [];    	 // limb presigns in the columns [NCOLM]

//nudes variables -
var axis;            // axis of next rotation
var bnums;           // TRUE if bar numbers to be displayed 
var comand;          // counter through all commands.
var df;              // interactive frame increment 
var ecount;          // number of entries in 'elist'
var ell1;            // ellipsoid to touch something  
var ell2;            // ellipsoid to be touched 
var ellpsd;          // active ellipsoid
var f;               // counter through frames 
var fast;            // multiplier of frame numbers
var fig;             // current figure 
var fnums;           // TRUE if frame numbers to be displayed 
var forward;         // TRUE for animation to go forwards 
var freeze;          // TRUE if animation frozen 
var fstop;           // last frame number of actions
var fslow;
var height = 512;    // height  of window in pixels
var hstart;          // frame at start of hold
var hend;            // frame at end of hold
var inmain;          // TRUE if still in main NUDES program 
var intersect;
var jcount ;
var join;            // joint for current bend command  
var k;
var length;          // length of next input string 
var lline;           // length of next input line 
var maxint;          // largest representable integer
var more;            // if > 0 means more actions for which stp>=fr 
var ne;              // number of ellipsoids in current frame
var nesave;
var nfaces;          // number of faces on sphere 
var nfigs;           // number of figures
var nfiles;          // number of texture map files
var njts;            // number of joints 
var nline;           // number of current nudes file line
var npfs;            // number of actions 
var nsph;            // number of chords around sphere 
var nsubs;           // number of subroutines 
var nvals;           // number of values in 'val' 
var nvars;           // number of variables in array val 
var ok;              // ok = 0 if ok, else problem reference number 
var p;               // counter through actions
var pause;           // TRUE if pausing on 1st and last frames 
var pok;             // true if positive integer read 
var prdone;          // TRUE if diagnostic printing already done
var ptype;           // code of current action 
var pp;
var donesurf;        // TRUE if 'surf' called from 'dotouch'
var refell;          // ellipsoid used as angular reference 
var shadow;          // TRUE if shadows wanted
var single;          // either TODO or DONE when frozen
var slow;            // number of pause calls between animating frames 
var start;           // pointer to next character on current input line 
var t;               // type of current action 
var var0;
var var1;
var var2;
var vstart;          // first frame from view command
var vstop;           // last frame from view command
var width = 512;     // height  of window 
var xw = 10;
var yw = 10;         // lower left corner of window 
var newcol = [];	 // [3]
var axlen = [];      // lengths of names [EMAX]
var ellen = [];		 // [EMAX]
var figlen = [];	 // [EMAX]
var jntlen = [];	 // [EMAX]
var fillen = [];	 // [EMAX]
var keylen = [];	 // [NKEYS]
var sublen = [];	 // [PMAX]
var varlen = [];	 // [PMAX]
var called = [];     // true if subroutine is called [PMAX]
var cline = [];      // line numbers in input file of each action [PMAX]
var coel = get2DArray(EMAX);		 // the 2 ellipsoids joined at a joint [EMAX][2]
var defined = [];    // TRUE if subroutine is defined [PMAX]
var distrn = [];     // how actions are distributed over frames [PMAX]
var ellfig = [];     // number of the figure containing each ellipsoid [EMAX]
var elist = [];      // array for lists of ellipsoids in current action [EMAX]
var figell = [];     // figell[i] - first ellipsoid in ith figure [EMAX]
var frames = [];     // original NUDES frame numbers [FMAX]
var frstart = [];    // frstart[i] - frame number of start of ith action [PMAX]
var frstop= [];      // frstop[i] - frame number of end of ith action  [PMAX]
var jlist = [];      // array for lists of joints in current action [EMAX]
var knee = [];       // knee[j] - true if jth joint is a knee i.e. flexes backwards [EMAX]
var nels = [];       // number of ellipsoids in each frame [FMAX]
var type = [];       // type of  action [PMAX]
var pf = get2DArray(PMAX);         // pf[i][j] - jth parameter of ith action-  +ve: itself, -ve: index into array val [PMAX][6]
var subact = get2DArray(PMAX);     // subact[i][] - action numbers of start and end of ith subroutine [PMAX][2]
var usevar = [];     // 0 if variable not used [PMAX]
var order = [		 // [3][3][3]
      [ [2,1,1],[1,3,4],[1,5,3] ],
      [ [3,1,5],[1,2,1],[4,1,3] ],
      [ [3,4,1],[5,3,1],[1,1,2] ]];
var perm = [		 // [3][3][3]
      [ [2,1,1],[1,3,4],[1,5,3] ],
      [ [3,1,5],[1,2,1],[4,1,3] ],
      [ [3,4,1],[5,3,1],[1,1,2] ]];
	  
/*
   keyword codes -
*/
var figure_keyword_code=  1;
var ellips_keyword_code=  2;
var joint_keyword_code=  3;
var accele_keyword_code=  5;
var subrou_keyword_code=  6;
var balanc_keyword_code=  7;
var attach_keyword_code=  8;
var detach_keyword_code=  9;
var decele_keyword_code= 10;
var grofig_keyword_code= 11;
var spinto_keyword_code= 12;
var moveby_keyword_code= 13;
var add_keyword_code= 14;
var touch_keyword_code= 15;
var stop_keyword_code= 16;
var spinby_keyword_code= 17;
var ground_keyword_code= 18;
var bendby_keyword_code= 19;
var set_keyword_code= 20;
var bendto_keyword_code= 21;
var dodebug_keyword_code= 22;
var repeat_keyword_code= 23;
var quadra_keyword_code= 24;
var linear_keyword_code= 25;
var observ_keyword_code= 26;
var moveto_keyword_code= 27;
var call_keyword_code= 28;
var endsub_keyword_code= 29;
var speed_keyword_code= 30;
var invert_keyword_code= 31;
var variable_keyword_code = 32;
var view_keyword_code= 33;
var groell_keyword_code= 34;
var grojnt_keyword_code= 35;
var angles_keyword_code= 36;
var centre_keyword_code= 37;
var flex_keyword_code= 38;
var rotate_keyword_code= 39;
var abduct_keyword_code= 40;
var negate_keyword_code= 41;
var subtra_keyword_code= 42;
var divide_keyword_code= 43;
var multip_keyword_code= 44;
var cubic_keyword_code= 46;
var place_keyword_code= 47;
var axes_keyword_code= 48;
var linkx_keyword_code= 49;
var colour_keyword_code= 50;
var print_keyword_code= 51;
var textur_keyword_code= 52;
var drag_keyword_code= 53;
var limits_keyword_code= 54;
var abut_keyword_code= 55;
var movjnt_keyword_code= 56;
var growto_keyword_code= 57;
var color_keyword_code= 58;
var center_keyword_code= 59;
var opacty_keyword_code= 60;
var lghtng_keyword_code= 61;
var allow_keyword_code= 62;
var forbid_keyword_code= 63;

var infile;
var nudesfile;
var figsfile;

var junk = [];

var buf;            	  // input buffer
var line;           	  // compl input buffer 
var lbnline = [];   	  // lbn file lines [LMAX]
var string;         	  // next set of non-blank characters from data file */
var name;           	  // name of input file
var finname;          	  // name of input file
var figsname;       	  // name of lintel nudes figures, declarations, and subroutines file */
var nudesname;      	  // name of intermediate nudes file
var ptitle;         	  // program title
var risesub = ['flow',
			   'fmed',
			   'fhigh'];
var xyz = ['mx my mz',
		   'wx wy wz'];
var aline = [];    // nudes input lines [PMAX]
var tname = [];    // name of texture map file [EMAX]
var jname = [];    // joint names [EMAX]
var sname = [];    // subroutine names [EMAX]
var vname = [];    // variable names [EMAX]
var axnam = [];    // first entry is the set of axis names 'x','y','z'. The rest are null [EMAX]
var ename = [];    // ellipsoid names [EMAX]
var fname = [];    // figure names [EMAX]
var tn3 = get2DArray(FMAX);// names of reduced texture map files [FMAX][EMAX]
var nullChar  = '\0'; // !!! null
var blank = ' ';
var dig = ['0','1','2','3','4','5','6','7','8','9','*'];
var dummy = 'dummy';
var every = 'every';
var nudes = 'nudes';
var world = 'world';
var variab = 'variab';
var expect = [' ',
			  'value',
			  'ellipsoid',
			  'joint',
			  'figure',
			  'axis',
			  'subroutine',
			  'variables',
			  'string'];

/*
   par[p,k] - the type of the kth parameter of the pth action -
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
*/
var par =  [ // [NKEYS][6]
      [0,0,0,0,0,0],//   0
      [0,0,0,0,0,0],//   1  figure
      [0,0,0,0,0,0],//   2 ellips
      [0,0,0,0,0,0],//   3 joint
      [0,0,0,0,0,0],//   4 
      [1,0,0,0,0,0],//   5 accele
      [0,0,0,0,0,0],//   6 subrou 
      [2,3,2,5,0,0],//   7 balanc
      [2,3,2,1,1,1],//   8 attach
      [2,3,4,0,0,0],//   9 detach
      [1,0,0,0,0,0],//  10 decele
      [4,2,1,1,1,0],//  11 grofig
      [4,2,2,1,1,1],//  12 spinto
      [4,2,1,1,1,0],//  13 moveby
      [7,1,1,0,0,0],//  14 add
      [2,2,2,2,3,5],//  15 touch
      [0,0,0,0,0,0],//  16 stop
      [4,2,2,1,5,0],//  17 spinby
      [4,0,0,0,0,0],//  18 ground
      [2,3,2,1,5,0],//  19 bendby
      [7,8,0,0,0,0],//  20 set
      [2,3,2,1,1,1],//  21 bendto
      [1,0,0,0,0,0],//  22 dodebug
      [1,0,0,0,0,0],//  23 repeat
      [1,0,0,0,0,0],//  24 quadra
      [1,0,0,0,0,0],//  25 linear
      [1,1,1,0,0,0],//  26 observ
      [4,2,1,1,1,0],//  27 moveto
      [8,0,0,0,0,0],//  28 call
      [0,0,0,0,0,0],//  29 endsub
      [0,0,0,0,0,0],//  30 speed
      [7,0,0,0,0,0],//  31 invert
      [0,0,0,0,0,0],//  32 variable
      [0,0,0,0,0,0],//  33 view
      [2,1,1,1,0,0],//  34 groell
      [2,3,1,1,1,0],//  35 grojnt
      [2,2,7,7,7,0],//  36 angles
      [2,7,7,7,0,0],//  37 centre
      [2,3,1,0,0,0],//  38 flex
      [2,3,1,0,0,0],//  39 rotate
      [2,3,1,0,0,0],//  40 abduct
      [7,0,0,0,0,0],//  41 negate
      [7,1,1,0,0,0],//  42 subtra
      [7,1,1,0,0,0],//  43 divide
      [7,1,1,0,0,0],//  44 multiply
      [0,0,0,0,0,0],//  45
      [0,0,0,0,0,0],//  46 cubic
      [1,1,1,0,0,0],//  47 place
      [2,7,7,7,0,0],//  48 axes
      [3,7,7,7,0,0],//  49 linkx
      [2,1,1,1,0,0],//  50 colour
      [7,0,0,0,0,0],//  51 print
      [2,9,1,1,0,0],//  52 texture
      [2,2,3,2,5,0],//  53 drag
      [0,0,0,0,0,0],//  54 limits
      [2,2,2,5,0,0],//  55 abut
      [3,2,1,1,1,0],//  56 movjnt
	  [4,2,1,1,1,0],//  57;
	  [2,1,1,1,0,0],//  58;
	  [2,7,7,7,0,0],//  59;	
	  [2,1,0,0,0,0],//  60;
	  [1,1,1,0,0,0],//  61;
	  [2,2,0,0,0,0],//  62;
	  [2,2,0,0,0,0]];// 63;

var keynam = [ // [NKEYS][BMAX]
      'keyword',            // 0
      'figure',             // 1
      'ellipsoid',     		// 2
      'joint',              // 3
      'copy',               // 4
      'accelerate', 		// 5
      'subroutine', 		// 6
      'balance',            // 7
      'attach',             // 8
      'detach',             // 9
      'decelerate', 		// 10
      'grofig',             // 11
      'spinto',             // 12
      'moveby',             // 13
      'add',                // 14
      'touch',              // 15
      'stop',               // 16
      'spinby',             // 17
      'ground',             // 18
      'bendby',             // 19
      'set',                // 20
      'bendto',             // 21
      'debug',              // 22
      'repeat',             // 23
      'quadratic',     		// 24
      'linear',             // 25
      'observe',            // 26
      'moveto',             // 27
      'call',               // 28
      'endsub',             // 29
      'speed',              // 30
      'invert',             // 31
      'variables',     		// 32
      'view',               // 33
      'groell',             // 34
      'grojnt',             // 35
      'angles',             // 36
      'centre',             // 37
      'flex',               // 38
      'rotate',             // 39
      'abduct',             // 40
      'negate',             // 41
      'subtract',         	// 42
      'divide',             // 43
      'multiply',         	// 44
      'read',               // 45
      'cubic',              // 46
      'place',              // 47
      'axes',               // 48
      'linkx',              // 49
      'colour',             // 50
      'print',              // 51
      'texture',            // 52
      'drag',               // 53
      'limit',              // 54
      'abut',               // 55
      'movjnt',				// 56
	  'growto',             // 57
	  'color',              // 58
	  'center',             // 59
	  'opacty',             // 60
	  'lghtng',             // 61
	  'allow',              // 62
	  'forbid' ];           // 63
	  /*
     code[p,k] - type of kth parameter of pth action using -

     0-illegal
     1-x coordinate
     2-y coordinate
     3-z coordinate
     4-angle 1
     5-angle 2
     6-angle 3
     7-x scaling factor
     8-y scaling factor
     9-z scaling factor
    10-value for a variable
    11,12,13-red green and blue colour coords, respectively,
               or image texture file reference ,xoffset and yoffset
    14-debug parameter

    21-axis
    22-joint
    23-reference ellipsoid
    24-moving or central ellipsoid
    25-figure
    27,28,29-names of variables
    30-touching or dragged ellipsoid (ell1)
    31-touched ellipsoid (ell2)
*/


var code = [              // [NKEYS][6]
	[0,0,0,0,0,0],      // 0
	[0,0,0,0,0,0],      // 1
	[0,0,0,0,0,0],      // 2
	[0,0,0,0,0,0],      // 3
	[0,0,0,0,0,0],      // 4
	[0,0,0,0,0,0],      // 5
	[0,0,0,0,0,0],      // 6
	[24,22,23,21,0,0],  // 7
	[24,22,23,1,2,3],   // 8
	[24,22,25,0,0,0],   // 9
	[0,0,0,0,0,0],      // 10
	[25,24,7,8,9,0],    // 11
	[25,24,23,4,5,6],   // 12
	[25,23,1,2,3,0],    // 13
	[27,1,2,0,0,0],     // 14
	[30,31,24,23,22,21],// 15
	[0,0,0,0,0,0],      // 16
	[25,24,23,4,21,0],  // 17
	[25,0,0,0,0,0],     // 18
	[24,22,23,4,21,0],  // 19
	[27,10,0,0,0,0],    // 20
	[24,22,23,4,5,6],   // 21
	[14,0,0,0,0,0],     // 22
	[0,0,0,0,0,0],      // 23
	[0,0,0,0,0,0],      // 24
	[0,0,0,0,0,0],      // 25
	[4,5,6,0,0,0],      // 26
	[25,24,1,2,3,0],    // 27
	[0,0,0,0,0,0],      // 28
	[0,0,0,0,0,0],      // 29
	[0,0,0,0,0,0],      // 30
	[27,0,0,0,0,0],     // 31
	[0,0,0,0,0,0],      // 32
	[0,0,0,0,0,0],      // 33
	[24,7,8,9,0,0],     // 34
	[24,22,7,8,9,0],    // 35
	[24,23,27,28,29,0], // 36
	[24,27,28,29,0,0],  // 37
	[24,22,4,0,0,0],    // 38
	[24,22,4,0,0,0],    // 39
	[24,22,4,0,0,0],    // 40
	[27,0,0,0,0,0],     // 41
	[27,1,2,0,0,0],     // 42
	[27,1,2,0,0,0],     // 43
	[27,1,2,0,0,0],     // 44
	[0,0,0,0,0,0],      // 45
	[0,0,0,0,0,0],      // 46
	[1,2,3,0,0,0],      // 47
	[24,27,28,29,0,0],  // 48
	[22,27,28,29,0,0],  // 49
	[24,11,12,13,0,0],  // 50
	[27,0,0,0,0,0],     // 51
	[24,11,12,13,0,0],  // 52
	[30,24,22,23,21,0], // 53
	[0,0,0,0,0,0],      // 54
	[30,31,23,21,0,0],  // 55
	[22,24,1,2,3,0],    // 56
	[0,0,0,0,0,0],      // growto 57;
	[24,11,12,13,0,0],  // color 58;
	[24,27,28,29,0,0],  // center 59;
	[30,31,0,0,0,0],    // opacity 60;
	[30,31,0,0,0,0],    // lighting 61;
	[30,31,0,0,0,0],    // allow 62;
	[30,31,0,0,0,0]];   // forbid 63;


/****************************************/

var menutext = [	// [NSYMS][4]
     'Bars',
     'Dirn',
     'Pins',
     'Face',
     'Limb',
     'Volm',
     'Area',
     'Rotn',
     'Keys',
     'Misc',
     'Ways',
     'ZZZZ',
     'Stav'];

var leg = 			//[12][3]   /* quaternion angles of 11 direction symbols */
                   	/* for walking */
         [[  0,  0,  0],
          [  0,  0, 30],
          [  0,315, 30],
          [  0,270, 30],
          [  0,225, 30],
          [  0,180, 30],
          [  0,180, 30],
          [  0,135, 30],
          [  0, 90, 30],
          [  0, 45, 30],
          [  0,  0, 30],
          [  0,  0,  0]];

var opp			  // [12]       /* opposite direction to a movement */
                  /* for the leg that is left behind */
                  =    [  0,  5,  7,  8,  9,  1,  1,  2,  3,  4,  5,  0];

var   stt           // [3][12][3] /* quaternion angles of 11 direction symbols */
                    /* for straight limbs */
                    /* at hip or shoulder */
     =  [[[  0,  0,  0], // null
			 [  0,  0, 35], //   1 R forward low
			 [  0,315, 35], //   2 R diagonally forward low
			 [ 68, 90,-35], //   3 R side low
			 [  0,225, 35], //   4 R diagonally back low
			 [  0,180, 35], //   5 R back low
			 [  0,180, 35], //   6 L back low
			 [  0,135, 35], //   7 L diagonally back low
			 [ 68, 90, 35], //   8 L side low
			 [  0, 45, 35], //   9 L diagonally forward low
			 [  0,  0, 35], //  10 L forward low
			 [  0,  0,  0]],//  11 in place low
			[[  0,  0,  0], // null
			 [  0,  0, 90], //   1 R forward middle
			 [337,339, 98], //   2 R diag forward middle
			 [315,270, 90], //   3 R side middle
			 [158,339, 98],
			 [180,  0, 90],
			 [180,  0, 90],
			 [202, 21, 98],
			 [114, 90, 90],
			 [ 22, 21, 98],
			 [  0,  0, 90],
			 [  0,  0,  0]],
			[[  0,  0,  0], // null
			 [  0,  0,135], // R forward high
			 [  0,315,135],
			 [  0,270,135],
			 [  0,225,135],
			 [  0,180,135],
			 [  0,180,135],
			 [  0,135,135],
			 [  0, 90,135],
			 [  0, 45,135],
			 [  0,  0,135],
			 [  0,  0,180]]];
var   trlx            // [3][12][3]  /* quaternion angles of 11 direction symbols */
                      /* for relaxed thighs */
     =  [[[  0,  0,  0], // null
			 [  0,  0, 55], // R low forward
			 [321,342, 70],
			 [180,270, 44], // R low side
			 [257,310, 60],
			 [180,  0, 30], // R low back
			 [180,  0, 30], // L low back
			 [ 99, 46, 58],
			 [ 87, 90, 44], // L low side
			 [ 42, 17, 66],
			 [  0,  0, 50],
			 [  0,  0,  0]],
			[[  0,  0,  0],
			 [  0,  0,110],
			 [341,339,107],
			 [ 45,277, 90],
			 [225,294, 95],
			 [180,  0, 80],
			 [180,  0, 80],
			 [135, 66, 95],
			 [315, 83, 90],
			 [ 19, 21,107],
			 [  0,  0,100],
			 [  0,  0,  0]],
			[[  0,  0,  0],
			 [  0,  0,145],
			 [352,338,148],
			 [ 68,275,135],
			 [192,291,135],
			 [180,  0,125],
			 [180,  0,125],
			 [167, 69,134],
			 [292, 85,135],
			 [  7, 22,148],
			 [  0,  0,145],
			 [  0,  0,180]]];
var   arlx                 // [3][12][3]  /* quaternion angles of 11 direction symbols */
     =  [[[  0,   0,   0], // null
          [  0,   0,  35], // R low forward
          [307, 346,  56],
          [287, 335, 102], // R low side
          [141, 342,  70],
          [180,   0,  50],
          [180,   0,  50],
          [222,  17,  66],
          [ 72,  25, 102], // L low side
          [ 57,  12,  54], // L low diag forward
          [  0,   0,  35], // L low forward
          [  0,   0,  0]], // L low centre
         [[  0,   0,   0], // null
          [  0,   0,  67],
          [328, 341,  79],
          [303, 331, 108], // R middle side
          [285, 328, 143],
          [270, 327, 179], // R middle back
          [ 90,  33, 179], // L middle back
          [ 75,  33, 143],
          [ 57,  29, 108], // L middle side
          [ 32,  19,  79],
          [  0,   0,  67], // L middle front
          [315,  60,  98]],// L middle centre
         [[  0,  0,  0],
          [  0,  0,145],
          [348,338,129],
          [247,275,135],
          [209,298,142],
          [ 90, 72,180],
          [ 90, 72,180],
          [150, 62,143],
          [112, 85,135],
          [ 12, 22,129],
          [  0,  0,125],
          [350,-45,166]]];

var   abnt            // [3][12][3]  /* quaternion angles of 11 direction symbols */
                      /* for 90 degree bent arms */
     =  [[[252 , 17,  95],
          [  0,   0,   0],
          [270,   0,  45],
          [270,   0,  90], // R low side
          [270,   0, 135],
          [270,   0, 180], // R low back 
          [ 90,   0, 180], // L low back
          [ 90,   0, 135],
          [ 90,   0,  90],
          [ 90,   0,  45],
          [  0,   0,   0], // L low forward
          [  0,   0,   0]],
         [[  0,   0,   0],
          [ 45, 300,  98],
          [ 67,  21,  98], // R middle side
          [225, 300,  98],
          [  0, 225, 135],
          [  0, 180, 135],
          [  0, 180, 135],
          [  0, 135, 135],
          [135,  60,  98],
          [292, 339,  98], // L middle side
          [315,  60,  98],
          [  0,   0,   0]],
         [[  0,   0,   0],
          [  0, 315, 135],
          [  0, 270, 135],
          [ 45, 225, 135],
          [ 90, 225, 135],
          [135, 225, 135],
          [135, 225, 135],
          [ 90,  45, 135],
          [ 45, 225, 225],
          [  0,  90, 135],
          [  0,  45, 135],
          [ 90,  67, 180]]];