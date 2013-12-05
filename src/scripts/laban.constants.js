var NULL = null;
var TRUE = 0;
var FALSE = 1;
var DONE = 0;
var TODO = 1;
var MAN = 0;
var WOMAN = 1;
var MAXINT = 1073741824;
var WINDOW_MODE;
var GLUT_KEY_ESCAPE = 27;
var BMAX = 256;         // size of character buffer
var EMAX = 1024;        // maximum number of ellipsoids
//var FACTOR    2         // number of y pixels per frame
var FMAX = 2048;        // maximum number of frames
var NKEYS = 64;         // number of keywords
var PMAX = 6000;        // maximum number of actions to perform
var SMAX = 100;         // maximum number of chords around sphere
var SMIN = 2;           // minimum number of chords around sphere
var SSTART = 20;        // initial number of chords around sphere

var LMAX = 5000;        // max number of laban score entries 
var TMAX = 30;          // max number of staff lines
var VMAX = 2048;        // max number of constants + variables
var NCOLM = 18;         // number of columns around staff
var STEP = 12;          // spacing of symbols 
var WIDTH = 1024;       // width of the score 
var NSYMS = 25;         // max number of items in each menu

var RELAX = 1;          // item number of 'relaxed' symbol
var BENT = 3;           // item number of 'bent' symbol
var STRAIGHT = 2;       // item number of 'straight' symbol
var STRETCH = 4;        // item number of 'stretched' symbol      
var FRONT = 100;        // front symbol found
var BACK = 200;         // back symbol found
var MLHAND = 1;         // man's left hand symbol found
var MRHAND = 2;         // man's right hand symbol found
var WLHAND = 10;        // woman's left hand symbol found
var WRHAND = 20;        // woman's right hand symbol found
var ARM = 'a';          // arm found in colm[]
var CHEST = 'c';        // chest found in colm[]

var LOW = 0;
var MED = 1;
var HIGH = 2;
var BLANK = 3;

var NO = 0;             // no hold
var CL = 1;             // closed hold: normal ballroom dancing position.
var PR = 2;             // promenade position: facing partner, bodies touching,
                        // but both prepared to travel to man's L.
var CP = 3;             // counter promenade position: facing partner, bodies touching,
                        // but both prepared to travel to man's R.
var DB = 4;             // double hold: facing partner, bodies apart,
                        // L hand to R, R hand to L.
var OP = 5;             // open hold: facing partner, bodies apart,
                        // man's L to lady's R, other hands free.
var CR = 6;             // crossed open hold: facing partner, bodies apart,
                        // man's R to lady's R, other hands free.
var OE = 7;             // open extended hold: both facing same way, bodies apart,
                        // man's L hand to lady's R, other hands free.
var CO = 8;             // counter open extended hold: both facing same way, bodies apart,
                        // man's R hand to lady's L, other hands free.
var SH = 9;             // shadow hold: both facing same way, bodies touching,
                        // L hand to L, R hand to R.
var SS = 10;            // semi-shadow hold: both facing same way, bodies touching, 
                        // man's L hand to lady's L,
                        // man's R hand on lady's R hip, lady's R hand free.
