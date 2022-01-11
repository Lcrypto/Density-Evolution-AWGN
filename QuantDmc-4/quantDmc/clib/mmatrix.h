
/* MWORD is set to be a double */
#define MMAXFLOAT    DBL_MAX   /* previously used this const: 3.40282347E+38F */
#define MMINFLOAT    DBL_MIN   /* 1.17549435E-38F */
typedef double MWORD;

#define EPS          2.7105053E-19F   /*guardband by 1E+1*/

#define MEX 1

#if MEX==1
#include "mex.h"
#endif

#define VARNAMELEN 30           /*maximum length for variable names*/


typedef struct {
    int nr;             /* number of rows */
    int nc;             /* number of colums */
    MWORD **elt;        /* pointer to rows */
    MWORD *cols;        /* pointer to contents of matrix */
    char name[VARNAMELEN+1];
    int type;           /* 0 for MWORDs, 1 for MWORD interpretted as strings */
    int init;           /* 0 if rows/contents not yet malloc'd, 1 otherwise */
} MATRIX;

/************************************************
    Initialize a matrix.  Unless a matrix in input from matlab,
    the declaration should look like this:
    MATRIX *x = matrix("x",0,0);    //for outputs
    MATRIX *x = matrix("x",nr,nc);   //for inputs
    x->elt[col][row] = 10;
    
    matrix_free(x);
    
    An alternate usage is:
    MATRIX dd = {1000,1,NULL,NULL,"d",0};
    MATRIX *d = &dd;
    
    An exception is when the matrix is an input from Matlab.
    The function matlab_to_matrix() invokes the 
    function matrix():
    
    MATRIX *a;

    a = matlab_to_matrix(prhs[0]);
    
*************************************************/

/************************************************
   *matrix() as a minimum returns a pointer to a that
   has MATRIX malloc'd is named.  If nc and nc are both
   positive, then space for the matrix contents will
   be malloc'd and initialized to zero

   matrix_init_ptr() takes a matrix that has already
   been MATRIX malloc'd (by *matrix("",0,0)) and will then
   malloc the contents and initialize them to zero.

   *matrix() is called when the matrix is first
   initialized, and matrix_init_ptr() is called from
   within functions that must produce a matrix as
   an output.

   *matrix can be invoked to initialize data in the following manner:
   MATRIX *data=matrix("*i data",4,1, 1,1,0,1);  OR
   MATRIX *data=matrix("*m data",4,1, 1.0,1.0,0,1.0);
   where the true "name" begins with the third character if
   the given name begins with "*i " for integers or "*m " for
   MWORDs.
*************************************************/

/*MATRIX *matrix(char *name,int nr,int nc);*/
MATRIX *matrix(char *name,int nr,int nc,...);
void matrix_init_ptr(MATRIX *m, int nr,int nc);
void matrix_free(MATRIX *m);
void matrix_null(MATRIX *m);

/************************************************
    returns 1 if matrix has been initialized,
   0 otherwise.
*************************************************/
int matrix_isinit(MATRIX *m);

/************************************************
    Resize an existing matrix.  Assumes that *m points
   to a MATRIX, and that MATRIX has a m->name.

   If *m has same size as requested size, then the funciton
   does nothing.  If not, then memory of the matrix, if
   initialized, is freed.  And then *m is initialized
   to requested size.
*************************************************/
void matrix_resize(MATRIX *m,int nr,int nc);

/************************************************
    Map a subset of elements from the input matrix to
   the output matrix.  Similar to matlab:
   out = in(r1:r2,c1:c2);
*************************************************/
void matrix_reindex(MATRIX *out,MATRIX *in,int c1,int c2,int r1,int r2);

/************************************************
    Repmat command.  Similar to out = repmat(in,m,n)
*************************************************/
void matrix_repmat(MATRIX *out,MATRIX *in,int m,int n);

/************************************************
    Length of matrix.  Returns an int = max(in->nc,in->nr)
*************************************************/
int matrix_length(MATRIX *in);

/************************************************
    Sum of all elements in matrix.
*************************************************/
MWORD matrix_sumsum(MATRIX *in);

/************************************************
    Sum,max min of elements in matrix in one dimension
*************************************************/
void matrix_sum(MATRIX *out,MATRIX *in1,int dim);
void matrix_max(MATRIX *out,MATRIX *in1,int dim);
void matrix_min(MATRIX *out,MATRIX *in1,int dim);

/************************************************
    returns a one if all elements of two matrices
   are equal, zero otherwise
*************************************************/
int matrix_equal(MATRIX *in1,MATRIX *in2);

/************************************************
    Copy a matrix from src to dest.
*************************************************/
void matrix_copy(MATRIX *dest,MATRIX *src);

/************************************************
    Sort entries of a matrix in ascending order.
   Each row is sorted independently.
*************************************************/
void matrix_sort(MATRIX *m,MATRIX *in);

/************************************************
    Add two matrices.  *in1 and *in2 must be the
   same size.
*************************************************/
void matrix_add(MATRIX *out,MATRIX *in1,MATRIX *in2);


/************************************************
    Perform a function like matlab
    OUT = (IN < IN2)

   If IN2 is a matrix of the same size as IN, then:
   "<" may be any of "<",">","<=",">=","=="

   If IN2 is an MWORD, then prefix the operator with "s":
   any of "s<","s>","s<=","s>=","s=="
   "s" is for scalar
*************************************************/
void matrix_op(MATRIX *out,MATRIX *in1,char *cs, ...);


/************************************************
    Convert a matrix entires to integer.  OP should be:
    "fix"
*************************************************/
void matrix_int(char *op,MATRIX *in,MATRIX *out);


/************************************************
    Convert a matrix entires to integer.  OP may be one of:
    "floor"
*************************************************/
void matrix_labsort(MATRIX *ilab, MATRIX *lab);

/*usage: MATRIX m; matrix_name(&m,"mymat");*/
void matrix_name(MATRIX *m, char *name);


/************************************************
    Set all entries of an existing matrix to
   a single value f.
*************************************************/
void matrix_clear(MATRIX *m,MWORD f);

/************************************************
   Finds the *row and *col of the first value that
   tests equal to f in matrix *m.  Starts at entry 0,0.
   Scans proceeds across increasing columns,
   then increasing rows.  Returns *row,*col=-1 if no
   match is found.

   Skips entires that are NaN.
*************************************************/
void matrix_find_equal_rc(MATRIX *m, MWORD f, int *row, int *col);

void matrix_write_quick(FILE *fp,MATRIX *m);
void matrix_disp(MATRIX *m);
void matrix_copy(MATRIX *dest,MATRIX *src);
MWORD matrix_find_min(MATRIX *m,MATRIX *out);
MWORD matrix_minmin(MATRIX *m,MWORD *f);
void matrix_clear_index(MATRIX *m,MWORD f,MATRIX *index);




#define NTAB 32
typedef struct
{
    long idum2;
    long iy;
    long iv[NTAB];
   long idum;
    int iset;     /* required for gasdevs */
    float gset;   /* required for gasdevs */
} RANDSTATE;
#undef NTAB

/************************************************
   Generates a matrix with random numbers between 0 and 1.
   Typical use:
      MATRIX *x;
      RANDSTATE *rs; 
      x = matrix(1,1000,"x");
      rs = randstate_init(0);
      matrix_rand(m,rs);
   Or,
      matrix_randq(m);
   "q" suffix funciton names are "quick" that use an internal
   seed of -1.
   
   randstate_init() performs a malloc for the returned
   RANDSTATE, which should ulitmately be free'd.
*************************************************/
RANDSTATE *randstate_init(long seed);
void matrix_rand(MATRIX *m,RANDSTATE *rs);
void matrix_randn(MATRIX *m,RANDSTATE *rs);

void matrix_randq(MATRIX *m);
void matrix_randnq(MATRIX *m);


/************************************************
   Displays an error message and forces an exit with
   return code 1 (if MEX==0) or uses a mex error message
   function (if MEX==1).  Example:

        error("matrix_add:must have same number of columns");
*************************************************/
void error(char *message);


/*
    MWORD es;
    es = findeps();
   fprintf(stdout,"EPS=%E\n",es);
*/
MWORD findeps();

/*NaNf makes *x into NaN (floats only)*/
void NaNf (float *x);
void NaNd (double *x);

/*returns 1 if x is NaN, 0 otherwise (floats only)*/
int IsNaNf (float x);
int IsNaNd (double x);

#if MEX==1
MATRIX *matlab_to_matrix(const mxArray *matlab);
void matrix_to_matlab_old(MATRIX *m,mxArray **matlab);
mxArray *matrix_to_matlab(MATRIX *m);
#endif


/*
    Philosophy. This is how I decided to initialize matrices:
   MATRIX *x;
   x = matrix();
   ...
   matrix_add(x,a,b);

    If I had done this instead:

   MATRIX *x;
   ...
   x = matrix_add(a,b);

   Then x=matrix(); isn't necessary.  But, it isn't well
   suited for repetitive processing, because the output
   is returned, it must be malloc'd every time--I cannot
   reuse already malloc'd matrices.

   I could also do this:
   MATRIX *x;
   x = matrix_init(3,1000,"x");
   ...
   matrix_add(x,a,b);
   which is OK, but that means I have to keep track of the
   size of outputs, which is a little bit tedius--why not
   have the program do it for me??

   Some problems might be mitigated by using &x, but I
   want to stay away from pointers to pointers. Yuck!
*/
