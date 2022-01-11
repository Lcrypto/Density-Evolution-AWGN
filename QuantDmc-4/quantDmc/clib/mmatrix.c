
#include <stdio.h>
#include <stdlib.h>  /* for malloc */
#include <string.h>  /* for strcpy */
#include <math.h>    /* for random number generation */
#include <time.h>    /* for seeding random number generators */
#include <stdarg.h>  /* for matrix_op which uses variable args */
#include <float.h>   /* for DBL_MAX and DBL_MIN, etc. */
#include "mmatrix.h"
#include "matrix.h"

#if MEX==1
	#include "mex.h"
#endif


MATRIX *matrix(char *name,int nr,int nc,...)
{
	va_list ap;
	int ii,jj,f=0;
   MATRIX *m;
   MWORD mw;

	va_start(ap,nc);
    if(*name == '*') {
		if(*(name+1) == 'i') {
	   	f = 1;
      } else {
      	if(*(name+1) == 'm') {
	   	   f = 2;
         } else {
         	error("matrix(): unrecognized option\n");
         }
      }
      name = name+3;
   }

   #if MEX==1
	   m = mxCalloc(1,sizeof(MATRIX));
   #else
	   m = malloc(sizeof(MATRIX));
   #endif

   if(m == NULL) {
	   printf("matrix(): matrix %s failed to malloc.\n",name);
   	error("Error in matrix()\n");
   }

   matrix_name(m,name);
   m->type = 0;

   if( (nc > 0) && (nr > 0) ) {
   	m->nr = nr;
	   m->nc = nc;
   	m->cols = malloc(sizeof(MWORD) * (m->nr) * (m->nc) );
	   m->elt = malloc(sizeof(m->cols) * (m->nc) );
	   for(jj=0;jj<m->nc;jj++)
   		m->elt[jj] = m->cols + (jj * (m->nr));
      m->init = 1;
		matrix_clear(m,0);
   } else {
	   m->cols = NULL;
   	m->elt = NULL;
	   m->nc = 0;
   	m->nr = 0;
      m->init = 0;
   }

	if(f==1) {
   	for(ii=0;ii<nr;ii++) {
      	for(jj=0;jj<nc;jj++) {
         	mw = (MWORD)va_arg(ap, int);
            m->elt[jj][ii] = mw;
         }
      }
   }

	if(f==2) {
   	for(ii=0;ii<nr;ii++) {
      	for(jj=0;jj<nc;jj++) {
         	mw = va_arg(ap, MWORD);
            m->elt[jj][ii] = mw;
         }
      }
   }

   return m;
}

/* BACKUP
MATRIX *matrix(char *name,int nr,int nc)
{
	int jj,f=0;
   MATRIX *m;

	m = malloc(sizeof(MATRIX));
   if(m == NULL) {
	   printf("matrix(): matrix %s failed to malloc.\n",name);
   	error("Error in matrix()\n");
   }

   matrix_name(m,name);
   m->type = 0;

   if( (nc > 0) && (nr > 0) ) {
   	m->nr = nr;
	   m->nc = nc;
   	m->cols = malloc(sizeof(MWORD) * (m->nr) * (m->nc) );
	   m->elt = malloc(sizeof(m->cols) * (m->nc) );
	   for(jj=0;jj<m->nc;jj++)
   		m->elt[jj] = m->cols + (jj * (m->nr));
		matrix_clear(m,0);
   } else {
	   m->cols = NULL;
   	m->elt = NULL;
	   m->nc = 0;
   	m->nr = 0;
   }

   return m;
}
*/

void matrix_init_ptr(MATRIX *m,int nr,int nc)
{
	int jj;

	
	if (matrix_isinit(m)) {
      if ((m->nc == nc) && (m->nr == nr) ) {
	      matrix_clear(m,0);
   	   return;
   	} else {
      	matrix_free(m);
      }
   }

   m->nr = nr;
   m->nc = nc;
   m->cols = malloc(sizeof(MWORD) * (m->nr) * (m->nc) );
   m->elt = malloc(sizeof(m->cols) * (m->nc) );

   for(jj=0;jj<m->nc;jj++)
   	m->elt[jj] = m->cols + (jj * (m->nr));

   m->type = 0;
	m->init = 1;

	matrix_clear(m,0);

   return;
}


void matrix_resize(MATRIX *m,int nr,int nc)
{
	int jj,init;

   init = matrix_isinit(m);

	if( init && (m->nc == nc) && (m->nr == nr) )
   	return;

   if(init)
   	matrix_free(m);

   m->nr = nr;
  	m->nc = nc;
   m->cols = malloc(sizeof(MWORD) * (m->nr) * (m->nc) );
  	m->elt = malloc(sizeof(m->cols) * (m->nc) );

   for(jj=0;jj<m->nc;jj++)
  		m->elt[jj] = m->cols + (jj * (m->nr));

	matrix_clear(m,0);
   m->type = 0;

   return;
}

void matrix_free(MATRIX *m)
{
	if(matrix_isinit(m)) {
		free(m->elt);
   	free(m->cols);

		m->cols = NULL;
   	m->elt = NULL;
	   m->nc = 0;
   	m->nr = 0;
   } else {
   	error("Tried to matrix_free() a matrix that was not initialized\n");
   }

	return;
}

void matrix_copy(MATRIX *dest,MATRIX *src)
{
   int ii,jj;

   /* if (!matrix_isinit(dest)) {
   	matrix_init_ptr(dest,src->nr,src->nc);
   } else {
		if( (dest->nc != src->nc) || (dest->nr != src->nr) ) {
	   	matrix_free(dest);
	   	matrix_init_ptr(dest,src->nr,src->nc);
      }
	} */
	
	matrix_init_ptr(dest,src->nr,src->nc);

   for(ii=0;ii < src->nr;ii++)
	   for(jj=0;jj < src->nc;jj++)
			dest->elt[jj][ii] = src->elt[jj][ii] ;
   return;
}

void matrix_reindex(MATRIX *out,MATRIX *in,int c1,int c2,int r1,int r2)
{
	int nr = r2-r1+1;
   int nc = c2-c1+1;
   int ci,ri;

   if (!matrix_isinit(out)) {
		matrix_init_ptr(out,nr,nc);
   } else {
		if( (out->nc != nc) || (out->nr != nr) ) {
  			matrix_free(out);
  			matrix_init_ptr(out,nr,nc);
      }
	}

   for(ci=0;ci<nc;ci++) {
		for(ri=0;ri<nr;ri++) {
    		out->elt[ci][ri] = in->elt[ci+c1][ri+r1];
      }
   }

	return;
}

void matrix_repmat(MATRIX *out,MATRIX *in,int m,int n)
{
   int nr = in->nr;
   int nc = in->nc;
   int mi,ni;
   int ici,iri;
   MATRIX *newin = matrix("newin",0,0);

   matrix_copy(newin,in);
   
   if (!matrix_isinit(out)) {
		matrix_init_ptr(out,nr*m,nc*n);
   } else {
		if( (out->nc != nc*n) || (out->nr != nr*m) ) {
  			matrix_free(out);
  			matrix_init_ptr(out,nr*m,nc*n);
      }
	}
   
   for(mi=0;mi<m;mi++) {
      for(ni=0;ni<n;ni++) {
         for(ici=0;ici<nc;ici++) {
		      for(iri=0;iri<nr;iri++) {
               out->elt[ni*nc + ici][mi*nr + iri] = newin->elt[ici][iri];
            }
         }
      }
   }
   
   matrix_free(newin);
   return;
}

void matrix_null(MATRIX *m)
{
   m->cols = NULL;
   m->elt = NULL;
   m->nc = 0;
   m->nr = 0;

   return;
}

void pakm(char *comment)
{
	printf("%s",comment);
	printf("Press Enter.");
   getchar();

   return;
}

void matrix_sort(MATRIX *m,MATRIX *in)
{
 	int r,c,ii;
   MWORD f;

   matrix_copy(m,in);

   for (r=0;r<in->nr;r++) {  /*do each row independently*/

 		for (ii=0;ii<in->nc;ii++) {
	   	for (c=0;c<(in->nc)-1;c++) {
				if (m->elt[c][r] > m->elt[c+1][r]) {
            	f = m->elt[c][r];
               m->elt[c][r] = m->elt[c+1][r];
               m->elt[c+1][r] = f;
            }
         }
      }

   } /*for r*/

	return;
}

void matrix_labsort(MATRIX *ilab, MATRIX *lab)
{
   MATRIX *temp = matrix("temp",1,lab->nc * lab->nr);
   MATRIX *temp2 = matrix("temp2",0,0);
	int szt=0,szt2=0;

   int r,c;

   MWORD *f = temp->cols;


	/* copy non-NaN entires to temp matrix */

   for (r=0;r<lab->nr;r++) {

   	for (c=0;c<lab->nc;c++) {

       	if( !IsNaNd(lab->elt[c][r])) {

         	*(f++) = lab->elt[c][r];

            szt++;

         }

   	}

   }


   if (szt == 0)

   	error("matrix_labsort: empty or all-NaN matrix\n");


   /* sort */

   temp->nc = szt;

	matrix_sort(temp2,temp);


   /* Eliminate duplicate entires.  Use temp for storage */

	szt = 1;

   temp->elt[0][0] = temp2->elt[0][0]; /* one must match */

   for (c=1;c<temp2->nc;c++) {

   	if(  temp2->elt[c][0] != temp->elt[szt-1][0] ) {

			temp->elt[szt][0] = temp2->elt[c][0] ;

         szt++;

      }

   }

	temp->nc = szt;

   matrix_copy(ilab,temp);  /*ilab is just the right size */


   matrix_free(temp2);

   matrix_free(temp);


   return;

}


int strhash(char *str)

{

	int rv=0;

   long sc=1;


	while(*str != 0) {

   	if(sc == 0)

      	error("str2num: well, that did not work\n");

   	rv = rv + (*str - 'a') * sc;

      printf("%ld\n",sc);

      str++;

      sc = sc << 5;

	}

   return rv;

}


void matrix_int(char *op,MATRIX *in,MATRIX *out)

{

	int ii;


	if(in != out) {

	   if (!matrix_isinit(out)) {
  			matrix_init_ptr(out,in->nr,in->nc);
	   } else {
			if( (out->nc != in->nc) || (out->nr != in->nr) ) {
   			matrix_free(out);
   			matrix_init_ptr(out,in->nr,in->nc);
	      }
		}
   }
   if (!strcmp(op,"fix"))
	{
     	for(ii=0;ii<(in->nr)*(in->nc);ii++)
     		out->cols[ii] = (int) in->cols[ii] ;
   } else {
     	printf("error matrix_int:%s",op);
      error("matrix_int:did not recognize oper\n");
   }

   return;
}


int matrix_length(MATRIX *in)

{

	if (!matrix_isinit(in))

   	error("matrix_length: input not initialized\n");


   if (in->nc > in->nr)

   	return in->nc;

   else

   	return in->nr;

}


MWORD matrix_sumsum(MATRIX *in)
{
	int ii;
   MWORD out=0;

	if (!matrix_isinit(in))
   	error("matrix_sumsum: input not initialized\n");

   for(ii=0;ii<(in->nr)*(in->nc);ii++)
    	out = out + in->cols[ii] ;
  	return out;
}


void matrix_op(MATRIX *out,MATRIX *in1,char *cs, ...)
{
	va_list ap;

   int ii;
   char *oper;
   MWORD in2d;
   MATRIX *in2;

   MWORD f,g;
   
	va_start(ap,cs);

	if(in1 != out) {
	   if (!matrix_isinit(out)) {
  			matrix_init_ptr(out,in1->nr,in1->nc);
	   } else {
			if( (out->nc != in1->nc) || (out->nr != in1->nr) ) {
   			matrix_free(out);
   			matrix_init_ptr(out,in1->nr,in1->nc);
	      }
		}
   }

	if( *cs == 's' ) {
	   oper = cs+1;

      in2d = va_arg(ap, MWORD);
		/*ii = *oper + (*oper+1)*256 ;*/
		switch ( *oper + 256* *(oper+1) ) {  /* kind of a hack for switching two chars */

	      case '<':
			   for(ii=0;ii<(in1->nr)*(in1->nc);ii++)
      	   	out->cols[ii] = (in1->cols[ii] < in2d) ;
         	break;
	      case '>':
			   for(ii=0;ii<(in1->nr)*(in1->nc);ii++)
      	   	out->cols[ii] = (in1->cols[ii] > in2d) ;
         	break;
	      case '<' + '='*256 :
			   for(ii=0;ii<(in1->nr)*(in1->nc);ii++)
      	   	out->cols[ii] = (in1->cols[ii] <= in2d) ;
         	break;
	      case '>' + '='*256 :
			   for(ii=0;ii<(in1->nr)*(in1->nc);ii++)
      	   	out->cols[ii] = (in1->cols[ii] >= in2d) ;
         	break;
	      case '=' + '='*256 :
			   for(ii=0;ii<(in1->nr)*(in1->nc);ii++)
      	   	out->cols[ii] = (in1->cols[ii] == in2d);
         	break;
	      case '.' + '^'*256 :
			   for(ii=0;ii<(in1->nr)*(in1->nc);ii++)
      	   	out->cols[ii] = pow(in1->cols[ii],in2d);
         	break;
	      case '*':
			   for(ii=0;ii<(in1->nr)*(in1->nc);ii++)
      	   	out->cols[ii] = (in1->cols[ii] * in2d) ;
         	break;
	      case '/':
			   for(ii=0;ii<(in1->nr)*(in1->nc);ii++)
      	   	out->cols[ii] = (in1->cols[ii] / in2d) ;
         	break;
	      case '-':
			   for(ii=0;ii<(in1->nr)*(in1->nc);ii++)
      	   	out->cols[ii] = (in1->cols[ii] - in2d) ;
         	break;
	      case '+':
			   for(ii=0;ii<(in1->nr)*(in1->nc);ii++)
      	   	out->cols[ii] = (in1->cols[ii] + in2d) ;
         	break;
	      default:
   	   	error("matrix_op: did not recognize operator\n");
      }

   } else {
      in2 = va_arg(ap,MATRIX *);
		if(in1->nc != in2->nc)
   		error("matrix_op: number of columns must be equal\n");
		if(in1->nr != in2->nr)
   		error("matrix_op: number of row must be equal\n");


	   oper = cs;
		switch ( *oper + 256* *(oper+1) ) {  /* kind of a hack for switching two chars */

	      case '.' + '*'*256 :
			   for(ii=0;ii<(in1->nr)*(in1->nc);ii++)
      	   	out->cols[ii] = in1->cols[ii] * in2->cols[ii] ;
         	break;
	      case '+':
			   for(ii=0;ii<(in1->nr)*(in1->nc);ii++)
      	   	out->cols[ii] = (in1->cols[ii] + in2->cols[ii]) ;
         	break;
	      case '-':
			   for(ii=0;ii<(in1->nr)*(in1->nc);ii++)
      	   	out->cols[ii] = (in1->cols[ii] - in2->cols[ii]) ;
         	break;
	      default:
   	   	error("matrix_op: did not recognize operator\n");
      }


/*
	   for(ii=0;ii<(in1->nr)*(in1->nc);ii++)
				out->cols[ii] = in1->cols[ii] ;
*/
   }
   return;
}
 void matrix_clear(MATRIX *m,MWORD f)
{
 	int ii;

	if(matrix_isinit(m)) {
	   for(ii=0;ii < m->nr * m->nc ;ii++)
				m->cols[ii] = f;
   } else {
      printf("name = %s\n",m->name);
   	error("matrix_clear() must have a defined matrix\n");
   }

   return;
}

int matrix_isinit(MATRIX *m)
{
 	return m->init;
}

void matrix_name(MATRIX *m, char *name)
{
 	strcpy(m->name,name);
	return;
}

void matrix_clear_index(MATRIX *m,MWORD f,MATRIX *index)
{
 	int ii;

   for(ii=0;ii < index->nr;ii++)
			m->elt[ (int) index->elt[ii][0] ][ (int) index->elt[ii][1] ] = f;
   return;
}

void matrix_add(MATRIX *out,MATRIX *in1,MATRIX *in2)
{
	int r,c;

	if(in1->nc != in2->nc)
   	error("matrix_add: number of columns must be equal\n");
	if(in1->nr != in2->nr)
   	error("matrix_add: number of row must be equal\n");

	matrix_init_ptr(out,in1->nr,in1->nc);

   for(r=0;r < in1->nr;r++)
	   for(c=0;c < in1->nc;c++)
			out->elt[c][r] = in1->elt[c][r] + in2->elt[c][r];

 	return;
}

void matrix_sum(MATRIX *out,MATRIX *in,int dim)
{
	int ii,jj;
	int nc = in->nc;
	int nr = in->nr;
	MATRIX *in1 =  matrix("in1",0,0);
	
	if (out == in) {
	   matrix_copy(in1,in);
	} else {
	   in1 = in;
	}
   
	if (dim==1) {
      if (!matrix_isinit(out)) {
   		matrix_init_ptr(out,1,nc);
      } else {
   		if((out->nc != nc) || (out->nr != 1) ) {
  	   		matrix_free(out);
     			matrix_init_ptr(out,1,nc);
         }
   	}
   	
   	for (ii=0;ii<nc;ii++) {
   	   for (jj=0;jj<nr;jj++) {
   	      out->elt[ii][0] = out->elt[ii][0] + in1->elt[ii][jj] ;
   	   }
   	}
	
	} else if (dim==2) {
      if (!matrix_isinit(out)) {
   		matrix_init_ptr(out,nr,1);
      } else {
   		if((out->nc != 1) || (out->nr != nr) ) {
  	   		matrix_free(out);
     			matrix_init_ptr(out,nr,1);
         }
   	}
   	
   	for (ii=0;ii<nc;ii++) {
   	   for (jj=0;jj<nr;jj++) {
   	      out->elt[0][jj] = out->elt[0][jj] + in1->elt[ii][jj] ;
   	   }
   	}
	
	
	
	} else {
	   error("matrix_sum: dim must be 1 or 2\n");
	}


	matrix_free(in1);
	
 	return;
}

void matrix_max(MATRIX *out,MATRIX *in,int dim)
{
	int ii,jj;
	int nc = in->nc;
	int nr = in->nr;
	MATRIX *in1 =  matrix("in1",0,0);
	MWORD f;
	
	if (out == in) {
	   matrix_copy(in1,in);
	} else {
	   in1 = in;
	}

	if (dim==1) {
      if (!matrix_isinit(out)) {
   		matrix_init_ptr(out,1,nc);
      } else {
   		if((out->nc != nc) || (out->nr != 1) ) {
  	   		matrix_free(out);
     			matrix_init_ptr(out,1,nc);
         }
   	}
   	
   	for (ii=0;ii<nc;ii++) {
   	   f = -MMAXFLOAT;
   	   for (jj=0;jj<nr;jj++) {
   	      f = (in1->elt[ii][jj] > f) ? in1->elt[ii][jj] : f;
   	   }
   	   out->elt[ii][0] = f;
   	}
	
	} else if (dim==2) {
      if (!matrix_isinit(out)) {
   		matrix_init_ptr(out,nr,1);
      } else {
   		if((out->nc != 1) || (out->nr != nr) ) {
  	   		matrix_free(out);
     			matrix_init_ptr(out,nr,1);
         }
   	}

   	for (jj=0;jj<nr;jj++) {
         f = -MMAXFLOAT;
   	   for (ii=0;ii<nc;ii++) {
   	      f = (in1->elt[ii][jj] > f) ? in1->elt[ii][jj] : f;
   	   }
 	      out->elt[0][jj] = f ;
   	}
	
	
	
	} else {
	   error("matrix_sum: dim must be 1 or 2\n");
	}


	matrix_free(in1);
	
 	return;
}

void matrix_min(MATRIX *out,MATRIX *in,int dim)
{
	int ii,jj;
	int nc = in->nc;
	int nr = in->nr;
	MATRIX *in1 =  matrix("in1",0,0);
	MWORD f;
	
	if (out == in) {
	   matrix_copy(in1,in);
	} else {
	   in1 = in;
	}

	if (dim==1) {
      if (!matrix_isinit(out)) {
   		matrix_init_ptr(out,1,nc);
      } else {
   		if((out->nc != nc) || (out->nr != 1) ) {
  	   		matrix_free(out);
     			matrix_init_ptr(out,1,nc);
         }
   	}
   	
   	for (ii=0;ii<nc;ii++) {
   	   f = MMAXFLOAT;
   	   for (jj=0;jj<nr;jj++) {
   	      f = (in1->elt[ii][jj] < f) ? in1->elt[ii][jj] : f;
   	   }
   	   out->elt[ii][0] = f;
   	}
	
	} else if (dim==2) {
      if (!matrix_isinit(out)) {
   		matrix_init_ptr(out,nr,1);
      } else {
   		if((out->nc != 1) || (out->nr != nr) ) {
  	   		matrix_free(out);
     			matrix_init_ptr(out,nr,1);
         }
   	}

   	for (jj=0;jj<nr;jj++) {
         f = MMAXFLOAT;
   	   for (ii=0;ii<nc;ii++) {
   	      f = (in1->elt[ii][jj] < f) ? in1->elt[ii][jj] : f;
   	   }
 	      out->elt[0][jj] = f ;
   	}
	
	
	
	} else {
	   error("matrix_sum: dim must be 1 or 2\n");
	}


	matrix_free(in1);
	
 	return;
}

int matrix_equal(MATRIX *in1,MATRIX *in2)
{
	int ii;

	if(in1->nc != in2->nc)
   	error("matrix_equal: number of columns must be equal\n");
	if(in1->nr != in2->nr)
   	error("matrix_equal: number of row must be equal\n");

   for(ii=0;ii < in1->nr * in1->nc;ii++) {
    	if(in1->cols[ii] != in2->cols[ii])
      	return 0;
   }

 	return 1;
}

void error(char *message)
{
	#if MEX==0
	printf("ERROR:\n");
	printf("%s",message);
   exit(1);
   #else
   mexErrMsgTxt(message);
	#endif
}

/*
MWORD matrix_minmin(MATRIX *m,MWORD *f)
{
 	int ii,jj,row,col;

   pakm("matrix_min():this code may contain an error: use of row before definition");
	pakm("PLease correct the rrror.");

   *f = MMAXFLOAT;
   for(ii=0;ii < m->nr;ii++)
	   for(jj=0;jj < m->nc;jj++)
      {
			if( m->elt[row][col] < *f )
         	*f = m->elt[row][col];
      }
   return *f;
}
*/

/*
MWORD matrix_find_min(MATRIX *m,MATRIX *out)
{
 	int ii,jj;
	MWORD f,t,*findrow,*findcol;
	int listlen=0;

	findrow = malloc((m->nr) * (m->nc));
   findcol = malloc((m->nr) * (m->nc));

   f = MMAXFLOAT;
   for(ii=0;ii < m->nr;ii++)
   {
	   for(jj=0;jj < m->nc;jj++)
      {
			t = m->elt[ii][jj];
         if (!IsNaNf(t)) {
            if( t < f ) {
               f = t;
         	   listlen=0;
            	findrow[listlen] = ii;
	            findcol[listlen++] = jj;
   	      } else if( t == f ) {
               findrow[listlen] = ii;
               findcol[listlen++] = jj;
   	      }
         }
      }
   }

   out->nr = listlen;
   out->nc = 2;
	if(listlen != 0)
   {
		matrix_init(out);
   	for(ii=0;ii<listlen;ii++)
	   {
			out->elt[ii][0] = findrow[ii];
			out->elt[ii][1] = findcol[ii];
	   }
   } else {
	   out->nc = 0;
      NaNd(&f);
   }

   free(findrow);
   free(findcol);

   return f;
}

*/

void matrix_find_equal_rc(MATRIX *m,MWORD f, int *row, int *col)
{
	int r,c;

   *row = -1;
   *col = -1;

   for(r=0;r< m->nr;row++)
   	for(c=0;c< m->nc;c++)
      	if( !IsNaNf(m->elt[c][r]) && (m->elt[c][r] == f) ) {
          	*row = r;
            *col = c;
            return;
         }

	return;
}


void matrix_write_quick(FILE *fp,MATRIX *m)
{
	int ii,jj;

	if(m->type == 0)
   {
		fprintf(fp,"%s = ...\n[",m->name);
   	for(ii=0;ii<(m->nr)-1;ii++)
		{
   	 	for(jj=0;jj<m->nc;jj++)
      		fprintf(fp,"%6.3f\t",m->elt[jj][ii]);
	      fprintf(fp,"; ...\n");
   	}
	   ii=m->nr-1;
 		for(jj=0;jj<m->nc;jj++)
    		fprintf(fp,"%6.3f\t",m->elt[jj][ii]);
	   fprintf(fp," ];\n");
   }

	if(m->type == 1)
   {
		fprintf(fp,"%s = [ '",m->name);
   	for(ii=0;ii<(m->nr)-1;ii++)
		{
   	 	for(jj=0;jj<m->nc;jj++)
      		fprintf(fp,"%c",(int) m->elt[ii][jj]);
	      fprintf(fp,"; ...\n");
   	}
	   ii=m->nr-1;
 		for(jj=0;jj<m->nc;jj++)
      		fprintf(fp,"%c",(int) m->elt[ii][jj]);
	   fprintf(fp,"' ];\n");
   }

	return;
}

void matrix_disp(MATRIX *m)
{
	int ii,jj;

	if(m->type == 0)
   {
		printf("%s = ...\n[",m->name);
   	for(ii=0;ii<(m->nr)-1;ii++)
		{
   	 	for(jj=0;jj<m->nc;jj++)
      		printf("%6.3f\t",m->elt[jj][ii]);
	      printf("; ...\n");
   	}
	   ii=m->nr-1;
 		for(jj=0;jj<m->nc;jj++)
    		printf("%6.3f\t",m->elt[jj][ii]);
	   printf(" ];\n");
   }

	if(m->type == 1)
   {
		printf("%s = [ '",m->name);
   	for(ii=0;ii<(m->nr)-1;ii++)
		{
   	 	for(jj=0;jj<m->nc;jj++)
      		printf("%c",(int) m->elt[ii][jj]);
	      printf("; ...\n");
   	}
	   ii=m->nr-1;
 		for(jj=0;jj<m->nc;jj++)
      		printf("%c",(int) m->elt[ii][jj]);
	   printf("' ];\n");
   }

	return;
}

#if MEX==1

void matrix_to_matlab_old(MATRIX *m,mxArray **matlab)
{
   double *fp;
   int r,c;

   *matlab = mxCreateDoubleMatrix(m->nc,m->nr,mxREAL);
   fp     = mxGetPr(*matlab);
	for(c=0;c < m->nc;c++)
      for(r=0;r < m->nr;r++)
			*(fp++) = m->elt[r][c] ;
	return;
} 

/*
   MATRIX test;
   test.nr=3;
   test.nc=4;
   matrix_malloc("test",&test);
   test.matrix[0][0] = 1;test.matrix[0][1] = 11;test.matrix[0][2] = 21;
   test.matrix[1][0] = 2;test.matrix[1][1] = 12;test.matrix[1][2] = 22;
   test.matrix[2][0] = 3;test.matrix[2][1] = 13;test.matrix[2][2] = 23;

	fl=fopen("mwq","wb");
   matrix_write_quick(fl,&test);
   fclose(fl);
   matrix_to_matlab(&test,&(plhs[0]));
*/
#endif


/*if I was smart, I could have used the NaN routines to find eps*/
MWORD findeps()
{
 	MWORD up=1,low=0,nlow,nup,delta,t;
   int ii,jj;

   for(ii=1;ii<30;ii++) {
   	delta = (up - low) / 10.0;
      t = low;
      for(jj=0;jj<10;jj++) {
			printf("testing=%30.30E\n",t);
         if ((1- (1 - t)) == 0 ) {
         	nlow = t;
			} else {
         	nup = t;
            break;
         }
       	t = t + delta;
      }
      printf("=============== %d\n",ii);
      up = nup;
      low = nlow;
      if( 1-(1-up) == 0)
      	break;
   }
   return up;
}



/*------------------------------------------------------*/
/*  NaN routines ---------------------------------------*/
/*------------------------------------------------------*/

  /***************************************************************
		http://www.psc.edu/general/software/packages/ieee/ieee.html

	The IEEE single precision floating point standard representation
   requires a 32 bit word, which may be represented as numbered from 0 to 31,
   left to right. The first bit is the sign bit, S, the next eight
   bits are the exponent bits, 'E', and the final 23 bits are the fraction 'F':

	  S EEEEEEEE FFFFFFFFFFFFFFFFFFFFFFF
	  0 1      8 9                    31

	The value V represented by the word may be determined as follows:


	If E=255 and F is nonzero, then V=NaN ("Not a number")
	If E=255 and F is zero and S is 1, then V=-Infinity
	If E=255 and F is zero and S is 0, then V=Infinity
	If 0<E<255 then V=(-1)**S * 2 ** (E-127) * (1.F) where "1.F" is intended to represent the binary number created by prefixing F with an implicit leading 1 and a binary point.

	  0 00000000 00000000000000000000000 = 0
	  1 00000000 00000000000000000000000 = -0

	  0 11111111 00000000000000000000000 = Infinity
	  1 11111111 00000000000000000000000 = -Infinity

	  0 11111111 00000100000000000000000 = NaN
	  1 11111111 00100010001001010101010 = NaN

  ***************************************************************/
struct Float
{
	unsigned char byte[4];
};

union RealFloat
/* This union allows us to write in a float and read out a */
/* struct Float; */
{
	struct Float sld;
   float f;
};

void NaNf (float *x)
{
	union RealFloat real;

	real.sld.byte[0]= 0xFF;
	real.sld.byte[1]= 0xFF;
	real.sld.byte[2]= 0xFF;
	real.sld.byte[3]= 0x7F;

	*x = real.f;
   return;
}

int parsef (float x)
{
	union RealFloat real;
   int ii;
   real.f = x;
   printf("%f=",x);
   for(ii=0;ii<4;ii++)
   {
    	printf("%2X ",(int) real.sld.byte[ii]);
   }
   printf("\n");
   return 0;
}

int IsNaNf (float x)
{
	union RealFloat real;
   real.f = x;

   if ( ((real.sld.byte[3] & 0x7F) == 0x7F) && ((real.sld.byte[2] & 0x80) == 0x80) )
      return 1;
   return 0;
}

/* NaN double stuff */
struct  Double
{
  unsigned char byte[8];
} ;

union RealDouble
/* This union allows us to write in a double and read out a */
/* struct Double; */
{
	struct Double sld;
   double d;
};

void NaNd (double *x)
{
	union RealDouble real;

	real.sld.byte[0]= 0xFF;
	real.sld.byte[1]= 0xFF;
	real.sld.byte[2]= 0xFF;
	real.sld.byte[3]= 0xFF;
	real.sld.byte[4]= 0xFF;
	real.sld.byte[5]= 0xFF;
	real.sld.byte[6]= 0xFF;
	real.sld.byte[7]= 0x7F;

	*x = real.d;
   return;
}

int parsed (double x,FILE *fl)
{
	union RealDouble real;
   int ii;
   real.d = x;
   fprintf(fl,"%lf=",x);
   for(ii=0;ii<8;ii++)
   {
    	fprintf(fl,"%2X ",(int) real.sld.byte[ii]);
   }
   fprintf(fl,"\n");
   return 0;
}

int IsNaNd (double x)
{
	union RealDouble real;
   real.d = x;

   if ( ((real.sld.byte[7] & 0x7F) == 0x7F) && ((real.sld.byte[6] & 0xF0) == 0xF0) )
      return 1;
   return 0;
}


/*------------------------------------------------------*/
/*  Random number generation routines ------------------*/
/*------------------------------------------------------*/

/* Numerical Recipes ran2: generates uniform deviates */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPSX 1.2e-7
#define RNMX (1.0-EPSX)

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

void matrix_randq(MATRIX *m)
{
   static long seed=-1;
	int ii;

   for(ii=0;ii < m->nr * m->nc ;ii++)
			m->cols[ii] = ran2(&seed) ;
   return;
}

float ran2s(RANDSTATE *rs)
{
	int j;
	long k;
	float temp;

	long idum2;
	long iy;
	long iv[NTAB];
	long idum;

   idum2 = rs->idum2;
   iy = rs->iy;
   for(j=0;j<NTAB;j++)
   	iv[j] = rs->iv[j];
   idum = rs->idum;

	if (idum <= 0) {
		if (-(idum) < 1) idum=1;
		else idum = -(idum);
		idum2=(idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(idum)/IQ1;
			idum=IA1*(idum-k*IQ1)-k*IR1;
			if (idum < 0) idum += IM1;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=(idum)/IQ1;
	idum=IA1*(idum-k*IQ1)-k*IR1;
	if (idum < 0) idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = idum;
	if (iy < 1) iy += IMM1;

   rs->idum2 = idum2;
   rs->iy = iy;
   for(j=0;j<NTAB;j++)
   	rs->iv[j] = iv[j];
   rs->idum = idum;

	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

void matrix_rand(MATRIX *m,RANDSTATE *rs)
{
	int ii;
   for(ii=0;ii < m->nr * m->nc ;ii++)
			m->cols[ii] = ran2s(rs) ;
   return;
}

/*
RANDSTATE randstate_init(long seed)
{
 	RANDSTATE rs;
   time_t t0;

   rs.idum2 = 123456789;
   rs.iy = 0;
   rs.idum = -1;
  	if(seed == 0L)
   {
		time(&t0);
   	rs.idum = -t0;
   } else {
   	rs.idum = -labs(seed);
   }

	rs.iset=0;

   return rs;
}
*/


RANDSTATE *randstate_init(long seed)
{
   time_t t0;
   RANDSTATE *rs;

   rs = malloc(sizeof(RANDSTATE));

   rs->idum2 = 123456789;
   rs->iy = 0;
   rs->idum = -1;
  	if(seed == 0L)
   {
		time(&t0);
   	rs->idum = -t0;
   } else {
   	rs->idum = -labs(seed);
   }

	rs->iset=0;

   return rs;
}


#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPSX
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software *;++71,&. */


/* Numerical Recipes gasdev: generates normal deviates */
float gasdev(long *idum)
{
	float ran2(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran2(idum)-1.0;
			v2=2.0*ran2(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}


void matrix_randnq(MATRIX *m)
{
   static long seed=-1;
	int ii;

   for(ii=0;ii < m->nr * m->nc ;ii++)
			m->cols[ii] = gasdev(&seed) ;
   return;
}

float gasdevs(RANDSTATE *rs)
{
	float fac,rsq,v1,v2;

   int iset;
	float gset;

   iset = rs->iset;
   gset = rs->gset;

	if  (iset == 0) {
		do {
			v1=2.0*ran2s(rs)-1.0;
			v2=2.0*ran2s(rs)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
      rs->iset = iset;
      rs->gset = gset;
		return v2*fac;
	} else {
		iset=0;
      rs->iset = iset;
      rs->gset = gset;
		return gset;
	}
}


void matrix_randn(MATRIX *m,RANDSTATE *rs)
{
	int ii;

   for(ii=0;ii < m->nr * m->nc ;ii++)
			m->cols[ii] = gasdevs(rs) ;
   return;
}


/*------------------------------------------------------*/

#if MEX==1
MATRIX *matlab_to_matrix(const mxArray *matlab)
{
   double *fp;
   int r,c;
   FILE *fl;
   MATRIX *m = matrix("matlab",0,0);

   /* m = malloc(sizeof(MATRIX)); */
   
   matrix_init_ptr(m,mxGetM(matlab),mxGetN(matlab));
   fp = mxGetPr(matlab);
   m->type = 0;
   /* matrix_name(m,mxGetName(matlab));  not compatible with V6.5 */

   for(c=0;c < m->nc;c++)
  		for(r=0;r < m->nr;r++) {
			m->elt[c][r] = *(fp++) ;
   }
    
   /* keep this to remind me how to debug
   fl=fopen("lut2","wb");
   matrix_write_quick(fl,m);
   fclose(fl);
   */
   
   return m;
}
#endif

#if MEX==1

mxArray *matrix_to_matlab(MATRIX *m)
{
   double *fp;
   int r,c;
   mxArray *matlab;

   matlab = mxCreateDoubleMatrix(m->nr,m->nc,mxREAL);
   fp     = mxGetPr(matlab);
	for(c=0;c < m->nc;c++)
      for(r=0;r < m->nr;r++)
			*(fp++) = m->elt[c][r] ;
	return matlab;
} 

/*
   MATRIX test;
   test.nr=3;
   test.nc=4;
   matrix_malloc("test",&test);
   test.matrix[0][0] = 1;test.matrix[0][1] = 11;test.matrix[0][2] = 21;
   test.matrix[1][0] = 2;test.matrix[1][1] = 12;test.matrix[1][2] = 22;
   test.matrix[2][0] = 3;test.matrix[2][1] = 13;test.matrix[2][2] = 23;

	fl=fopen("mwq","wb");
   matrix_write_quick(fl,&test);
   fclose(fl);
   matrix_to_matlab(&test,&(plhs[0]));
*/
#endif
