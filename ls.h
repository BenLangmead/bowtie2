/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

/* Code in this file is ultimately based on:

   qsufsort.c
   Copyright 1999, N. Jesper Larsson, all rights reserved.

   This file contains an implementation of the algorithm presented in "Faster
   Suffix Sorting" by N. Jesper Larsson (jesper@cs.lth.se) and Kunihiko
   Sadakane (sada@is.s.u-tokyo.ac.jp).

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.*/

#ifndef LS_H_
#define LS_H_

#include <iostream>
#include <limits>
#include <stdint.h>

template<typename T>
class LarssonSadakane {
	T *I, /* group array, ultimately suffix array.*/
	*V,   /* inverse array, ultimately inverse of I.*/
	r,    /* number of symbols aggregated by transform.*/
	h;    /* length of already-sorted prefixes.*/

	#define LS_KEY(p)          (V[*(p)+(h)])
	#define LS_SWAP(p, q)      (tmp=*(p), *(p)=*(q), *(q)=tmp)
	#define LS_SMED3(a, b, c)  (LS_KEY(a)<LS_KEY(b) ?                        \
			  (LS_KEY(b)<LS_KEY(c) ? (b) : LS_KEY(a)<LS_KEY(c) ? (c) : (a))  \
			: (LS_KEY(b)>LS_KEY(c) ? (b) : LS_KEY(a)>LS_KEY(c) ? (c) : (a)))

	/* Subroutine for select_sort_split and sort_split. Sets group numbers for a
	   group whose lowest position in I is pl and highest position is pm.*/

	inline void update_group(T *pl, T *pm) {
	   T g;
	   g=(T)(pm-I);                 /* group number.*/
	   V[*pl]=g;                    /* update group number of first position.*/
	   if (pl==pm)
		  *pl=-1;                   /* one element, sorted group.*/
	   else
		  do                        /* more than one element, unsorted group.*/
			 V[*++pl]=g;            /* update group numbers.*/
		  while (pl<pm);
	}

	/* Quadratic sorting method to use for small subarrays. To be able to update
	   group numbers consistently, a variant of selection sorting is used.*/

	inline void select_sort_split(T *p, T n) {
	   T *pa, *pb, *pi, *pn;
	   T f, v, tmp;

	   pa=p;                        /* pa is start of group being picked out.*/
	   pn=p+n-1;                    /* pn is last position of subarray.*/
	   while (pa<pn) {
		  for (pi=pb=pa+1, f=LS_KEY(pa); pi<=pn; ++pi)
			 if ((v=LS_KEY(pi))<f) {
				f=v;                /* f is smallest key found.*/
				LS_SWAP(pi, pa);       /* place smallest element at beginning.*/
				pb=pa+1;            /* pb is position for elements equal to f.*/
			 } else if (v==f) {     /* if equal to smallest key.*/
				LS_SWAP(pi, pb);       /* place next to other smallest elements.*/
				++pb;
			 }
		  update_group(pa, pb-1);   /* update group values for new group.*/
		  pa=pb;                    /* continue sorting rest of the subarray.*/
	   }
	   if (pa==pn) {                /* check if last part is single element.*/
		  V[*pa]=(T)(pa-I);
		  *pa=-1;                   /* sorted group.*/
	   }
	}

	/* Subroutine for sort_split, algorithm by Bentley & McIlroy.*/

	inline T choose_pivot(T *p, T n) {
	   T *pl, *pm, *pn;
	   T s;

	   pm=p+(n>>1);                 /* small arrays, middle element.*/
	   if (n>7) {
		  pl=p;
		  pn=p+n-1;
		  if (n>40) {               /* big arrays, pseudomedian of 9.*/
			 s=n>>3;
			 pl=LS_SMED3(pl, pl+s, pl+s+s);
			 pm=LS_SMED3(pm-s, pm, pm+s);
			 pn=LS_SMED3(pn-s-s, pn-s, pn);
		  }
		  pm=LS_SMED3(pl, pm, pn);      /* midsize arrays, median of 3.*/
	   }
	   return LS_KEY(pm);
	}

	/* Sorting routine called for each unsorted group. Sorts the array of integers
	   (suffix numbers) of length n starting at p. The algorithm is a ternary-split
	   quicksort taken from Bentley & McIlroy, "Engineering a Sort Function",
	   Software -- Practice and Experience 23(11), 1249-1265 (November 1993). This
	   function is based on Program 7.*/

	inline void sort_split(T *p, T n)
	{
	   T *pa, *pb, *pc, *pd, *pl, *pm, *pn;
	   T f, v, s, t, tmp;

	   if (n<7) {                   /* multi-selection sort smallest arrays.*/
		  select_sort_split(p, n);
		  return;
	   }

	   v=choose_pivot(p, n);
	   pa=pb=p;
	   pc=pd=p+n-1;
	   while (1) {                  /* split-end partition.*/
		  while (pb<=pc && (f=LS_KEY(pb))<=v) {
			 if (f==v) {
				LS_SWAP(pa, pb);
				++pa;
			 }
			 ++pb;
		  }
		  while (pc>=pb && (f=LS_KEY(pc))>=v) {
			 if (f==v) {
				LS_SWAP(pc, pd);
				--pd;
			 }
			 --pc;
		  }
		  if (pb>pc)
			 break;
		  LS_SWAP(pb, pc);
		  ++pb;
		  --pc;
	   }
	   pn=p+n;
	   if ((s=(T)(pa-p))>(t=(T)(pb-pa)))
		  s=t;
	   for (pl=p, pm=pb-s; s; --s, ++pl, ++pm)
		  LS_SWAP(pl, pm);
	   if ((s=(T)(pd-pc))>(t=(T)(pn-pd-1)))
		  s=t;
	   for (pl=pb, pm=pn-s; s; --s, ++pl, ++pm)
		  LS_SWAP(pl, pm);

	   s=(T)(pb-pa);
	   t=(T)(pd-pc);
	   if (s>0)
		  sort_split(p, s);
	   update_group(p+s, p+n-t-1);
	   if (t>0)
		  sort_split(p+n-t, t);
	}

	/* Bucketsort for first iteration.

	   Input: x[0...n-1] holds integers in the range 1...k-1, all of which appear
	   at least once. x[n] is 0. (This is the corresponding output of transform.) k
	   must be at most n+1. p is array of size n+1 whose contents are disregarded.

	   Output: x is V and p is I after the initial sorting stage of the refined
	   suffix sorting algorithm.*/

	inline void bucketsort(T *x, T *p, T n, T k)
	{
	   T *pi, i, c, d, g;

	   for (pi=p; pi<p+k; ++pi)
		  *pi=-1;                   /* mark linked lists empty.*/
	   for (i=0; i<=n; ++i) {
		  x[i]=p[c=x[i]];           /* insert in linked list.*/
		  p[c]=i;
	   }
	   for (pi=p+k-1, i=n; pi>=p; --pi) {
		  d=x[c=*pi];               /* c is position, d is next in list.*/
		  x[c]=g=i;                 /* last position equals group number.*/
		  if (d == 0 || d > 0) {    /* if more than one element in group.*/
			 p[i--]=c;              /* p is permutation for the sorted x.*/
			 do {
				d=x[c=d];           /* next in linked list.*/
				x[c]=g;             /* group number in x.*/
				p[i--]=c;           /* permutation in p.*/
			 } while (d == 0 || d > 0);
		  } else
			 p[i--]=-1;             /* one element, sorted group.*/
	   }
	}

	/* Transforms the alphabet of x by attempting to aggregate several symbols into
	   one, while preserving the suffix order of x. The alphabet may also be
	   compacted, so that x on output comprises all integers of the new alphabet
	   with no skipped numbers.

	   Input: x is an array of size n+1 whose first n elements are positive
	   integers in the range l...k-1. p is array of size n+1, used for temporary
	   storage. q controls aggregation and compaction by defining the maximum value
	   for any symbol during transformation: q must be at least k-l; if q<=n,
	   compaction is guaranteed; if k-l>n, compaction is never done; if q is
	   INT_MAX, the maximum number of symbols are aggregated into one.

	   Output: Returns an integer j in the range 1...q representing the size of the
	   new alphabet. If j<=n+1, the alphabet is compacted. The global variable r is
	   set to the number of old symbols grouped into one. Only x[n] is 0.*/

	inline T transform(T *x, T *p, T n, T k, T l, T q)
	{
	   T b, c, d, e, i, j, m, s;
	   T *pi, *pj;

	   for (s=0, i=k-l; i; i>>=1)
		  ++s;                      /* s is number of bits in old symbol.*/
	   e=std::numeric_limits<T>::max()>>s; /* e is for overflow checking.*/
	   for (b=d=r=0; r<n && d<=e && (c=d<<s|(k-l))<=q; ++r) {
		  b=b<<s|(x[r]-l+1);        /* b is start of x in chunk alphabet.*/
		  d=c;                      /* d is max symbol in chunk alphabet.*/
	   }
	   m=(((T)1)<<(r-1)*s)-1;            /* m masks off top old symbol from chunk.*/
	   x[n]=l-1;                    /* emulate zero terminator.*/
	   if (d<=n) {                  /* if bucketing possible, compact alphabet.*/
		  for (pi=p; pi<=p+d; ++pi)
			 *pi=0;                 /* zero transformation table.*/
		  for (pi=x+r, c=b; pi<=x+n; ++pi) {
			 p[c]=1;                /* mark used chunk symbol.*/
			 c=(c&m)<<s|(*pi-l+1);  /* shift in next old symbol in chunk.*/
		  }
		  for (i=1; i<r; ++i) {     /* handle last r-1 positions.*/
			 p[c]=1;                /* mark used chunk symbol.*/
			 c=(c&m)<<s;            /* shift in next old symbol in chunk.*/
		  }
		  for (pi=p, j=1; pi<=p+d; ++pi)
			 if (*pi)
				*pi=j++;            /* j is new alphabet size.*/
		  for (pi=x, pj=x+r, c=b; pj<=x+n; ++pi, ++pj) {
			 *pi=p[c];              /* transform to new alphabet.*/
			 c=(c&m)<<s|(*pj-l+1);  /* shift in next old symbol in chunk.*/
		  }
		  while (pi<x+n) {          /* handle last r-1 positions.*/
			 *pi++=p[c];            /* transform to new alphabet.*/
			 c=(c&m)<<s;            /* shift right-end zero in chunk.*/
		  }
	   } else {                     /* bucketing not possible, don't compact.*/
		  for (pi=x, pj=x+r, c=b; pj<=x+n; ++pi, ++pj) {
			 *pi=c;                 /* transform to new alphabet.*/
			 c=(c&m)<<s|(*pj-l+1);  /* shift in next old symbol in chunk.*/
		  }
		  while (pi<x+n) {          /* handle last r-1 positions.*/
			 *pi++=c;               /* transform to new alphabet.*/
			 c=(c&m)<<s;            /* shift right-end zero in chunk.*/
		  }
		  j=d+1;                    /* new alphabet size.*/
	   }
	   x[n]=0;                      /* end-of-string symbol is zero.*/
	   return j;                    /* return new alphabet size.*/
	}
	
	public:

	/* Makes suffix array p of x. x becomes inverse of p. p and x are both of size
	   n+1. Contents of x[0...n-1] are integers in the range l...k-1. Original
	   contents of x[n] is disregarded, the n-th symbol being regarded as
	   end-of-string smaller than all other symbols.*/

	void suffixsort(T *x, T *p, T n, T k, T l)
	{
	   T *pi, *pk;
	   T i, j, s, sl;

	   V=x;                         /* set global values.*/
	   I=p;

	   if (n>=k-l) {                /* if bucketing possible,*/
		  j=transform(V, I, n, k, l, n);
		  bucketsort(V, I, n, j);   /* bucketsort on first r positions.*/
	   } else {
		  transform(V, I, n, k, l, std::numeric_limits<T>::max());
		  for (i=0; i<=n; ++i)
			 I[i]=i;                /* initialize I with suffix numbers.*/
		  h=0;
		  sort_split(I, n+1);       /* quicksort on first r positions.*/
	   }
	   h=r;                         /* number of symbols aggregated by transform.*/

	   while (*I>=-n) {
		  pi=I;                     /* pi is first position of group.*/
		  sl=0;                     /* sl is negated length of sorted groups.*/
		  do {
			 if ((s=*pi) <= 0 && (s=*pi) != 0) {
				pi-=s;              /* skip over sorted group.*/
				sl+=s;              /* add negated length to sl.*/
			 } else {
				if (sl) {
				   *(pi+sl)=sl;     /* combine sorted groups before pi.*/
				   sl=0;
				}
				pk=I+V[s]+1;        /* pk-1 is last position of unsorted group.*/
				sort_split(pi, (T)(pk-pi));
				pi=pk;              /* next group.*/
			 }
		  } while (pi<=I+n);
		  if (sl)                   /* if the array ends with a sorted group.*/
			 *(pi+sl)=sl;           /* combine sorted groups at end of I.*/
		  h=2*h;                    /* double sorted-depth.*/
	   }

	   for (i=0; i<=n; ++i)         /* reconstruct suffix array from inverse.*/
		  I[V[i]]=i;
	}
};

#endif /*def LS_H_*/
