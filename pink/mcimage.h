/*
Copyright ESIEE (2009) 

m.couprie@esiee.fr

This software is an image processing library whose purpose is to be
used primarily for research and teaching.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software. You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/
/* $Id: mcimage.h,v 1.9 2006/02/28 07:49:12 michel Exp $ */
/* ============== */
/* prototypes for mcimage.c    */
/* ============== */

extern struct xvimage *allocimage(char * name, int32_t rs, int32_t cs, int32_t d, int32_t t);
extern void razimage(struct xvimage *f);
extern struct xvimage *allocheader(char * name, int32_t rs, int32_t cs, int32_t d, int32_t t);
extern int32_t showheader(char * name);
extern void freeimage(struct xvimage *image);
extern struct xvimage *copyimage(struct xvimage *f);
extern int32_t copy2image(struct xvimage *dest, struct xvimage *source);
extern int32_t equalimages(struct xvimage *im1, struct xvimage *im2);
extern void list2image(struct xvimage * image, double *P, int32_t n);
extern double * image2list(struct xvimage * image, int32_t *n);

extern void writeimage(
  struct xvimage * image,
  char *filename
);

extern void writese(
  struct xvimage * image,
  char *filename,
  int32_t x, int32_t y, int32_t z
);

extern void writelongimage(
  struct xvimage * image,
  char *filename
);

extern void writerawimage(
  struct xvimage * image,
  char *filename
);

extern void writeascimage(
  struct xvimage * image,
  char *filename
);

extern void printimage(
  struct xvimage * image
);

extern void writergbimage(
  struct xvimage * redimage,
  struct xvimage * greenimage,
  struct xvimage * blueimage,
  char *filename
);

extern struct xvimage * readimage(
  char *filename
);

extern struct xvimage * readheader(
  char *filename
);

extern struct xvimage * readse(char *filename, int32_t *x, int32_t *y, int32_t*z);

extern struct xvimage * readlongimage(
  char *filename
);

extern int32_t readrgbimage(
  char *filename,
  struct xvimage ** r,
  struct xvimage ** g,
  struct xvimage ** b
);

extern int32_t readbmp(
  char *filename, 
  struct xvimage ** r, 
  struct xvimage ** g, 
  struct xvimage ** b
);

extern void writebmp(
  struct xvimage * redimage,
  struct xvimage * greenimage,
  struct xvimage * blueimage,
  char *filename
);

extern int32_t readrgb(
  char *filename, 
  struct xvimage ** r, 
  struct xvimage ** g, 
  struct xvimage ** b
);

extern int32_t convertgen(struct xvimage **f1, struct xvimage **f2);
extern int32_t convertlong(struct xvimage **f1);
extern int32_t convertfloat(struct xvimage **f1);
