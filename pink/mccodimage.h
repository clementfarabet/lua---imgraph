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
/* $Id: mccodimage.h,v 1.10 2006/04/24 15:07:28 michel Exp $ */
#ifndef SHRT_MIN
#define SHRT_MIN -32767 
#define SHRT_MAX +32767 
#define USHRT_MAX 65535 
#define INT_MIN -32767 
#define INT_MAX +32767 
#define UINT_MAX 65535 
#define LONG_MIN -2147483647 
#define LONG_MAX +2147483647
#define ULONG_MAX 4294967295
#endif
#define NDG_MAX 255            /* niveau de gris max */
#define NDG_MIN 0              /* niveau de gris min */

/* definitions for data storage type,
   uint32_t data_storage_type; */
#define	VFF_TYP_BIT		0	/* pixels are on or off (binary image)*/
                                        /* Note: This is an X11 XBitmap 
					   with bits packed into a byte and
					   padded to a byte */
#define	VFF_TYP_1_BYTE		1	/* pixels are byte (uint8_t) */
#define	VFF_TYP_2_BYTE		2	/* pixels are two byte (int16_t) */
#define	VFF_TYP_4_BYTE		4	/* pixels are four byte (integer) */
#define	VFF_TYP_FLOAT		5	/* pixels are float (single precision)*/
#define VFF_TYP_DOUBLE		9	/* pixels are float (double precision)*/
#define	VFF_TYP_COMPLEX		10      /* pixels are complex float */
#define VFF_TYP_DCOMPLEX	11	/* double complex */

#ifndef KHOROS
struct xvimage {
  char *name;
  uint32_t row_size;
  uint32_t col_size;
  uint32_t num_data_bands;	    /* Number of bands per data pixel,
					   or number of bands per image, or
					   dimension of vector data, or
					   number of elements in a vector */

  uint32_t data_storage_type;  /* storage type for disk data */
  double xdim, ydim, zdim;          /* voxel dimensions in real world */
  uint8_t imagedata[1];
};
#endif

#define UCHARDATA(I) ((uint8_t *)((I)->imagedata))
#define USHORTDATA(I) ((uint16_t *)((I)->imagedata))
#define ULONGDATA(I) ((uint32_t *)((I)->imagedata))
#define SLONGDATA(I) ((int32_t *)((I)->imagedata))
#define FLOATDATA(I) ((float *)((I)->imagedata))
#define DOUBLEDATA(I) ((double *)((I)->imagedata))
#define colsize(I)  ((I)->col_size)
#define rowsize(I)  ((I)->row_size)
#define depth(I)    ((I)->num_data_bands)
#define datatype(I) ((I)->data_storage_type)

/*            
		Codage du voisinage


		3	2	1			
		4	X	0
		5	6	7
*/
#define EST 0
#define NORD 2
#define OUEST 4
#define SUD 6
#define NORD_EST 1
#define NORD_OUEST 3
#define SUD_OUEST 5
#define SUD_EST 7
#define DEVANT 8
#define DERRIERE 10

#define nonbord(p,rs,N) ((p%rs!=rs-1)&&(p>=rs)&&(p%rs!=0)&&(p<N-rs))
#define nonbord3d(p,rs,ps,N) ((p>=ps)&&(p<N-ps)&&(p%ps>=rs)&&(p%ps<ps-rs)&&(p%rs!=0)&&(p%rs!=rs-1))

/* ============== */
/* prototypes     */
/* ============== */

extern int32_t voisin(int32_t i, int32_t k, int32_t rs, int32_t nb);
extern int32_t voisin2(int32_t i, int32_t k, int32_t rs, int32_t nb);
extern int32_t voisin6(int32_t i, int32_t k, int32_t rs, int32_t n, int32_t nb);
extern int32_t bord(int32_t i, int32_t rs, int32_t nb);
extern int32_t bord3d(int32_t i, int32_t rs, int32_t ps, int32_t nb);
extern int32_t voisin26(int32_t i, int32_t k, int32_t rs, int32_t n, int32_t nb);
extern int32_t voisin18(int32_t i, int32_t k, int32_t rs, int32_t n, int32_t nb);
extern int32_t voisins4(int32_t i, int32_t j, int32_t rs);
extern int32_t voisins8(int32_t i, int32_t j, int32_t rs);   
extern int32_t voisins6(int32_t i, int32_t j, int32_t rs, int32_t ps);   
extern int32_t voisins26(int32_t i, int32_t j, int32_t rs, int32_t ps);   
extern int32_t voisin5(int32_t i, int32_t k, int32_t rs, int32_t nb);
extern int32_t voisin6b(int32_t i, int32_t k, int32_t rs, int32_t nb, int32_t par);
extern int32_t voisinNESO(int32_t i, int32_t k, int32_t rs, int32_t nb);
extern int32_t voisinNOSE(int32_t i, int32_t k, int32_t rs, int32_t nb);
extern int32_t voisin14b(int32_t i, int32_t k, int32_t rs, int32_t ps, int32_t N);
extern int32_t voisinONAV(int32_t i, int32_t k, int32_t rs, int32_t ps, int32_t N );
extern int32_t voisinENAR(int32_t i, int32_t k, int32_t rs, int32_t ps, int32_t N );
extern int32_t voisinENAV(int32_t i, int32_t k, int32_t rs, int32_t ps, int32_t N );
extern int32_t voisinONAR(int32_t i, int32_t k, int32_t rs, int32_t ps, int32_t N );
extern uint32_t maskvois26(uint8_t *F, uint32_t bitmask, int32_t i, int32_t rs, int32_t ps, int32_t N);
