/* /proj/stdevs/stdev0f/SLIBS/alr.dev/SCCS/s.chanmatfuns.h */
/* alr support chanmatfuns.h 3.2 97/03/23 */


extern MATRIX *create_matrix(),
     *matcopy(),
     *extract_rows(),
     *matadd(),
     *matsub(),
     *matmult(),
     *transp(),
     *Cchol(),
     *sweep(),
     *col_1s(),
     *matabs(),
     *matexp(),
     *px1_times_pxq(),
     *pxq_divby_px1(),
     *scalar_times_matrix(),
     *ident(),
     *form_diag(),
     *diaginv(),
     *diagsqrt(),
     *corner(),
     *covlag(),
     *toeplitz(),
     *band(),
     *matcumnorm(),
     *extract_cols(),
     /* following two functions added by pj catalano  */
     *matnpdf(),
     *matncdf();

extern double matmax(), elsum();

extern void matdump(), plug(), destroy_matrix(), fmatdump();

extern MATRIX *luinv();


