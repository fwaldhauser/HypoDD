/*
 * $Header: /home1/crhet/julian/hypoDD/RCS/compat.h,v 1.1 2001/02/13 21:21:10 julian Exp $
 * ANSI-traditional_C compatibility macros
 */
#ifndef COMPAT_H
#define COMPAT_H

#ifdef __STDC__
typedef void* GENPTR;
#define PARMS(x) x
#else /* __STDC__ */
typedef char* FENPTR;
#define PARMS(x) ()
#define const	/* Nothing */
#endif /* __STDC__ */
#endif /* COMPAT_H */
