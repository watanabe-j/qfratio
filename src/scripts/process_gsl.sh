#! /bin/sh

## Shell script to process code files copied from GSL


## Replace '#include <gsl/...>' with '#include "..."' using relative path,
## so that include option in Makevars can just be 'PKG_CPPFLAGS = -I.'
localize_header_inclusion() {
    INDIR=$1

    ## First process headers in GSL root dir
    GSL_ROOT_HEADERS=`cd $INDIR && ls *.h`
    for HEADER in $GSL_ROOT_HEADERS
    do
        sed -i -E "s|#include *<gsl/$HEADER>|#include \"$HEADER\" // edited for qfratio|"    $INDIR/*.h
        sed -i -E "s|#include *<gsl/$HEADER>|#include \"../$HEADER\" // edited for qfratio|" $INDIR/*/*
    done
    
    ## Then loop across subdirectories
    SUBDIRS=`cd $INDIR && ls -d */`
    for DIR_ in $SUBDIRS;
    do
        ## Remove trailing slash from DIR_ (though not strictly necessary)
        DIR=${DIR_%/}
        SUBDIR_HEADERS=`cd $INDIR/$DIR && ls *.h`
        for HEADER in $SUBDIR_HEADERS
        do
            sed -i -E "s|#include *<gsl/$HEADER>|#include \"$DIR/$HEADER\" // edited for qfratio|"    $INDIR/*.h
            sed -i -E "s|#include *<gsl/$HEADER>|#include \"$HEADER\" // edited for qfratio|"         $INDIR/$DIR/*
            sed -i -E "s|#include *<gsl/$HEADER>|#include \"../$DIR/$HEADER\" // edited for qfratio|" $INDIR/*/*
        done
    done
    
    # ## These will replace every GSL header inclusions, but all subdirs should
    # ## be specified in Makevars with 'PKG_CPPFLAGS = -I. -I./gsl/err ...'
    # sed -i -E 's/#include *<gsl\/(.+)>/#include "\1"/' $INDIR/*.h
    # sed -i -E 's/#include *<gsl\/(.+)>/#include "\1"/' $INDIR/*/*
}

## Manually comment out unnecessary lines
comment_out_lines(){
    INDIR=$1
    DIR_ERR=$INDIR/err
    DIR_INT=$INDIR/integration
    DIR_PLY=$INDIR/poly
    DIR_RTS=$INDIR/roots
    DIR_SPF=$INDIR/specfunc

    ## err/error.c: remove '#include <gsl/gsl_message.h>'
    sed -i -E  '26s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_ERR/error.c

    ## err/gsl_errno.h: remove stream-related types and functions
    sed -i -E   '80,81s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_ERR/gsl_errno.h
    sed -i -E   '88,89s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_ERR/gsl_errno.h
    sed -i -E  '97,100s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_ERR/gsl_errno.h

    ## integration/gsl_integration.h:
    ## leave only *_workspace*, *_qk15, *_qk, *_qagil
    sed -i -E  '87,117s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_INT/gsl_integration.h
    sed -i -E '131,152s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_INT/gsl_integration.h
    sed -i -E '176,192s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_INT/gsl_integration.h
    sed -i -E '200,405s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_INT/gsl_integration.h

    ## integration/qags.c:
    ## leave only qags, gsl_integration_qagiu, iu_transform
    sed -i -E '39,141s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_INT/qags.c

    ## poly/gsl_poly.h: leave only gsl_poly_eval
    sed -i -E      '25s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_PLY/gsl_poly.h
    sed -i -E   '50,56s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_PLY/gsl_poly.h
    sed -i -E  '69,102s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_PLY/gsl_poly.h
    sed -i -E '105,179s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_PLY/gsl_poly.h
    # sed -i -E '144,146s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_PLY/gsl_poly.h
    # sed -i -E '158,161s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_PLY/gsl_poly.h
    # sed -i -E '164,179s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_PLY/gsl_poly.h

    ## roots/convergence.c: leave only gsl_root_test_delta
    sed -i -E   '25,57s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_RTS/convergence.c
    sed -i -E   '76,86s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_RTS/convergence.c

    ## roots/gsl_roots.h: remove fdfsolver-related parts
    sed -i -E   '59,75s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_RTS/gsl_roots.h
    sed -i -E  '93,107s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_RTS/gsl_roots.h
    sed -i -E '121,123s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_RTS/gsl_roots.h

    ## specfunc/bessel.c: remove #include gsl_sf_airy.h
    sed -i -E '26s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_SPF/bessel.c

    ## specfunc/exp.c, gsl_sf_exp.h: remove gsl_sf_exp_mult_e10_e()
    sed -i -E '184,227s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_SPF/exp.c
    sed -i -E   '64,68s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_SPF/gsl_sf_exp.h

    ## specfunc/gamma.c, gsl_sf_gamma.h: remove gsl_sf_choose_e()
    sed -i -E '1584,1625s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_SPF/gamma.c
    sed -i -E '1677,1680s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_SPF/gamma.c
    sed -i -E   '156,161s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_SPF/gsl_sf_gamma.h

    ## specfunc/gsl_sf_hyperg.h: remove *_0F1 and *_2F0
    sed -i -E   '40,48s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_SPF/gsl_sf_hyperg.h
    sed -i -E '142,149s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_SPF/gsl_sf_hyperg.h

    ## specfunc/psi.c: remove complex-related parts
    sed -i -E      '30s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_SPF/psi.c
    sed -i -E '458,552s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_SPF/psi.c
    sed -i -E '797,830s/^(.+)/\/\/ \1 \/\/ edited for qfratio/' $DIR_SPF/psi.c
}

## Use R's error handler
edit_error_handler(){
    INDIR=$1
    DIR_ERR=$INDIR/err

    sed -i -E  '41,47s/^( *)(.+)/\1\/\/ \2 \/\/ edited for qfratio/' $DIR_ERR/error.c
    sed -i -E  '48i\  error("Problem in %s: %d: %s\\n  This is unexpected; contact qfratio maintainer", file, line, reason); // added for qfratio' $DIR_ERR/error.c
}

## Add necessary preprocessor directives
add_directives(){
    INDIR=$1
    DIR_ERR=$INDIR/err
    DIR_INT=$INDIR/integration
    DIR_SPF=$INDIR/specfunc

    ## err/error.c: use R's error hander
    sed -i -s '24i#include <R_ext/Error.h> // added for qfratio' $DIR_ERR/error.c

    ## integration/qelg.c: requires several include directives
    sed -i -s '20i#include <stddef.h> // added for qfratio\n#include <math.h> // added for qfratio\n#include "../gsl_machine.h" // added for qfratio\n#include "../gsl_minmax.h" // added for qfratio\n' $DIR_INT/qelg.c

    ## integration/..: requires "integration.h"
    sed -i -s '20i#include "gsl_integration.h" // added for qfratio\n' $DIR_INT/initialise.c $DIR_INT/qpsrt.c $DIR_INT/util.c
    sed -i -s '24i#include "gsl_integration.h" // added for qfratio\n' $DIR_INT/qpsrt2.c
    sed -i -s  '1i#include "gsl_integration.h" // added for qfratio\n' $DIR_INT/reset.c $DIR_INT/set_initial.c

    ## integration/util.c: requires <math.h>
    sed -i -s '20i#include <math.h> // added for qfratio' $DIR_INT/util.c

    ## integration/util.c,qags.c: move #include "qpsrt.c" to avoid warning for implicit declaration
    sed -i -s '22i#include "qpsrt.c" // added for qfratio' $DIR_INT/util.c
    sed -i -E 's|#include "qpsrt.c"|// #include "qpsrt.c" // edited for qfratio|' $DIR_INT/qags.c

    ## integration/positivity.c: requires <math.h> and "../gsl_machine.h"
    sed -i -s  '4i#include <math.h> // added for qfratio\n#include "../gsl_machine.h" // added for qfratio\n'  $DIR_INT/positivity.c

    ## specfunc/chebyshev.h: requires #ifndef ... #endif directives
    sed -i -s '22i#ifndef __GSL__CHEBYSHEV_H__ // added for qfratio\n#define __GSL__CHEBYSHEV_H__ // added for qfratio\n' $DIR_SPF/chebyshev.h
    sed -i -s '$i\\n#endif // added for qfratio\n' $DIR_SPF/chebyshev.h

    ## specfunc/cheb_eval.c: requires several include diriectives
    sed -i -s '1i#include <math.h> // added for qfratio\n#include "chebyshev.h" // added for qfratio\n#include "gsl_sf_result.h" // added for qfratio\n#include "../gsl_machine.h" // added for qfratio\n#include "../err/gsl_errno.h" // added for qfratio\n' $DIR_SPF/cheb_eval.c
}

## Add license statements
add_license_statement(){
    INDIR=$1
    sed -i -s '1i// This file is taken from GSL version 2.8 and distributed as part of qfratio\n// with modification, in accordance with the GNU General Public License\n// version 3.  All modified lines are marked with comments.\n// - 2023 Junya Watanabe\n' $INDIR/*\.h $INDIR/*/*\.c $INDIR/*/*\.h
}


localize_header_inclusion $*

comment_out_lines $*

edit_error_handler $*

add_directives $*

add_license_statement $*
