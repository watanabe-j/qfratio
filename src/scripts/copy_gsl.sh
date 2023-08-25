#! /bin/sh
copy_from_gsl_to_src(){
    ORIGIN_ROOT=$1
    DEST_ROOT=$2
    FILE_LIST=$3
    for FILE in `cat $FILE_LIST`
    do
        # PAR_DIR=${FILE%/*}
        # if [ -d $ORIGIN_ROOT/$PAR_DIR ]
        # then
        #     DEST_DIR=$DEST_ROOT/$PAR_DIR
        #     if ! [ -d $DEST_DIR ]
        #     then
        #         mkdir -p $DEST_DIR
        #     fi
        # fi
        # cp $ORIGIN_ROOT/$FILE $DEST_ROOT/$FILE
        install -D $ORIGIN_ROOT/$FILE $DEST_ROOT/$FILE
    done
}

copy_from_gsl_to_src $*
