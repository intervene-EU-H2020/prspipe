
mounts() {
    echo "/dhc/ /dhc/" 
}

hooks() {
    echo "wd=${PWD}" >> ${ENROOT_ENVIRON}
}

rc(){
    cd $wd && bash
}
