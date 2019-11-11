# if no arguments
if [[ $# -eq 0 ]] ; then
    echo -e 'Provide exactly one argument <lab_number> eg:\n./archive.sh 1'
    exit 1
fi

tar -czvf IMN_Rajchel_Tomasz_lab_$1.tar.gz lab0$1