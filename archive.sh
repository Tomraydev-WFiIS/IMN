# if no arguments
if [[ $# -eq 0 ]] ; then
    echo 'Provide exactly one argument <lab_number>'
    exit 1
fi

tar -czvf IMN_Rajchel_Tomasz_lab_$1.tar.gz lab0$1