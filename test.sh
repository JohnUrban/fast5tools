for dir in 01 02 03 04 05 06 07 08 09 10 11 11b 12 13 14 15 16; do
 echo rundata/$dir
 python testf5class.py rundata/$dir/example.fast5
 echo
 echo
done
