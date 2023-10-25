#get last n lines of matrix file
sed '$d' $1 | tail -n $2
