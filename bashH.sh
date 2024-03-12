input=$1
file="/home/mgray7/output3/fasta/$input.txt";
file2="/home/mgray7/output3/structuresH/$input.txt";
exec 5<$file
exec 6<$file2
while read line1 <&5 && read line2 <&6; do

if ((i % 2 == 1))  

then
# /usr/bin/time -o out.txt -f "%e\t%M" ./build/Spark -P "/home/mgray7/Spark/src/params/parameters_DP09_Vienna.txt" -p -d1 -r $line2  $line1;
# ./build/Spark -P "/home/mgray7/Spark/src/params/parameters_DP09_Vienna.txt" -d1 -r $line2  $line1 > "/home/mgray7/Spark/out.txt";
/home/mgray7/HFold/HFold -s $line1 -r $line2 > "/home/mgray7/Spark/out.txt";


cat "/home/mgray7/Spark/out.txt" >> "/home/mgray7/output3/proof/HFold/$input.txt"
# cat "/home/mgray7/Spark/out.txt" >> "/home/mgray7/output3/proof/Spark/$input.txt"

# echo /usr/bin/time -o out.txt -f "%e\t%M" ./build/Spark -P "/home/mgray7/Spark/src/params/parameters_DP09_Vienna.txt" -p -d1 -r \"$line2\"  $line1;
# echo $line2
fi
i=$((i+1));
done
echo "first is $i";
#   echo "${line}";
exec 5</dev/null
exec 6</dev/null
