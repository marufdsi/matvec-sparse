#for i in 36 32 18 16 9 8 4 2 1
for i in 35 34 33 31 30 3 5 6 7 10 11 12 13 14 15 17 19 20 21 22 23 24 25 26 27 28 29 
do 
	./runScatter.sh $i $1 $2 $3 0 $4
	./runScatter.sh $i $1 $2 $3 1 $4
	#./runCompact.sh $i $1 $2 $3 0 $4
	#./runCompact.sh $i $1 $2 $3 1 $4
	#./runBalance.sh $i $1 $2 $3 0 $4
	#./runBalance.sh $i $1 $2 $3 1 $4
done

