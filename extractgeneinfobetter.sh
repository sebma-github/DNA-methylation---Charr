#!/bin/bash

## This script allows me to extract the methylation information from all 48 sample files, at a specific 1000bp DMR in a specific chromosome
## It was written to test my theory that global methylation over a 1000bp region could still be interpreted as 
## semi-valid proof of the differences between morphs at different times. Or at least construct graphs.
##Needs to be run in the directory where the .cov files are

#for i in {2..389}; #For DMRs between stages 
for i in {2..141}; #For DMRs between morphs
#for i in {2..15};   #For DMRs-qPCR genes

do
	#gene=((i-1)) doesn't matter
        name=$(awk '{ if(NR=='$i') print $1 }' /home/sebastien/tablesforextraction/DMRsbetweenmorphsforextraction.tsv)
        scaff=$(awk '{ if(NR=='$i') print $2 }' /home/sebastien/tablesforextraction/DMRsbetweenmorphsforextraction.tsv)
        postemp=$(awk '{ if(NR=='$i') print $3 }' /home/sebastien/tablesforextraction/DMRsbetweenmorphsforextraction.tsv)
	pos=$(awk '{ if(NR=='$i') print $4 }' /home/sebastien/tablesforextraction/DMRsbetweenmorphsforextraction.tsv)
		
		#Enter the scaffold number related to the DMR of interest
		for f in ./*.cov; 
		do
			grep -n "$scaff" "$f" > ""$f"temp" 
		done

		#Enter the lines which match the coordinates indicated by the DMR of interest       
		#for f in ./*.covtemp;     							#OLD VERSION, DO NOT USE
		#do
		#	awk '{if($2~/^'$pos'.../){print $0} else {}}' "$f" > ""$f"final"
		#wait
		#	mv ""$f"final" "$f"
		#done
		#wait

		#Keep all lines that are above the start of the windowtile
                for f in ./*.covtemp;
                do
                        awk '{if($2>='$postemp'){print $0} else {}}' "$f" > ""$f"final"
                wait
                        mv ""$f"final" "$f"
                done
                wait

		#Keep all lines that are below the end of the windowtile
                for f in ./*.covtemp;
                do
                        awk '{if($2<='$pos'){print $0} else {}}' "$f" > ""$f"final"
                wait
                        mv ""$f"final" "$f"
                done
                wait

		#Removes from those files all the columns that I am not interested in (1, 2 and 4)
		for f in ./*.covtemp; 
		do
			cut -f1,2,4 --complement "$f" > ""$f"cut"
		wait
			mv ""$f"cut" "$f"
		done
		wait

		#Add headers so I can create a comprehensible df (will change the names later when I am in R)
		for f in ./*.covtemp; 
		do
			x=${f%best.covtemp}
			y=${x#./*}
			head=( POS ""$y"meth" ""$y"unmeth" )
			( IFS=$'\t'; echo "${head[*]}"; cat "$f" ) > ""$f"head"
		wait
			mv ""$f"head" "$f"
		done
		wait

		#Should add a function to remove all lines where meth + unmeth <10. This prevents me from having to deal with cytosines covered <10X when making the graphs.
		#Considering that the DMRs were calculated based on cytosines covered >10X, it makes sense to only have those in the graphs as well.
		#However, that means that I need to work around NAs 
		for f in ./*.covtemp; 
                do
			awk '$2+$3<10{next}1' "$f" > "$f""cov"
		wait
			mv ""$f"cov" "$f"
		done
	

		#Merge all the files together based on their first column. This is a bad way of doing things because I should't have one file as a base, as I will only join positions that
		#are only present in the first file.
		#UNLESS I create a dummy file with a vector of every one of the thousand numbers. Then I am not losing any information.
		
		#Put all files one after the other
		cat *.covtemp > temp.covmodel
		#Keep only the position column
		awk '{print $1}' temp.covmodel > temp2.covmodel
		#Sort the positions and remove duplicates
		sort -k 1,1 -u temp2.covmodel > temp3.covmodel
		#Add 		
		#cp temp3.covmodel temp3.test #for testing


		#BEfore joining all the files together, I meed to replace the empty values with some string eg "NA"
		#Otherwise if I join all the files, some of the numbers might get shifted to the left.
		#To do so, need to join temp3.covmodel to each file individually.

		for f in ./*.covtemp;
		do
			join <(sort temp3.covmodel) <(sort "$f") -a 1 | column -t | tac > individual.cov
		wait
			mv individual.cov "$f"
		done

		#Still that doesn't solve the problem because there ISNT any blank space in the first place !
		#So  first let's retransform the data with 3 columns. 

		for f in ./*.covtemp;
		do
			awk 'BEGIN{ OFS="\t"}{ print $1, $2, $3}' "$f" > "$f""goodcolumns"
		wait
			mv "$f""goodcolumns" "$f"
		done			

		#Now that there are 3 columns I can replace the blank space with missing values.
		#The followimg line replaces the blank with 0 but I would like to replace it by a string.

		for f in ./*.covtemp;
		do
			awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }; 1' "$f" > "$f""NA"
		wait
			mv "$f""NA" "$f"
		done


                # Now I can join them together. 
                for f in ./*.covtemp; 
                do
                       join <(sort temp3.covmodel) <(sort "$f") -a 1 | column -t | tac > temp.cov
                       mv temp.cov temp3.covmodel
                done
                wait
		
		#Now in order for me to work on the data in R, I need to be sure that everything is tab sep.
		#So I recreate all the columns manually

		awk 'BEGIN{ OFS="\t"}{ print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44, $45, $46, $47, $48, $49, $50, $51, $52, $53, $54, $55, $56, $57, $58, $59, $60, $61, $62, $63, $64, $65, $66, $67, $68, $69, $70, $71, $72, $73, $74, $75, $76, $77, $78, $79, $80, $81, $82, $83, $84, $85, $86, $87, $88, $89, $90, $91, $92, $93, $94, $95, $96, $97}' temp3.covmodel > ./DMRtest/"$name".covdf


		#I think that works
		#mv temp3.covmodel ./DMRtest/"$name".covdf

                wait
		
		rm *.covtemp temp.covmodel temp2.covmodel temp3.covmodel
		echo "Genes extracted: $i-1"
done

echo "Information about all genes has been extracted"
