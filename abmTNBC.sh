#!/bin/bash

for i in {0..9999}
do
   echo $i> NumRun.txt
   #sed -r 's/\s+//g' NumRun.txt
   #cat NumRun.txt | tr -d " \t\n\r"
   #tr -d "[:blank:]" < NumRun.txt | cat -n
   awk '{$1=$1}1' OFS= NumRun.txt
   awk '{gsub(" ","",$0);print }' NumRun.txt
   cp ABMparamsBefore.txt ABMparamsBefore$i.txt 
   cp ABMparamsAfter.txt ABMparamsAfter$i.txt 
   java -jar HAL-freq.jar
done
