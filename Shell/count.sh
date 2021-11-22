#!/bin/bash

for module in {Berghei,Chabaudi,Falciparum,Knowlesi,Malariae,Ovale,Vivax}
do
	cut -f1 rls_${module}/Final_${module}.tsv | sort | uniq -c | sed 's/^ //g' | awk '{print $1,$2}' > rls_${module}/Final_${module}.count
done
