## Haybaler Taxonomy

#### 1. Try to extract genus or species name from the reference sequence names.

reference |1_CP002843_1_Streptococcus_parasanguinis_ATCC_15912__complete_genome_BAC
----------|-------------------------------------------------------------------------
genus	  |Streptococcus
species	  |Streptococcus parasanguinis

Filter in functions find_genus and find_species work okay. Most of the time it is possible to identify genus and species for 90-95% of all chromosomes. Use --test_reference True to test how good a reference genome works

#### 2. Get the taxonomic ID from the names via name2taxid
Extracted name as input

reference name                  | filtered  name | name2taxid name |	name2taxid ID
--------------------------------|----------------|-----------------|------------------
CR382121.1_Kluyveromyces_lactis |Kluyveromyces   |Kluyveromyces    |	4910
FM992695.1_Candida_dubliniensis |Candida	 |Candida	   |	5475
FM992688.1_Candida_dubliniensis	|Candida	 |Candida   	   |	1535326
AE016814.2_Ashbya_gossypii	|Ashbya	         |Candida 	   |	5475
AE016815.5_Ashbya_gossypii	|Ashbya	         |Candida 	   |	1535326

The first row is how it should be. The next show a Problem:

Some chromosomes (e.g. Candida) are actually two organisms sharing the same genus name. Name2taxid knows both organisms and gives both taxIDs as output. That causes all following chromosome outputs to be shifted and we do not know which candida with which TaxID is the one from our reference. 

#### 3. Get the lineage via pytaxonkit.lineage
TaxID as input

reference name 			|filtered name |name2taxid name |name2taxid ID |lineage name  |lineage ID
--------------------------------|--------------|----------------|--------------|--------------|-----------
CR382121.1_Kluyveromyces_lactis |Kluyveromyces |Kluyveromyces   |4910 	       |Kluyveromyces |4910
AE016820.4_Ashbya_gossypii 	|Ashbya        |Ashbya 		|33170         |Eremothecium  |33170
1_CP017599_1_Moorea_producens 	|Moorea        |Moorea	        |619344        |Ianmoorea     |619344
x 				| x 	       |Moorea	        |1155738       |Moorea        |1155738

The first row shows how it should be. The next show a problem:

Sometimes lineage outputs a different name than what was the input for name2taxid. That is because names change for various reasons and are updated. So name2taxid accepts and knows the old names but lineage gives the new ones as output. If the name just changes, it does not matter. But if this is combined with the issue from before, it makes it harder to determine the right TaxID. 

#### 4. Filter for double TaxID for the same name

I. Filter for multiple rows with the same name but different TaxID’s

II. Filter for multiple lines in the input with the same name but different names in the output

Check if exact one of these double (sometimes even triple) chromosomes is a bacterium by using the lineage output.
If so, take the chromosome which is the bacteria, else it is not possible to identify which TaxID is the right one.

#### Info
Uncomment the last lines in the report function to get these tables with all different steps as columns


