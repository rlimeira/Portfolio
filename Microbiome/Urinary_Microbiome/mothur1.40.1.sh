#****************START MOTHUR VERSION 1.40.1 SCRIPT*****************
# Author: Roberto Limeira, derived from Scripts by Gina Kuffel, Krystal White, and Megan Pierce
# Last updated: 04/26/2018

#BEFORE YOU BEGIN
#1. Log into mongo1 server ssomsmp1.luhs.org
#2. Using your command line terminal of choice navigate to your working directory using the cd command. 
#		a. Type "cd" space and then the directory of choice
#			i. TIP: capitalization matters, spaces matter, do not include spaces in your directory names or file names, use instead the underscore (_)
#			ii. To navigate back to your home directory just type "cd" and hit enter
#			iii. To navigate back a file use "cd ../"
#			iv. To see the files in the directory you are in use the "ls" command 
#3. Transfer your .fastq.gz or .fastq files to your working directory by using the cp command. 
#		a. when you are in your destination directory type "cp" space then the entire path of the files you want to copy with a "/*.gz" or "/*.fastq" at the end of the path, followed by another space and "./". The * will transfer all files ending with those characters. The ./ will transfer those files into the directory you are in. 
#			i. EXAMPLE: I am in my /data/rlimeira/Studies directory and I want to copy all fastq.gz files from the /scratch/microbiome directory ignoring all other files. I then enter "cp /scratch/microbiome/*.gz ./"
#4. Transfer your core mothur files over to your working directory using the cp command. 
#		a. Type "cd /scratch/microbiome/Core_Mothur_files/mothur1.40/* ./" This will copy all files within that directory (denoted by * with a blank after it), to your current directory (denoted by ./)
#5. Run screen (screen will allow you to close your terminal window and let the program run in the background). This may be done at any time prior to starting your mothur run (step 7)
#		a. Check to make sure you don't have an existing screen session running by typing "screen -ls" 
#			i. every time you type "screen" a new screen session will begin, if you do not type "exit" in a screen session and just 
#			detatch by pressing control a, then control d, the screen session will still be active. Don't have multiple screen versions 
#			open without a need. This could cause unecessary problems later on.  
#		b. Run screen by typing "screen", or continue your previous screen session by typing "screen -r"
#6. Create the sability files by typing "bash Make_Stability_Files.sh"
#		a. The command should show on screen the contents of the stability.files (i.e. results of this script). double check that your stability files worked. They should have your sample names on the left column followed by the forward R1 and backwards R2 fastq file names. This is important for the first command of the script below. 
#7. Run mothur by typing "mothur mothur1.40.1.sh" 
#		a. Check to make sure that your stability.files have worked by making sure no errors appear on the screen for a few seconds. 
#8. Detatch screen by pressing control a and control d
#9. Exit the server by typing "exit" 

#AFTER SCRIP HAS FINISHED (A few minutes to hours to a day or so depending on how many fastq files you have and how many sequences are in these files.)
#1. Log into mongo1 server
#2. type screen -r to resume your previous screen session
#3. check that mothur ran ok by making sure that there are no errors on the summary that you see. There may be warnings. That is ok. 
#		a. If you have any issues you can open the latest mothur log file. 
#			i. Type ls and enter you will get a list of all the files in the directory. 
#			ii. Find the file starting with mothur and ending with logfile with a series of numbers in between. These numbers are linux time which is a cronological count forward from a data back in the 70s. The file with the largest number is the most recent. mothur creates a new file every time you open the aplicaiton by typing in "mothur" on command line. 
#4. Exit screen by typing "exit"
#5. Exit command line by typing "exit"
#6. Run the R script OTUTaxRelab_RL041918.R for relative abundance file with OTU designations. 

#Things to keep in mind. 
#1. The purpose of the majority of the commands is to decrease PCR errors, chimeras, and sequencing errors while paring down numbers of sequences so that the entire process is run relatively quickly. Some commands below are specifically placed so that the total number of sequences analyzed are decreased. At most steps, sequences removed are saved as a separate file. So it is expected that you will have working files and files containing your discarded sequences so you can troubleshoot. 
#2. Some settings were chosen by authors above to fit our amplicon length, region sequenced, and to fall within the lowest error rate possible. P. Schloss' group worked on determining some of the other settings through experimentation.
#3. Throughout your commands there are two main files as outputs from mothur: the fasta file and the count file. Fasta will have your sequences with their names, the count files will have the sequence names and your samples and the count of those sequences within the samples. Other files of interest include the .taxa files that the taxonomic identification of your OTUs and the summary files, that will summarize a particular step. The output from a particular command is listed at the end of the command in the logfile. Whether or not you can open a file generated from mothur is dependent on your computer and program's limitations to file size. So, check the file size before trying to open a file. 

#STATS: According to P. Schloss the starting error rate is %0.06

make.contigs(file=stability.files, processors=8)

#This will read each fastq file in the stability.files and make the contigs. Contigs are contiguous sequences in one direction derived from the forward (R1) and backwards (R2) fastq files you receive from the sequencer. 

summary.seqs(fasta=stability.trim.contigs.fasta)

#This gives you the summary of what is in your contigs files. Here you are able to see the percentage of seuqences that fall within the range of length, ambiguity or polymers in bases. 

screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, summary=stability.trim.contigs.summary, maxambig=0, minlength=275, maxlength=300)

#This eliminates or screens the sequences based on the parameters entered. These parameters are longer than P. Schloss' MiseqSOP parameters (maxlength=250) because he is performing a single PCR step whereas we are performing two: amplicon and library prep PCRs. 

unique.seqs(fasta=stability.trim.contigs.good.fasta)

#This further pares down the sequences to only unique sequences increasing speed of downstream processing.

count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)

#Counts the number of sequences in order to check if the unique.seqs command previously worked in these files and so you can get an idea how many repeated sequences you had. This also produces the count table. 

summary.seqs(count=stability.trim.contigs.good.count_table)

#This gives you the summary of your sequences again so that you may check against the last summary.seqs command and see how many seuqences you've eliminated. You can also check that you've removed sequences with ambiguities and 

pcr.seqs(fasta=silva.nr_v123.align, start=11894, end=25319, keepdots=F, processors=8)

#Takes the silva database and narrows it down (trims it) to the sequence region you want. You can think of it as performing an amplification PCR on the silva database wihtout the multiple copies of identical sequences. 

rename.file(input=silva.nr_v123.pcr.align, new= silva.v4.fasta)

#Renames a file to shorten the commands a bit. 

summary.seqs(fasta=silva.v4.fasta)

#Gives you the summary of the database file you just created so you can check the length. 

align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.v4.fasta)

#Aligns your sequences to the database region file you created above. 

summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table)

#Gives you the summary of the sequences once more to show you how long the sequences are once you aligned them to the database file. The sequences are very long because some sequences aligned at one end in the region and others to another end. They are aligned with the bases in columns and with periods marking blank spaces or spaces needing a base. Sequences at this step will look somewhat like .......AGTCGAATC...... (only much larger). Most of the sequences should have aligned to the correct location. Meaning that the dots will all be in the same region, same as the nucleotides, when ligned up one on top of each other (This is what the ...unique.align file looks like if you could open it. Don't try it. It's too large.). You can further narrow down to your sequences to the region of interest by using the screen.seqs and filter.seqs comands below. 

screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary=stability.trim.contigs.good.unique.summary, start=1968, end=11550, maxhomop=8)

#Removes outliers from your sequences. 

summary.seqs(fasta=current, count=current)

#Summarizes your sequences so you can see how many were eliminated with the screen.seqs commands above. 

filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)

#Filters out the blanks in the sequences (trump=.). Drastically reducing the size. Shows the output after the command.

unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)

#Keeps only the unique sequences. There might have been a base outside the aligment when we ran align.seqs that when you run filter.seqs it removes those bases. This then could make that sequence identical to another sequence. This removes that redundance. 

pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=2)

#Further narrows down sequence numbers by clustering sequences together based on number of diffs (nucleotides). In this command the sequence with the most abundant number of reads is trusted, the second sequence with the most abundant reads is compared to the first one. If it is within 2 nucleotides of difference they are clustered into one sequence. This goes down the line repeating with consecutive lower abundance reads. In this comand you lose a bit of resolution in sequences in return for drastic reduction in sequence number. This 2 diffs account for 1.6% of the 250nt reads. It is good to keep this percentage below 3% if your OTU definition will be at 3%. 

#STATS: Theoretical Error Rate at this step is: 0.02% from 0.06%

chimera.vsearch(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)

#Searches for chimeras within the sequences and eliminate them. Dereplicate=t means that if the command finds a sequence to be chimeric within a group it will not flag the sequence within all groups as chimeric. 

#What are chimeras and why are they a big bioinformatic problem? 
#As polymerase is amplifying your region of interest it might fall off before it finishes (this could be because of short PCR times, temperature variations between different sequences, sequences that make the polymerase fall off for some weird reason etc.). The polymerase doesn't just randomly fall off there are chimeras that occur often meaning that the polymerase are often falling off at the same spot. So you can have chimeras that make up a large amount of your total sequences. It is not uncommon, at the end of remove.seqs below for your total number of sequences to drop off by ~20%. However, the number of unique sequences removed should be much less ~<1%. 
#This incomplete amplicon DNA will become primer for another piece of DNA. If it is for a different region of the DNA we cath it at steps above. However, it could be from the amplicon region that you are trying to amplify, only from another species. So now you have a piece of amplicon which part of it belongs to one species and the other part belongs to another species. You could even have tripple chimeras or quadruple. These are tough to discern completely and is the biggest source of errors or mistakes with this whole process. Chimeras are also the reason that it is important to use good polymerases and reagents. 

#Simple explanation of how this function works.
#There are many aproaches to finding chimeras. The way this particular function works is thorugh the denovo method. The most aboundant sequence is most likely not a chimera, neither is the second. You then take the third most abundant sequence, cut it in half and compare the two ends to the first two. If it is flagged as chimera it is tossed, if not, it is added to the reference pool and the process is repeated.

#STATS: Sensitiviy is 88% and specificity is 94% so our false discover rate (FDR)= 2.6% (2.6% of thigs we call chimeric arent chimeric), Final percent of unique sequences that are chymeric is about 25% (but if you factor in the counts it is ~1%) 10 to 20% of all your dataset will be removed as falgged for chimera.

remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

#Removes the sequences flagged as chimeras. 

summary.seqs(fasta=current, count=current)

#Gives you the summary of your sequences so that you can see how many sequences flagged as chimeras were removed. the "current" function returns the current file mothur has been working with. If you are running a script you dont need to wory about this but if you are entering the comands through experimentation or troubleshooting you can check what your current files are with the get.current() command. Also, if you want to set a current file or processors you can use the set.current() command. 

classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=trainset16_022016.pds.fasta, taxonomy= trainset16_022016.pds.tax, cutoff=80)

#Classifies your sequences by matching them against the trainsets you supplied earlier (part of the core mothur files). This is important for the final classify OTU function below and so that we can remove non-bacteria sequences in the next step below. 

#How are sequences classified in this command? 
#There are many ways to classify sequences. This command uses the classification technique in the Wang paper with a confidence cutoff of 80%. The Wang method uses the kmer search method, with default kmer size of 8, and a Bayesian approach. Kmer of 8 means it randomly splits the database and your unknown sequences in 8 nucleotide long strands and builds a profile of them. Bayesian part comes when it then calculates probability of getting your unknown sequence kmer profile based on the kmer profile for each genus. Highest probability is the match right? Not quite, it will then calculate a confidence score by randomly selecting a subset of kmers for your sequence, randomly generate a new profile and based on that will calculate a probability for each genera, (does this 100 times "bootstrapping procedure"). Then based on the results of the bootstraping it counts the number of time the genus comes up. It then gives the percent of bootstraps that result in that genus and requires it to be >80% in order to classify it as that genus. If it falls below it goes up to the next level.

remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

#Since we are looking at only bacteria it removes any sequences that might have been matched to the described taxons. 

summary.tax(taxonomy=current, count=current)

#Gives you the summary of the sequences for you to check how many have been removed since last time you checked. 

cluster.split(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, splitmethod=opti, cutoff=0.03)

#This command will cluster sequences into OTUs based on a 3% cuttoff. This basically means that your sequences of ~300bp will cluster as part of an OTU if they have differences of 9 basepairs between each other (well, not quite, see below). 

#How to think about OTUs when looking at the data
#These OTUS could be thought of as species or even strain differences depending on the species, the size of the DNA you are analyzing. The depth in taxa with which you can classify your OTUs will depend on those factors as well as your classification method and database quality. Schloss' point of view is to talk about his findings as OTUs and only bring out taxonomic names to please clinicians. This is because taxonomy is not necessarily a mathematically objective art (sometimes major differences in genomes make the difference between taxonometric groups, sometimes its minor differences). Taxonomy is different depending on who you ask. And it is fluid; it is constantly changing. 

#How are the sequences clustered into OTUs with this command? 
#	Lets say we will be classifying OTUS at 3% (which is what we are doing here. Keep this in mind.). 
#	There are two aproaches to cluster a dataset: heuristic and higherarchical. Higherarchical is slower but more acurate and that is what this uses. In the higherarchical aproach one could go by three methods; furthest neighbor, nearest neighbor, and average. Furthest neighbor is when all the sequences in the OTU (or bin) is within 3% of all other sequences in the OTU. This leads to small groups and if you look at it graphically it appears as tight clusters with things on the outside. Nearest neighbor is when each sequence in OTU is within 3% of at leas one other sequence in OTU. This leads to long chains as groups. Average is when on average all sequences are within 3% of the other sequences in the OTU. Average appears ameboid in graphical representations. Now, which of these do we use? 
#	Well, looking at the quality of bin clustering there could be 4 types of occurences: true positive= two sequences are in the same OTU and are <3% different from each other; true negative= two sequences in different OTUs and are >3% different from each other; false postiives= two sequences that are more than 3% different but were placed in the same OTU; false negatives= two sequences that are <3% different but are in different OTUs. Using the Matthews correlation coefficient (MCC) you can determine which of the three (furthest, nearest, and average) is the best one to use. Turns out furthest neighbor has a lot of false negative, but no false positives, and nearest neighbor has a lot of false positives, but no false negatives, average seems to be the best with the MCC. Mothur uses the average method but with a few futher steps in the calculation. Why did you tell me this convoluted story, dude? 
# Well, the new split method (OptiClust), that came out in version 1.39, uses the MCC. It randomly joins sequences or OTUs with each other and looks at the change in the MCC between them. It keeps doing this until the change in MCC is zero. Schloss found that this method is better than average neighbor alone (used in previous versions) and it uses less RAM. 

make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)

#This makes the shared file that has the OTU listed in the first column and your samples as the top row and the count of each OTU per sample fillin in the matrix. 

classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=0.03)

#Takes your OTUs and classifies them against the sequence taxonomy file derived in the classify.seqs command above. 

#STATS: Roberto ran the seq.error command below which calculates the sequencing error rate based on known sequences from a purchased mock community. These were from sequence run 2018-1 and were whole cell mock community standards run through our entire extraction protocol and DNA mock community standard run through either the extraction method or added at the PCR amplification step. The results were great! Error rates were very low for all mock community variations ~ 6x10^-6. 

#*****************END MOTHUR VERSION 1.40.1 SCRIPT**********************

# ADDITIONAL COMMANDS IF NEEDED 
#(Uncomment if a mock community is added) get.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, groups=Mock)
#(Uncomment if a mock community is added) seq.error(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table, reference=HMP_MOCK.v35.fasta, aligned=F)
#(Uncoment if a mock community is added) remove.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, groups=Mock)
#Sub-sample data using a size suitable for the microbiome being studied.
#sub.sample(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared, size=5000)
