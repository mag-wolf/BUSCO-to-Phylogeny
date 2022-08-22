
#!/bin/bash
#
# A pipeline created to automatize all my steps for phylogenetic analyses using the BUSCO toolkit.
# In the end you'll have a concatenated supermatrix tree and a consensus supertree.
#
#
# Does require specific conda environments!
#
# by Magnus Wolf 2021 (magnus.wolf@senckenberg.de)
#
# In the following part you need to provide paths and desired settings for one or all of the different scripts as well as some general
# tree constructing steps:
#
#################################general dependencies and parameters#####################################################################
source ~/anaconda3/etc/profile.d/conda.sh               #make sure to be able to switch between conda envs
WORKINGDIR=/path/to/your/working/dir/                   #path to the Working directory
ulimit -n 2048                                          #get more doable jobs-temp files, I suggest not to change that.
THREADS=N                                               #number of threads
TASKS=N                                                 #number of parallel task to run alignment or trimming processes
MAFFTenv=MAFFTenv                                       #name of a conda environment with the newest MAFFT installation
CLIPKITenv=CLIPKITenv                                   #name of a conda environment with the newest CLIPKIT installation
FASCONCAT=$WORKINGDIR"/FASconCAT-G_v1.04.pl"            #path to the FASconCAT-G perl script by Patrick Kueck from Koenig Museum Bonn, no change necessary when the script is in the working directory
IQTREEenv=IQTREEenv                                     #name of a conda environment with the newest IQTREE installation
ASTRAL=$WORKINGDIR"/Astral/astral.5.7.3.jar"            #path to the Astral java script for supermatrix trees, no change necessary when the directory is in the working directory
#
##################################BUSCO-to-Phylogeny dependencies and parameters#########################################################
ODB10=whatever_odb10                                    #name of the OrthoDB database that should be used for all species
SPECIES=human                                           #name of the species that augustus is trained for and should be used in BUSCO
export BUSCO_CONFIG_FILE=$WORKINGDIR"/config.ini"       #provide the path to a BUSCO config file after the #=# sign, not necessary when config is in Working Dir
EXTRACTSCOS=$WORKINGDIR"/extract_scos_BUSCO.py"         #path to my python script extract_scos_BUSCO.py, no change necessary when the script is in the working directory
BUSCOenv=BUSCOenv                                       #name of BUSCO5 env, should have the newest BUSCO5 installed.
#                                                       #use this command to install BUSCO5: conda install -n BUSCO5env -c bioconda -c conda-forge busco=5.1.2
#
#
# Now that you provided the dependencies, you can choose which script you wanne run, everything you'll call #TRUE# will run!
# You can find a detailed discription for each step in a makeshift square down below. Additionaly, this discription will pop up if you run the respective
# step.
#
##################################Options for BUSCO-to-PHYLOGENY#########################################################################
runbusco=TRUE                                           #this will take a loooooong time!
findsharedscos=TRUE
makealignments=TRUE
trimmgenealignments=TRUE
constructgenetrees=TRUE
filtergenetrees=TRUE
concatgenealignments=TRUE
trimmsupermatrix=TRUE
constructsupermtree=TRUE                                #this will take a loooooong time!
constructsuperttree=TRUE
#
#
#########################################################################################################################################
# Now you finished. Go into the Working Directory and hit #bash BUSCO-to-Phylogeny.sh#. Then everything you chosed will run one by one.
# I recommend using some sort of log pipeing since long pipelines like this one tend to produce a lot of hidden errors. Try e.g.
# #2>&1 | tee error.log# behind the bash command.
#
# After this line I wouldn't touch anything...
#########################################################################################################################################

if [[ "$runbusco" = TRUE ]]
        then
        echo "########################start runbusco#############################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#this is runbusco of BUSCO-to-PHYLOGENY, a script for fast        #"
        echo "#phylogenetic studies based on the BUSCO pipline. This part will  #"
        echo "#run BUSCO and keeps all single copy orthologs. For space saving  #"
        echo "#purposes every other results that come with a BUSCO run will be  #"
        echo "#removed. Requires a working conda env with BUSCO and a regarding #"
        echo "#config file available from path. Most functions in this config   #"
        echo "#file will be overwritten by my command but make sure that the    #"
        echo "#config file states that it uses the #short# training.            #"
        echo "#please provide also a lineage database from OrthoDB10 and an     #"
        echo "#augustus training species above.                                 #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        conda activate $BUSCOenv
        mkdir BUSCO
        cd BUSCO
        for file in ./../*.fasta;
        do
                bn=`basename $file .fasta`
                busco --cpu $THREADS --mode genome --in $file --out $bn"_BUSCO" --lineage $ODB10 --augustus_species $SPECIES
                cd $bn"_BUSCO"
                mv "run_"$ODB10 "./../run_"$bn
                cd ..
                rm -r $bn"_BUSCO"
                cd "run_"$bn
                rm -r augustus_output #apparently BUSCO5 does not produce these output directories anoymore so might be redundant
                rm -r blast_output
                rm -r hmmer_output
                cd ./../
        done
        rm busco_downloads
        cd ./../
        conda deactivate
        wait
        fi

if [[ "$findsharedscos" = TRUE ]]
        then
        echo "########################start findsharedscos#######################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#this is findsharedscos of BUSCO-to-PHYLOGENY, a script for       #"
        echo "#fast phylogenetic studies based on the BUSCO pipline. This part  #"
        echo "#will screen through all BUSCO output directories to find         #"
        echo "#orthologs shared between all species resulted from the step above#"
        echo "#Make sure to provide the path to my script extract_scos_BUSCO.py #"
        echo "#Outputs a directory in the BUSCO folder filled with concatenated #"
        echo "#fasta files. Per default, 25% missing data in SCOS is allowed    #"
        echo "#in the extract_scos_BUSCO.py script. If you want to change it,   #"
        echo "#change the line ...exludecount = counter // 4... to whatever you #"
        echo "#desire. 4 means 1/4 are allowed to be absent.                    #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        cd BUSCO
        mkdir SCOS
        python $EXTRACTSCOS $WORKINGDIR"/BUSCO/"
        mv *.fasta ./SCOS
        cd ./../
        wait
        fi

if [[ "$makealignments" = TRUE ]]
        then
        echo "########################start makealignments#######################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#this is makealignments of BUSCO-to-PHYLOGENY, a script for       #"
        echo "#fast phylogenetic studies based on the BUSCO pipline. This part  #"
        echo "#will make alignments for each shared single copy ortholog found  #"
        echo "#in the steps above and located in the BUSCO/SCOS/ directory.     #"
        echo "#Requires a working conda env with mafft > v7 installed. Since    #"
        echo "#gene sequences are short, mafft is using the most accurate       #"
        echo "#option with 1000 iterative refinements incorporating local       #"
        echo "#pairwise alignment information. Afterwards it will also remove   #"
        echo "#every character of the IUPAC code that represents more then one  #"
        echo "#amino acid and will remove every MSA with unequal length.        #"
        echo "#FASTCONCAT can't handle thouse...This was the easiest method to  #"
        echo "#implement and there might be better ones!                        #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        conda activate $MAFFTenv
        mkdir ALIGNMENTS
        cd ALIGNMENTS
        for file in ./../BUSCO/SCOS/*.fasta;
        do
                ln -s $file
        done
        task1(){
                 genename=$(echo $1 | cut -d "_" -f1)
                 mafft --maxiterate 1000 --localpair $1 > $genename"_alig.fasta";
        }
        N=$TASKS
        (
        for file in *complete.fasta;
        do
                ((i=i%N)); ((i++==0)) && wait
                task1 "$file" &
        done
        )
        wait
        for file in *_alig.fasta;
        do
                filename=$(echo $file | cut -d "_" -f1)
                while IFS= read -r line;
                do
                        if [[ $line == ">"* ]];
                        then
                                echo $line >> $filename"_alig_noIUPAC.fasta"
                        else
                                echo $line | sed "s/J/-/g" | sed "s/Z/-/g" | sed "s/B/-/g" | sed "s/X/-/g" >> $filename"_alig_noIUPAC.fasta"
                        fi
                done < "$file"
        done
        for file in *_alig_noIUPAC.fasta;
        do
                filename=$(echo $file | cut -d "_" -f1)
                seqlength=0
                adder=0
                while IFS= read -r line;
                do
                        if [[ $line == ">"* ]];
                        then
                                echo $seqlength >> $filename"_noIUPAC_temp.txt"
                                seqlength=0
                        else
                                adder=$(echo $line | wc -c)
                                seqlength=$(($seqlength + $adder))
                        fi
                done < "$file"
                echo $seqlength >> $filename"_noIUPAC_temp.txt"
                uniqlengths=$(cat $filename"_noIUPAC_temp.txt" | uniq | wc -l)
                if [ $uniqlengths -gt 2 ];
                then
                        #echo $file
                        #echo $uniqlengths
                        rm $file
                fi
        done
        rm *_noIUPAC_temp.txt
        rm *_alig.fasta
        rm *complete.fasta
        cd ./../
        conda deactivate
        wait
        fi

if [[ "$trimmgenealignments" = TRUE ]]
        then
        echo "########################start trimmgenealignments##################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#this is trimmgenealignments of BUSCO-to-PHYLOGENY, a script      #"
        echo "#for fast phylogenetic studies based on the BUSCO pipline. This   #"
        echo "#part will trimm all individual SCOS alignments in the ALIGNMENTS #"
        echo "#directory. Make sure to provide a working conda env with ClipKIT #"
        echo "#installed. The program will only keep parsimony-informative      #"
        echo "#sites and completely conserved sites. It also incorporates a     #"
        echo "#dynamic gap filtering. Will output timmed.fasta files.           #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        cd ALIGNMENTS
        conda activate $CLIPKITenv
        task2(){
                 genename=$(echo $1 | cut -d "_" -f1)
                 clipkit $1 -m kpic-smart-gap -o $genename"_alig_trimmed.fasta"
                 clipkit $1 -m kpi-smart-gap -o $genename"_parinfosite.fasta"
        }
        N=$TASKS
        (
        for file in *alig_noIUPAC.fasta;
        do
                ((i=i%N)); ((i++==0)) && wait
                task2 "$file" &
        done
        )
        wait
        cd ./../
        wait
        fi

if [[ "$constructgenetrees" = TRUE ]]
        then
        echo "########################start constructgenetrees###################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#this is constructgenetrees of BUSCO-to-PHYLOGENY a script        #"
        echo "#for fast phylogenetic studies based on the BUSCO pipline. This   #"
        echo "#part will construct individual gene trees for each of the        #"
        echo "#alignments in the ALIGNMENTS directory using IQtree. This        #"
        echo "#requires a working conda env with IQtree installed. Make sure    #"
        echo "#to provide the name of the env above. Will make 1000 bootstrap   #"
        echo "#replications and will use 5 threads since more will not increase #"
        echo "#the speed with so short sequences. You'll find constructed trees #"
        echo "#in the GENETREES directory.                                      #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        conda activate $IQTREEenv
        mkdir GENETREES
        cd GENETREES
        for file in ./../ALIGNMENTS/*.fasta;
        do
                ln -s $file
        done
        for file in *trimmed.fasta;
        do
                iqtree -s $file -bb 1000 -nt 5
        done
        cd ./../
        conda deactivate
        wait
        fi

if [[ "$filtergenetrees" = TRUE ]]
        then
        echo "########################start filtergenetrees######################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#this is filtergenetrees of BUSCO-to-PHYLOGENY a script           #"
        echo "#for fast phylogenetic studies based on the BUSCO pipline. This   #"
        echo "#part uses the maximum likelihood distance of the alignment that  #"
        echo "#was calculated by iqtree to filter out genetrees and alignments  #"
        echo "#that are either too conserved to produce a meaningful tree       #"
        echo "#or had a too high genetic distance that might speak for a        #"
        echo "#miss-alignment artefact. To do so, it will calculate the         #"
        echo "#distribution of maximum likelihood distances of all alignments   #"
        echo "#and will then infere the 0.05 and 0.95 quantiles. All alignments #"
        echo "#and genetrees within these quantiles will be removed. A report,  #"
        echo "#which files were remove is printed to the console.               #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        cd GENETREES
        numofspec=$(cat *_alig_trimmed.fasta | grep ">" | sort | uniq | wc -l)
        numofspec2=$(($numofspec + 1))
        for file in *.mldist;
        do
                sum=0
                column=1
                filename=$(echo $file | cut -d "." -f1)
                while [[ $column -le $numofspec2 ]];
                do
                        ((column = column + 1))
                        adder=$(cat $file | awk FNR!=1 | sed "s/  / /g" | sed "s/  / /g" | sed "s/  / /g" | sed "s/  / /g" | sed "s/  / /g" | cut -d " " -f$column | awk '{s+=$1} END {print s}')
                        sum=$(echo "scale=6; $sum + $adder" | bc)
                        adder=0
                done
                filelength1=$(cat $file | wc -l)
                filelength2=$(echo "scale=6; $filelength1 - 1" | bc)
                relativesum=$(echo "scale=6; $sum / $filelength2" | bc)
                sumint=$(echo $relativesum | cut -d "." -f1)
                echo $sumint $filename >> mldist_dist.txt
                relativesum=0
                sumint=0
                sum=0
        done
        for file in mldist_dist.txt;
        do
                numofgenetrees=$(cat $file | wc -l)
                lowquantfloat=$(bc<<<$numofgenetrees*0.05)
                lowquantint=$(echo $lowquantfloat | cut -d "." -f1)
                cat $file | sort -k1 -h | head -$lowquantint >> mldist_toremove.txt
                cat $file | sort -k1 -h | tail -$lowquantint >> mldist_toremove.txt
        done
        for file in mldist_toremove.txt;
        do
                while IFS= read -r line;
                do
                        toremove=$(echo $line | cut -d " " -f2)
                        rm $toremove*
                done < "$file"
        done
        numofremoved=$(cat mldist_toremove.txt | wc -l)
        echo "Removed "$numofremoved" files because their maximum likelihood distance was either in the low 5percent quantile or the high 5percent quantile."
        echo "Therefore, it might be that these files were outlier and either to conserved for a meaningful tree or had an alignment artefact that"
        echo "resulted in an exeptional high genetic distance. Here a list of the files removed:"
        echo "#int of total relative genetic distance #filename"
        cat mldist_toremove.txt
        rm mldist_dist.txt
        rm mldist_toremove.txt
        cd ./../
        wait
        fi

if [[ "$concatgenealignments" = TRUE ]]
        then
        echo "########################start concatgenealignments#################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#this is concatgenealignments of BUSCO-to-PHYLOGENY a script      #"
        echo "#for fast phylogenetic studies based on the BUSCO pipline. This   #"
        echo "#part will concatenate all gene alignments in the ALIGNMENTS      #"
        echo "#directory using the FASconCAT_v1.11.pl perl script written by    #"
        echo "#Patrick Kueck from the Museum Koenig (https://github             #"
        echo "#.com/PatrickKueck/FASconCAT) with default options. Make sure to  #"
        echo "#provide the path to the script above. Will output a so called    #"
        echo "#FcC_supermatrix.fas in a dedicated SUPERMATRIX directory.        #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir SUPERMATRIX
        cd SUPERMATRIX
        for file in ./../GENETREES/*.treefile;
        do
                bn=`basename $file _trimmed.fasta.treefile`
                ln -s "./../GENETREES/"$bn"_noIUPAC.fasta"
        done
        perl $FASCONCAT -s
        rm *.fasta
        cd ./../
        wait
        fi

if [[ "$trimmsupermatrix" = TRUE ]]
        then
        echo "########################start trimmsupermatrix#####################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#this is trimmsupermatrix of BUSCO-to-PHYLOGENY a script          #"
        echo "#for fast phylogenetic studies based on the BUSCO pipline. This   #"
        echo "#part will trimm the supermatrix made in the step above using the #"
        echo "#same clipkit command as used to trimm the genealignments. This is#"
        echo "#part is only necessary when no trimming was applied to the       #"
        echo "#genealignments. Otherwise it might be redundant. Use if desired. #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        cd SUPERMATRIX
        conda activate $CLIPKITenv
        clipkit FcC_supermatrix.fas -m kpic-smart-gap -o FcC_supermatrix_trimmed.fas
        conda deactivate
        cd ./../
        wait
        fi

if [[ "$constructsupermtree" = TRUE ]]
        then
        echo "########################start constructsupermtree##################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#this is constructgenetrees of BUSCO-to-PHYLOGENY a script        #"
        echo "#for fast phylogenetic studies based on the BUSCO pipline. This   #"
        echo "#part will construct a supermatrix tree using IQtree and the      #"
        echo "#supermatrix alignment constructed above. Make sure that the      #"
        echo "#alignment is trimmed and that you provided the name of the       #"
        echo "#working conda env above. Will make 1000 bootstrap replications.  #"
        echo "#Constructed tree is lated placed in the SUPERMATRIX directory.   #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        conda activate $IQTREEenv
        cd SUPERMATRIX
        iqtree -s FcC_supermatrix_trimmed.fas -bb 1000 -nt $THREADS
        cd ./../
        conda deactivate
        wait
        fi

if [[ "$constructsuperttree" = TRUE ]]
        then
        echo "########################start constructsuperttree##################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#this is constructsuperttree of BUSCO-to-PHYLOGENY a script       #"
        echo "#for fast phylogenetic studies based on the BUSCO pipline. This   #"
        echo "#part will construct a consensus supertree tree from the gene     #"
        echo "#trees constructed in the steps above. This will be done by using #"
        echo "#Astral-III, a java script that takes a concatenated newick file  #"
        echo "#and can be downloaded at https://github.com/smirarab/ASTRAL.     #"
        echo "#Make sure to provide the path to the java script above. If your  #"
        echo "#java version isn't up to date create a conda env with a new      #"
        echo "#java installation and then run the overall script.               #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir SUPERTREE
        cd SUPERTREE
        for file in ./../GENETREES/*.treefile;
        do
                ln -s $file
        done
        for file in *.treefile;
        do
                filename=$(echo $file | cut -d "_" -f1)
                mv $file $filename"_1000bt.nex"
        done
        cat *_1000bt.nex > GENETREES_concat.nex
        java -jar $ASTRAL -i GENETREES_concat.nex -o SUPERTREE.nex
        rm *1000bt.nex
        cd ./../
        wait
        fi

