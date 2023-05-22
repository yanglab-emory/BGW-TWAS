#!/usr/bin/bash

################################################################
# Step 2: Prune blocks
# Select a subset of genome blocks (up to ${max_blocks}) for joint model training by BGW-TWAS
# Cis blocks are always selected
# Rank trans blocks with minimum single-variant eQTL p-value < ${p_thresh} by the smallest p-value within block (from the smallest to the largest), and select top ranked trans blocks up to ${max_blocks}.
################################################################

# Variable needed for pruning genome blocks
###
# --wkdir : Specify a working directory
# --gene_name : Specify the gene name/id that should be the same used in `GeneExpFile`
# --GeneExpFilehead : Directory of the file containing a list of fileheads of segmented genotype files
# --Score_dir : Specify summary score statistic file directory
# --p_thresh : Specify p-value threshold
# --max_blocks : Specify maximum genome block number

#################################
VARS=`getopt -o "" -a -l \
wkdir:,gene_name:,GeneExpFile:,Genome_Seg_Filehead:,Score_dir:,p_thresh:,max_blocks:,clean_output: \
-- "$@"`

if [ $? != 0 ]
then
    echo "Please provide required input arguments. Terminating....." >&2
    exit 1
fi

eval set -- "$VARS"

while true
do
    case "$1" in
        --wkdir|-wkdir) wkdir=$2; shift 2;;
        --gene_name|-gene_name) gene_name=$2; shift 2;;
        --GeneExpFile|-GeneExpFile) GeneExpFile=$2; shift 2;;
        --Genome_Seg_Filehead|-Genome_Seg_Filehead) Genome_Seg_Filehead=$2; shift 2;;
        --Score_dir|-Score_dir) Score_dir=$2; shift 2;;
        --p_thresh|-p_thresh) p_thresh=$2; shift 2;;
        --max_blocks|-max_blocks) max_blocks=$2; shift 2;;
        --clean_output|-clean_output) clean_output=$2; shift 2;;
        --) shift;break;;
        *) echo "Wrong input arguments!";exit 1;;
        esac
done


echo Summary Score statistic directory:  $Score_dir

##########################################
# Setting Default Input Argument Values
##########################################
p_thresh=${p_thresh:-0.00001}
max_blocks=${max_blocks:-100}
clean_output=${clean_output:-0}

echo ${gene_name} with up to ${max_blocks} genome blocks and p-value threshold ${p_thresh}

## directory with eQTL summary statistics
# Score_dir=${wkdir}/${gene_name}_scores
cd ${wkdir}

if [ -z ${Score_dir} ] ; then
    Score_dir=${wkdir}/${gene_name}_scores
    echo Default summary score statistics file directory: $Score_dir
else
    echo Summary Score statistics file directory is provided as ${Score_dir} ;
fi

if [ -s ${Genome_Seg_Filehead} ]; then
> ${wkdir}/${gene_name}_ranked_segments.txt
cat ${Genome_Seg_Filehead} | while read filehead ; do
if [ -s  ${Score_dir}/${filehead}.score.txt.gz ] ; then
    nsnp=$(zcat ${Score_dir}/${filehead}.score.txt | wc -l)
    if [ "$nsnp" -gt 1 ] ; then
        zcat ${Score_dir}/${filehead}.score.txt | awk -v var=$filehead 'NR == 2 {line = $0; min = $13}; NR >2 && $13 < min {line = $0; min = $13}; END{print var, min}' >> ${wkdir}/${gene_name}_ranked_segments.txt
    else
        echo ${Score_dir}/${filehead}.score.txt.gz do not have any SNPs...
    fi
else
    echo A non-empty ${Score_dir}/${filehead}.score.txt.gz file dose not exist !
fi
done
else
    echo ${Genome_Seg_Filehead} is empty. Please check.
    exit
fi

# Grep gene info from ${GeneExpFile}
gene_info=$(grep -w ${gene_name} ${GeneExpFile} | head -n1)
target_chr=$( echo ${gene_info} | awk 'FS {print $1}');
start_pos=$(echo ${gene_info} | awk 'FS {print $2}');
start_pos=$((start_pos - 1000000))
end_pos=$(echo ${gene_info} | awk 'FS {print $3}');
end_pos=$((end_pos + 1000000))
echo ${gene_name} from CHR $target_chr ranging from $start_pos to $end_pos

## for gene_ranked_segs, FS=_.,
> ${wkdir}/${gene_name}_cis_segments.txt
> ${wkdir}/${gene_name}_trans_segments.txt

cat ${wkdir}/${gene_name}_ranked_segments.txt | while read line ; do

    filehead=$(echo $line | awk -F " " '{print $1}' )
    pval=$(echo $line | awk -F " " '{if(NF>1) print $2}' )
    comp_pval=$(awk -v pval=$pval -v p_thresh=$p_thresh 'BEGIN{ print (pval+0)<(p_thresh+0) }')
    # echo Compare if $pval less than $p_thresh gives comp_pval = $comp_pval ...

    row_1=$(zcat ${Score_dir}/${filehead}.score.txt | grep -v "#CHROM" | head -n 1)
    chr=$(echo ${row_1} | awk '{print $1}' )
    start=$(echo ${row_1} | awk '{print $2}')

    if [ "$chr" -eq "${target_chr}" ] ; then

        end=$( zcat ${Score_dir}/${filehead}.score.txt | tail -n1 | awk '{print $2}' )

        if [ "$end" -gt "$start_pos" ] && [ "$start" -lt "$start_pos" ] ; then
            echo -e "${filehead}\t${pval}" >> ${gene_name}_cis_segments.txt
        elif [ "$start" -lt "$end_pos" ]  && [ "$end" -gt "$end_pos" ] ; then
            echo -e  "${filehead}\t${pval}" >> ${wkdir}/${gene_name}_cis_segments.txt
        elif [ "$comp_pval" -eq 1 ] ; then
            echo -e  "${filehead}\t${pval}" >> ${wkdir}/${gene_name}_trans_segments.txt
        fi

        elif [ "$comp_pval" -eq 1 ] ; then
            echo -e  "${filehead}\t${pval}" >> ${wkdir}/${gene_name}_trans_segments.txt
        else
        continue;
    fi

done

n_cis=$(wc -l ${wkdir}/${gene_name}_cis_segments.txt | awk '{print $1}')
n_trans=$(wc -l ${wkdir}/${gene_name}_trans_segments.txt | awk '{print $1}')
max_trans_blocks=$((max_blocks - n_cis))

cat ${wkdir}/${gene_name}_cis_segments.txt > ${wkdir}/${gene_name}_select_segments.txt
if [ "$n_trans" -gt "$max_trans_blocks" ] ; then
        sort -g -k 2 -u ${wkdir}/${gene_name}_trans_segments.txt | head -${max_trans_blocks} >> ${wkdir}/${gene_name}_select_segments.txt
    else
        sort -g -k 2 -u ${wkdir}/${gene_name}_trans_segments.txt >> ${wkdir}/${gene_name}_select_segments.txt
fi

### Fileheads of segmented genome blocks
if [ -s ${wkdir}/${gene_name}_select_segments.txt ] ; then
    cut -f1 ${wkdir}/${gene_name}_select_segments.txt > ${wkdir}/${gene_name}_select_filehead.txt
    filehead=${wkdir}/${gene_name}_select_filehead.txt
else
    echo "Selected genome blocks for BGW-TWAS is empty. Please check gene annotation."
    exit
fi

if [ $clean_output -eq 1  ]; then
   rm -f ${wkdir}/${gene_name}_cis_segments.txt  ${wkdir}/${gene_name}_trans_segments.txt ${wkdir}/${gene_name}_ranked_segments.txt
fi

echo Complete Step2 for pruning genome blocks!

exit
