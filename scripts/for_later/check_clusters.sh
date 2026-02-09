echo RepeatExplorer
while read -r fasta; do 
  name=$(echo $fasta | sed 's/.fasta$//' | xargs basename)
  echo $name
  echo raw: $(grep ">" -c $fasta)
  echo mmseqs_out: $(grep -ic $name /home/agorbunov/bat_bchrom/analysis/annotation/repeatmasker/custom_database/cluster/mmseqs_out_rep_seq.fasta)
done < <(find  /home/agorbunov/bat_bchrom/analysis/annotation/repeatexplorer/fastas -type f)
  
echo TideCluster
echo raw: $(grep ">" -c /home/agorbunov/bat_bchrom/analysis/annotation/fastas/tidecluster_tarean.fasta)
  echo mmseqs_out: $(grep -ic tidecluster /home/agorbunov/bat_bchrom/analysis/annotation/repeatmasker/custom_database/cluster/mmseqs_out_rep_seq.fasta)

echo RepeatModeler
echo raw: $(grep ">" -c /home/agorbunov/bat_bchrom/analysis/annotation/fastas/repeatmodeler.fasta)
  echo mmseqs_out: $(grep -ic rnd /home/agorbunov/bat_bchrom/analysis/annotation/repeatmasker/custom_database/cluster/mmseqs_out_rep_seq.fasta)
