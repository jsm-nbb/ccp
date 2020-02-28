#!/bin/bash

if [ "$1" = "" ]; then
 echo "This tool searches the Pfam domain in the input FAST file and annotates species names."
 echo "$0 <MAG fasta file>"
 exit 1
fi

set -eux
set -o pipefail

i="$1"
echo "##Search Pfam domain"

sdir=$(dirname `readlink -f $0 || echo $0`)
export PATH=$sdir/hmmer-3.1b2/bin:"$PATH"
if [ "${PERL5LIB:-}" = "" ]; then
 export PERL5LIB=$sdir/PfamScan
else
 export PERL5LIB=$sdir/PfamScan:"$PERL5LIB"
fi

transeq -frame 6 $i -outseq $i.aa;
$sdir/PfamScan/pfam_scan.pl -fasta $i.aa -dir $sdir/PfamScan -cpu 8 > $i.aa.pfam;
echo -e "id\t$i" > $i.aa.pfam.txt;
cat $i.aa.pfam|grep -v "^#"|sed '/^$/d'|awk '{cnt[$6]++; name[$6]=$7} END{for(i in cnt){print i"\t"cnt[i]"\t"name[i]}}'|sort -k 2,2nr|cut -f 1,2 >> $i.aa.pfam.txt

$sdir/merge_table.pl -k $i.aa.pfam.txt $sdir/db.tsv |sed 's/\t\t/\t0\t/g; s/\t\t/\t0\t/g; s/\t$/\t0/g' > $i.input.tsv
perl -F'\t' -nlae 'BEGIN{@val=()}; print STDERR "read $. records"; for($i=0;$i<=$#F;$i++){$val[$i][$. - 1]=$F[$i]}; $n=$. - 1; $m=$#F;
                   END{print STDERR "Calculating correlation coefficient..."; @x=(); @y=();
                       for($i=1;$i<=$n;$i++){$x[$i]=$val[1][$i];};
                       for($j=2;$j<=$m;$j++){
                               for($i=1;$i<=$n;$i++){$y[$i]=$val[$j][$i]};
                               $sum_x=0; $sum_y=0; $sum_dx2=0; $sum_dy2=0; $sum_dxdy=0;
                               for($i=1;$i<=$n;$i++){$sum_x += $x[i]; $sum_y += $y[i]};
                               $m_x = $sum_x / $n; $m_y = $sum_y / $n;
                               for($i=1;$i<=$n;$i++){$sum_dx2 += ($x[$i] - $m_x) ** 2; $sum_dy2 += ($y[$i] - $m_y) ** 2; $sum_dxdy += ($x[$i] - $m_x) * ($y[$i] - $m_y)};
                               if($sum_dx2 == 0.0 || $sum_dy2 == 0.0){}else{$cor=$sum_dxdy / sqrt($sum_dx2) / sqrt($sum_dy2); print $val[$j][0]."\t".$cor}
                       }
                   }' $i.input.tsv |sort -k2,2gr |head -n 10| tee ccp.output.txt
