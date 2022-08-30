# ==================
# = call genotypes =
# ==================

# 1. jgil to call genotypes
#    performed on zenyatta
#    but copied the vcf to hyperion
# ============================================================

~/software/jgil-1.6/jgil -q 10 -Q 10 -m 2 -n 20 -p 20 -f ~/storage/dgrpFreeze2/dmel.fa -t 30 -b freeze2SamplePlusRelated -o Freeze2_JGIL_ > jgil.log 2>&1 &

# 2. parse genotypes, need to get all genotypes for use later
#    can extract genotypes for the 50 lines
# ============================================================

perl jgil/parseVCFall.pl | awk -F " " '($5 + $6 >= 40 && $6/($5+$6) >= 0.01) || $1 == "chr" ' > jgil/all.tgeno &

# 3. prepare genotyeps for plink
# ============================================================

tail -n+2 jgil/all.tgeno | perl -wne 'chomp $_; @line = split / /, $_;
  if ($line[0] eq "2L") { print "1 "; } elsif ($line[0] eq "2R") { print "2 "; } elsif ($line[0] eq "3L") { print "3 "; }
  elsif ($line[0] eq "3R") { print "4 "; } elsif ($line[0] eq "X") { print "5 "; } elsif ($line[0] eq "4") { print "6 "; } else { next; }
  print $line[0], "_", $line[1], " 0 ", $line[1];
  for (my $i = 6; $i <= $#line; $i++) { print " ", $line[$i], " ", $line[$i]; }
  print "\n";' | sed 's/-/0/g' > jgil/all.tped &
head -n 1 jgil/all.tgeno | cut -d " " -f 7- | sed 's/ /\n/g' | awk '{print $1" "$1" 0 0 2 NA"}' > jgil/all.tfam
~/software/plink-v1.90b3w/plink --silent --tped jgil/all.tped --tfam jgil/all.tfam --make-bed --out jgil/all &

awk '{ print $1" "$1 }' jgil/diallelLineIDs > jgil/diallel.line.txt
~/software/plink-v1.90b3w/plink --silent --bfile jgil/all --keep jgil/diallel.line.txt --maf 0.1 --geno 0.2 --make-bed --out jgil/diallel &

# 4. calculate relationship matrix
# ============================================================

~/software/plink-v1.90b3w/plink --silent --bfile jgil/diallel --recode12 --transpose --out gwas/diallel &
