#1/bin/bash
#Author: Joseph Sevigny
#Purpose: Generic blast for blobtools create json command.


blast_flavor=megablast
database_path="/export/data/bio/ncbi/blast/db/v5/nt"
outname=$(basename $1)

blastn \
-task $blast_flavor \
-query $1 \
-db $database_path \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-culling_limit 5 \
-num_threads 8 \
-evalue 1e-25 \
-out $outname.vs.nt.cul5.1e25.megablast.out
