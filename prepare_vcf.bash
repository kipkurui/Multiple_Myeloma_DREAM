 # inside each of the folders, we run these commands
 
 parallel --gnu gunzip  ::: *gz
 parallel --gnu bgzip  ::: *vcf
 parallel --gnu tabix -p vcf ::: *gz
 
 #After this we can now be able to run use the Pysam command
