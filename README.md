# Minimax approach for gene-set enrichment analysis of rare copy number variants

This is the code base to perform gene-set enrichment analyses using the MORST approach to rare copy number variant analysis. 

To set up the environment you can use the `setup_conda_env.sh`, which can be run in any linux system and Window's WSL. 

## Data simulation

Simulated data for gene-set enrichment analyses was run using the `GItest_demo.r` from <https://www4.stat.ncsu.edu/~jytzeng/Software/CCRET/software_ccret.php>.

The CCRET result for this dataset is `CCRET p-value of GI effect adjusting for length and dosage = 0 6.109567e-09`.

The data is in the `data/` folder:
* `mycnv2_ds.txt` contains the gene dosage matrix (with entries 0,1,2,3,4).
* `mycnv2_gi.txt` contains the gene overlaps matrix (0,1,2)
* `mycnv2_ln.txt` contains the event length matrix. 
* `mycnv2_yy.txt` contains the phenotypes (binary 0,1). 


Using MORST on the same data set with the code in `MORST_CNV.R` we get a `p = 2.634e-09, at a 0.05 alpha level`. 

## TODO:
 * Perform simulations using different amounts of signal sparsity
 * How does this vary with LD/genotype correlation
 * Apply region based analyses