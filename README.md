<<<<<<< HEAD

# SCJ Vote scaffolder: A multiple reference-based scaffolder

  ## Algorithm
 This program is based on [this paper](https://www.researchgate.net/publication/49853612_SCJ_A_Breakpoint-Like_Distance_that_Simplifies_Several_Rearrangement_Problems). For any adjacency, if more than half reference genomes agree with it, we just try to use it to scaffold our target genome. 

## Input data
One target genome: 
* Please name your target genome as **"target.all"**.
* Every line in "target.all" should look like **"<gene_order> <gene_id> <contig_name> <contig_orientation>"**.

k reference genomes:
* Please name your reference genomes from **"ref1.all"** to **"ref\<k>.all"**.
* Your "ref.all" should has the same form as "target.all". 

We show an example input file in "Input".

## Output
The output of **SCJ Vote scaffolder.exe** will be placed in the folder you specified.

The scaffolding result will be named as **ScaffoldResult**. Our scaffolder also outputs two simple dotplots, **DotplotBefore.png** and **DotplotAfterSCF.png**.
  

## Relevant Paper

[SCJ: A Breakpoint-Like Distance that Simplifies Several Rearrangement Problems](https://www.researchgate.net/publication/49853612_SCJ_A_Breakpoint-Like_Distance_that_Simplifies_Several_Rearrangement_Problems)
=======
# SCJ-Vote-scaffolder
A multiple reference-based scaffolder based on SCJ distance.
>>>>>>> origin/main
