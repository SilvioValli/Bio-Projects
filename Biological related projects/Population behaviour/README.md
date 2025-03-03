<h1>Population's Allele Behaviour</h1>

<h2>Description</h2>
This mini-project was inspired by the end-of-chapter exercise on object-oriented programming in Dr. Martin Jones' book “Advanced Python for Biologists”.

</br>This problem is designed in that way:
- **Population** of 100 individuals
- Each individual has **3 Locus** 
- Each locus is defined by one between **2 possible alleles** , ex: 'A' or 'a'
- Each allele has its own **fitness**

Example:
- Population of 50 individuals = [ind0,ind1, ind2, ind3, ind4, ind5, ..., ind49]
- ind1 = 'locus1, locus2, locus3' = 'AbC'
- Locus1 = 'A' (in that case) but it could be 'A' or 'a'
- fitness('A') = 1
- fitness('a') = 0.78

Objectives:
- Single Allele, Locus and Individuals need to have their own custom classes, so 3 in total -> class's names("Allele", "Locus", "Individual")
- Simulate a RANDOM SELECTION event -> at each generation the individuals are removed if their fitness is lower than a random number between 0 and 1
- Simulate a REPRODUCTION event -> at each generation the remaining individuals from the random selection mate to give life to new offspring. Each new individual will have a random genotype whose alleles will be randomly chosen from those present in the population.
- Combine the RANDOM SELECTION event with the REPRODUCTION -> in order to mantain a population of 100 individuals but with a different proportion of the alleles with the respect of the previous generation population
- Write a "csv" file all the frequencies for each locus in the population and VISUALIZE the changes with a graphical representation
- Create a Pipeline 
<br />



<h2>Languages Used</h2>

- <b>Python</b> (Pandas, Matplotlib, Seaborn, Random and Built-in Python functions/classes)

<h2>Environments Used </h2>

- <b>Jupyter Notebook</b>

<h2>Final results</h2>
<img src="https://i.imgur.com/QFpFLtb.png" height="50%" width="50%" alt=""/>
