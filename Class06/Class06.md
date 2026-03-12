# Class 6: R Functions
Jolie Adair (PID A17337497)

- [Background](#background)
- [Our first functions](#our-first-functions)
- [A second function](#a-second-function)
- [A Protein generating function](#a-protein-generating-function)

## Background

ALL FUNCTIONS IN R HAVE AT LEAST 3 THINGS:

- A **name** that we use to call the cfunction.
- One or more input **arguments**
- The **body** the lines of R code that do the work

## Our first functions

Let’s write a silly wee function called `add()` to add some numbers (the
input arguments)

``` r
add <- function(x, y) {
  x + y
}
```

Now we can use this function

``` r
add(x=100, y=1)
```

    [1] 101

``` r
add(x=c(100, 1, 100), y=1)
```

    [1] 101   2 101

> Q. What if I give a multiple element vector to `x` and `y`

``` r
add(c(100, 1), y=c(100,1))
```

    [1] 200   2

> Q. What if I give three inputs to the functions?

``` r
#add(x=c(100,1), y=1, z=1 )
```

> Q. What if I give only one input to the add function?

``` r
addnew <-  function(x, y=1) {
  x + y
}
```

``` r
addnew(x=100)
```

    [1] 101

``` r
addnew(c(100,1), 100)
```

    [1] 200 101

If we write out functions with input arguments having no default value
then the user will be required to set them when they use the function.
We can give our input argument “default” values by setting them equal to
some sensible value - e.g. y=1 in the `add()`

## A second function

Let’s try something more interesting: Make a sequence generating tool..

The `sample()` function can be a useful starting point here:

``` r
sample(1:10, size=4)
```

    [1] 1 3 4 8

> Q. Generate 9 random numbers taken from the input vector x=1:10?

``` r
sample(1:10, size=9)
```

    [1]  2  7  9  5  3  6  1 10  8

> Q. Generate 12 random numbers taken from the input vector x=1:10?

``` r
sample(1:10, size=12, replace = TRUE)
```

     [1]  7 10  2  4  1  2  7  4  7  1 10  4

> Q. Write code for the `sample()` functions that generates nucleotide
> sequences of length 6?

``` r
sample(x=c("A", "G", "T", "C") , size =6, replace = TRUE, ) 
```

    [1] "A" "G" "T" "G" "T" "A"

> Q. Write a first function `generate_dna()` that returns a **user
> specified length** DNA sequence:

``` r
generate_dna <- function(len=6) {
  sample(x=c("A", "G", "T", "C") , size =len, replace = TRUE, )
  }
```

``` r
generate_dna(len = 100)
```

      [1] "T" "A" "T" "T" "T" "T" "A" "T" "G" "A" "G" "G" "A" "A" "T" "T" "C" "A"
     [19] "A" "A" "T" "A" "A" "T" "A" "T" "C" "C" "G" "G" "C" "C" "T" "A" "C" "T"
     [37] "T" "G" "G" "A" "G" "C" "G" "G" "T" "T" "T" "C" "T" "A" "A" "C" "T" "C"
     [55] "G" "G" "T" "T" "A" "T" "G" "A" "G" "C" "T" "T" "T" "G" "G" "T" "T" "T"
     [73] "A" "C" "G" "C" "T" "C" "G" "A" "G" "T" "A" "C" "A" "G" "G" "C" "T" "T"
     [91] "G" "G" "A" "A" "G" "A" "C" "C" "C" "T"

> **Key Points**

Every function in R looks fundamentally the same in the terms of its
structure. Basically 3 things: name, input, and body

    name <- function(input) {
      body
    }

> Functions can have multiple inputs. These can be **required**
> arguments or **optional** arguments. With optional arguments having a
> set value.

> Q. Modify and improve our `generate_dna()` function to return it’s
> generated sequence in a more standar format like “AGTAGTA” rather than
> the vector “A”, “G”, “T”, “C”

``` r
generate_dna <- function(len=6, fasta =TRUE) {
  
  
 ans <- sample(x=c("A", "G", "T", "C") , 
         size =len, 
         replace = TRUE, )
 if(fasta) {
   cat("Single element vector output")
  ans <- paste(ans, collapse = "")
 } else {
   cat("Multi-element vector output")
 }
  return(ans)
 }

generate_dna(fasta=FALSE)
```

    Multi-element vector output

    [1] "T" "A" "C" "G" "A" "A"

``` r
generate_dna(fasta=TRUE)
```

    Single element vector output

    [1] "GCATTG"

The `paste()` function - it’s job is to join up or stick together(a.k.a.
paste) input strings together

``` r
paste("diego", "loves R", sep=" does not ")
```

    [1] "diego does not loves R"

Flow control means where the R brain goes in your code

``` r
good_mood <- FALSE

if(good_mood) {
  cat("Great!")
} else {
  cat("Bummer!")
}
```

    Bummer!

## A Protein generating function

> Q. Write a function that generates a user specified length protein
> sequence.

> Q. Use that function to generate random protein sequences between
> length 6 to 12.

> Q. Are any of your sequnces unique i.e. not found anywhere in nature?

``` r
 aa <- c("A", "R", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "N")
```

``` r
generate_protein <- function(len)  {
  # The amino-acids to sample from
  
 aa <- c("A", "R", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "N")
 #Draw n=len amino acids to make our sequence
 
 ans <- sample(aa, size = len, replace = T)
 ans <- paste(ans, collapse= "")
 return(ans)
}
```

``` r
myseq <- generate_protein(42)
myseq
```

    [1] "WIMKDQTRVKNFYYMIVCKTKEYKVHLLSCASAVHNGERPKM"

> Q.2 answer

``` r
generate_protein(6)
```

    [1] "TMDGFR"

``` r
generate_protein(7)
```

    [1] "KLECIDS"

> Q. 3 answer

``` r
for(i in 6:12) {
  # FASTA ID line ">id"
  cat(">", i, sep="", "\n")
  # Protein Sequence
  cat( generate_protein(i), "\n" )
}
```

    >6
    PPTNEH 
    >7
    GNQKVQF 
    >8
    TPTSTWIN 
    >9
    VPCKFCPPF 
    >10
    QHRFEPNAHC 
    >11
    KAVPLAVWDYQ 
    >12
    DQDTHKFKSMQL 
