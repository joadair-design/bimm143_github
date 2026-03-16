# Class 18: Pertussis and the CMI-PB project
Jolie Adair (PID: A17337497)

- [Background](#background)
- [CDC Tracking Data](#cdc-tracking-data)
  - [Exploring CMI-PB data](#exploring-cmi-pb-data)
  - [Side-Note: Working with dates](#side-note-working-with-dates)
  - [Joining multiple tables](#joining-multiple-tables)

## Background

Pertussis (whooping cough) is a common lung infection caused by the
bacteria B. Pertussis. It can infect anyone but is most deadly for
infants (under 1 year of age).

# CDC Tracking Data

The CDC data

> Q1. With the help of the R “addin” package datapasta assign the CDC
> pertussis case number data to a data frame called cdc and use ggplot
> to make a plot of cases numbers over time.

``` r
library(datapasta)

cdc <- data.frame(
                                 Year = c(1922L,1923L,1924L,1925L,
                                          1926L,1927L,1928L,1929L,1930L,1931L,
                                          1932L,1933L,1934L,1935L,1936L,
                                          1937L,1938L,1939L,1940L,1941L,1942L,
                                          1943L,1944L,1945L,1946L,1947L,
                                          1948L,1949L,1950L,1951L,1952L,
                                          1953L,1954L,1955L,1956L,1957L,1958L,
                                          1959L,1960L,1961L,1962L,1963L,
                                          1964L,1965L,1966L,1967L,1968L,1969L,
                                          1970L,1971L,1972L,1973L,1974L,
                                          1975L,1976L,1977L,1978L,1979L,1980L,
                                          1981L,1982L,1983L,1984L,1985L,
                                          1986L,1987L,1988L,1989L,1990L,
                                          1991L,1992L,1993L,1994L,1995L,1996L,
                                          1997L,1998L,1999L,2000L,2001L,
                                          2002L,2003L,2004L,2005L,2006L,2007L,
                                          2008L,2009L,2010L,2011L,2012L,
                                          2013L,2014L,2015L,2016L,2017L,2018L,
                                          2019L,2020L,2021L,2022L,2023L,2024L, 2025L),
         Cases = c(107473,164191,165418,152003,
                                          202210,181411,161799,197371,
                                          166914,172559,215343,179135,265269,
                                          180518,147237,214652,227319,103188,
                                          183866,222202,191383,191890,109873,
                                          133792,109860,156517,74715,69479,
                                          120718,68687,45030,37129,60886,
                                          62786,31732,28295,32148,40005,
                                          14809,11468,17749,17135,13005,6799,
                                          7717,9718,4810,3285,4249,3036,
                                          3287,1759,2402,1738,1010,2177,2063,
                                          1623,1730,1248,1895,2463,2276,
                                          3589,4195,2823,3450,4157,4570,
                                          2719,4083,6586,4617,5137,7796,6564,
                                          7405,7298,7867,7580,9771,11647,
                                          25827,25616,15632,10454,13278,
                                          16858,27550,18719,48277,28639,32971,
                                          20762,17972,18975,15609,18617,
                                          6124,2116,3044,7063,21996,22538 )
       )
```

> Q. I want a plot of year vs cases

``` r
library(ggplot2)

ggplot(cdc) +
  aes(x = Year, y = Cases) +
  geom_point() +
  geom_line() +
  labs(title = "CDC Pertussis Cases Over Time",
       x = "Year",
       y = "Number of Cases")
```

![](Class18_files/figure-commonmark/unnamed-chunk-2-1.png)

> Q2. Using the ggplot geom_vline() function add lines to your previous
> plot for the 1946 introduction of the wP vaccine and the 1996 switch
> to aP vaccine (see example in the hint below). What do you notice?

After the introduction of the wP vaccine in 1946, pertussis cases
dropped dramatically, but cases began to increase again after the switch
to the aP vaccine in 1996.

``` r
ggplot(cdc) +
  aes(x = Year, y = Cases) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 1946, col="orange", lty=2) +
  geom_vline(xintercept = 1996, col="red", lty =2) +
  geom_vline(xintercept = 2020, col="blue", lty =2) +
  labs(title = "CDC Pertussis Cases Over Time",
       x = "Year",
       y = "Number of Cases")
```

![](Class18_files/figure-commonmark/unnamed-chunk-3-1.png)

> Q3. Describe what happened after the introduction of the aP vaccine?
> Do you have a possible explanation for the observed trend?

After the aP vaccine was introduced, pertussis cases gradually
increased, possibly because immunity from the acellular vaccine weakens
faster than the original whole-cell vaccine.

## Exploring CMI-PB data

The CMI-PB project \< https://www.cmi-pb.org/ \>: The mission of CMI-PB
is to provide the scientific community with a comprehensive,
high-quality and freely accessible resource of Pertussis booster
vaccination.

Basically, make available a large dataset on the immune response to
Pertussis. They use a “booster” vaccination as a proxy for Pertussis
infection.

They make their data available as JSON format API. We can read this into
R with the `read_json()` function from th **jsonlite** package.

``` r
library(jsonlite)
subject <- read_json("https://www.cmi-pb.org/api/subject", simplifyVector = TRUE) 

head(subject)
```

      subject_id infancy_vac biological_sex              ethnicity  race
    1          1          wP         Female Not Hispanic or Latino White
    2          2          wP         Female Not Hispanic or Latino White
    3          3          wP         Female                Unknown White
    4          4          wP           Male Not Hispanic or Latino Asian
    5          5          wP           Male Not Hispanic or Latino Asian
    6          6          wP         Female Not Hispanic or Latino White
      year_of_birth date_of_boost      dataset
    1    1986-01-01    2016-09-12 2020_dataset
    2    1968-01-01    2019-01-28 2020_dataset
    3    1983-01-01    2016-10-10 2020_dataset
    4    1988-01-01    2016-08-29 2020_dataset
    5    1991-01-01    2016-08-29 2020_dataset
    6    1988-01-01    2016-10-10 2020_dataset

> Q4. How many aP and wP infancy vaccinated subjects are in the dataset?

87

``` r
table(subject$infancy_vac)
```


    aP wP 
    87 85 

> Q5. How many Male and Female subjects/patients are in the dataset?

way more females at 112 than males at 60

``` r
table(subject$biological_sex)
```


    Female   Male 
       112     60 

> Q6. What is the breakdown of race and biological sex (e.g. number of
> Asian females, White males etc…)?

The dataset is mostly composed of White and Asian participants, with 48
White females and 32 White males, and 32 Asian females and 12 Asian
males, while other racial groups are represented in much smaller
numbers.

> Q. is this repersentative of the US population?

No, it is not representative of the U.S. population because the dataset
is dominated by White and Asian participants while other racial groups
are underrepresented.

## Side-Note: Working with dates

``` r
library(lubridate)
```


    Attaching package: 'lubridate'

    The following objects are masked from 'package:base':

        date, intersect, setdiff, union

``` r
today()
```

    [1] "2026-03-15"

``` r
today() - ymd("2000-01-01")
```

    Time difference of 9570 days

``` r
time_length(today() - ymd("2000-01-01"), "years")
```

    [1] 26.20123

> Q7. Using this approach determine (i) the average age of wP
> individuals, (ii) the average age of aP individuals; and (iii) are
> they significantly different?

``` r
subject$age <- today() - ymd(subject$year_of_birth)
```

``` r
time_length(mean(subject$age[subject$infancy_vac == "wP"]), "years")
```

    [1] 36.84728

``` r
time_length(mean(subject$age[subject$infancy_vac == "aP"]), "years")
```

    [1] 28.09658

The average age of wP individuals is about 36.8 years while the average
age of aP individuals is about 28.1 years, showing that wP individuals
are older because the whole-cell vaccine was used earlier.

> Q8. Determine the age of all individuals at time of boost?

The ages at the time of booster vaccination vary across individuals but
generally range from about 18 to about 50 years old.

``` r
subject$age_at_boost <- ymd(subject$date_of_boost) - ymd(subject$year_of_birth)
```

``` r
time_length(subject$age_at_boost, "years")
```

      [1] 30.69678 51.07461 33.77413 28.65982 25.65914 28.77481 35.84942 34.14921
      [9] 20.56400 34.56263 30.65845 34.56263 19.56194 23.61944 27.61944 29.56331
     [17] 36.69815 19.65777 22.73511 35.65777 33.65914 31.65777 25.73580 24.70089
     [25] 28.70089 33.73580 19.73443 34.73511 19.73443 28.73648 27.73443 19.81109
     [33] 26.77344 33.81246 25.77413 19.81109 18.85010 19.81109 31.81109 22.81177
     [41] 31.84942 19.84942 18.85010 18.85010 19.90691 18.85010 20.90897 19.04449
     [49] 20.04381 19.90691 19.90691 19.00616 19.00616 20.04381 20.04381 20.07940
     [57] 21.08145 20.07940 20.07940 20.07940 32.26557 25.90007 23.90144 25.90007
     [65] 28.91992 42.92129 47.07461 47.07461 29.07324 21.07324 21.07324 28.15058
     [73] 24.15058 24.15058 21.14990 21.14990 31.20876 26.20671 32.20808 27.20876
     [81] 26.20671 21.20739 20.26557 22.26420 19.32375 21.32238 19.32375 19.32375
     [89] 22.41752 20.41889 21.41821 19.47707 23.47707 20.47639 21.47570 19.47707
     [97] 35.90965 28.73648 22.68309 20.83231 18.83368 18.83368 27.68241 32.68172
    [105] 27.68241 25.68378 23.68241 26.73785 32.73648 24.73648 25.79603 25.79603
    [113] 25.79603 31.79466 19.83299 21.91102 27.90965 24.06297 23.90965 27.12115
    [121] 22.12183 23.12115 26.17933 22.17933 29.17728 29.23477 26.23682 28.29295
    [129] 31.29363 26.29432 24.35044 27.35113 25.40999 32.41068 27.56194 27.41136
    [137] 24.50650 22.56263 29.56057 21.69473 26.69678 31.90691 19.90691 23.90691
    [145] 20.90623 31.00616 23.00616 35.00616 32.00548 32.00548 31.04449 28.12047
    [153] 25.11978 26.11910 26.19302 22.19302 26.19302 23.19507 29.19370 27.32923
    [161] 30.32717 24.55852 30.55715 32.55852 30.55715 22.67488 26.67488 32.67625
    [169] 20.67625 31.75086 20.86516 36.06297

> Q9. With the help of a faceted boxplot or histogram (see below), do
> you think these two groups are significantly different?

Yes, the two groups appear significantly different because the wP group
is generally older while the aP group is mostly younger.

``` r
ggplot(subject) +
  aes(time_length(age, "year"),
      fill = as.factor(infancy_vac)) +
  geom_histogram(show.legend = FALSE) +
  facet_wrap(vars(infancy_vac), nrow = 2) +
  xlab("Age in years")
```

    `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

![](Class18_files/figure-commonmark/unnamed-chunk-16-1.png)

## Joining multiple tables

``` r
table(subject$race, subject$biological_sex)
```

                                               
                                                Female Male
      American Indian/Alaska Native                  0    1
      Asian                                         32   12
      Black or African American                      2    3
      More Than One Race                            15    4
      Native Hawaiian or Other Pacific Islander      1    1
      Unknown or Not Reported                       14    7
      White                                         48   32

We can read more tables from the CMI-PB database

``` r
speciman <- read_json("http://cmi-pb.org/api/v5_1/specimen", simplifyVector = T)
ab_titer <- read_json("http://cmi-pb.org/api/v5_1/plasma_ab_titer", simplifyVector = T)
```

``` r
head(speciman)
```

      specimen_id subject_id actual_day_relative_to_boost
    1           1          1                           -3
    2           2          1                            1
    3           3          1                            3
    4           4          1                            7
    5           5          1                           11
    6           6          1                           32
      planned_day_relative_to_boost specimen_type visit
    1                             0         Blood     1
    2                             1         Blood     2
    3                             3         Blood     3
    4                             7         Blood     4
    5                            14         Blood     5
    6                            30         Blood     6

``` r
head(ab_titer)
```

      specimen_id isotype is_antigen_specific antigen        MFI MFI_normalised
    1           1     IgE               FALSE   Total 1110.21154       2.493425
    2           1     IgE               FALSE   Total 2708.91616       2.493425
    3           1     IgG                TRUE      PT   68.56614       3.736992
    4           1     IgG                TRUE     PRN  332.12718       2.602350
    5           1     IgG                TRUE     FHA 1887.12263      34.050956
    6           1     IgE                TRUE     ACT    0.10000       1.000000
       unit lower_limit_of_detection
    1 UG/ML                 2.096133
    2 IU/ML                29.170000
    3 IU/ML                 0.530000
    4 IU/ML                 6.205949
    5 IU/ML                 4.679535
    6 IU/ML                 2.816431

To make sense of all this data we need to “join” (a.k.a.”merge” or
“link”) all these tables together. only then will uou know that a give
Ab measurment was collected on a ceretin data form a certion data from a
certin wP or aP subject (from the `subject` table)

we can use **dplyer** and the `*_join()` family of functions to do this.

> Q9. Complete the code to join specimen and subject tables to make a
> new merged data frame containing all specimen records along with their
> associated subject details:

Yes, the two groups appear significantly different because the wP group
is generally older while the aP group is mostly younger.

``` r
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
meta <- inner_join(subject, speciman)
```

    Joining with `by = join_by(subject_id)`

``` r
head(meta)
```

      subject_id infancy_vac biological_sex              ethnicity  race
    1          1          wP         Female Not Hispanic or Latino White
    2          1          wP         Female Not Hispanic or Latino White
    3          1          wP         Female Not Hispanic or Latino White
    4          1          wP         Female Not Hispanic or Latino White
    5          1          wP         Female Not Hispanic or Latino White
    6          1          wP         Female Not Hispanic or Latino White
      year_of_birth date_of_boost      dataset        age age_at_boost specimen_id
    1    1986-01-01    2016-09-12 2020_dataset 14683 days   11212 days           1
    2    1986-01-01    2016-09-12 2020_dataset 14683 days   11212 days           2
    3    1986-01-01    2016-09-12 2020_dataset 14683 days   11212 days           3
    4    1986-01-01    2016-09-12 2020_dataset 14683 days   11212 days           4
    5    1986-01-01    2016-09-12 2020_dataset 14683 days   11212 days           5
    6    1986-01-01    2016-09-12 2020_dataset 14683 days   11212 days           6
      actual_day_relative_to_boost planned_day_relative_to_boost specimen_type
    1                           -3                             0         Blood
    2                            1                             1         Blood
    3                            3                             3         Blood
    4                            7                             7         Blood
    5                           11                            14         Blood
    6                           32                            30         Blood
      visit
    1     1
    2     2
    3     3
    4     4
    5     5
    6     6

let’s do one more `inner_join()` to join the `ab_titer` which all out
`meta` data

> Q10. Now using the same procedure join meta with titer data so we can
> further analyze this data in terms of time of visit aP/wP, male/female
> etc.

``` r
abdata <- inner_join(ab_titer, meta)
```

    Joining with `by = join_by(specimen_id)`

``` r
head(abdata)
```

      specimen_id isotype is_antigen_specific antigen        MFI MFI_normalised
    1           1     IgE               FALSE   Total 1110.21154       2.493425
    2           1     IgE               FALSE   Total 2708.91616       2.493425
    3           1     IgG                TRUE      PT   68.56614       3.736992
    4           1     IgG                TRUE     PRN  332.12718       2.602350
    5           1     IgG                TRUE     FHA 1887.12263      34.050956
    6           1     IgE                TRUE     ACT    0.10000       1.000000
       unit lower_limit_of_detection subject_id infancy_vac biological_sex
    1 UG/ML                 2.096133          1          wP         Female
    2 IU/ML                29.170000          1          wP         Female
    3 IU/ML                 0.530000          1          wP         Female
    4 IU/ML                 6.205949          1          wP         Female
    5 IU/ML                 4.679535          1          wP         Female
    6 IU/ML                 2.816431          1          wP         Female
                   ethnicity  race year_of_birth date_of_boost      dataset
    1 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    2 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    3 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    4 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    5 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    6 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
             age age_at_boost actual_day_relative_to_boost
    1 14683 days   11212 days                           -3
    2 14683 days   11212 days                           -3
    3 14683 days   11212 days                           -3
    4 14683 days   11212 days                           -3
    5 14683 days   11212 days                           -3
    6 14683 days   11212 days                           -3
      planned_day_relative_to_boost specimen_type visit
    1                             0         Blood     1
    2                             0         Blood     1
    3                             0         Blood     1
    4                             0         Blood     1
    5                             0         Blood     1
    6                             0         Blood     1

> Q11. How many specimens (i.e. entries in abdata) do we have for each
> isotype?

Six antibody isotypes are measured in the dataset: IgE (6698), IgG
(7265), IgG1 (11993), IgG2 (12000), IgG3 (12000), and IgG4 (12000).

``` r
table(abdata$isotype)
```


      IgE   IgG  IgG1  IgG2  IgG3  IgG4 
     6698  7265 11993 12000 12000 12000 

> Q hoe many diffrent \`antigen’ values are measured?

There are 16 different antigen values measured in the dataset.

``` r
table(abdata$antigen)
```


        ACT   BETV1      DT   FELD1     FHA  FIM2/3   LOLP1     LOS Measles     OVA 
       1970    1970    6318    1970    6712    6318    1970    1970    1970    6318 
        PD1     PRN      PT     PTM   Total      TT 
       1970    6712    6712    1970     788    6318 

> Q12. What are the different \$dataset values in abdata and what do you
> notice about the number of rows for the most “recent” dataset?

There are four dataset values: 2020_dataset, 2021_dataset, 2022_dataset,
and 2023_dataset. The most recent datasets contain fewer rows than the
earlier datasets, likely because fewer samples have been collected or
processed so far.

``` r
table(abdata$dataset)
```


    2020_dataset 2021_dataset 2022_dataset 2023_dataset 
           31520         8085         7301        15050 

Lets focus on IgG istype

``` r
igg <- abdata |>
  filter(isotype=="IgG")
head(igg)
```

      specimen_id isotype is_antigen_specific antigen        MFI MFI_normalised
    1           1     IgG                TRUE      PT   68.56614       3.736992
    2           1     IgG                TRUE     PRN  332.12718       2.602350
    3           1     IgG                TRUE     FHA 1887.12263      34.050956
    4          19     IgG                TRUE      PT   20.11607       1.096366
    5          19     IgG                TRUE     PRN  976.67419       7.652635
    6          19     IgG                TRUE     FHA   60.76626       1.096457
       unit lower_limit_of_detection subject_id infancy_vac biological_sex
    1 IU/ML                 0.530000          1          wP         Female
    2 IU/ML                 6.205949          1          wP         Female
    3 IU/ML                 4.679535          1          wP         Female
    4 IU/ML                 0.530000          3          wP         Female
    5 IU/ML                 6.205949          3          wP         Female
    6 IU/ML                 4.679535          3          wP         Female
                   ethnicity  race year_of_birth date_of_boost      dataset
    1 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    2 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    3 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    4                Unknown White    1983-01-01    2016-10-10 2020_dataset
    5                Unknown White    1983-01-01    2016-10-10 2020_dataset
    6                Unknown White    1983-01-01    2016-10-10 2020_dataset
             age age_at_boost actual_day_relative_to_boost
    1 14683 days   11212 days                           -3
    2 14683 days   11212 days                           -3
    3 14683 days   11212 days                           -3
    4 15779 days   12336 days                           -3
    5 15779 days   12336 days                           -3
    6 15779 days   12336 days                           -3
      planned_day_relative_to_boost specimen_type visit
    1                             0         Blood     1
    2                             0         Blood     1
    3                             0         Blood     1
    4                             0         Blood     1
    5                             0         Blood     1
    6                             0         Blood     1

We make a plot of `MFI_normalised` values fpr all `antigen` vlaues

> Q13. Complete the following code to make a summary boxplot of Ab titer
> levels (MFI) for all antigens:

The antigens PT, FIM2/3, and FHA appear to have the widest range of
antibody titer values.

``` r
ggplot(igg) +
  aes(MFI_normalised, antigen) +
  geom_boxplot()
```

![](Class18_files/figure-commonmark/unnamed-chunk-27-1.png)

The antigens “PT”, “FIM2/3” and “FHA” apeaer to have the widest range of
values

> Q14. What antigens show differences in the level of IgG antibody
> titers recognizing them over time? Why these and not others?

The antigens PT, FHA, FIM2/3, and PRN show the largest differences in
IgG antibody titers. These antigens are components of the pertussis
vaccine, so they trigger stronger immune responses, while unrelated
antigens such as OVA or measles show much smaller changes.

> Q15. Filter to pull out only two specific antigens for analysis and
> create a boxplot for each. You can chose any you like. Below I picked
> a “control” antigen (“OVA”, that is not in our vaccines) and a clear
> antigen of interest (“PT”, Pertussis Toxin, one of the key virulence
> factors produced by the bacterium B. pertussis).

Two antigens were selected for analysis: PT, a pertussis toxin antigen
found in the vaccine, and OVA, which serves as a control antigen not
related to pertussis.

``` r
pt_ova <- igg |>
  filter(antigen %in% c("PT", "OVA"))
```

``` r
ggplot(pt_ova) +
  aes(MFI_normalised, antigen) +
  geom_boxplot()
```

![](Class18_files/figure-commonmark/unnamed-chunk-29-1.png)

For this comparison I selected the control antigen OVA and the pertussis
antigen PT. OVA is not part of the pertussis vaccine and serves as a
baseline control, while PT (Pertussis toxin) is a major antigen in the
vaccine and should show a stronger immune response.

> Q16. What do you notice about these two antigens time courses and the
> PT data in particular?

The PT antigen shows a clear increase in IgG antibody levels after the
booster vaccination, while the control antigen OVA shows little change
over time because it is not related to the pertussis vaccine.

> Q17. Do you see any clear difference in aP vs. wP responses?

Yes, there appears to be some differences between aP and wP responses,
with wP individuals generally showing slightly stronger antibody
responses to several pertussis antigens compared to aP individuals.

``` r
ggplot(igg) +
  aes(MFI_normalised, antigen, col=infancy_vac) +
  geom_boxplot()
```

![](Class18_files/figure-commonmark/unnamed-chunk-30-1.png)

``` r
ggplot(igg) +
  aes(MFI_normalised, antigen,) +
  geom_boxplot() +
  facet_wrap(~infancy_vac)
```

![](Class18_files/figure-commonmark/unnamed-chunk-31-1.png)

Is there a difference with time (I.e. before booster shot vs. after
booster shot)

``` r
ggplot(igg) +
  aes(MFI_normalised, antigen, col=infancy_vac) +
  geom_boxplot() +
  facet_wrap(~visit)
```

![](Class18_files/figure-commonmark/unnamed-chunk-32-1.png)

``` r
abdata.21 <- abdata %>% filter(dataset == "2021_dataset")

abdata.21 %>%
  filter(isotype == "IgG",  antigen == "PT") %>%
  ggplot() +
    aes(x=planned_day_relative_to_boost,
        y=MFI_normalised,
        col=infancy_vac,
        group=subject_id) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_vline(xintercept=14, linetype="dashed") +
  labs(title="2021 dataset IgG PT",
       subtitle = "Dashed lines indicate day 0 (pre-boost) and 14 (apparent peak levels)")
```

![](Class18_files/figure-commonmark/unnamed-chunk-33-1.png)

> Q18. Does this trend look similar for the 2020 dataset?

Yes, the trend looks similar in the 2020 dataset. IgG antibody levels
against PT increase after the booster vaccination, peak around day 14,
and then gradually decline, which is consistent with a typical immune
response.
