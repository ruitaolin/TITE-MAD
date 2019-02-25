# Time-to-Event Model-Assisted Design
R codes to implement time-to-event keyboard (TITE-Keyboard) and time-to-event mTPI (TITE-mTPI) designs in phase I dose-finding trials with late-onset toxicities.

# Description
Two useful strategies to speed up drug development are to increase the patient accrual rate and use novel adaptive designs. Unfortunately, these two strategies often conflict when the evaluation of the outcome cannot keep pace with the patient accrual rate and thus the interim data cannot be observed in time to make adaptive decisions. A similar logistic difficulty arises when the outcome is late onset. Based on a novel formulation and approximation of the likelihood of the observed data, the time-to-event model-assisted designs can handle toxicity data that are pending due to fast accrual or late-onset toxicity, and facilitate seamless decision making in phase I dose-finding trials. The time-to-event model-assisted designs consider each dose separately and the dose escalation/de-escalation rules can be tabulated before the trial begins, which greatly simplifies trial conduct in practice compared to that under existing methods.

# Functions
The repository includes two functions:
* get.boundary.tite.R: The R code that includes the function ```get.boundary.tite()``` to select the next dose level for the new patients by the NOC design.
```rscript
get.boundary.tite(target, ncohort, cohortsize, n.earlystop, prior.p, p.saf, p.tox, cutoff.eli, extrasafe, offset, print, output, design)
```
* get.oc.tite.R: The R code that includes the function ```get.oc.tite()``` to obtain the operating characteristics of the time-to-event model-assisted design for single agent trials with late-onset toxicities by simulating trials.
```rscipt
get.oc.tite(target, p.true, ncohort, cohortsize, maxt, prior.p, accrual, maxpen, dist1, dist2, alpha, n.earlystop, startdose, p.saf, p.tox, cutoff.eli, extrasafe, offset, ntrial, seed, design)
```


# Inputs
* ```target```: The target toxicity probability, e.g., ```target<-0.33```.
* ```dlt```: A vector of length *n* that stores the toxicity outcome for each patient, where *n* is the total number of patients so far.
* ```dose.level```: A vector of length *n* that stores the dose level assigned to each patient.
* ```ndose```: Number of prespecified dose levels of the new drug.
* ```epi```: A small positive value that defines the neighbourhood of the target toxicity probability.
* ```a```: The feasibility bound for overdose control, as default, ```a<-0.35```. 
* ```eta```: The dose-switching cutoff, as default, ```eta<-0.60```.
* ```lambda```: The dose-elimination cutoff, as default, ```lambda<-0.85```.
* ```enter.time```: A vector of length *n* that stores the day of arrival of each patient under late-onset cases.
* ```dlt.time```: A vector of length *n* that stores the time-to-toxicity outcome of each patient; If the subject has not experienced the DLT by the decision-making time, his time-to-toxcity outcome is *0*.
* ```current.time```: The arrival time of the new patient, or the decision-making time to decide the next dose level.
* ```tau```: The length of follow-up period.


#Example
We apply the NOC and fNOC designs to the sonidegib trial.
* Based on the accumulated data, two DLTs were observed at dose level 3 when patient 13 arrived on day 130. At this moment, patients 6, 8, 9, 11, and 12 were still under the follow-up of evaluation without experiencing any DLT, which led to a total of five missing toxicity outcomes. We utilize the following code to decide the dose level for patient 13.
```rscript
target <- 0.33
ndose <- 5
enter.time <- c(4,7,19,29,31,50,58,67,78,91,100,118)
dlt.time <- c(0,0,0,0,0,0,65,0,0,29,0,0)
dose.level <- c(1,1,1,2,2,2,3,3,3,3,3,3)
tau <- 90
current.time <- 130
get.next.fnoc(target, enter.time, dlt.time, current.time,tau, dose.level, ndose)
```
The output is given by 
```rscript
The feasibility bound for overdose control rule is 0.35 
The dose-switching cutoff is 0.6 
The dose-elimination cutoff is 0.85 
The posterior model probabilities are 0.02 0.16 0.55 0.2 0.07 
The posterior probability that the current dose level is overly toxic is 0.4749329 
The next dose level is 2 
      Patient No. Dose level Day of arrival Observed DLT Time to DLT Fractional DLT
 [1,]           1          1              4            0          91      0.0000000
 [2,]           2          1              7            0          91      0.0000000
 [3,]           3          1             19            0          91      0.0000000
 [4,]           4          2             29            0          91      0.0000000
 [5,]           5          2             31            0          91      0.0000000
 [6,]           6          2             50            0          91      0.0000000
 [7,]           7          3             58            1          65      1.0000000
 [8,]           8          3             67            0          91      0.1428571
 [9,]           9          3             78            0          91      0.1428571
[10,]          10          3             91            1          29      1.0000000
[11,]          11          3            100            0          91      0.1428571
[12,]          12          3            118            0          91      0.2207792
```
* On the other hand, if we had treated these missing data as no DLTs, we can utilize the ```get.next.noc``` function to generate the next dose level, as given by
```rscript 
target <- 0.33
ndose <- 5
dlt <- c(0,0,0,0,0,0,1,0,0,1,0,0)
dose.level <- c(1,1,1,2,2,2,3,3,3,3,3,3)
get.next.noc(target, dlt, dose.level, ndose)

-----------------------output------------------------
The feasibility bound for overdose control rule is 0.35 
The dose-switching cutoff is 0.6 
The dose-elimination cutoff is 0.85 
The posterior model probabilities are 0.01 0.08 0.49 0.28 0.14 
The posterior probability that the current dose level is overly toxic is 0.3341323 
The next dose level is 3 
      Patient No. Dose level DLT
 [1,]           1          1   0
 [2,]           2          1   0
 [3,]           3          1   0
 [4,]           4          2   0
 [5,]           5          2   0
 [6,]           6          2   0
 [7,]           7          3   1
 [8,]           8          3   0
 [9,]           9          3   0
[10,]          10          3   1
[11,]          11          3   0
[12,]          12          3   0
```
* In the end, three dose levels had been explored with four DLTs out of 18 patients at dose level 2 and four DLTs out of 9 patients at dose level 3 being observed. The MTD can be estimated by the ```select.mtd.noc``` function. 
```rscript 
target <- 0.33
ndose <- 5
dlt <- c(0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,1,1,1,1,0,1,0,0,0,1,0)
dose.level <- c(1,1,1,2,2,2,3,3,3,3,3,3,2,2,2,2,2,2,3,3,3,2,2,2,2,2,2,2,2,2)
select.mtd.noc(target, dlt, dose.level, ndose)

-----------------------output------------------------
The posterior model probabilities are 0.03 0.55 0.36 0.05 0.01 
The MTD is the dose level 2 
```
# Authors and Reference
* Ruitao Lin and Ying Yuan
* Lin, R. and Yuan, Y. (2019) “Time-to-event model-assisted designs for dose-finding trials with delayed toxicity”, Biostatistics, in press.

