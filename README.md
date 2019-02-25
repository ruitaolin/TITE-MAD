# Time-to-Event Model-Assisted Design
R codes to implement time-to-event keyboard (TITE-Keyboard) and time-to-event mTPI (TITE-mTPI) designs in phase I dose-finding trials with late-onset toxicities.

# Description
Two useful strategies to speed up drug development are to increase the patient accrual rate and use novel adaptive designs. Unfortunately, these two strategies often conflict when the evaluation of the outcome cannot keep pace with the patient accrual rate and thus the interim data cannot be observed in time to make adaptive decisions. A similar logistic difficulty arises when the outcome is late onset. Based on a novel formulation and approximation of the likelihood of the observed data, the time-to-event model-assisted designs can handle toxicity data that are pending due to fast accrual or late-onset toxicity, and facilitate seamless decision making in phase I dose-finding trials. The time-to-event model-assisted designs consider each dose separately and the dose escalation/de-escalation rules can be tabulated before the trial begins, which greatly simplifies trial conduct in practice compared to that under existing methods.

# Functions
The repository includes two functions:
* get.boundary.tite.R: The R code that includes the function ```get.boundary.tite()``` to select the next dose level for the new patients by the NOC design.
```rscript
get.boundary.tite(target, ncohort, cohortsize, n.earlystop, prior.p, p.saf, p.tox, cutoff.eli, 
extrasafe, offset, print, output, design)
```
* get.oc.tite.R: The R code that includes the function ```get.oc.tite()``` to obtain the operating characteristics of the time-to-event model-assisted design for single agent trials with late-onset toxicities by simulating trials.
```rscipt
get.oc.tite(target, p.true, ncohort, cohortsize, maxt, prior.p, accrual, maxpen, dist1, dist2, alpha, 
n.earlystop, startdose, p.saf, p.tox, cutoff.eli, extrasafe, offset, ntrial, seed, design)
```


# Inputs
* ```target```: The target toxicity rate, e.g., ```target<-0.33```.
* ```p.true```: A vector containing the true toxicity probabilities of the investigational dose levels.
* ```ncohort```: The total number of cohorts.
* ```cohortsize```: The cohort size.
* ```maxt```: The maximum follow-up time.
* ```prior.p```: A vector of length 3, which specifies the prior probability that the time to toxicity lies inside the time interval ```(0,maxt/3)```, ```(maxt/3,2maxt/3)```, ```(2maxt/3,1)```. The default value is ```prior.p=c(1/3,1/3,1/3)```. 
* ```accrual```: The accrual rate, i.e., the number of patients accrued in 1 unit of time.
* ```dist1```: The underlying distribution of the time to toxicity outcomes; ```dist1=1``` is the uniform distribution, ```dist1=2``` corresponds to the Weibull distribution, ```dist1=3``` is the log-logistic distribution.
* ```dist2```: The underlying distribution of patient arrival time; ```dist2=1``` is the uniform distribution, ```dist2=2``` is the exponential distribution
* ```alpha``` A number from (0,1) that controls ```alpha*100%``` events in ```(0, maxt/2)```. The default is ```alpha=0.5```.             
* ```n.earlystop```: The early stopping parameter. If the number of patients treated at the current dose reaches ```n.earlystop```, stop the trial and select the MTD based on the observed data. The default value ```n.earlystop=100``` essentially turns off this type of early stopping.
* ```startdose```: The starting dose level for the trial, e.g., ```startdose<-1```.
* ```p.saf```: The lower bound of the proper dosing interval, e.g., ```p.saf=target-0.05```.
* ```p.tox```: The upper bound of the proper dosing interval, e.g., ```p.saf=target-0.05```.
#' @param cutoff.eli the cutoff to eliminate an overly toxic dose for safety.
#'                  We recommend the default value of (\code{cutoff.eli=0.95}) for general use.
#' @param extrasafe set \code{extrasafe=TRUE} to impose a more stringent stopping rule
#' @param offset a small positive number (between 0 and 0.5) to control how strict the
#'               stopping rule is when \code{extrasafe=TRUE}. A larger value leads to a more
#'               strict stopping rule. The default value \code{offset=0.05} generally works well.
#' @param ntrial the total number of trials to be simulated.

# Example
We consider use the TITE-keyboard design as an illustration. 
* The following code can reproduce Table 1 given in the paper with a target toxicity rate of 30%. 
```rscript
get.boundary.tite(target=0.3,ncohort=4,cohortsize=3, design=1) 
```
The output is given by 
```rscript
$boundary_tab
   #Patients #DLTs observed #Pending patients Ecalation      Stay               De-escalation          
1  3         0              < 2               Yes            No                 No                     
2  3         0              > 1               Suspend        Suspend            Suspend                
3  3         1              0                 No             Yes                No                     
4  3         1              1 ~ 2             No             ESS>= 2.88         ESS< 2.88              
5  3         2              < 2               No             No                 Yes                    
6  3         3              0                 No             No                 Yes & Eliminate        
7  6         0              < 7               Yes            No                 No                     
8  6         1              < 2               Yes            No                 No                     
9  6         1              2 ~ 3             ESS> 4.07      ESS<= 4.07         No                     
10 6         1              4 ~ 5             ESS> 4.07      ESS in [2.88,4.07] ESS< 2.88              
11 6         2              0                 No             Yes                No                     
12 6         2              1 ~ 4             No             ESS>= 5.75         ESS< 5.75              
13 6         3              < 4               No             No                 Yes                    
14 6         > 3            < 3               No             No                 Yes & Eliminate        
15 9         0              < 10              Yes            No                 No                     
16 9         1              < 5               Yes            No                 No                     
17 9         1              5 ~ 6             ESS> 4.07      ESS<= 4.07         No                     
18 9         1              7 ~ 8             ESS> 4.07      ESS in [2.88,4.07] ESS< 2.88              
19 9         2              0                 Yes            No                 No                     
20 9         2              1 ~ 3             ESS> 8.15      ESS<= 8.15         No                     
21 9         2              4 ~ 7             ESS> 8.15      ESS in [5.75,8.15] ESS< 5.75              
22 9         3              0                 No             Yes                No                     
23 9         3              1 ~ 6             No             ESS>= 8.63         ESS< 8.63              
24 9         4              < 6               No             No                 Yes                    
25 9         > 4            < 5               No             No                 Yes & Eliminate        
26 12        0              < 13              Yes            No                 No                     
27 12        1              < 8               Yes            No                 No                     
28 12        1              8 ~ 9             ESS> 4.07      ESS<= 4.07         No                     
29 12        1              10 ~ 11           ESS> 4.07      ESS in [2.88,4.07] ESS< 2.88              
30 12        2              < 4               Yes            No                 No                     
31 12        2              4 ~ 6             ESS> 8.15      ESS<= 8.15         No                     
32 12        2              7 ~ 10            ESS> 8.15      ESS in [5.75,8.15] ESS< 5.75              
33 12        3              < 4               No             Yes                No                     
34 12        3              4 ~ 9             No             ESS>= 8.63         ESS< 8.63              
35 12        4              0                 No             Yes                No                     
36 12        4              1 ~ 8             No             ESS>= 11.5         ESS< 11.5              
37 12        5, 6           < 8               No             No                 Yes                    
38 12        > 6            < 6               No             No                 Yes & Eliminate   
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

