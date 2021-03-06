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
* ```cutoff.eli```: The cutoff to eliminate an overly toxic dose for safety. We recommend the default value of (```cutoff.eli=0.95```) for general use.
* ```extrasafe```: Set ```extrasafe=TRUE``` to impose a more stringent stopping rule
* ```offset```: A small positive number (between 0 and 0.5) to control how strict the stopping rule is when ```extrasafe=TRUE```. A larger value leads to a more strict stopping rule. The default value ```offset=0.05``` generally works well.
* ```seed```: The seed for random number generation. 
* ```ntrial```: The total number of trials to be simulated.
* ```print```: To print out the boundary results.
* ```design```: the design indicator: ```design=1``` means the TITE-keyboard design, ```desing=2``` means the TITE-mTPI design. 

# Outputs
* ```get.boundary.tite()``` will generate a decision table that includes optimal dose escalation and deescalation boundaries.
* ```get.oc.tite()``` will return the operating characteristics of the time-to-event model-assisted design as a data frame,
including: (1) selection percentage at each dose level (```selpercent```);
           (2) the number of patients treated at each dose level (```nptsdose```);
           (3) the number of toxicities observed at each dose level (```ntoxdose```);
           (4) the average number of toxicities (```totaltox```);
           (5) the average number of patients (```totaln```);
           (6) the percentage of early stopping without selecting the MTD (```pctearlystop```); 
           (7) the average number of suspending times (```npend```);
           (8) the average trial duration needed for the trial based on the design (```duration```).

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
* Suppose the target toxicity rate is 30%, assessment window is 3 months, and patients are enrolled at the rate of 2 patients per month. A maximum of 36 patients will be recruited in cohorts of 3. Suppose 6 dose levels are considered and the true toxicity rates are (0.13, 0.28, 0.41, 0.50, 0.60, 0.70). The time-to-toxicity ourcomes are simulated from Weibull distributions by controlling that 50% of the toxicity events occur in the latter half of the assessment window and the patient accrual time is uniformly distributed. We can use the following code to reproduce the simulation results in Table 2 (scenario 1). 
```rscript 
target <- 0.3
p.true <- c(0.13,0.28,0.41,0.50,0.60,0.70)
get.oc.tite(target,p.true,ncohort=12,cohortsize=3,accrual=2,dist1=2,dist2=1,alpha=0.5,ntrial=5000,maxt=3)

-----------------------output------------------------
$selpercent
[1] 14.18 58.32 22.76  3.88  0.36  0.02

$npatients
[1] 11.9652 15.1296  6.9450  1.5702  0.2346  0.0168

$ntox
[1] 1.5652 4.2772 2.8722 0.7742 0.1386 0.0126

$totaltox
[1] 9.64

$totaln
[1] 35.8614

$percentstop
[1] 0.48

$duration
[1] 25.2315

$sdduration
[1] 3.220733
```

# Authors and Reference
* Ruitao Lin and Ying Yuan
* Lin, R. and Yuan, Y. (2019) “Time-to-event model-assisted designs for dose-finding trials with delayed toxicity”, Biostatistics, in press.

