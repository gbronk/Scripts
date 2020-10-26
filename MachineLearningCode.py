# The Code that Predicts the Future! The Future!
# Gabriel Bronk, PhD
# 2020

# This Python code loads patient data and uses it to predict whether patients
# will get sepsis. This code is inteded for teaching college or advanced
# high school students how to do machine learning (using the specific machine
# learning example of a generalized linear model). 


#1 - Load Python Modules:
import numpy
import pandas 
import statsmodels.api as StatsModelsModule
from scipy import stats
from matplotlib import pyplot 


#2 - Load Data, and split it into two tables called SepsisOutcomes and BodyMeasurements:
DataTable = pandas.read_csv("PatientData.txt",sep='\t')
SepsisOutcomes = DataTable[['SepsisOrNot']]
BodyMeasurements = DataTable.iloc[:,1:673:24]
#=================================================================================
# Student Challenge (to be done later):
# Write a "for loop" to loop through data from every time point.
#=================================================================================


#3 - Adjust Data, which will make analysis easier:
BodyMeasurementAverages = BodyMeasurements.mean(axis = 0)
BodyMeasurementStandardDeviations = BodyMeasurements.std(axis = 0)
AdjustedBodyMeasurements = (BodyMeasurements - BodyMeasurementAverages)/BodyMeasurementStandardDeviations


#4 - Add A Column of Ones To the AdjustedBodyMeasurements Table:
NumberOfPatients = len(AdjustedBodyMeasurements)
OnesArray = [1]*NumberOfPatients
AdjustedBodyMeasurements['Constant'] = OnesArray


#5 - Split the Data into Training Data and Test Data: 
N = round(0.8*NumberOfPatients)
BodyMeasurementsForTraining = AdjustedBodyMeasurements.iloc[0:N,:]
SepsisOutcomesForTraining = SepsisOutcomes.iloc[0:N,:]

BodyMeasurementsForTesting = AdjustedBodyMeasurements.iloc[N:(NumberOfPatients+1),:]
SepsisOutcomesForTesting = SepsisOutcomes.iloc[N:(NumberOfPatients+1),:]


#6 - Fit the Model to the Data:
OutputFromStatsModels = StatsModelsModule.GLM(SepsisOutcomesForTraining, BodyMeasurementsForTraining, family=StatsModelsModule.families.Binomial())
FittedModel = OutputFromStatsModels.fit()
print(FittedModel.summary())


###########
#7 - Use the Fitted Model to Predict the Outcome for Other Patients:
PredictedSepsisOutcomes = FittedModel.predict(BodyMeasurementsForTesting)


#8 - Determine How Accurate the Model's Predictions Are:
CorrectlyPredictedSickPeople = 0
CorrectlyPredictedHealthyPeople = 0
TotalNumberOfSickPeople = 0
TotalNumberOfHealthyPeople = 0

LengthOfTestingData = len(SepsisOutcomesForTesting)

for PatientNumber in numpy.arange(0,LengthOfTestingData):

    RealOutcome = SepsisOutcomesForTesting.iloc[PatientNumber][0]
    PredictedOutcome = PredictedSepsisOutcomes.iloc[PatientNumber]
     
    if RealOutcome == 1:
        TotalNumberOfSickPeople = TotalNumberOfSickPeople + 1
        if PredictedOutcome > 0.5:
            CorrectlyPredictedSickPeople = CorrectlyPredictedSickPeople + 1
    #=================================================================================
    # Student Challenge (to be done now):
    # Do the same kind of thing to count up the TotalNumberOfHealthyPeople
    # and the CorrectlyPredictedHealthyPeople
    #=================================================================================
           
Sensitivity = (CorrectlyPredictedSickPeople/TotalNumberOfSickPeople)*100
print(Sensitivity)
##Specificity = (CorrectlyPredictedHealthyPeople/TotalNumberOfHealthyPeople)*100
##print(Specificity)




