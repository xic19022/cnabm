# Introduction

As the COVID-19 pandemic impacts global health in an unprecedented way, how to effectively control the spread of the disease at school and then make safe reopening plans has become a prime concern of local governments and school officials. We thus propose a contact network agent-based model (CN-ABM) to simulate different disease transmission scenarios on campus at the micro-scale. The model establishes a contact network for each agent, and then evaluates the change of the agent's health status via interpersonal contacts in different activity environments. Based on the CN-ABM, we identify how different community risk levels, teaching modalities, and vaccination rates could shape the epidemic curve.The proposed CN-ABM model and its applications to a real-world campus scenario lay the methodological foundation for risk assessments of COVID-19 at the micro-scale and can provide proactive intervention strategies to safely reopen schools. 

# CN-ABM Data and Codes 

**CN-ABM** represents contact network agent-based model. The CN-ABM model is developed to simulate the cases of COVID-19 infection at the micro-scale level. 

The data and codes can be accessed on this Github repo.

## Data

[*Personal_attribute.csv*] contains the personal attribute variable for the model, including: residential building id, room id, school, major, vaccination, incubation period, time of S→E, time of Ip→I, infection period, time of I→R, and health_level.
The data of personal information is generated according to campus portal data. And data of health is generated based on lognormal distribution and normal distribution.

[*Building_attribute.csv*] contains the building attribute variable for the model, including: building id, building type, loacaltion of X and Y.
The coordinates are derived from the Baidu Map application programming interface (API) [https://lbsyun.baidu.com/] and have been projected.

[*Trajectory_location.mat*] contains the weekly trajectory of each agents in minute in the form of locations.
The data is generated based on agent-based model according to agent attributes, agent activity environment, and agent activity patterns.

[*Trajectory_building_id.mat*] is corresponding to locations ("X" and "Y") in [*Trajectory_localtion.mat*] and "building id" in [*Building_attribute.csv*].

## Codes

One MATLAB code is developed for the CN-ABM model.

[*simulation.m*] is the main code for model simulation. It may take a few hours to derive the full results.

**Inputs:** Input data for running [*simulation.m*] include the personal attributes(file: [*Personal_attribute.csv*]), building attributes(file: [*Personal_attribute.csv*]), trajectory in the form of location (file: [*Trajectory_location.mat*]) and building id (file:[*Trajectory_building_id.mat*]). Also, serval parameters are needed for model initialization, including number of simulation weeks, initial infetion rate, vaccine efficacy, vaccination rate, student composition, infection rate and  average number of close contacts in each activity environments.

**Outputs:** The outputs of [*simulation.m*] is a csv file of active infectious cases(%). 

We have also genereted an example of active infectious cases(%) with low community risk, students mostly in-class (25% of distance-learning),and vaccination rate of 45%, as shown in [*result_S25_V45.csv*].
