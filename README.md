# Introduction

As the COVID-19 pandemic impacts global health in an unprecedented way, how to effectively control the spread of the disease at school and then make safe reopening plans has become a prime concern of local governments and school officials. We thus propose a contact network agent-based model (CN-ABM) to simulate different disease transmission scenarios on campus at the micro-scale. The model establishes a contact network for each agent, and then evaluates the change of the agent's health status via interpersonal contacts in different activity environments. Based on the CN-ABM, we identify how different community risk levels, teaching modalities, and vaccination rates could shape the epidemic curve.The proposed CN-ABM model and its applications to a real-world campus scenario lay the methodological foundation for risk assessments of COVID-19 at the micro-scale and can provide proactive intervention strategies to safely reopen schools. 

**CN-ABM** is short for **contact network agent-based model**. The CN-ABM model is developed to simulate the cases and epi curves of COVID-19 infection at the micro-scale. The data and codes can be accessed on this Github repo.

# Methodology

The transmission process of the SARS-CoV-2 virus can be regarded as the transition of an individual's health status via interactions with infectious agents. This process can be illustrated by [Figure 1](https://github.com/xic19022/cnabm/blob/figures/contact_network.png), where an agent’s health status could change if the agent’s contact network entails an infectious agent. Additionally, the internal transmission process inside a school can be influenced by the community spread beyond the school environment. We thus incorporate the community risk as an external node into the contact network, meaning that each agent can become infectious at a given probability even without contacting other internal agents.

![Figure 1](https://github.com/xic19022/cnabm/blob/figures/contact_network.png)

<em>Figure 1. Schematic illustration of the contact network</em> 

To construct the contact network, all agents’ daily activity patterns, including their movement trajectories, locations of stay in an activity environment (e.g., dining hall, residential building), and periods of stay, must be acquired. 

After the contact network is established, we use a modified SEIR model with six heath statuses to illustrate a complete infection cycle for an agent, including susceptible (S), exposed (E), pre-symptomatic (Ip), infected (I), recovered/removed (R), and vaccinated (V). These six health statuses constitute five infection phases, as shown in [Figure 2](https://github.com/xic19022/cnabm/blob/figures/SEIR.png).

![Figure 2](https://github.com/xic19022/cnabm/blob/figures/SEIR.png)

<em>Figure 2. An infection cycle with the change of an agent’s health status</em>

The workflow of the simulation is shown in [Figure 3](https://github.com/xic19022/cnabm/blob/figures/workflow.png). For each scenario, we have performed the simulation for twenty-five weeks at a one-minute timestamp; we also repeated each simulation five times to account for the stochastic nature of the disease transmission. 

![Figure 3](https://github.com/xic19022/cnabm/blob/figures/workflow.png)

<em>Figure 3. Workflow of the CN-ABM model</em>

# Study Area

Our study area is a university campus located in Southern China, as shown in [Figure 4](https://github.com/xic19022/cnabm/blob/figures/study_area.png). We consider four different activity environments: dining hall, lecture hall (\*including library), residential building, and outdoors. These activity environments have different infection rates. Other buildings (e.g., administrative buildings) are excluded from the simulation.

![Figure 4](https://github.com/xic19022/cnabm/blob/figures/study_area.png)

<em>Figure 4. Building footprints in the study area</em>

# Data

* <em>Personal_attribute.csv</em> contains the personal attributes, including: residential building id, room id, school, major, vaccination, incubation period, time of S→E, time of Ip→I, infection period, time of I→R, and health level. Personal information is generated based on the campus portal data. Health level is generated based on the lognormal distribution and normal distribution.

* <em>Building_attribute.csv</em> contains the building attributes, including: building id, building type, and X-Y coordinates. The coordinates are derived from the [Baidu Map application programming interface (API)](https://lbsyun.baidu.com/) and have been projected.

* <em>Trajectory_location.mat</em> (in five zipped files) contains the trajectory of each agent at each timestamp (one minute).

* <em>Trajectory_building_id.mat</em> contains activity locations along the agents' trajectories, including (X-Y) in the <em>Trajectory_location.mat</em> and "building id" in <em>Building_attribute.csv</em>.

# Codes

One MATLAB code is developed for the CN-ABM model.

* <em>simulation.m*</em> contains the main codes for the model simulation. It may take a few hours to run the model and generate the full set of results.

**Inputs:** Input data for running em>simulation.m*</em> include the personal attributes <em>Personal_attribute.csv</em>, building attributes <em>Building_attribute.csv</em>, agents' trajectories <em>Trajectory_location.mat</em> and the agents' actitivty building id <em>Trajectory_building_id.mat</em>. Also, serval parameters are needed for model initialization, including the number of simulation weeks, initial infetion rates, vaccine efficacy, vaccination rate, student composition, infection rate, and the average number of close contacts in each activity environment.

**Outputs:** The outputs of <em>simulation.m</em> is a csv file of active infectious cases(%).

We have also genereted an example of active infectious cases (%) with a low community risk, mostly in-class modality (25% of distance-learning), and a vaccination rate of 45%, as shown in <em>result_S25_V45.csv</em>.
