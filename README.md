# Predict Transparent Semiconductors using Machine Learning
Transparent semiconductors may be the future for flat panel display!<br>

So what determine the transparency of the semiconductors? <br>
Short answer: The band gap. <br>
Then what determine the stability of the transparent semiconductors? <br>
Short answer: The formation energy. <br>

This [kaggle competition](https://www.kaggle.com/c/nomad2018-predict-transparent-conductors/) predicts the band gap and the formation energy of 600 semiconductors, given structural properties of 2400 semiconductors as the training set.

Key Results/Highlights:
1. Predicted bandgap energy is used as a feature to predict the formation energy, as they are highly correlated. <br>
2. CV score on formation energy prediction is increased by a few percentages when switched over to gradient boosting regression, due to the high bias (train score > CV score) proned by random forest regression.
3. Percentage of In and Al are important to determine the bandgap, while percentage of Ga is important to the formation energy.

Band gap and formation energy prediction on randomly selected test samples. <br>
Bandgap prediction accuracy = 95.0 % <br>
<img src=bandgap.png> <br>
Formation prediction accuracy = 89.4 % <br>
<img src=formation.png> <br>

For exploratory data analysis, head over to [exploratory_data_analysis.ipynb](exploratory_data_analysis.ipynb). <br>
For feature selection based on correlation, head over to [feature_correlation.ipynb](feature_correlation.ipynb). <br>
For machine learning application, head over to [machine_learning.ipynb](machine_learning.ipynb). <br>

Python source codes are in the [source folder](source). Inside the folder, simply run "python machine_learning.py". <br>
For features generation, run "python get_feature.py" will generate train5.csv and test5.csv. Both have been generated in the folder.
