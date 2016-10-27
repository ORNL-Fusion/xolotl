This file gives instructions on how to use modelError scripts.

This series of python scripts requires an installation of the Uncertainty Quantification Toolkit (UQTk),
which is a library and tool collection developed in Sandia National Laboratory for UQ in numerical model predictions.
The download link is provided in http://www.sandia.gov/UQToolkit/.
	
There are five scripts here including a pre-processing and a post-processing steps.
The usage of all the scripts are discussed below.


----------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------
Step0_Preprocessing
----------------------------------------------------------------------------------------------------
This step takes any raw uncertain data, and prepares it for parametric fitting. 
Simpson's 1/3 model is used here to obtain a normalized data, by taking the integral of
y data, I, and divide all y values to I. It must be noted that depth (xdata.txt) 
must be in nanometers to be used in other steps.
----------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------
 
 
----------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------- 
Step1_ParametricFit 
---------------------------------------------------------------------------------------------------- 
There are three folders under Step1_ParametricFit for three different parametric equations; Gumbel, 
Weibull, and LogNormal. Under each of these folders, there are two python, one cpp scripts, and a Makefile. 
The first one to use is a cpp file, which calls UQTk functions to infer parameters of the fit equation; 
gumbel_infer.cpp, weibull_infer.cpp, and log_normal_infer.cpp. To build and run these cpp files, the first 
two lines of the Makefile must be arranged to locate the UQTk installation directory, and Fortran libraries', 
respectively (the fortran version must be the same as the one that you used when you installed UQTk). 
When, you type "make", an executable file will be generated; gumbel_infer, weibull_infer, or log_normal_infer.
If you execute this file (e.g. ./weibull_infer), two data files; xdata.txt and Normalized_ydata.txt 
will be read and fit into a parametric equation. Hence, the coefficients of the selected parametric 
equation will be inferred. The code also gives the designated standard deviations; pushed-forward 
posterior, posterior predictive, and model-error.

The first one of the python scripts is main.py. It has two jobs. First, it extracts the inferred parameters 
of the equations, and also the related errors. Then, it plots these errors as colored areas with the mean 
prediction, which is a red line. Second, the script sends the generated files under Step2 and 
Step 3 folders.

The second python script is plotchains.py. As its name implies, the code generates MCMC chain plots for 
the inferred parameters with the value of the parameter on vertical axis and MCMC step on horizontal axis.
----------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------


----------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------- 
Step2_SurrogateConstruction 
---------------------------------------------------------------------------------------------------- 
A surrogate model is developed in this step. A parameter domain file is generated, using the results of 
the previous step. The code calls Xolol using the generated domain. Then, it plots data-versus-model 
for surrogate accuracy assessment. Finally, the needed files are sent to the next step for forward 
propagation.

For this step, a few directories must be added to the bash, and exported (.bashrc (Linux) or .profile (MAC)). 
If you are using eclipse, you should configure your run by adding the following into the project's 
environment: PYTHONPATH, UQTK_INS, UQPCDIR, UQBINDIR, PATH, LD_LIBRARY_PATH, PAPI_PREFIX, MESA_PREFIX, 
EAVL_PREFIX.

To do this, click on Run > Run Configurations > Environment > New, and add each item. Select 
"append environment to native environment" option. 
----------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------


----------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------- 
Step3_ForwardPropagation 
---------------------------------------------------------------------------------------------------- 
Forward propagation of the initial uncertainties takes place in this step, and the results are saved in 
full_propagation.dat file. UQTk's model inference code for a non-inference, propagation only regime runs, 
and the result file is sent to Step 4 folder.

PYTHONPATH, UQTK_INS, UQPCDIR, and UQBINDIR, must be added to the environment as similar as in Step 2.
----------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------


----------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------- 
Step4_Postprocessing 
---------------------------------------------------------------------------------------------------- 
The results obtained in the previous steps are plotted. The code generates a plot with a mean prediction 
result and pushed-forward posterior, posterior predictive, and model-error standard deviations.
----------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------
