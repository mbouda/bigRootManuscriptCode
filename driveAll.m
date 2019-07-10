%This script holds code to drive all calculations for Big Root and Ohm's  
%law models using site data
%
%It was used on data from the US-Wrc Ameriflux site's 2010 drought to
%produce the case study in the Big Root manuscript. The underlying data are
%available from the site's PI. All results obtained with the present code
%are contained in this repository.
%
%Each function called in this script runs simulations for one of the
%sections of the manuscript, independently of the others.
%
%All code is provided on an as-is basis with no guarantee of producing any
%results and will require adaptation before running correctly on
%any particular system other than the one originally used.
%
% Martin Bouda, July 2019

%% Section 3.1 Model Inversion
%synthetic, single root datasets, inversions 
singleRootLibrary(); %produces a dataset for a single, synthetic root

%inversions of big-root model to this dataset:
fitCase1B(pwd);  %using xylem water potentials, flows, and collar water potential
fitCase1BG(pwd);  %using xylem water potentials, flows, and collar flow (gradient)
fitCase1BQR(pwd);  %using layer flows and collar water potential
fitCase1BGQR(pwd);  %using layer flows and collar flow (gradient).

%results collated and analysed by:
run postProcCase1b.m

% Reported and underlying results contained in:
load('errorTable.mat');
load('moments1b.mat');
    %full dataset of analytical solutions and fitted model predictions
    %(file named case1BFits.mat) not included in gitHub due to size (>8GB),
    %but similar raw results (different random seed) can be reproduced
    %using code provided.


%% Section 3.2.4 Model comparison

% Underlying dataset, constructed from US-Wrc Ameriflux data:
load('uswrcDroughtDataSet.mat');
    %The dataset does not include original soil moisture values, to which
    %the author does not have rights and which are available from the
    %site's PI 
    
 % Inversions of each model to larger, well-distributed dataset
[qrModPar,qrModSer,qrModBR]=baseLineOpt();

load('baseLineRes'); 
    %best-fit values of model parameters and resulting uptakes for all
    %layers and times


 % Inversions of each model using 1000, randomly selected subsets of 3
 % 1-Day periods each
[eBR,ePar,eSer,thModBR,thModPar,thModSer]=bootStrapOpt(); 

load('bootsX'); %contains all sets of model parameters 
                %fit by inversion to randomly chosen 3-day subsets of data
                %predicted flows, soil moistures and error distributions
                %can be generated from these with the code provided. 

%% SI 1: Sensitivity Analysis

sensitivity();
load('errorSens'); %Not yet in repo: too large?

%% SI 2: Fits of unconstrained RSA Stencils

unconstrained();
load('unconstrainedResults');


