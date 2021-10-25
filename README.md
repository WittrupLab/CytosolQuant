# CytoQuant

This software package can be used to estimate cytosolic release amounts of fluorescently labelled macromolecules and estimate single-cell dose-response relationsships. By combining quantifications with continuous measurements of a fluorescent reporter gene, dose-response relationships can be visualised and a mathematical model for cellular responses can be estimated. The software runs on MATLAB 2018B or later (requires Statistics and Machine Learning Toolbox). To run the examples below, set the “current folder” in MATLAB to the folder containing the script, then run the script. Exemplary data of cells expressing a destabilised eGFP treated with siGFP-AF647 lipoplexes are included. 

Example use:

code/siRNA/inputEstimationRNA.m

Fits three siRNA quantification models to cytosolic siRNA measurements. Press button during run for next cell. Set ‘buttonPressFlag = 0’ for uninterrupted run. Total run time is 10 min on standard computer, ~10 s for a single cell.

code/QuantilePlots/plotReleaseSortedeGFP.m

Plots eGFP responses for quartiles of cells treated with siGFP-1 and with quantification R2>0.75

code/eGFP/ploteGFPmeans.m

Plots eGFP downregulation model estimates relative to individual traces of eGFP expression with groups divided by release magnitude quintiles

# Acknowledgement

This work was carried out in the labs of Anders Wittrup and Jonas Wallin, Lund University, and supported by funding from the Swedish Society for Medical Research (SSMF). 


# Licence

Copyright 2021 Jonas Wallin, Anders Wittrup

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
