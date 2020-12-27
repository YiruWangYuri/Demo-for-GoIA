# Practical Globally Optimal Consensus Maximization by Branch-and-Bound based on Interval Arithmetic

About 
===============
 
Consensus maximization is widely used in robust model fitting. Here, we achieve globally optimal consensus maximization by Branch-and-Bound framework and draw the idea of interval arithmetic-based bound calculation back on the map. We provide the detailed derivation of interval arithmetic-based bound calculation for consensus maximization problems with both linear and quasi-convex residuals. Extensive experiments show that the proposed method can better deal with larger number of data points and higher outlier ratios than existing global methods.

-----------------------------------------

Problem list in the demo
===============

**Linear**:
For a linear model, two kinds of constraints are employed to fix the unknown scale of model. (Linear regression and Unit Norm constriant) Hence, you can choose "GoIA-LR" and "GoIA-UN" to slove the linear problem. 
1. Synthetic linear model. 
2. Plane fitting.
3. Translation estimation. 
4. Decomposed homography estimation.

**Nonlinear**: 
5. Affine registration. 
6. Triangulation. 

-----------------------------------------
 
Getting Started
=============== 

1. Clone this repository.  
2. Run function "demo()" in MATLAB. 
 
------------------------

Contact
------------------------
 
Email: yiruwang18@fudan.edu.cn
 
 


 

 
 
