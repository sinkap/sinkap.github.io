---
layout: post
title:  Gaussian Jordan Elimnation
date:   2012-10-22 11:34:20
categories: Mathematics
highlight: true
author: Kapileshwar Singh
image: /img/2012-10-22-gaussian-jordan-elimination/three-planes.png
---

This post is going to explain one of basic building blocks for solving "Linear Algebra Equations". Consider a set of equations:

$$
	\begin{align}
		a_{0,0}x_1 + a_{0,1}x_2 + \cdots + a_{0,M} = b_0 \\\\

		a_{1,0}x_1 + a_{1,1}x_2 + \cdots + a_{1,M} = b_1 \\\\

		\vdots \hspace{25mm} \vdots \\\\

		a_{N,0}x_1 + a_{N,1}x_2 + \cdots + a_{N,M} = b_N \\\\
	
	\end{align}

$$

This is a system on M unknowns \\(x_0, x_1 \cdots x_M \\) and N equations. Each variable can be thought of a degree of freedom and each equation can be thought of as a constraint. Think about a three variable situation, like a position of a person in a 3-D coordinate. Without any constraints, he has three degrees of freedom in the x, y and z direction. If we are given three equations describing his position(each equation in x, y and z represents a plane in 3-D), we can pin point his co-ordinates in the 3-D space.

![three-planes]

#### Validation

* If \\(M > N\\), the number of unknowns is greater than the number of equations, the system is said to be undetermined and has infinitely many solutions. The solution space can be restricted by Compressed Sensing.
* If \\(M < N\\), the number of equations are greater than the number of variables, the system is said to be overdetermined. Here the general approach is to find the best fit solution (i.e R.M.S error values are a minimum for all equations)
* If \\(M = N\\), the system is consistent if the following caveats are satisfied:

	* No row should be a linear combination of the other row, this leads to row degeneracy
	* If all the equations have a certain variable in the exact same linear combination, the system is afflicted by column degeneracy
* Both these equations effective result in the removal of  a constraint and thus the system becomes indeterminable.

#### Pivoting

In order to obtain more accurate results and reduce round-off errors, a technique called Pivoting is used. Pivoting is done to convert a matrix to its row echelon form.

#### What is row echelon form?

A matrix is said to be in row echelon form if:

* All non zero rows are above the zero rows.
* The first non zero number in a row from the left called the Leading coefficient or Pivot should be strictly to the right of the leading coefficient of row above it.
* All entries in a column below the leading coefficient must be zero

Here is an example of a matrix in row echelon form:

$$ \left[ \begin{array}{ccccc} 1 & a_0 & a_1 & a_2 & a_3 \\ 0 & 0 & 1 & a_4 & a_5 \\ 0 & 0 & 0 & 1 & a_6 \end{array} \right] $$

Pivoting can be done in two ways:

##### Partial Pivoting 

In this the algorithm selects element the largest absolute value and shuffles the rows in such a way that it lies along the diagonal

##### Complete Pivoting

The algorithm scans the whole matrix for the largest element and shuffles both columns and rows to place the pivot along a diagonal //(a_{ii} //)

#### The Algorithm

We will be using an example matrix to illustrate this Algorithm (which is given in the text-book **Numerical Recipes in C++**:

$$ A = \left[ \begin{array}{ccccc} 1 & 2 & 3 & 4 & 5 \\ 2 & 3 & 4 & 5 & 1 \\ 3 & 4 & 5 & 1 & 2 \\ 4 & 5 & 1 & 2 & 3 \\ 5 & 4 & 3 & 2 & 1\end{array} \right] $$

The equations we aim at solving are:

$$ 
	A \cdot Y = I \\
	A \cdot X_1 = B_1 \\
	A \cdot X_2 = B_2 \\
$$

The algorithm takes two inputs, matrix A (coefficient matrix) and B (solution vector). The inverse of the matrix is returned in A and the variable vector is returned in B.

##### Step 1: Finding the Pivot Element

In the first step the algorithm iterates through the matrix and finds the largest element, in the first iteration the pivot element is the largest element of the last row. In our case it comes out to be five and is in the fist column, so there is no need for a column swap, it only needs to be swapped with t he first row. This swap is maintained in a two book-keeping arrays storing the actual position of pivot, so that the result can be restored.

The next time the algorithm searches for a Pivot element, it excludes \\(R_1\\) and \\(C_1\\) from the search.

![pivot11]

![pivot21]

##### Step 2: Normalizing the row

Before we understand the first step we need to understand why this actually works. Using our transformations we are basically converting the matrix into the identity matrix I. Therefore,

$$if A = I; I \cdot X_1 = B_1^\prime \implies X_1 = B_1^\prime$$

Where \\(B_1^\prime\\) is the transformed solution vector

As we are using the equation \\(A \cdot Y = I\\) to determine the inverse of the matrix we store the result back in A.

This step can be further subdivided into two sub-steps:


* Normalize a given row by the Pivot element, So now our matrix equations looks like:

  $$ A = \left[ \begin{array}{ccccc} \frac{5}{5} & \frac{1}{5} & \frac{2}{5} & \frac{3}{5} & \frac{4}{5} \\ 2 & 3 & 4 & 5 & 1 \\ 3 & 4 & 5 & 1 & 2 \\ 4 & 5 & 1 & 2 & 3 \\ 1 & 2 & 3 & 4 & 5\end{array} \right] \cdot Y = \left[ \begin{array}{ccccc} \frac{1}{5} & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 1\end{array} \right] $$

  The solution vector also gets transformed as:

  $$ \left[ \begin{array}{c} \frac{b_0}{5} \\ b_1 \\ b_2 \\ b_3 \\ b_4\end{array} \right]$$

* The next step is to reduce each element below the Pivot by subtracting the right amount of first row:

  $$ A = \left[ \begin{array}{ccccc} 1 & 0.2 & 0.4 & 0.6 & 0.8 \\ 0 & 2.6 & 3.2 & 3.8 & -0.6 \\ 0 & 3.4 & 3.8 & -0.8 & -0.4 \\ 0 & 4.2 & -0.6 & -0.4 & 0.2 \\ 0 & 1.8 & 2.6 & 3.4 & 4.2\end{array} \right] \cdot Y = \left[ \begin{array}{ccccc} 0.2 & 0 & 0 & 0 & 0 \\ -0.4 & 1 & 0 & 0 & 0 \\ -0.6 & 0 & 1 & 0 & 0 \\ -0.8 & 0 & 0 & 1 & 0 \\ -0.2 & 0 & 0 & 0 & 1\end{array} \right]$$

and similar transforms are performed on the solution vector

We, will discuss certain parts of the second iteration as they are slightly different from the first:

![pivot]

Now while iterating for the second column, the largest element found its at \\( R_4,C_4\\)

Here, there is no need for swapping as the pivot is found along the diagonal itself.
At the end we have done pivoting for all columns and have reduced our matrix, but we need to accommodate for the shuffling that we have done. Let us say that our book-keeping arrays are:

![arrays]

Let us take the first case:

As the row and column number was not the same, there is an initial swap that needs to restored back. So, we swap \\(C_4\ with\ C_0\\). A row operation in the input appears as a column operation in its inverse a (explains the shuffling of columns instead of rows)

[pivot]: /img/2012-10-22-gaussian-jordan-elimination/pivot.png
{: .image-quarter }
[arrays]: /img/2012-10-22-gaussian-jordan-elimination/arrays.png
{: .image-quarter }
[three-planes]: /img/2012-10-22-gaussian-jordan-elimination/three-planes.png
{: .image-quarter }
[pivot11]: /img/2012-10-22-gaussian-jordan-elimination/pivot11.png
{: .image-quarter }
[pivot21]: /img/2012-10-22-gaussian-jordan-elimination/pivot21.png
{: .image-quarter }
