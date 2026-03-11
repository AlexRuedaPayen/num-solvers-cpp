<h1>Numerical Solvers – Conjugate Gradient and Schwarz Preconditioning</h1>

<p>
<strong>Project written in April–May 2019</strong><br>
Course project carried out under the supervision of 
<strong>Prof. Xavier Claeys</strong> and <strong>Prof. Bertrand Thierry</strong>, 
as part of the course <em>UE 4M053</em> at <strong>Sorbonne University (Université Paris VI)</strong>.
</p>

<p>
This project was developed jointly with <strong>Nathan Cohen</strong>.
</p>

<hr>

<h2>Project Overview</h2>

<p>
This program aims to compare the efficiency of the <strong>Schwarz preconditioner</strong> 
when applied to the <strong>Conjugate Gradient algorithm</strong>.
The Schwarz preconditioners used in this project are selected randomly and are not based on a strong theoretical framework.
</p>

<p>
To perform the required inversions on small matrices, we compute their 
<strong>LU decomposition</strong> and then solve the associated triangular systems.
</p>

<p>
Matrices are stored in three different formats:
</p>

<ul>
<li><strong>Dense format</strong> (class <code>DenseMatrix</code>)</li>
<li><strong>Sparse Dictionary format</strong> (class <code>SparseMatrix</code>)</li>
<li><strong>CSR – Compressed Sparse Row format</strong> (class <code>CSR</code>)</li>
</ul>

<p>
However, only the dense format is used for the large matrices included in the examples.
One Laplacian matrix of size 100 is used in the <code>SparseMatrix</code> format.
</p>

<p>
The implementations of the Sparse Dictionary and CSR classes were written with the goal of optimizing
their use in the contexts where they are most appropriate.
</p>

<hr>

<h2>Project Structure</h2>

<p>
We created a virtual class <code>SolveMatrix</code>, which serves as the parent class of the three matrix format classes.
This class contains implementations of the iterative solvers. Since these solvers share identical code
across matrix formats, this approach avoids code duplication.
</p>

<p>
Block matrices are defined through template classes that take one of the three matrix formats as an argument.
Mutual friend relationships were defined between the matrix classes.
</p>

<hr>

<h2>Class: SolveMatrix</h2>

<p>
<code>SolveMatrix</code> is an abstract (virtual) class designed to group together a set of methods
that can be shared by its derived classes.
</p>

<p>
This class also defines several pure virtual methods that are used inside the shared algorithms.
Their implementations are provided in each derived class.
This allows the algorithms to be written once while still operating on the specific matrix implementation.
</p>

<p>
The solver methods <code>cg</code> and <code>pcg</code> are marked as <strong>private</strong>
to prevent improper usage. They can only be accessed through the method <code>Question</code>.
</p>

<p>
The <code>MinRes</code> solver is public and is used in the <code>Block</code> class.
</p>

<p>
The implementations derived from the pure virtual method <code>pcg</code> are located in 
<code>block_code.hpp</code>, which allows access to functions defined in the template class 
<code>BlockMatrix&lt;Matrix&gt;</code>.
</p>

<p>
Note: Ideally, the <code>SolveMatrix</code> class should have a virtual destructor.
Since memory usage is limited in this project and due to limited experience with such implementations,
this was left unimplemented.
</p>

<hr>

<h2>Matrix Formats</h2>

<p>
Three matrix storage formats were implemented. In all cases, matrix dimensions are fixed before construction.
A private mutator was added to avoid invalidating previously computed LU decompositions.
</p>

<h3>DenseMatrix</h3>

<p>
The most intuitive format. Matrix coefficients are stored in a static array.
</p>

<h3>SparseMatrix (Dictionary format)</h3>

<p>
Only non-zero coefficients are stored in a dictionary that associates each pair
<code>(row, column)</code> with its value.
</p>

<h3>CSR (Compressed Sparse Row)</h3>

<p>
Non-zero coefficients are stored along with their corresponding column index in arrays
associated with each row. This format is efficient when each row contains only a small
number of non-zero elements.
</p>

<p>
Each matrix class also stores pointers to matrices representing the LU decomposition.
These pointers are initialized by a private method <code>LU1()</code> and default to <code>nullptr</code>.
This allows reuse of the decomposition without recomputing it.
</p>

<p>
The LU decomposition assumes that the matrix has no zero pivots.
Handling general pivoting would require implementing a permutation class,
which was outside the scope of this project.
</p>

<p>
Note: LU decomposition for the CSR format currently does not work correctly.
The issue seems related to the use of <code>erase</code> on dynamically allocated lists.
The code was left as a possible starting point for further investigation.
</p>

<hr>

<h2>Block and BlockMatrix Classes</h2>

<p>
These classes are templates that take a matrix class as an argument.
They rely on the matrix classes sharing the same interface (e.g., <code>MvProd</code>).
</p>

<p>
The LU decomposition defined in <code>SolveMatrix</code> is used when solving linear systems
inside the method <code>MvProdInv</code>.
</p>

<p>
Since the convergence of the <strong>MinRes</strong> method depends on the eigenvalues of the matrix
rather than its dimension, it is reasonable to use it for large matrices,
while LU-based solvers are more efficient for small matrices.
</p>

<p>
Small block matrices are used as preconditioners.
This provides reasonable efficiency while keeping inversion costs low,
thanks to the reuse of stored LU decompositions.
</p>

<hr>

<h2>Main Function</h2>

<p>
A constructor was implemented to read matrices from input files using an <code>istream</code>.
The list of matrix files is stored in the file <code>MatrixFiles</code>.
</p>

<p>
Adding or removing matrices from the input folder automatically updates the execution workflow.
</p>

<p>
The program then calls the method <code>Question</code> on each matrix,
which compares the performance of the standard Conjugate Gradient method
and the preconditioned Conjugate Gradient method.
</p>

<p>
The method <code>Question2</code> compares the performance of the LU and MinRes solvers
(measured as execution time required to reach an error of approximately 10<sup>-6</sup>).
Results are written to <code>./output/LUvsMinRes.json</code>.
</p>

<hr>

<h2>ExecuteThisFile Script</h2>

<p>
The script <code>ExecuteThisFile</code> asks the user whether the JSON outputs should be recomputed
(<code>Y</code> or <code>N</code>).
</p>

<p>
It then:
</p>

<ul>
<li>Compiles the project using <code>make</code> / <code>Makefile</code></li>
<li>Runs <code>./bin/main</code></li>
<li>Generates plots using <code>matplotlib</code></li>
<li>Saves results in the <code>Output</code> folder</li>
</ul>

<p>
Running this script produces the results required for the final question of exercise 2.
</p>

<p>
Total computation time is approximately <strong>5 minutes</strong>,
with matrices 8, 14, and 15 being particularly slow.
</p>

<p>
To speed up execution, the call to <code>Question2</code> in <code>main.cpp</code> can be commented out.
</p>

<hr>

<h2>Testing</h2>

<p>
Most tests were performed by temporarily making methods public,
adding debugging outputs such as <code>std::cout &lt;&lt; "ok"</code>,
and checking for segmentation faults or suspicious execution times.
</p>

<p>
We also used the debugger <strong>gdb</strong> when necessary.
</p>

<hr>

<h2>Preconditioner Choice</h2>

<p>
Different parameters were tested for the block lists (chosen randomly).
</p>

<p>
Using a large number of small lists appeared to give the best performance.
To reduce both the number of vectors and matrix sizes, we chose:
</p>

<p>
<code>l = log₂(n)</code> and <code>r = l − 1</code>
</p>

<p>
This choice was made empirically and does not rely on a formal mathematical justification.
</p>
