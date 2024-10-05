
# 🔧 FEM Code in Fortran

This repository contains a **Finite Element Method (FEM)** solver written in **Fortran**, designed for classical FEM analysis of 2D frames. The main program `main.py` reads an input file that defines the nodes, elements, materials, loads, and boundary conditions, and outputs the results of the FEM analysis.

## 📋 Input File Format

The input file must follow this structured format with the following sections:

1. **2D FRAME Header**
   Specifies the problem type and units:
   ```
   2D FRAME [m], [kN]
   ```

2. **NODES Section**
   Defines the number of nodes and their coordinates:
   ```
   NODES <NNODES> <IDNODE, X, Y>
   ```

3. **SECTIONS Section**
   Defines the cross-sectional properties for each section:
   ```
   SECTIONS <NSEC> <IDSEC, AREA, INERTIA, DEPTH, SHEAR>
   ```

4. **MATERIALS Section**
   Describes the material properties:
   ```
   MATERIALS <NMAT> <IDMAT, YOUNG, POISSON, THERMAL, WEIGHT>
   ```

5. **ELEMENTS Section**
   Specifies the connectivity of the elements:
   ```
   ELEMENTS <NELE> <IDELE, NODE1, NODE2, IDSEC, IDMAT>
   ```

6. **RESTRAINTS Section**
   Defines the boundary conditions (restraints) at specific nodes:
   ```
   RESTRAINTS <NRES> <NODE, DIRECTION>
   ```

7. **ELEMENT LOADS Section**
   Specifies the loads applied to elements:
   ```
   ELEMENT LOADS <NLOADS> <ELEMENT, PX, PY1, PY2, DTTOP, DTBOT>
   ```

## 🚀 How to Compile and Run the Code

1. **Compile the Fortran Code**:  
   Ensure you have a Fortran compiler installed. You can use `gfortran` to compile the code:
   ```bash
   gfortran -o fem_solver main.f90
   ```

2. **Run the Code**:  
   After compiling, you can run the program with your input file:
   ```bash
   ./fem_solver <input_file_path>
   ```

The program will process the input, perform the FEM analysis, and output the results.

## 📚 Future Work

- Add support for 3D frame analysis.
- Implement advanced material models.
- Improve visualization of results.

## 🧑‍💻 About

This project was developed during my studies in mechanical engineering, using Fortran to solve complex structural problems with FEM. Feel free to explore and contribute to the project!
