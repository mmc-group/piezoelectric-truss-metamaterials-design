### Piezoelectric truss metamaterials: data-driven design and additive manufacturing

This repository contains the codes used for homogenization, ML-based desgin, and print path planning for the truss metamaterials presented in the manuscript:

**"Piezoelectric truss metamaterials: data-driven design and additive manufacturing"**

The corresponding dataset used in this work can be found under the following link: [Dataset underlying the publication "Piezoelectric truss metamaterials: data-driven design and additive manufacturing"](https://doi.org/10.5281/zenodo.17041483)

--- 

**Structure**
<pre>
│   LICENSE
│   README.md
│
├───data
│       Readme.md
│
├───deterministic-inverse-design
│   │   optimization.py
│   │   parameters.py
│   │   train.py
│   │
│   ├───data
│   │       Readme.md
│   │
│   └───src
│           errorAnalysis.py
│           lattice_utils.py
│           loadDataset.py
│           model_utils.py
│           normalization.py
│           voigt_rotation.py
│           __init__.py
│
├───generative-design
│   ├───data
│   │       Readme.md
│   │
│   ├───dataGeneration
│   │       generation.py
│   │       main.py
│   │       moveNodes.py
│   │
│   ├───fwdTrain
│   │       main.py
│   │       validation.py
│   │       __init__.py
│   │
│   ├───invOpt
│   │       invOpt_serial.py
│   │
│   ├───models
│   │       model.py
│   │       parameters.py
│   │       utils.py
│   │       __init__.py
│   │
│   └───results
│           best_e_model.pt
│           best_model.pt
│           e_model.pt
│           model.pt
│
├───print-path-opt
│   │   adjacency_table.m
│   │   Angle_bw_lines.m
│   │   bridge_check.m
│   │   checkPointOnSegment.m
│   │   directed_edges.m
│   │   DistBetween2Segment.m
│   │   edge_validation.m
│   │   endpointsofAonB.m
│   │   Fleurys_algorithm.m
│   │   gcode_write.m
│   │   GenEulerPath.m
│   │   intersection_point.m
│   │   is_valid.m
│   │   lattice_plot.m
│   │   linexline.m
│   │   Main.m
│   │   path_plot.m
│   │   RearrangeEdges.m
│   │   rotate_arbt_axis.m
│   │   sort_vertically.m
│   │   split_current_group.m
│   │   tessellate.m
│   │
│   ├───outputs
│   └───Unit_Cell_Coordinates
│           BaseUC7_UC.mat
│           BCC_loop_UC.mat
│           BCC_UC.mat
│           Bow_Tie_3D.mat
│           Cube_UC.mat
│           Octahedron_UC.mat
│           Octet_UC.mat
│
└───truss-homogenization-piezo
        ang_ax.m
        bluewhitered.m
        calc_rad.m
        corner_nodes.m
        dof.m
        edge_nodes.m
        effective_piezo_surface.m
        Effective_stiffness.m
        effective_surface_plot.m
        electric_field.m
        e_rotation_tensor.m
        face_nodes.m
        fill_vec_map.m
        fnalize_vec_map.m
        Generate_Lattice.m
        homogenization.m
        identify_face_nodes.m
        lattice_data.m
        lattice_vectors.m
        main.m
        nodesInit.csv
        PBC_mat.m
        plotcube.m
        ptb_mask.csv
        remove_redundancy.m
        rgsearch.m
        rotate_arbt_axis.m
        sanity_check_fem.m
        ScaleToUnit.m
        StiffnessMatrix.m
        truss_plotting.m
</pre>
#### Usage and purpose of the different codes:

# 1. Machine learning based optimization: 
   
## Design Space 1: 

### Overview

This repository provides a framework to inversely design truss metamaterials with specific piezoelectric properties using deep learning. It extends the approach introduced in ['Inverting the structure-property map of truss metamaterials via deep learning'](https://www.pnas.org/content/119/1/e2111505119) by enabling targeted optimization of piezoelectric truss structures.

### Key Scripts:
- `train.py`: Trains a neural network to learn the relationship between truss structures and their material properties.
- `optimization.py`: Uses the trained network to optimize truss metamaterials for specific piezoelectric properties.

### Training the Network (`train.py`)

The `train.py` script trains a deep learning model to predict the mechanical properties of truss metamaterials based on their structure. This script should be used if you need to retrain the network from scratch, e.g., to adapt the framework to new datasets or refine performance.

#### Running the Training Script
```sh
python train.py
```

### Optimizing Truss Metamaterials (`optimization.py`)

The `optimization.py` script utilizes the trained network to inversely design trusses that exhibit desired piezoelectric properties.

#### Running the Optimization Script
```sh
python optimization.py <objective>
```

#### Optimization Parameters:
- `<objective>`: Specifies the target property to optimize (e.g., maximizing or minimizing specific piezoelectric components). Refer to the code to use/add new objectives

#### Output
The optimized truss designs will be stored in the `results/` directory in `.csv` format. These files contain the optimized structural parameters required to achieve the specified properties.

### Requirements

- Python (tested on version 3.7+)
- Python packages:
  - PyTorch (CUDA-supported version recommended)
  - Pandas
  - NumPy
  - Matplotlib (optional, for visualization)



All codes were adapted from the original publication:

[https://github.com/jhbastek/InvertibleTrussDesign](https://github.com/jhbastek/InvertibleTrussDesign)


---


## Design Space 2: 

### Overview

This repository provides a generative modeling framework for constructing a unified design space for inverse design with truss lattices, as described in ['Unifying the design space and optimizing linear and nonlinear truss metamaterials by generative modeling'](https://www.nature.com/articles/s41467-023-42068-x). The framework is adapted from ['Zheng et al. (2023)'](https://www.nature.com/articles/s41467-023-42068-x) and enables both forward property prediction and inverse design optimization of truss structures with tailored piezoelectric properties.

#### Training the Network (`main.py`)

The `main.py` script trains the generative modeling framework and the property predictor using the provided dataset. This script should be run if you need to retrain the models from scratch.

#### Running the Training Script
```sh
python main.py
```

### Validating the Trained Model (`validation.py`)

The `validation.py` script loads the pre-trained models and performs the following tasks:
- Predicts the reconstructed truss structures and their corresponding effective piezelectric matrices for a test dataset.
- Randomly samples 1000 points from a multivariate Gaussian distribution and calculates the validity score.

#### Running the Validation Script
```sh
python validation.py
```

### Inverse Design (`invOpt.py`)

The `invOpt.py` script enables inverse design of truss structures with desired piezoelectric properties using gradient-based optimization.

#### Running the Inverse Design Script
```sh
python invOpt.py <objective> <seed>
```

#### Optimization Parameters:
- `<objective>`: Specifies the target piezoelectric matrix component to optimize (e.g., maximizing or minimizing specific components such as $e_{15}$, $e_{31}$, etc.)

#### Example Usage:
To optimize for maximizing the piezoelectric component \(e_{15}\) with random seed=1, run the following:
```python
python invOpt.py +e15 1
```

### Requirements

- Python (tested on version 3.10.4)
- Python packages:
  - PyTorch (tested on CUDA version 11.0)
  - Pandas
  - NumPy
  - SciPy


All codes were adapted from the original publication:
[https://github.com/li-zhengz/UnifyingTrussDesignSpace](https://github.com/li-zhengz/UnifyingTrussDesignSpace)



# 2. Homogenization
The codes used for homogenizing the truss unit cells for thier effective properties are contained in `truss-homogenization-piezo/`
   
- `main.m` --- Generates the dataset for ML training by homogenizing trusses of dataset I and dataset II for electromechanical properties.
- For dataset I, the unit cell geometries and topologies are generated by processing the input features strored in `data/input_DS1.mat`. 
- For dataset II, the geometry and topologies are stored in `data/nodes_DS2.mat`
and `data/connect_DS2.mat`.

---


# 3. Print Path Planning

Codes utilized for the graph-based print path planning of the lattices are contained in `print-path-opt/` (based on the graph theory-based approach of Weeks et al.)  

**Ref:** Weeks, R. D., Truby, R. L., Uzel, S. G. & Lewis, J. A. Embedded 3D printing of multimaterial polymer lattices via graph-based print path planning. *Adv. Mater.* **35**, 2206958 (2023).  

## Files and Their Functions  

- **`main.m`**  
  - Loads coordinates and the connectivity of the lattice to be printed from `print-path-opt/Unit_Cell_Coordinates/`.  
  - Generates `path` for the lattice with the sequence of nodes to be followed (includes both travel and printing motions).  

- **`gcode_write.m`**  
  - Can be used as a function with `main.m` to generate g-code from the variables `path` and `lattice_coords` or from a separate file containing path and coordinates.  
  - Creates a `.txt` file with g-code coordinates for travel and printing motions.  
  - Printing and extrusion parameters can be tuned to adjust for the printer specifications used.  
  - Saves g-codes to `print-path-opt/outputs/`.  

## Example Usage  

To generate the path for an octet lattice with a **2×2×2** tessellation:  

1. Load `Unit_Cell_Coordinates/Octet_UC.mat`.  
2. Set `tess = [2,2,2]`.  
3. Run `Main.m`.  
4. The coordinates of the points in the G-code format will be written to a `.txt` file and saved to `'print-path-opt/outputs/'`.  
   The **"Start G-code"** and **"End G-code"** lines, specific to the printer being used, can be added to the `.txt` file.  
   Finally, rename the file with a `.gcode` extension to use it for printing.  