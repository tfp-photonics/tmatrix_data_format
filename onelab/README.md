## Installation

### Prerequisites

1. **OneLab Package**
   - Download the OneLab package from [here](https://onelab.info/#Download).

2. **Python Environment**
   - Ensure you have a Python environment set up. You can use `virtualenv` or `conda` to create an isolated environment.

3. **Treams Package**
   - Install the `treams` package in your Python environment. You can do this via pip:
     ```bash
     pip install treams
     ```

## Usage

### Initial Setup

1. **Launch the App**
   - Open the Gmsh application in interactive mode.

2. **Open a Model**
   - Go to the `File` menu, select `Open`, and choose a GetDP `.pro` file, for example, `models/Magnetometer/magnetometer.pro`.

3. **Run the Model**
   - Press `Run`. This step will fix the configuration files on your computer and allow you to verify that everything is functioning correctly with the results of this calculation.

### Running the Computation

After verifying this setup, you can run the computation using the command line.

1. **Execute the Driving Python File**
   - Run the compute file in the bash shell from its location directory:
     ```bash
     <path_to_onelab_installation>/./gmsh sphere_onelab.py -
     ```
