## Setting Up the Environment

To set up the conda environment for this project, follow these steps:


1. Navigate to the project directory:
    ```sh
    cd meep
    ```

3. Create the conda environment from the `environment.yml` file:
    ```sh
    conda env create --name meep -f environment.yml
    ```

4. Activate the newly created environment:
    ```sh
    conda activate meep
    ```

## Running the Code

Once the environment is set up, you can run the main.py() script. It is possible to parallelize the simulations:
    ```
    mpirun -np N_of_processes python main.py
    ```