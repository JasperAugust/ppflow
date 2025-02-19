# Project Plan: Hybrid Wet Lab and AI Peptide Design Using PPFlow

## Overview

Integrate AI-driven peptide design with experimental wet lab validation to optimize and characterize peptide-receptor interactions. Utilize the PPFlow model to first validate the model's performance against published results and subsequently design and analyze novel peptides targeting specific receptors.

## Objectives

1. **Validation Phase**
   - Reproduce and verify PPFlow model results as presented in the research paper.
   
2. **Optimization and Characterization Phase**
   - Design optimized peptides targeting a chosen receptor.
   - Synthesize and experimentally validate the designed peptides.
   - Characterize the binding interactions and structural properties of the peptide-receptor complex.

## Project Phases

### Phase 1: Validation of PPFlow Model

**Duration:** 3 Months

**Tasks:**

1. **Setup Computational Environment**
   - Create and activate the `ppflow` Conda environment.
   - Install necessary packages as per the [Installation](#installation) section.

2. **Data Preparation**
   - Download the processed datasets from [Google Drive](https://drive.google.com/drive/folders/1ce5DVmZz0c-p3PKrGDQoU_C9MD3cWLNq).
   - Verify data integrity and organize directories as specified.

3. **Reproduce Published Results**
   - Run training scripts to replicate model configurations:
     ```bash
     python train_ppf.py
     python train_diffpp.py
     ```
   - Generate peptide samples using:
     ```bash
     python codesign_diffpp.py 
     python codesign_ppflow.py
     ```

4. **Compare Results**
   - Use evaluation scripts to assess model performance:
     ```bash
     python evaluation/eval_struct_seq.py
     python evaluation/eval_sidechain_acc.py
     python evaluation/eval_bind.py
     ```
   - Ensure alignment with published metrics (e.g., FoldX Energy, TM Score, Validity).

5. **Troubleshooting**
   - Address any discrepancies or errors by reviewing and adjusting code as necessary.
   - Utilize logs and metrics for debugging.

**Deliverables:**

- Successfully reproduced results matching the published paper.
- Documentation of any deviations or issues encountered during reproduction.

### Phase 2: Optimization and Characterization

**Duration:** 6 Months

**Tasks:**

1. **Select Target Receptor-Peptide Pair**
   - Choose a receptor of interest based on research goals.
   - Ensure availability of structural data for both receptor and native peptide.

2. **Design Optimized Peptides Using PPFlow**
   - Configure PPFlow with the selected receptor:
     ```bash
     python codesign_ppflow.py --config ./configs/test/optimize_ppflow.yml --tag target1
     ```
   - Generate multiple peptide candidates optimizing for binding affinity and specificity.
   
3. **Data Integration and Iteration**
   - Analyze experimental data to validate computational predictions.
   - Refine the PPFlow model based on experimental feedback to enhance design accuracy.

4. **Documentation and Reporting**
   - Compile comprehensive reports detailing methodologies, results, and interpretations.
   - Prepare data for publication or internal review.

**Deliverables:**

- A set of optimized peptides with validated binding properties.
- Structural data of peptide-receptor complexes.
- Comprehensive reports and potential publications.

## Timeline

| Phase                               | Duration   | Key Milestones                              |
|-------------------------------------|------------|---------------------------------------------|
| **Phase 1: Validation**             | 3 Months   | - Environmental setup<br>- Data preparation<br>- Result reproduction<br>- Comparison and documentation |
| **Phase 2: Optimization & Characterization** | 6 Months   | - Target selection<br>- Peptide design<br>- Synthesis<br>- Experimental validation<br>- Data integration<br>- Reporting |


## Risk Management

| Risk                                 | Mitigation Strategy                          |
|--------------------------------------|-----------------------------------------------|
| **Computational Errors**            | Regular code reviews and debugging sessions. |
| **Peptide Synthesis Delays**        | Establish multiple synthesis partnerships.   |
| **Discrepancies Between AI and Wet Lab Results** | Iterate model based on feedback; enhance data integration. |
| **Resource Limitations**            | Secure necessary funding and resources early.|
| **Data Integrity Issues**           | Implement strict data validation protocols.  |

## Installation

### Create and Activate Conda Environment

```bash
conda create -n ppflow python==3.9
conda activate ppflow
```

### Install Basic Packages

```bash
# Install requirements
pip install -r requirements.txt

pip install easydict
pip install biopython
# MMseqs2
conda install bioconda::mmseqs2

# Alternative: Open Babel and RDKit
conda install -c openbabel openbabel
conda install -c conda-forge rdkit

# Alternative for Visualization: py3dmol
conda install -c conda-forge py3dmol
```

### Packages for Training and Generating

#### Install PyTorch 1.13.1 with Compatible CUDA Version

```bash
# Install geomstats
conda install -c conda-forge geomstats

# Install torch-scatter
pip install torch-scatter -f https://data.pyg.org/whl/torch-1.13.0+cu117.html  

# OR: Stable torch-scatter
pip install ./temp/torch_scatter-2.1.1+pt113cu117-cp39-cp39-linux_x86_64.whl 
```

## Dataset

We provide the processed dataset of `PPBench2024` through Google Drive: [PPBench2024](https://drive.google.com/drive/folders/1ce5DVmZz0c-p3PKrGDQoU_C9MD3cWLNq), together with processed `PPDBench`.

### Structure

```plaintext
- data
    - processed
        - cluster_result_all_seqs.fasta
        - cluster_result_cluster.tsv
        - cluster_result_rep_seq.fasta
        - parsed_pair.pt
        - receptor_sequences.fasta
        - split.pt
    - processed_bench
        - cluster_result_all_seqs.fasta
        - cluster_result_cluster.tsv
        - cluster_result_rep_seq.fasta
        - parsed_pair.pt
        - receptor_sequences.fasta
        - split.pt
    - pdb_benchmark.pt
    - pdb_filtered.pt
```

To obtain raw datasets for preprocessing:

1. Download `datasets_raw.zip` from [Google Drive](https://drive.google.com/drive/folders/1ce5DVmZz0c-p3PKrGDQoU_C9MD3cWLNq).
2. Unzip to get:
    ```plaintext
    - dataset
        - PPDbench
            - 1cjr
                - peptide.pdb
                - receptor.pdb
            - 1cka
                - peptide.pdb
                - receptor.pdb
            ...      
        - ppbench2024
            - 1a0m_A
                - peptide.pdb
                - receptor.pdb
    ```

## Training and Generating

### Train PPFlow

```bash
python train_ppf.py
```

### Train DiffPP

```bash
python train_diffpp.py
```

### Generate Peptides After Training

```bash
python codesign_diffpp.py 
python codesign_ppflow.py
```

## Packages and Scripts for Evaluation

### Docking and Other Evaluation

#### Vina Docking

Install required packages:

```bash
conda install -c conda-forge vina
pip install meeko
pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3.git@aee55d50d5bdcfdbcd80220499df8cde2a8f4b2a
pip install pdb2pqr
```

Refer to `./tools/dock/vinadock.py` for the Python interface example.

#### HDock

1. Install `libfftw3`:
    ```bash
    apt-get install -y libfftw3-3
    ```

2. Download HDock from [HDock Lite](http://huanglab.phys.hust.edu.cn/software/hdocklite/).

3. Install/unzip to `./bin`:
    ```plaintext
    - bin
        - hdock
            - 1CGl_l_b.pdb
            - 1CGl_r_b.pdb
            - createpl
            - hdock
    ```

Refer to `./tools/dock/hdock.py` for the Python interface example.

#### Pyrosetta

1. Sign up at [Pyrosetta Downloads](https://www.pyrosetta.org/downloads).
2. After license authorization, install via:

    ```bash
    conda config --add channels https://yourauthorizedid:password@conda.rosettacommons.org 
    conda install pyrosetta   
    ```

#### FoldX

1. Register and download FoldX from [FoldX Suite](https://foldxsuite.crg.eu/foldx4-academic-licence).
2. Copy to `./bin` and unzip:
    ```plaintext
    - bin
        - FoldX
            - foldx
    ```

Refer to `./tools/score/foldx_energy.py` for the Python interface example.

#### ADCP

1. Obtain ADFRsuite from [ADFRsuite Downloads](https://ccsb.scripps.edu/adcp/downloads/).
2. Copy `ADFRsuite_x86_64Linux_1.0.tar` to `./bin`.

    ```plaintext
    - bin
        - ADFRsuite_x86_64Linux_1.0
            - Toos
                - CCSBpckgs.tar.gz
                - ...
            - ADFRsuite_Linux-x86_64_1.0_install.run
            - uninstall
    ```

3. Add to `PATH`:

    ```bash
    export PATH=/absolute/path/to/ppflow/bin/ADFRsuite_x86_64Linux_1.0/bin:$PATH
    ```

Refer to `./tools/dock/adcpdock.py` for the Python interface example.

#### TMscore

Install TMscore:

```plaintext
- bin
    - TMscore
        - TMscore 
        - TMscore.cpp
```

#### PLIP for Interaction Analysis

1. Install Open Babel:

    ```bash
    conda install -c conda-forge openbabel
    ```

2. Clone PLIP:

    ```bash
    cd ./bin
    git clone https://github.com/pharmai/plip.git
    cd plip
    python setup.py install
    ```

3. Create alias:

    ```bash
    alias plip='python3 /absolute/path/to/ppflow/bin/plip/plip/plipcmd.py'
    ```

Refer to `./tools/interaction/interaction_analysis.py` for the Python interface example.

## Citation

If our paper or the code in the repository is helpful to you, please cite the following:

```bibtex
@inproceedings{lin2024ppflow,
	author = {Lin, Haitao and Zhang, Odin and Zhao, Huifeng and Jiang, Dejun and Wu, Lirong and Liu, Zicheng and Huang, Yufei and Li, Stan Z.},
	title = {PPFlow: Target-Aware Peptide Design with Torsional Flow Matching},
	year = {2024},
	booktitle={International Conference on Machine Learning},
	URL = {https://www.biorxiv.org/content/early/2024/03/08/2024.03.07.583831},
	eprint = {https://www.biorxiv.org/content/early/2024/03/08/2024.03.07.583831.full.pdf},
}
```

## Additional Tips

- **Version Control:** Use Git to track changes and collaborate effectively.
- **Documentation:** Maintain clear documentation for both computational scripts and wet lab protocols.
- **Regular Meetings:** Schedule periodic check-ins to ensure alignment between AI and wet lab teams.
- **Data Backup:** Implement robust data backup solutions to prevent loss of critical data.
- **Ethical Considerations:** Ensure all experiments comply with ethical guidelines and institutional regulations.

---

Good luck with your peptide design project! Reach out if you encounter any issues or need further assistance.