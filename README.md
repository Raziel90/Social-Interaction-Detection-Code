### Social-Interaction-Detection-Code
Code Related to the Roman paper: "Automatic detection of human interactions from RGB-D data for social activity classification"

If you use the dataset or the code for your research, please cite our RO-MAN 2017 that describes the data collection in detail
```
@article{coppola2017automatic,
  title={Automatic detection of human interactions from RGB-D data for social activity classification},
  author={Coppola, Claudio and Cosar, Serhan and Faria, Diego R and Bellotto, Nicola and others},
  year={2017},
  publisher={IEEE}
}
```


The code relies on the usage of the **UoL 3D Social Interaction Dataset** described in 
https://lcas.lincoln.ac.uk/wp/research/data-sets-software/uol-3d-social-interaction-dataset/


The code is written in Matlab.
To be able to run it, the dataset must be downloaded and the code configured.

# Configuration

The **config** in the K2_Code folder has to be edited with the path of four elements:
1. **The folder where this repository is placed**
2. **Path to the data folder** ( in the dataset: [...]/uol_social_interaction_dataset/social_interaction_segmentation/extrated_data)
3. **Path to the csv annotation file**
4. **Path to where the output files are going to be placed**


# Run

To run the experiments, after the configuration, run the following code files:

```
k2_init 
k2_hmm_full_crossvalidation
```

Take note that the parameters of the experiments can be configured in **k2_hmm_full_crossvalidation.m**

