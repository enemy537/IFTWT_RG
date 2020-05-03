# IFTWT_RG

## Description

Point Cloud segmentation techinique based on Watershed Transform implemented with IFT (seeds from Region Growing)

## Support

Build Platform           | Status
------------------------ | ------------------------------------------------------------------------------------------------- |
Ubuntu                   | [![Status][ci-ubuntu-16.04]][ci-latest-build] <br> [![Status][ci-ubuntu-18.04]][ci-latest-build]  |


## Requirements

GCC	 7.1 (supports OpenMP 4.5) *needed for parallelism
CMake	 3.8
PCL	 1.8
OpenCV	 3.2

## Citation
Please cite these papers in your publications if it helps your research.

@article{paiva2020historical,
  title={Historical Building Point Cloud Segmentation Combining Hierarchical Watershed Transform and Curvature Analysis},
  author={Paiva, Pedro VV and Cogima, Camila K and Dezen-Kempter, Eloisa and Carvalho, Marco AG},
  journal={Pattern Recognition Letters},
  year={2020},
  publisher={Elsevier}
}

Most of the point clouds used to test our method are also availabe on CHAS dataset. Check it out, is public!

@dataset{paiva_pedro_victor_vieira_de_2019_2609498,
  author       = {Paiva, Pedro Victor Vieira de and
                  Cogima, Camila Kimi and
                  Dezen-Kempter, Eloisa and
                  Carvalho, Marco Antonio Garcia de},
  title        = {{CHAS - Cultural Heritage Architectural 
                   Segmentation dataset}},
  month        = mar,
  year         = 2019,
  note         = {{This dataset was funded by São Paulo Research 
                   Foundation trough e-Science program grants:
                   \#2016/04991-0, \#2017/02787-9, \#2017/01237-5.}},
  publisher    = {Zenodo},
  version      = 1,
  doi          = {10.5281/zenodo.2609498},
  url          = {https://doi.org/10.5281/zenodo.2609498}
}
