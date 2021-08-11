.. ENIGMA documentation master file, created by
   sphinx-quickstart on Wed Aug 11 07:18:25 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ENIGMA -  A fast and accurate deconvolution algorithm based on regularized matrix completion algorithm
==================================

.. image:: https://github.com/WWXkenmo/ENIGMA/blob/main/figure_git.jpg

Cell type-specific gene expression (CSE) brings novel insights into physiological and pathological processes compared with bulk tissue gene expression. Although fluorescence-activated cell sorting and single-cell RNA sequencing (scRNA-seq) are two widely used techniques to detect CSE, the constraints of cost and labor force make it impractical as a routine on large patient cohorts. 

Here, we present **ENIGMA**, an algorithm that using matrix completion to accurately deconvolute bulk RNA-seq into CSE matrices and cell type fraction matrices without the need of physical sorting or sequencing of single cells. We demonstrated the superior performance of ENIGMA to existing algorithms while requiring much less running time on both simulated and realistic datasets. ENIGMA could also reveal a monocyte to macrophage transition in arthritis patients and a beta cell-specific module associated with senescence and apoptosis in type 2 diabetes. Together, ENIGMA improves the CSE estimation and extends our understandings for diverse biological processes.

For more information, please refer to a preprint in `bioRxiv <https://www.biorxiv.org/content/10.1101/2021.06.30.450493v1>`_.

Reference
---------

.. code-block:: latex

    @article{Wang2021.06.30.450493,
      title={Improved estimation of cell type-specific gene expression through deconvolution of bulk tissues with matrix completion},
      author={Wang, Weixu and Yao, Jun and Wang, Yi and Zhang, Chao and Tao, Wei and Zou, Jiahua and Ni, Ting},
      journal={bioRxiv},
      year={2021},
      publisher={Cold Spring Harbor Laboratory}
    }

.. toctree::
   :maxdepth: 1
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
