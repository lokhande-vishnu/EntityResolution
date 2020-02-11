# Accelerating Column Generation via Flexible Dual Optimal Inequalities with Application to Entity Resolution

## Abstract
In this paper, we introduce a new optimization approach to Entity Resolution. Traditional approaches tackle entity resolution with hierarchical clustering, which does not benefit from a formal optimization formulation. In contrast, we model entity resolution as correlation-clustering, which we treat as a weighted set-packing problem and write as an integer linear program (ILP). In this case sources in the input data correspond to elements and entities in output data correspond to sets/clusters. We tackle optimization of weighted set packing by relaxing integrality in our ILP formulation. The set of potential sets/clusters can not be explicitly enumerated, thus motivating optimization via column generation. In addition to the novel formulation, we also introduce new dual optimal inequalities (DOI), that we call flexible dual optimal inequalities, which tightly lower-bound dual variables during optimization and accelerate column generation. We apply our formulation to entity resolution (also called de-duplication of records), and achieve state-of-the-art accuracy on two popular benchmark datasets.

## Paper
Link to arxiv is https://arxiv.org/abs/1909.05460

## Other Project particulars
Slides - "slides_AAAI20.pdf" 
Poster - "poster_AAAI20.pdf"
AAAI Paper - "AAAI-LokhandeV.2787.pdf"
