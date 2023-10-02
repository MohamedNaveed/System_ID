# The Information-State Based Approach to Linear System Identification

This work considers the problem of system identification for linear systems. We propose a new system realization approach that uses an “information-state” as the state vector,
where the “information-state” is composed of a finite number
of past inputs and outputs. The system identification algorithm
uses input-output data to fit an autoregressive moving average
model (ARMA) to represent the current output in terms of finite
past inputs and outputs. This information-state-based approach
allows us to directly realize a state-space model using the
estimated ARMA or time-varying ARMA parameters for linear
time-invariant (LTI) or linear time-varying (LTV) systems,
respectively. The paper develops the theoretical foundation for
using ARMA parameters-based system representation using
only the concept of linear observability, details the reasoning
for exact output modeling using only the finite history, and
shows that there is no need to separate the free and the forced
response for identification. The proposed approach is tested on
various different systems, and the performance is compared
with state-of-the-art system identification techniques - Time-varying Observer Kalman Filter Identification (TV-OKID).

# Reference
M. N. G. Mohamed, R. Goyal, S. Chakravorty and R. Wang, "The Information-State Based Approach to Linear System Identification," 2023 American Control Conference (ACC), San Diego, CA, USA, 2023, pp. 301-306, doi: 10.23919/ACC55779.2023.10156137. Link: https://ieeexplore.ieee.org/abstract/document/10156137. Preprint: https://arxiv.org/abs/2211.10583 

Citation:
@INPROCEEDINGS{10156137,
  author={Mohamed, Mohamed Naveed Gul and Goyal, Raman and Chakravorty, Suman and Wang, Ran},
  booktitle={2023 American Control Conference (ACC)}, 
  title={The Information-State Based Approach to Linear System Identification}, 
  year={2023},
  volume={},
  number={},
  pages={301-306},
  doi={10.23919/ACC55779.2023.10156137}}


