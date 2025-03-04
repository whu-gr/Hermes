# Hermes
Official implementation for "Streamlining Data Transfer in Collaborative SLAM through Bandwidth-aware Map Distillation"

### Installation

- First you need to install [Gurobi](https://www.gurobi.com/) in your machine
- Then follow the [ORB-SLAM3](https://github.com/UZ-SLAMLab/ORB_SLAM3) installation guidelines to build Hermes front-end, remember to include the Gurobi path in your FindGUROBI.cmake.

### Running Hermes

- Currently, you can follow the workflow of [ORB-SLAM3](https://github.com/UZ-SLAMLab/ORB_SLAM3) to run the front-end odometry of Hermes.

### More Information

This code includes only the agent-side components of the paper, specifically KeyFrame Designation and Map Distillation. Due to the strong dependencies on communication protocols and hardware, we are currently refining this part of the code to ensure easier deployment across different machines. Additionally, we are working on streamlining the entire project to enhance code readability. We plan to release the full codebase along with detailed documentation by April, 2025.
=======
Official implementation for "Streamlining Data Transfer in Collaborative SLAM through Bandwidth-aware Map Distillation"
