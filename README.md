# bargaining-distortion
Code to measure the efficiency of various bargaining methods.

duke.py uses the Euclidean distance metric.
dukeL22.py uses the Euclidean squared distance metric.
duke1.py uses Manhattan/Taxicab distance .

The distortion function takes the average distortion of a set of 2D points the user can input. 
The grid function takes the worst-case distortion of a set of points (of user-chosen size) on a square grid, also of user-chosen size.

Used in [Fain, B., Goel, A., Munagala, K., Prabhu, N.  "Random Dictators with a Random Referee: Constant Sample Complexity Mechanisms for Social Choice." Proceedings of the AAAI Conference on Artificial Intelligence. Vol. 33. 2019.](https://arxiv.org/abs/1811.04786)
