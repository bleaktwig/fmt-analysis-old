# FMT Analysis
Code used for the FMT detector residual analysis.

### Useful Information
Reconstruction:
* Currently, the FVT engine handles FMT data.
* The engine grabs the DC tracks and reconstructs them, updating them with the FMT cluster data.
* Reconstruction works in a similar fashion to the DC's:
    * Clusters are made from the hits via a Cluster Finding algorithm.
    * Crosses are constructed from this clusters by grouping them.
    * The track is updated with these crosses via a Kalman Filter algorithm.

Plotting Residuals:
* Residuals are the difference between the DC track and the FMT clusters.
* Looking at the residuals gives us an idea of how to fix misalignments in the geometry.

### TODOs
* Separate the residual plots in the 4 FMT "regions" to improve the visibility of the data.
    * **Task assigned to Luis & Claudio, due 2020-05-01.**
* Plot the residuals as a function of energy.
    * **Task assigned to Luis & Claudio, due 2020-05-01.**
* Check if there's a better way to match tracks & clusters without biasing the distribution.
    * **Task assigned to Bruno & Jorge, due 2020-05-01.**
* Try small shifts to the FMT to improve alignment:
    * Shift the whole FMT detector by small nudges in x, y, z.
    * Shift each layer few millimeters.
    * **Task assigned to Bruno, due 2020-05-01.**
* Analyze discontinuities in the residual vs strip plots and try to find an explanation.
* Add timing info (Tmin) to the clusters.
    * **Task assigned to Bruno.**
