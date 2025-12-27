# MERMAIDS-Project

This project focuses on the development of algorithms to detect earthquakes within time series data derived from underwater hydrophones.
It was carried out as part of the research term at Mines de Paris and supervised by Karin Sigloch (sigloch@geoazur.unice.fr) and Sébastien Bonnieux (bonnieux@geoazur.unice.fr) of the GéoAzur team at the CNRS campus in Sophia-Antipolis.

## Repository Structure

This repository contains the various codes used to complete this project, organized as follows:

### Data Analysis & Visualization

* **`get_data.ipynb`**: Reads the different recordings and saves the timestamps of the various registered earthquakes.
* **`info_events.csv`**: Contains all available information regarding the known earthquakes in the dataset.
* **`spectro.ipynb`**: Generates spectrograms to visually identify seismic events.

### Artificial Intelligence

* **`random_forest.ipynb`**: Creates the initial random forest model utilizing 2880 variables.
* **`light_random_forest.ipynb`**: A lightweight version of the random forest model that achieves similar performance with reduced complexity.

### Statistics

* **`sta_lta.ipynb`**: Implementation of the Short Term Average over Long Term Average (STA/LTA) algorithm to detect earthquakes.

## Data Availability

The raw data required to execute these notebooks is not public; it was provided to us by Sébastien Bonnieux (bonnieux@geoazur.unice.fr).

## Documentation

See *[rapport.pdf](https://github.com/enzo-boulin/MERMAIDS-Project/blob/main/rapport.pdf)* to understand the project workflow and details.

## Contact

For more information, please send an email to enzo.boulin@etu.minesparis.psl.eu.