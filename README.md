# ChimeraX-EpitopeMapper

A ChimeraX bundle for epitope mapping, providing tools to identify and analyze potential epitopes on protein structures.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
  - [Command Description](#command-description)
  - [Example](#example)
- [Contributing](#contributing)
- [License](#license)

## Features

- **Epitope Scoring**: Calculates epitope scores based on Gaussian blur, protrusion scores, and B-cell epitope propensity scales.
- **Spatial Clustering**: Identifies conformational regions by clustering high-scoring residues in 3D space.
- **Robust Logging**: Provides detailed logging within the ChimeraX HTML log viewer.
- **Edge-Case Handling**: Gracefully handles missing protein atoms or density map data.

## Installation

To install the `ChimeraX-EpitopeMapper` bundle, follow these steps:

1. **Clone the repository**:
   ```sh
   git clone https://github.com/j-h-greene/ChimeraX-EpitopeMapper.git