# Correlation Plot Feature

## Overview
The correlation scatter plot feature allows users to visualize correlations between two properties that have been evaluated in the script editor. Each point on the scatter plot represents a snapshot/frame of the trajectory.

## Features
- **Property Selection**: Dropdown menus to select X and Y axis properties from evaluated script properties
- **Multiple Series Support**: Automatically creates multiple series when array properties are selected
- **Interactive Plotting**: 
  - Hover over points to see frame number and property values in a tooltip
  - Click on any point to jump to that trajectory frame
- **Color-coded Series**: Each series gets a different color for easy distinction
- **Real-time Updates**: Updates when new properties are evaluated in the script editor

## Usage
1. Evaluate properties in the script editor (e.g., calculate distances, angles, etc.)
2. Open the Correlation window from the main menu
3. Select X-axis property from the dropdown
4. Select Y-axis property from the dropdown  
5. Click "Generate Plot" to create the scatter plot
6. Hover over points to see details
7. Click on points to jump to that frame in the trajectory

## Example Use Cases
- Plot distance vs angle correlations
- Visualize how different structural properties relate to each other
- Identify outlier frames by clicking on extreme points
- Study temporal evolution patterns across multiple properties

## Implementation Details
- Component: `src/components/correlation/correlation.cpp`
- Uses ImPlot for scatter plot rendering
- Integrates with the existing script evaluation system
- Supports both scalar and array properties
- Maintains interactivity with the main application's animation system