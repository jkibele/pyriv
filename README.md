# PyRiv

Minimum aquatic distance (MAD), the shortest path between two points that does not cross land, is a useful metric in the study aquatic organisms. However, MAD can be difficult and time consuming to calculate, particularly for anadromous species where both coastal and riverine distances must be considered. Researchers are often forced to choose between imprecise straight line substitutes and laborious manual path tracing techniques that don’t scale well to large data sets. PyRiv is a free and open source Python library that was created to address this problem. It uses a novel network graph approach to find MAD paths around complex coastlines, and it employs existing hydrography datasets to navigate river networks. By combining these methods, PyRiv can determine MAD between any two points whether they’re on the same river, different rivers, offshore, or any combination. Given a point shapefile as input, PyRiv can return a line shapefile representing MAD paths between all the points as well as a distance matrix in CSV format. Once a network has been prepared for a given area, PyRiv is simple to use and the calculations are fast. 

## Current Status

PyRiv is still under development. The existing code is functional, but messy. It'll be getting cleaned up over the next few months.
